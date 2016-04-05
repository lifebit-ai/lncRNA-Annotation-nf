/*
 * Copyright (c) 2016, Centre for Genomic Regulation (CRG) and the authors.
 *
 *   This file is part of 'lncRNA-Annotation-NF'.
 *
 *   lncRNA-Annotation-NF is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   lncRNA-Annotation-NF is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with lncRNA-Annotation-NF.  If not, see <http://www.gnu.org/licenses/>.
 */

/* 
 * Main lncRNA-Annotation-NF pipeline script
 *
 * @authors
 * Sarah Djebali
 * Thomas Derrien
 * Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * Evan Floden <evanfloden@gmail.com> 
 */


params.name          ="lncRNA Annotation from Pig RNA-Seq"
params.genome        ="$baseDir/tutorial/genome/NEW_susScr102vega.fa"
params.annotation    ="$baseDir/tutorial/annotation/NEW_ensembl.83.vega.62.gtf"
params.reads         ="$baseDir/tutorial/reads/*_{1,2}.fastq.gz"
params.overhang      ='99'
params.threads       ='1'
params.output        ="results/"


log.info "lncRNA Annotation - N F  ~  version 0.1"
log.info "====================================="
log.info "name                   : ${params.name}"
log.info "genome                 : ${params.genome}"
log.info "annotation             : ${params.annotation}"
log.info "reads                  : ${params.reads}"
log.info "STAR overhang          : ${params.overhang}"
log.info "STAR threads           : ${params.threads}"
log.info "output                 : ${params.output}"
log.info "\n"


/*
 * Input parameters validation
 */

genomeFile             = file(params.genome)
annotationFile         = file(params.annotation) 

/*
 * validate input files/
 */
if( !annotationFile.exists() ) exit 1, "Missing annotation file: ${annotationFile}"

if( !genomeFile.exists() ) exit 1, "Missing genome directory: ${genomeFile}"

/*
 * Create a channel for read files 
 */
 
Channel
    .fromPath( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .map { path -> 
       def prefix = readPrefix(path, params.reads)
       tuple(prefix, path) 
    }
    .groupTuple(sort: true)
    .set { read_files } 


process index {
    input:
    file genomeFile
    file annotationFile
    
    output:
    file "STARgenome" into STARgenomeIndex
    file 'genome.length' into genomeLengths    

    script:
    //
    // STAR Generate Index and create genome length file
    //
    """
        mkdir STARgenome
        STAR --runThreadN ${params.threads} --runMode genomeGenerate --genomeDir STARgenome \
             --genomeFastaFiles ${genomeFile} --sjdbGTFfile ${annotationFile} \
             --sjdbOverhang ${params.overhang} --outFileNamePrefix STARgenome \

        fastalength ${genomeFile} > 'genome.length'
    """
}


process mapping {
    tag "reads: $name"

    input:
    file STARgenome from STARgenomeIndex.first()
    set val(name), file(reads:'*') from read_files

    output:
    set val(name), file("STAR_${name}") into STARmappedReads 

    script:
    //
    // STAR Mapper
    //
    """
        STAR --genomeDir ${STARgenome} --readFilesIn ${reads} --readFilesCommand zcat \
             --outFilterType BySJout --outSAMunmapped Within --outSAMtype BAM SortedByCoordinate \
             --outSAMattrIHstart 0 --outFilterIntronMotifs RemoveNoncanonical \
             --runThreadN ${params.threads} --quantMode TranscriptomeSAM --outWigType bedGraph --outWigStrand Stranded \
             --outFileNamePrefix ${name}
        mkdir STAR_${name}
        mv ${name}Aligned* STAR_${name}/.
        mv ${name}Signal* STAR_${name}/.
        mv ${name}SJ* STAR_${name}/.
        mv ${name}Log* STAR_${name}/.
    """
   
}


process cufflinks {

    input:
    file annotationFile
    set val(name), file(STAR_alignment) from STARmappedReads
    file genomeLength from genomeLengths.first()

    output:
    file "CUFF_${name}" into cufflinksTranscripts
     

    script:
    //
    // Cufflinks
    //
    """
        # Cufflinks
        mkdir CUFF_${name}
        cufflinks -p ${params.threads} -g ${annotationFile} -o CUFF_${name}  \
            --overlap-radius 5 --library-type fr-firststrand \
            ${STAR_alignment}/${name}Aligned.sortedByCoord.out.bam

        # Post-process: remove exons exceeding length of chromosome/scaffold
        cuff=CUFF_${name}/transcripts.gtf
        awk -v fileRef=${genomeLength} 'BEGIN{while (getline < fileRef >0){lg[\$2]=\$1}} \
            {nbex[\$12]++; line[\$12,nbex[\$12]]=\$0}\$4<1||\$5>lg[\$1]{ko[\$12]=1}END\
            {for(t in nbex){if(ko[t]!=1){for(k=1; k<=nbex[t]; k++){print line[t,k]}}}}'\
            \$cuff | awk -f ${baseDir}/bin/gff2gff.awk > cufflinks_ok.gtf
        cp cufflinks_ok.gtf CUFF_${name}/.
        """
}


// ===================== UTILITY FUNCTIONS ============================


/* 
 * Helper function, given a file Path 
 * returns the file name region matching a specified glob pattern
 * starting from the beginning of the name up to last matching group.
 * 
 * For example: 
 *   readPrefix('/some/data/file_alpha_1.fa', 'file*_1.fa' )
 * 
 * Returns: 
 *   'file_alpha'
 */
 
def readPrefix( Path actual, template ) {

    final fileName = actual.getFileName().toString()

    def filePattern = template.toString()
    int p = filePattern.lastIndexOf('/')
    if( p != -1 ) filePattern = filePattern.substring(p+1)
    if( !filePattern.contains('*') && !filePattern.contains('?') ) 
        filePattern = '*' + filePattern 
  
    def regex = filePattern
                    .replace('.','\\.')
                    .replace('*','(.*)')
                    .replace('?','(.?)')
                    .replace('{','(?:')
                    .replace('}',')')
                    .replace(',','|')

    def matcher = (fileName =~ /$regex/)
    if( matcher.matches() ) {  
        def end = matcher.end(matcher.groupCount() )      
        def prefix = fileName.substring(0,end)
        while(prefix.endsWith('-') || prefix.endsWith('_') || prefix.endsWith('.') ) 
          prefix=prefix[0..-2]
          
        return prefix
    }
    
    return fileName
}


