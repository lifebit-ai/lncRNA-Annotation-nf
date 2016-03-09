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


params.name          = "lncRNA Annotation from Pig RNA-Seq"
params.genome        = "$baseDir/tutorial/genome/Scrofa10.2.chr1.fa"
params.annotation    = "$baseDir/tutorial/annotation/Sus_scrofa.Sscrofa10.2.62.gtf"
params.reads         = "$baseDir/tutorial/reads/*_{1,2}.fastq.gz"
params.overhang      = '99'
params.threads       = '1'
params.output        = "results/"


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

genomeFile              = file(params.genome)
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
      
    script:
    //
    // STAR Generate Index
    //
    """
        mkdir STARgenome
        STAR --runThreadN ${params.threads} --runMode genomeGenerate --genomeDir STARgenome \
             --genomeFastaFiles ${genomeFile} --sjdbGTFfile ${annotationFile} \
             --sjdbOverhang ${params.overhang} --outFileNamePrefix STARgenome \
             --genomeSAindexNbases 3
    """
}


process mapping {
    tag "reads: $name"

    input:
    file STARgenome from STARgenomeIndex.first()
    set val(name), file(reads:'*') from read_files

    output:
    file "${name}" into STARmappedReads 

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
        mkdir ${name}
        mv ${name}Aligned* ${name}/.
        mv ${name}Signal* ${name}/.
        mv ${name}SJ* ${name}/.
        mv ${name}Log* ${name}/.
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


