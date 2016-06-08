############################################################
# Dockerfile for lncRNA-Annotation-nf
# Based on debian
############################################################

# Set the base image to Ubuntu
FROM debian:jessie

# File Author / Maintainer
MAINTAINER Evan Floden <evanfloden@gmail.com>

# Install compiler and perl stuff
RUN apt-get update && apt-get install --yes --no-install-recommends \
 wget \
 ed \
 less \
 locales \
 vim-tiny \
 git \
 cmake \
 build-essential \
 gcc-multilib \
 apt-utils \
 perl \
 python \
 expat \
 libexpat-dev \
 libarchive-zip-perl \
 libdbd-mysql \
 libdbd-mysql-perl \
 libdbd-pgsql \
 libgd-gd2-perl \
 libgd2-noxpm-dev \
 libpixman-1-0 \
 libpixman-1-dev \
 graphviz \
 libxml-parser-perl \
 libsoap-lite-perl \
 libxml-libxml-perl \
 libxml-dom-xpath-perl \
 libxml-libxml-simple-perl \
 libxml-dom-perl \
 cpanminus \
 && rm -rf /var/lib/apt/lists/*


# Install perl modules
RUN cpanm --force CPAN::Meta \
 XML::Parser \
 readline \ 
 Term::ReadKey \
 YAML \
 Digest::SHA \
 Module::Build \
 ExtUtils::MakeMaker \
 Test::More \
 Data::Stag \
 Config::Simple \
 Statistics::Lite \
 Statistics::Descriptive \
 Parallel::ForkManager \
 GD \
 GD::Graph \
 GD::Graph::smoothlines \
 Test::Most \
 Algorithm::Munkres \
 Array::Compare Clone \
 PostScript::TextBlock \
 SVG \
 SVG::Graph \
 Set::Scalar \
 Sort::Naturally \
 Graph \
 GraphViz \
 HTML::TableExtract \
 Convert::Binary::C \
 Math::Random \
 Error \
 Spreadsheet::ParseExcel \
 XML::Parser::PerlSAX \
 XML::SAX::Writer \
 XML::Twig XML::Writer \
 && rm -rf /root/.cpanm/work


# Install BioPerl last built
RUN cpanm -v  \
 CJFIELDS/BioPerl-1.6.924.tar.gz \
 && rm -rf /root/.cpanm/work


# Install R
RUN echo "deb http://cran.rstudio.com/bin/linux/debian jessie-cran3/" >>  /etc/apt/sources.list \
 && apt-key adv --keyserver keys.gnupg.net --recv-key 381BA480 \
 && apt-get update --fix-missing \
 && apt-get -y install r-base

# Install R libraries
RUN R -e 'install.packages("ROCR", repos="http://cloud.r-project.org/"); install.packages("randomForest",repos="http://cloud.r-project.org/")'

# Install Star Mapper
RUN wget -qO- https://github.com/alexdobin/STAR/archive/2.5.2a.tar.gz | tar -xz \ 
 && cd STAR-2.5.2a \
 && make STAR

# Install FEELnc
RUN wget -q https://github.com/tderrien/FEELnc/archive/a6146996e06f8a206a0ae6fd59f8ca635c7d9467.zip \
 && unzip a6146996e06f8a206a0ae6fd59f8ca635c7d9467.zip \ 
 && mv FEELnc-a6146996e06f8a206a0ae6fd59f8ca635c7d9467 /FEELnc \
 && rm a6146996e06f8a206a0ae6fd59f8ca635c7d9467.zip


ENV FEELNCPATH /FEELnc
ENV PERL5LIB $PERL5LIB:${FEELNCPATH}/lib/


# Install Cufflinks/Cuffmerge 
RUN wget -qO- http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz | tar xz 

ENV PATH $PATH:/FEELnc/bin/LINUX:/FEELnc/utils:/FEELnc/scripts/:/cufflinks-2.2.1.Linux_x86_64:/STAR-2.5.2a/bin/Linux_x86_64/
