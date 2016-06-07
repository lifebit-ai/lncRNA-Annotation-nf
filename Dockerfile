############################################################
# Dockerfile for lncRNA-Annotation-nf
# Based on debian
############################################################

# Set the base image to Ubuntu
FROM debian:jessie

# File Author / Maintainer
MAINTAINER Evan Floden <evanfloden@gmail.com>

# Update the repository sources list
RUN apt-get update

# Install compiler and perl stuff
RUN apt-get install --yes --no-install-recommends \
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
 libexpat-dev 

# Install cpanminus for perl modules 
RUN apt-get install -y cpanminus

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
 Parallel::ForkManager

RUN apt-get install --yes \
 libarchive-zip-perl

# Install related DB modules
RUN apt-get install --yes \
 libdbd-mysql \
 libdbd-mysql-perl \
 libdbd-pgsql

# Install GD
RUN apt-get remove --yes libgd-gd2-perl

RUN apt-get install --yes \
 libgd2-noxpm-dev

RUN cpanm GD \
 GD::Graph \
 GD::Graph::smoothlines 


# Install BioPerl dependancies, mostly from cpan
RUN apt-get install --yes \
 libpixman-1-0 \
 libpixman-1-dev \
 graphviz \
 libxml-parser-perl \
 libsoap-lite-perl 

RUN cpanm Test::Most \
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
 XML::Twig XML::Writer

RUN apt-get install -y \
 libxml-libxml-perl \
 libxml-dom-xpath-perl \
 libxml-libxml-simple-perl \
 libxml-dom-perl

# Install BioPerl last built
RUN cpanm -v  \
 CJFIELDS/BioPerl-1.6.924.tar.gz 


# Install R
RUN echo "deb http://cran.rstudio.com/bin/linux/debian jessie-cran3/" >>  /etc/apt/sources.list &&\
 apt-key adv --keyserver keys.gnupg.net --recv-key 381BA480 &&\
 apt-get update --fix-missing && \
 apt-get -y install r-base

# Install R libraries
RUN R -e 'install.packages("ROCR", repos="http://cloud.r-project.org/"); install.packages("randomForest",repos="http://cloud.r-project.org/")'

# Install Star Mapper
RUN wget https://github.com/alexdobin/STAR/archive/2.5.2a.tar.gz &&\
 tar -xzf 2.5.2a.tar.gz && \
 cd STAR-2.5.2a &&\
 make STAR

# Install FEELnc
RUN wget https://github.com/tderrien/FEELnc/archive/a6146996e06f8a206a0ae6fd59f8ca635c7d9467.zip &&\
 unzip a6146996e06f8a206a0ae6fd59f8ca635c7d9467.zip &&\ 
 mv FEELnc-a6146996e06f8a206a0ae6fd59f8ca635c7d9467 /FEELnc

ENV FEELNCPATH /FEELnc
ENV PERL5LIB $PERL5LIB:${FEELNCPATH}/lib/

RUN rm -rf /root/.cpanm/work/
