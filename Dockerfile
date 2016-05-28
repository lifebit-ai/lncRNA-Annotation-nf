############################################################
# Dockerfile for lncRNA-Annotation-nf
# Based on debian
############################################################

# Set the base image to Ubuntu
FROM debian

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
 expat \
 libexpat-dev 

# Install cpanminus for perl modules 
RUN apt-get install -y cpanminus

# Install perl modules
RUN cpanm --force CPAN::Meta \
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
 cd FEELnc-a6146996e06f8a206a0ae6fd59f8ca635c7d9467 &&\
 export FEELNCPATH=${PWD} &&\
 export PERL5LIB=$PERL5LIB:${FEELNCPATH}/lib/
