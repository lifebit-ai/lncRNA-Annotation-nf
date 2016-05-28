############################################################
# Dockerfile for lncRNA-Annotation-nf
# Based on debian
############################################################

# Set the base image to Ubuntu
FROM debian::jessie

# File Author / Maintainer
MAINTAINER Evan Floden <evanfloden@gmail.com>

# Update the repository sources list
RUN apt-get update

# Install compiler and perl stuff
RUN apt-get install --yes --no-install-recommends \
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
RUN cpanm CPAN::Meta \
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
RUN R -e 'install.packages('ROCR'); install.packages('randomForest')'

# Install Star Mapper
RUN wget https://github.com/alexdobin/STAR/archive/STAR_2.5.2a.tar.gz &&\
 tar -xzf STAR_2.5.2a.tar.gz &&\
 cd STAR_2.5.2a &&\
 make STAR

# Install FEELnc
RUN git clone git@github.com:tderrien/FEELnc.git &&\
 cd FEELnc &&\
 export FEELNCPATH=${PWD} &&\
 export PERL5LIB=$PERL5LIB:${FEELNCPATH}/lib/ &&\
