FROM debian:latest

MAINTAINER Evan Floden <evanfloden@gmail.com>

RUN  apt-get update \
  && apt-get install -y \
  curl \
  wget \
  libcurl4-gnutls-dev \
  build-essential \
  vim \
  libz-dev

#
# Install STAR
# https://github.com/alexdobin/STAR
#
RUN wget https://github.com/alexdobin/STAR/archive/STAR_2.4.2a.tar.gz && \
    tar -zxf STAR_2.4.2a.tar.gz && \
    cd STAR-STAR_2.4.2a/source && \
    make STAR && \
    export PATH=$PATH:/STAR-STAR_2.4.2a/source

ENTRYPOINT ["/bin/sh", "-c"]

