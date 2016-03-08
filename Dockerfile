FROM ubuntu:12.04

FROM ubuntu:14.04
RUN  apt-get update \
  && apt-get install -y wget \
  && rm -rf /var/lib/apt/lists/*


MAINTAINER Evan Floden <evanfloden@gmail.com>
		
#
# Install STAR
# https://github.com/alexdobin/STAR
#
RUN wget https://github.com/alexdobin/STAR/archive/STAR_2.4.2a.tar.gz && \
    tar -zxf STAR_2.4.2a.tar.gz && \
    cd STAR-STAR_2.4.2a && \
    export PATH=$PATH:~/STAR-STAR_2.4.2a/bin/Linux_x86_64
    
	
ENTRYPOINT ["/bin/sh", "-c"]	
    
