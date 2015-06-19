# Set the base image to Ubuntu
FROM debian:latest

# File Author / Maintainer
MAINTAINER Matthieu Foll <follm@iarc.fr>

# Update the repository sources list
RUN apt-get update -y

# Install dependences
RUN DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y \
	g++ \
	make \
	wget \
	git \
	zlib1g-dev \
	python \
	libncurses5-dev \
	ca-certificates \
	dialog \
	bzip2

# Install bedtools 
RUN git clone https://github.com/arq5x/bedtools2.git && \
	cd bedtools2 &&  \
	make && \
	make install && \
	cd .. && \
	rm -rf bedtools2

# Install samtools specific version manually
#RUN wget https://github.com/samtools/samtools/releases/download/1.2/samtools-1.2.tar.bz2 && \
#	tar -jxvf samtools-1.2.tar.bz2 && \
#	cd samtools-1.2 && \
#	make && \
#	make install && \
#	cd .. && \
#	rm -rf samtools-1.2 samtools-1.2.tar.bz2

# Install samtools from github repos (htslib needed first)
RUN git clone git://github.com/samtools/htslib.git && \
	cd htslib && \
	make && \
	make install && \
	cd .. && \
	git clone git://github.com/samtools/samtools.git && \
	cd samtools && \
	make && \
	make install && \
	cd .. && \
	rm -rf htslib samtools

# Install pileup2base 
RUN wget -P /usr/local/bin/ https://github.com/mfoll/pileup2base/raw/follm/master/pileup2baseindel.pl && \
	chmod +x /usr/local/bin/pileup2baseindel.pl

# Install R and package qvalue from bioconductor 
RUN DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y r-base && \
	Rscript -e 'source("http://bioconductor.org/biocLite.R");biocLite("qvalue")'
