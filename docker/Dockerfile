# Set the base image to Ubuntu
FROM ubuntu:latest

# File Author / Maintainer
MAINTAINER Matthieu Foll <follm@iarc.fr>

# Update the repository sources list
RUN apt-get update -y

# Install dependences
RUN apt-get install -y g++ make wget git zlib1g-dev python libncurses5-dev

# Install bedtools 
RUN git clone https://github.com/arq5x/bedtools2.git && cd bedtools2 &&  make && make install && cd .. && rm -rf bedtools2

# Install samtools 
RUN wget https://github.com/samtools/samtools/releases/download/1.2/samtools-1.2.tar.bz2 && tar -jxvf samtools-1.2.tar.bz2 && cd samtools-1.2 && make && make install && cd .. && rm -rf samtools-1.2 samtools-1.2.tar.bz2

# Install pileup2base 
RUN wget -P /usr/local/bin/ https://github.com/mfoll/pileup2base/raw/patch-1/pileup2baseindel.pl && chmod +x /usr/local/bin/pileup2baseindel.pl

# Install R and package qvalue from bioconductor 
RUN apt-get install -y r-base && Rscript -e 'source("http://bioconductor.org/biocLite.R");biocLite("qvalue")'
