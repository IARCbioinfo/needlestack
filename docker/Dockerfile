# Set the base image to nfcore base
FROM nfcore/base

# File Author / Maintainer
MAINTAINER Tiffany Delhomme <delhommet@students.iarc.fr>

RUN DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y r-base build-essential
RUN conda install -c bioconda bedtools=2.27.1
RUN conda install -c bioconda -c conda-forge samtools=1.9
RUN wget https://raw.githubusercontent.com/IARCbioinfo/mpileup2readcounts/master/mpileup2readcounts.cc && g++ -std=c++11 -O3 mpileup2readcounts.cc -o mpileup2readcounts && cp mpileup2readcounts /usr/local/bin && rm mpileup2readcounts.cc

