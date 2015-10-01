# Set the base image to Ubuntu
# Set the base image to Debian
FROM debian:latest

# File Author / Maintainer
MAINTAINER Matthieu Foll <follm@iarc.fr>

RUN apt-get clean && \
	apt-get update -y && \
	# Add R new version repos and update the repository sources list
	DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y software-properties-common && \
	DEBIAN_FRONTEND=noninteractive add-apt-repository "deb http://cran.rstudio.com//bin/linux/debian jessie-cran3/" && \
	apt-key adv --keyserver keys.gnupg.net --recv-key 381BA480 && \
	apt-get update -y && \

	# Install dependences
	DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y \
	g++ \
	make \
	git \
	zlib1g-dev \
	python \
	libncurses5-dev \
	ca-certificates \
	dialog \
	wget \
	bzip2 \
	apt-utils && \

	# Install bedtools specific version manually
	wget https://github.com/arq5x/bedtools2/releases/download/v2.24.0/bedtools-2.24.0.tar.gz && \
	tar -zxf bedtools-2.24.0.tar.gz && \
	cd bedtools2 &&  \
	make && \
	make install && \
	cd .. && \
	rm -rf bedtools2 bedtools-2.24.0.tar.gz && \
	
	# Install vcflib 
	git clone --recursive git://github.com/ekg/vcflib.git && \
	cd vcflib &&  \
	make && \
	cp bin/vcfoverlay /usr/local/bin/ && \
	cd .. && \
	rm -rf vcflib && \

	# Install samtools specific version manually
	wget https://github.com/samtools/samtools/releases/download/1.2/samtools-1.2.tar.bz2 && \
	tar -jxf samtools-1.2.tar.bz2 && \
	cd samtools-1.2 && \
	make && \
	make install && \
	cd .. && \
	rm -rf samtools-1.2 samtools-1.2.tar.bz2 && \

	# Install R and package qvalue from bioconductor 
	DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y r-base && \

	# Remove unnecessary dependences
	DEBIAN_FRONTEND=noninteractive apt-get remove -y \
	g++ \
	make \
	git \
	zlib1g-dev \
	python \
	libncurses5-dev \
	dialog \
	apt-utils \
	wget \
	bzip2 \
	software-properties-common && \

	# Clean
	DEBIAN_FRONTEND=noninteractive apt-get autoremove -y && \ 
	apt-get clean
