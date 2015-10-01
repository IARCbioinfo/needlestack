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
	apt-utils && \

	# Install bedtools 
	git clone https://github.com/arq5x/bedtools2.git && \
	cd bedtools2 &&  \
	make && \
	make install && \
	cd .. && \
	rm -rf bedtools2 && \
	
	# Install vcflib 
	git clone --recursive git://github.com/ekg/vcflib.git && \
	cd vcflib &&  \
	make && \
	cp bin/vcfoverlay /usr/local/bin/ && \
	cd .. && \
	rm -rf vcflib && \

	# Install samtools from github repos (htslib needed first)
	git clone git://github.com/samtools/htslib.git && \
	cd htslib && \
	make && \
	make install && \
	cd .. && \
	git clone git://github.com/samtools/samtools.git && \
	cd samtools && \
	make && \
	make install && \
	cd .. && \
	rm -rf htslib samtools && \

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
	software-properties-common && \

	# Clean
	DEBIAN_FRONTEND=noninteractive apt-get autoremove -y && \ 
	apt-get clean
