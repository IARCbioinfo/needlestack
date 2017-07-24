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

	# Install R 
	DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y r-base && \

	# Install Bioconductor packages
	DEBIAN_FRONTEND=noninteractive apt-get -y install gfortran && \
	DEBIAN_FRONTEND=noninteractive apt-get -y install libxml2-dev && \
	DEBIAN_FRONTEND=noninteractive apt-get -y install libssl-dev && \
	DEBIAN_FRONTEND=noninteractive apt-get -y install libcurl4-openssl-dev && \	
	
	Rscript -e "source('http://bioconductor.org/biocLite.R'); biocLite('Gviz'); biocLite('TxDb.Hsapiens.UCSC.hg19.knownGene'); biocLite('TxDb.Hsapiens.UCSC.hg18.knownGene'); biocLite('TxDb.Hsapiens.UCSC.hg38.knownGene'); biocLite('org.Hs.eg.db')" && \
	
	# mpileup2readcounts compilation
	wget https://raw.githubusercontent.com/IARCbioinfo/mpileup2readcounts/master/mpileup2readcounts.cc && \
	g++ -std=c++11 -O3 mpileup2readcounts.cc -o mpileup2readcounts && \
	cp mpileup2readcounts /usr/local/bin && \
	rm mpileup2readcounts.cc && \

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
	gfortran \
	libssl-dev \
	libcurl4-openssl-dev \
	software-properties-common && \

	# Clean
	DEBIAN_FRONTEND=noninteractive apt-get autoremove -y && \ 
	apt-get clean
