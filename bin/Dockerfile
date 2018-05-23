# Set the base image to Debian
FROM debian:jessie

# File Author / Maintainer
MAINTAINER Matthieu Foll <follm@iarc.fr>

RUN mkdir -p /var/cache/apt/archives/partial && \
	touch /var/cache/apt/archives/lock && \
	chmod 640 /var/cache/apt/archives/lock && \
	apt-key adv --keyserver keyserver.ubuntu.com --recv-keys F76221572C52609D && \
	apt-get clean && \
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
	
	# Install R 
	DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y r-base && \

	# Install Bioconductor packages
	DEBIAN_FRONTEND=noninteractive apt-get -y install gfortran && \
	DEBIAN_FRONTEND=noninteractive apt-get -y install libxml2-dev && \
	DEBIAN_FRONTEND=noninteractive apt-get -y install libssl-dev && \
	DEBIAN_FRONTEND=noninteractive apt-get -y install libcurl4-openssl-dev && \	
	
	Rscript -e "source('http://bioconductor.org/biocLite.R'); biocLite('Gviz'); biocLite('TxDb.Hsapiens.UCSC.hg19.knownGene'); biocLite('TxDb.Hsapiens.UCSC.hg18.knownGene'); biocLite('TxDb.Hsapiens.UCSC.hg38.knownGene'); biocLite('org.Hs.eg.db')" && \

	# Install bedtools specific version manually
	wget https://github.com/arq5x/bedtools2/releases/download/v2.25.0/bedtools-2.25.0.tar.gz && \
	tar -zxf bedtools-2.25.0.tar.gz && \
	cd bedtools2 &&  \
	make && \
	make install && \
	cd .. && \
	rm -rf bedtools2 bedtools-2.25.0.tar.gz && \
	
	# Install samtools specific version manually
	wget https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2 && \
	tar -jxf samtools-1.3.1.tar.bz2 && \
	cd samtools-1.3.1 && \
	make && \
	make install && \
	cd .. && \
	rm -rf samtools-1.3.1 samtools-1.3.1.tar.bz2 && \
	
	# mpileup2readcounts compilation
	wget https://raw.githubusercontent.com/IARCbioinfo/mpileup2readcounts/master/mpileup2readcounts.cc && \
	g++ -std=c++11 -O3 mpileup2readcounts.cc -o mpileup2readcounts && \
	cp mpileup2readcounts /usr/local/bin && \
	rm mpileup2readcounts.cc && \
	
  # retrieve files dependencies
  wget https://github.com/IARCbioinfo/needlestack/tree/v1.1b/bin/glm_rob_nb.r && \
	wget https://github.com/IARCbioinfo/needlestack/tree/v1.1b/bin/needlestack.r && \
  wget https://github.com/IARCbioinfo/needlestack/tree/v1.1b/bin/plot_alignments.r && \
  wget https://github.com/IARCbioinfo/needlestack/tree/v1.1b/bin/bed_cut.r && \
  wget https://github.com/IARCbioinfo/needlestack/tree/v1.1b/bin/plot_rob_nb.r && \
  wget https://github.com/IARCbioinfo/needlestack/tree/v1.1b/bin/hg19_chromosomeNames2UCSC.txt && \
  wget https://github.com/IARCbioinfo/needlestack/tree/v1.1b/bin/hg38_chromosomeNames2UCSC.txt && \
  
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
	gfortran \
	libssl-dev \
	libcurl4-openssl-dev \
	software-properties-common && \

	# Clean
	DEBIAN_FRONTEND=noninteractive apt-get autoremove -y && \ 
	apt-get clean

