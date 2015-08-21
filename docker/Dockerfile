# Set the base image to Ubuntu
FROM debian:latest

# File Author / Maintainer
MAINTAINER Matthieu Foll <follm@iarc.fr>

# Install pileup2base 
#RUN wget -P /usr/local/bin/ https://github.com/mfoll/pileup2base/raw/follm/master/pileup2baseindel.pl && \
ADD https://github.com/mfoll/robust-regression-caller/raw/master/bin/pileup2baseindel.pl /usr/local/bin/
ADD https://github.com/mfoll/robust-regression-caller/raw/master/bin/bed_large_cut.r /usr/local/bin/
ADD https://github.com/mfoll/robust-regression-caller/raw/master/bin/pileup_nbrr_caller_vcf.r /usr/local/bin/

RUN apt-get update -y && \
	# Add R new version repos and update the repository sources list
	DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y \
	software-properties-common && \
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

	# Make pileup2base executable
	chmod +x /usr/local/bin/pileup2baseindel.pl && \
	# Make R scripts executable
	chmod +x /usr/local/bin/bed_large_cut.r && \
	chmod +x /usr/local/bin/pileup_nbrr_caller_vcf.r && \

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

# Install samtools specific version manually
#RUN wget https://github.com/samtools/samtools/releases/download/1.2/samtools-1.2.tar.bz2 && \
#	tar -jxvf samtools-1.2.tar.bz2 && \
#	cd samtools-1.2 && \
#	make && \
#	make install && \
#	cd .. && \
#	rm -rf samtools-1.2 samtools-1.2.tar.bz2

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
#	Rscript -e 'source("http://bioconductor.org/biocLite.R");biocLite("qvalue")' && \

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
