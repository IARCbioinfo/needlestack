# Robust regression multi-sample somatic variant caller

Warning: development in progress, unreliable results warrantied. 

You will need a set of [BAM files](https://samtools.github.io/hts-specs/) (called `*.bam`) grouped in a single folder along with their [index files](http://www.htslib.org/doc/samtools.html) (called `*.bam.bai`), a [bed file](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) and a reference [fasta file](https://en.wikipedia.org/wiki/FASTA_format) (eventually compressed with [bgzip](http://www.htslib.org/doc/tabix.html)) along with its [faidx index](http://www.htslib.org/doc/faidx.html) (and `*.gzi` faidx index if compressed).

## Quick start

Works under most Linux distributions and Apple OS X.

1. Install [java](https://java.com/download/) JRE if you don't already have it.

2. Install [nextflow](http://www.nextflow.io/).

	```bash
	curl -fsSL get.nextflow.io | bash
	```
	And move it to a location in your `$PATH` (`/usr/local/bin` for example here):
	```bash
	sudo mv nextflow /usr/local/bin
	```
3. Install [docker](https://www.docker.com).
	
	This is very system specific (but quite easy in most cases), follow  [docker documentation](https://docs.docker.com/installation/). Also follow the optional configuration step called `Create a Docker group` in their documentation.

4. Optionally download a sample dataset.

	```bash
	git clone --depth=1 https://github.com/mfoll/NGS_data_test.git
	```
5. Run the pipeline.
	
	Here on the example dataset downloaded above:
	```sh
	cd NGS_data_test/1000G_CEU_TP53/
	nextflow run mfoll/robust-regression-caller -with-docker mfoll/robust-regression-caller \
			--bed TP53_all.bed --nsplit 10 --bam_folder BAM/ --fasta_ref 17.fasta.gz
	```
	
	You will find a [VCF file](https://samtools.github.io/hts-specs/) called `all_variants.vcf` in the `BAM/` folder once done.
	
	The first time it will take more time as the pipeline will be downloaded from github and the docker container from [dockerhub](https://hub.docker.com/r/mfoll/robust-regression-caller/).


## Detailed instructions

If you can't install [docker](https://www.docker.com) or don't want to use it, the pipeline will also work if you install [perl](https://www.perl.org),  [bedtools](http://bedtools.readthedocs.org/en/latest/), [samtools](http://www.htslib.org), vcfoverlay from [vcflib](https://github.com/ekg/vcflib) and Rscript from [R](https://www.r-project.org) and put them in your path (executables are assumed to be respectively called `perl`, `bedtools`, `samtools`, `vcflib` and `Rscript`). In this case, remove the `-with-docker` option from step 5 above.

The exact same pipeline can be run on your computer or on a HPC cluster, simply by adding a [nextflow configuration file](http://www.nextflow.io/docs/latest/config.html) to choose your scheduler. See the nextflow documentation [here](http://www.nextflow.io/docs/latest/config.html#scope-executor) how to do that.

`--bed`, `--bam_folder` and `--fasta_ref` are compulsary. The optional parameters with default values are:

| Parameter | Default value | Description |
|-----------|--------------:|-------------| 
| min_dp    |            50 | Minimum coverage in at least one sample to consider a site |
| min_ao | 5 | Minimum number of non-ref reads in at least one sample to consider a site|
| nsplit | 1 | Split the bed file in nsplit pieces and run in parallel |
| min_qval | 50 | qvalue in Phred scale to consider a variant |
| sor_snv | 4 | Strand bias SOR threshold for snv |
| sor_indel | 10 | Strand bias SOR threshold for indels |
| map_qual | 20 | Min mapping quality (passed to samtools) |
| base_qual | 20 | Min base quality (passed to samtools) |
| max_DP | 30000 | Downsample coverage per sample (passed to samtools) |
| sample_names | "BAM" | Put "FILE" to use the bam file names as sample names and "BAM" to use the sample name filed from the bam files |
| all_sites | "FALSE" | Output all sites, even when no variant is detected (but still affected by min_dp and min_ao) |
| do_plots | "TRUE" | Produce pdf plots of regressions |

Simply add the parameters you want in the command line like `--min_dp 1000` for exmaple to change the min coverage.

[![Join the chat at https://gitter.im/mfoll/robust-regression-caller](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/mfoll/robust-regression-caller?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

[![Circle CI](https://circleci.com/gh/mfoll/robust-regression-caller/tree/master.svg?style=shield)](https://circleci.com/gh/mfoll/robust-regression-caller/tree/master) 

[![Docker Hub](https://img.shields.io/badge/docker-ready-blue.svg)](https://hub.docker.com/r/mfoll/robust-regression-caller/)
