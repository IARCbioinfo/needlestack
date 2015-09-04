# Robust regression multi-sample variant caller

## Quick start

Works under most Linux distributions and MacOS.

Install [nextflow](http://www.nextflow.io) (see [here](http://www.nextflow.io/docs/latest/getstarted.html#installation) how to do this).

Install [docker](https://www.docker.com) (see [here](https://docs.docker.com/installation/) how to do this). 

If you can't install it or don't want to use docker, it will also work if you install [perl](https://www.perl.org),  [bedtools](http://bedtools.readthedocs.org/en/latest/), [samtools](http://www.htslib.org), vcfoverlay from [vcflib](https://github.com/ekg/vcflib) and Rscript from [R](https://www.r-project.org) and put them in your path (executables are assumed to be respectively called bedtools, samtools, vcflib and Rscript). 

The exact same pipeline can be run on your computer or on a HPC cluster, simply by adding a [nextflow configuration file](http://www.nextflow.io/docs/latest/config.html) to choose your scheduler. See the nextflow documentation [here](http://www.nextflow.io/docs/latest/config.html#scope-executor) how to do that.

You will need a set of BAM files (called `*.bam`) in a folder  (called `BAM` here) along with their index (called `*.bam.bai`), a `bed` file and a reference `fasta` file (eventually gz compressed) along with its index (`*.fai`, and `*.gzi` if compressed)  
Run it the variant calling simply using:
```sh
nextflow run mfoll/robust-regression-caller --bed your_bed_file.bed --nsplit 10 --bam_folder BAM/ --fasta_ref my_ref.fasta
```
You will find both a VCF file in the `BAM` folder once done.

`--bed`, `--bam_folder` and `--fasta_ref` are compulsary. The optional parameters with default values are:
```sh
min_dp = 50  # minimum coverage in at least one sample to consider a site
min_ao = 5 # minimum number of non-ref reads in at least one sample to consider a site
nsplit = 1 # split the bed file in nsplit pieces and run in parallel 
min_qval = 50 # qvalue in Phred scale to consider a variant 
sor_snv = 4 # strand bias SOR threshold for snv
sor_indel = 10 # strand bias SOR threshold for indels
map_qual = 20 # min mapping quality (passed to samtools)
base_qual = 20 # min base quality (passed to samtools)
max_DP = 30000 # downsample coverage per sample (passed to samtools)
sample_names = "BAM" # put FILE to use the bam file names as sample names and BAM to use the sample name filed from the bam files
all_sites = "FALSE" # output all sites, even when no variant is detected (but still affected by min_dp and min_ao)
do_plots = "TRUE" # produce pdf plots of regressions 
```
Simply add the parameters you want in the command line like `--min_dp 1000` for exmaple to change the min coverage.

[![Join the chat at https://gitter.im/mfoll/robust-regression-caller](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/mfoll/robust-regression-caller?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

[![Circle CI](https://circleci.com/gh/mfoll/robust-regression-caller/tree/master.svg?style=shield)](https://circleci.com/gh/mfoll/robust-regression-caller/tree/master) 

[![Docker Hub](https://img.shields.io/badge/docker-ready-blue.svg)](https://hub.docker.com/r/mfoll/robust-regression-caller/)
