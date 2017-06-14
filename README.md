<img align="center" src="https://cloud.githubusercontent.com/assets/3366818/12507528/99afd6fa-c0f6-11e5-993b-5ac925d178fa.png" width="700">

## A multi-sample somatic variant caller

<img align="center" src="https://cloud.githubusercontent.com/assets/3366818/10489240/7da79144-729c-11e5-8cb3-0225106d9b06.jpg" width="600">

[![Join the chat at https://gitter.im/iarcbioinfo/needlestack](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/iarcbioinfo/needlestack?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge) [![Circle CI](https://circleci.com/gh/IARCbioinfo/needlestack/tree/master.svg?style=shield&circle-token=402d456a7c50af352bb4e1a52425ce0fe645f78f)](https://circleci.com/gh/IARCbioinfo/needlestack/tree/master) [![Docker Hub](https://img.shields.io/badge/docker-ready-blue.svg)](https://hub.docker.com/r/iarcbioinfo/needlestack/)

Warning: development in progress, unreliable results warranted.

Please wait upcoming publication before using it in production.

Contact: follm@iarc.fr

Workflow representation

## Description

Needlestack is an ultra-sensitive multi-sample variant caller for Next Generation Sequencing (NGS) data. It is based on the idea that analysing several samples together can help estimate the distribution of sequencing errors to accurately identify variants. It has been initially developed for somatic variant calling using very deep NGS data from circulating free DNA, but is also applicable to lower coverage data like Whole Exome Sequencing (WES) or even Whole Genome Sequencing (WGS). It is a highly scalable and reproducible pipeline thanks to the use of [nextflow](http://www.nextflow.io/) and [docker](https://www.docker.com) technologies.

Here is a summary of the method:

- At each position and for each candidate variant, we model sequencing errors using a negative binomial regression with a linear link and a zero intercept. The data is extracted from the BAM files using [samtools](http://www.htslib.org).
- Genetic variants are detected as being outliers from the error model. To avoid these outliers biasing the regression we use a robust estimator for the negative binomial regression (published [here](http://www.ncbi.nlm.nih.gov/pubmed/25156188) with code available [here](https://github.com/williamaeberhard/glmrob.nb)).
- We calculate for each sample a p-value for being a variant (outlier from the regression) that we further transform into q-values to account for multiple testing.


## Dependencies

Needlestack works under most Linux distributions and Apple OS X.

1. This pipeline is based on [nextflow](https://www.nextflow.io). As we have several nextflow pipelines, we have centralized the common information in the [IARC-nf](https://github.com/IARCbioinfo/IARC-nf) repository. Please read it carefully as it contains essential information for the installation, basic usage and configuration of nextflow and our pipelines.

2. External software :
- Install [bedtools](http://bedtools.readthedocs.org/en/latest/)
- Install [samtools](http://www.htslib.org)
- Install Rscript from [R](https://www.r-project.org)
- Install g++ compiler
- Compile the file *mpileup2readcounts.cc* located [here](https://github.com/IARCbioinfo/mpileup2readcounts)
- Add the previous tools to your path (executables are assumed to be respectively called `samtools`, `Rscript` and `mpileup2readcounts`)

You can avoid installing all the external software by only installing Docker. See the [IARC-nf]((https://github.com/IARCbioinfo/IARC-nf) repository for more information.

To use needlestack without nextflow, in addition to the previous tools, download the files in this [bin](https://github.com/IARCbioinfo/needlestack/tree/gabriela_cpp/bin) directory and add them to your path.


## Input
| Type      | Description     |
|-----------|---------------|
| [BAM files](https://samtools.github.io/hts-specs/)  | BAM files (called `*.bam`) grouped in a single folder along with their [index files](http://www.htslib.org/doc/samtools.html) (called `*.bam.bai`). A minimum of 20 BAM files is recommended. |
| [fasta file](https://en.wikipedia.org/wiki/FASTA_format)  |A reference fasta file (eventually compressed with [bgzip](http://www.htslib.org/doc/tabix.html)) along with its [faidx index](http://www.htslib.org/doc/faidx.html) (and `*.gzi` faidx index if compressed). |
| [bed file](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) | Optional input file. Otherwise the variant calling is performed on the whole reference provided.
|

A sample dataset is available at this [address](https://github.com/mfoll/NGS_data_test.git) for testing.

## Parameters

Type `--help` to get the full list of options.

  * #### Mandatory
| Name      | Example value | Description     |
|-----------|---------------|-----------------|
| --bam_folder    |     BAM/ | Folder containing all the BAM files |
| --fasta_ref   |        fasta_name.fa   | A reference fasta file (eventually compressed with [bgzip](http://www.htslib.org/doc/tabix.html)) along with its [faidx index](http://www.htslib.org/doc/faidx.html) (and `*.gzi` faidx index if compressed) |
| --output_vcf   |         vcf_name.vcf | The name of the VCF file generated by needlestack |

  * #### Optional
| Name      | Default value | Description     |
|-----------|---------------|-----------------|
| min_dp    |            30 | Minimum median coverage to consider a site. In addition, at least 10 samples have to be covered by min_dp. |
| min_ao | 3 | Minimum number of non-ref reads in at least one sample to consider a site |
| nsplit | 1 | Split the bed file in nsplit pieces and run in parallel |
| min_qval | 50 | qvalue threshold in [Phred scale](https://en.wikipedia.org/wiki/Phred_quality_score) to consider a variant |
| sb_type | SOR | Strand bias measure, either SOR, RVSB or FS |
| sb_snv | 100 or 1000 | Strand bias threshold for SNVs (100 (1000 if FS) = no filter) |
| sb_indel | 100 or 1000 | Strand bias threshold for indels (100 (1000 if FS) = no filter)|
| map_qual | 20 | Min mapping quality (passed to samtools) |
| base_qual | 20 | Min base quality (passed to samtools) |
| max_dp | 30000 | Downsample coverage per sample (passed to samtools) |
| plots | SOMATIC if tn_pairs provided, ALL if not | To create pdf plots of regressions in the output. To remove pdf plots set --do_plot to NONE (See *Plot options* paragraph).|
| ref_genome |  - | Reference genome for the annotions on the alignments plots. Examples : *Hsapiens.UCSC.hg19*, *Hsapiens.UCSC.hg19*, *Hsapiens.UCSC.hg38*, *Mmusculus.UCSC.mm10*. The right terminology is described below (*Plot options* paragraph). WARNING : this option is mandatory if the --do_alignments flag is used|
| output_folder | --bam_folder | Output folder, by default equals to the input bam folder |
| bed |  - | BED file containing a list of regions (or positions) where needlestack should be run |
| region |  - | A region in format CHR:START-END where calling should be done |
| tn_pairs | - | A tab-delimited file containing two columns (normal and tumor sample names) for each sample in line. This enables matched tumor/normal pair calling features (see below) |
| sigma_normal | 0.1 | Sigma parameter for negative binomial modeling of expected germline allelic fraction. We strongly recommend not to change this parameter unless you really know what it means |

By default, if neither `--bed` nor `--region` are provided, needlestack would run on whole reference, building a bed file from fasta index inputed.
If `--bed` and `--region` are both provided, it should run on the region only.

  * #### Flags

Flags are special parameters without value.

| Name      | Description     |
|-----------|-----------------|
| --help    | Display help |
| all_SNVs | Put this argument to output all SNVs, even when no variant is detected. Note that positions with zero coverage for all samples might still be missing depending on how the region split is performed |
| input_vcf | A VCF file (from GATK) where calling should be done. Needlestack will extract DP and AO from this VCF (DP and AD fields) and annotate it with phred q-value score (`FORMAT/QVAL` field), error rate (`INFO/ERR`) and over-dispersion sigma parameter (`INFO/SIG`). WARNING: by default, only work with split (coreutils) version > 8.13 |
| power_min_af | Allelic fraction used to classify genotypes to 0/0 or ./. depending of the power to detect a variant at this fraction (see below) |
| extra_robust_gl | Add this argument to perform extra-robust regression (useful for common germline SNPs, see below) |
| no_labels | Put this argument for not labeling the outliers on regression plots |
| no_indels | Put this argument to do not perform the variant calling on insertions and deletions |
| no_contours | Put this argument to do not plot qvalues contours (for qvalue threshold={10,30,50,70,100} by default) and do not plot minimum detectable allelic fraction in function of coverage |
| use_file_name | Put this argument to use the bam file names as sample names. By default the sample name is extracted from the bam file SM tag |
| do_alignments | Put this argument to add alignments plots to the pdf plots of regressions |

Simply add the parameters you want in the command line like `--min_dp 1000` for example to change the min coverage or `--all_SNVs` to output all sites.

## Usage

1. Optionally download a sample dataset.

	```bash
	git clone --depth=1 https://github.com/IARCbioinfo/data_test
	```
	
2. Run the pipeline using nextflow.

	Here on the example dataset downloaded above:
	```bash
	cd data_test
	nextflow run iarcbioinfo/needlestack -with-docker  \
	         --bed BED/TP53_all.bed --bam_folder BAM/BAM_multiple/ --ref REF/17.fasta.gz --output_vcf all_variants.vcf
	```

	You will find a [VCF file](https://samtools.github.io/hts-specs/) called `all_variants.vcf` in the `BAM_multiple/` folder once done.

	The first time it will take more time as the pipeline will be downloaded from github and the docker container from [dockerhub](https://hub.docker.com/r/iarcbioinfo/needlestack/).

	Creating an alias for the long command above can be useful. For example:
	```bash
	alias needlestack='nextflow run iarcbioinfo/needlestack -with-docker'
	```

	If you want to permanantly add this alias (and not just for your current session), add the above  line to your `~/.bashrc` file (assuming you are using bash).

	It will allow you to do this:
	```bash
	needlestack --bed BED/TP53_all.bed --bam_folder BAM/BAM_multiple/ --ref REF/17.fasta.gz --output_vcf all_variants.vcf
	```

	Official releases can be found [here](https://github.com/iarcbioinfo/needlestack/releases/). There is a corresponding official [docker container](https://hub.docker.com/r/iarcbioinfo/needlestack/) for each release and one can run a particular version using (for example for v0.3):
		```bash
		nextflow run iarcbioinfo/needlestack -r v0.3 -with-docker \
		        --bed BED/TP53_all.bed --bam_folder BAM/BAM_multiple/ --ref REF/17.fasta.gz --output_vcf all_variants.vcf
		```
2.  Run the pipeline without nextflow.

	Here on the example dataset downloaded above:
	```bash
	cd data_test/
	needlestack.sh --region=17:7572814-7573814 --bam_folder=BAM/BAM_multiple --ref=REF/17.fasta.gz --output_vcf=all_variants.vcf
	```

	You will find a [VCF file](https://samtools.github.io/hts-specs/) called `all_variants.vcf` in the `BAM/BAM_multiple/` folder once done.



## Output
  | Type      | Description     |
  |-----------|---------------|
  | VCF file    | ...... |
  | PDF files    | ...... |


## Detailed description 

### Germline, somatic, matched Tumor-Normal pairs calling and contamination

When using matched tumor/normal, Needlestack can classify variants (VCF `FORMAT/STATUS` field) according to the following table:
![Status Table](STATUS_TABLE.png)

For this one need to provide a tab-delimited file containing two columns with normal and tumor sample names using the `--tn_pairs` option. The first line of this file is a header with `TUMOR` and `NORMAL` keywords. When one normal or one tumor is missing, one can write `NA`. In this mode, the parameter `power_min_af` defines the allelic fraction in the tumor one is trying to detect to classify genotypes as `./.` or `0/0` depending on the power to detect this allelic fraction. Variants found as somatic in a tumor, but germline in another sample of the series will be flagged as `POSSIBLE_CONTAMINATION`. We found this particularly important, as needlestack is very sensitive to low allelic fractions, to filter out contamination among samples for pooled exome capture.

In other cases (when there is no `--tn_pairs` parameter defined), genotypes are defined as `./.` or `0/0` assuming one is looking for allelic fractions expected for germline variants (negative binomial distribution centered at 0.5 with over-dispersion parameter sigma=`sigma_normal`, with `sigma_normal=0.1` by default). If you are looking for somatic variants without matched-normal and assuming you are interesting to correctly distinguish `./.` and `0/0`genotypes, you can set the `power_min_af` parameter to the lowest allelic fraction of somatic variants you are interested with (and your coverage allows you to find).  Note that this is by far not the most common situation, and that in most cases you don't have to worry about the `power_min_af` parameter.

### Plot options

#### Two options to consider :

1. --plots :

	* SOMATIC : To produce pdf regression plots only for somatic variants. Default value when using matched tumor/normal.
	* ALL : To produce pdf regression plots for all variants. Default value when not using matched tumor/normal.
	* NONE : To remove pdf regression plots from the output

2. --do_alignments : To add the alignments plots to the regression plots. If this option is set to "true", the name of the reference genome (--ref_genome option) needs to be provided to choose the correct annotation (See *--ref_genome option* below). False is the default value.

#### Bioconductor packages to install for plotting the alignments :

- The [Gviz](https://bioconductor.org/packages/release/bioc/html/Gviz.html) package
- An [Annotation package for TxDb objects](http://bioconductor.org/packages/release/BiocViews.html#___TxDb).
- A [Genome wide annotation](https://bioconductor.org/packages/release/BiocViews.html#___OrgDb), it contains mappings between Entrez Gene identifiers and GenBank accession numbers. Examples : the Genome wide annotation package for Human : *org.Hs.eg.db* and for the mouse : *org.Mm.eg.db*.

Examples :

If one used the UCSC version of the reference Human genome hg19, the fowolling packages should be installed :
- *Gviz*
- *TxDb.Hsapiens.UCSC.hg19.knownGene* (hg18 and hg38 UCSC version can also be used)
- *org.Hs.eg.db*

These packages exist for other organisms than Human but have not been tested.

One can for example generate the alignments plot for data issued from the mouse by installing :
- *TxDb.Mmusculus.UCSC.mm10.knownGene*(mm9 UCSC version can also be used)
- *org.Mm.eg.db*

For the other organisms the packages need to have the same nomenclature as the ones listed above.

#### --ref_genome option :

The argument corresponds to the TxDb annotation package name without its extrimities. For the UCSC version of the reference Human genome hg19, one needs to set --ref_genome to "*Hsapiens.UCSC.hg19*".
Note that the packages chosen for the annotations are compatible with the UCSC notations since most of the Gviz fonctionalities can handle these notations.

By default, when using matched tumor/normal (tn_pairs option), needlestack will produce pdf plots of regressions only for somatic variants and without the alignment plots; when not, needlestack will produce them for all variants and without the alignment plots.

![Example of an alignment plot](alignments.png "Example of an alignment plot")

#### The alignment plot from top to bottom :

- Chromosome representation : a red vertical line shows the position of the variant
- Genomic axis associated with the alignment.
- BAM alignment(s) (50 bases on both sides of the variant) : the variant position is highlighted in red.
- The reference genome : the variant position is highlighted in red.
- Genome annotation : the yellow blocks represent exons, the variant position is highlighted in red. The annotation is represented only if the variant is not in an intergenic region.
- Zoom-out on the genome annotation : representation of the whole gene whose name is on the right side of the gene annotation. The annotation is represented only if the variant is not in an intergenic region. A red vertical line shows the position of the variant.
- Genomic axis corresponding to the previous genome annotation. A red vertical line shows the position of the variant.

### Notes

#### Common variants

Needlestack is made to identify rare variants (i.e. only a few samples in your set of samples have a particular variant), because of the robust regression used. Therefore, common SNPs (>10%) or strong somatic hotspots will be missed. The optional `extra_robust_gl` can overcome partially this issue for common germline mutations: it first discard high allelic fraction (>20%, assuming these are likely true variants) before fitting the regression model when between 10% and 50% of sample have such a high allelic fraction. A flag is written in the VCF `INFO/WARN` field when this happened (`EXTRA_ROBUST_GL`). Additionally, when an other allele than the reference allele is the most common, it is taken as the reference allele and a flag is also written in the VCF (`INFO/WARN=INV_REF`).

#### Strand bias

For conventional variant callers, GATK WES/WGS [recommended values](http://gatkforums.broadinstitute.org/discussion/5533/strandoddsratio-computation) for SOR strand bias are SOR < 4 for SNVs and < 10 for indels. We haven't found this particularly useful, but a more detailed evaluation is necessary. For amplicon based targeted sequencing, RVSB>0.85 seems to reduce some erros. There is no hard filter by default as this is easy to do afterward using [bcftools filter](http://samtools.github.io/bcftools/bcftools.html#filter) command.

## Directed Acyclic Graph
<img align="center" src="https://cloud.githubusercontent.com/assets/3366818/15250619/eb6afcea-1925-11e6-9f91-e1d6ceecbe19.jpg" width="600">

## Contributions

  | Name      | Email | Description     |
  |-----------|---------------|-----------------|
  | Matthieu Foll*    |          FollM@iarc.fr | Developer to contact for support (link to specific gitter chatroom) |
  | Tiffany Delhomme*    |            delhommet@students.iarc.fr | Developer |
  |    |           | Tester |

