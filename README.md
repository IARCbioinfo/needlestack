<img align="center" src="https://cloud.githubusercontent.com/assets/3366818/12507528/99afd6fa-c0f6-11e5-993b-5ac925d178fa.png" width="700">
# A multi-sample somatic variant caller

<img align="center" src="https://cloud.githubusercontent.com/assets/3366818/10489240/7da79144-729c-11e5-8cb3-0225106d9b06.jpg" width="600">

Warning: development in progress, unreliable results warrantied. 

Please wait upcoming publication before using it in production. 

Contact: follm@iarc.fr

## Input

- A set of [BAM files](https://samtools.github.io/hts-specs/) (called `*.bam`) grouped in a single folder along with their [index files](http://www.htslib.org/doc/samtools.html) (called `*.bam.bai`). A minimum of 20 BAM files is recommended. 
- A [bed file](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) .
- A reference [fasta file](https://en.wikipedia.org/wiki/FASTA_format) (eventually compressed with [bgzip](http://www.htslib.org/doc/tabix.html)) along with its [faidx index](http://www.htslib.org/doc/faidx.html) (and `*.gzi` faidx index if compressed).

## Quick start

Needlestack works under most Linux distributions and Apple OS X.

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
	```bash
	cd NGS_data_test/1000G_CEU_TP53/
	nextflow run iarcbioinfo/needlestack -with-docker iarcbioinfo/needlestack \
	         --bed TP53_all.bed --bam_folder BAM/ --fasta_ref 17.fasta.gz
	```
	
	You will find a [VCF file](https://samtools.github.io/hts-specs/) called `all_variants.vcf` in the `BAM/` folder once done.
	
	The first time it will take more time as the pipeline will be downloaded from github and the docker container from [dockerhub](https://hub.docker.com/r/iarcbioinfo/needlestack/).

	Creating an alias for the long command above can be useful. For example:
	```bash
	alias needlestack='nextflow run iarcbioinfo/needlestack -with-docker iarcbioinfo/needlestack'
	```
	
	If you want to permanantly add this alias (and not just for your current session), add the above  line to your `~/.bashrc` file (assuming you are using bash).
	
	It will allow you to do this:
	```bash
	needlestack --bed TP53_all.bed --bam_folder BAM/ --fasta_ref 17.fasta.gz
	```
	
	See below for a more complete alias we recommend using.
	
6. Update the pipeline

	You can update the nextflow sofware and the pipeline itself simply using:
	```bash
	nextflow -self-update
	nextflow pull iarcbioinfo/needlestack
	```

	You can also automatically update the pipeline when you run it by adding the option `-latest` in the `nextflow run` command. Doing so you will always run the latest version from [Github](https://github.com/iarcbioinfo/needlestack).

	Official releases can be found [here](https://github.com/iarcbioinfo/needlestack/releases/). There is a corresponding official [docker container](https://hub.docker.com/r/iarcbioinfo/needlestack/) for each release and one can run a particular version using (for example for v0.2):
	```bash
	nextflow run iarcbioinfo/needlestack -r v0.2 -with-docker iarcbioinfo/needlestack:v0.2 \
	         --bed TP53_all.bed --bam_folder BAM/ --fasta_ref 17.fasta.gz
	```

## Detailed instructions

If you can't install [docker](https://www.docker.com) or don't want to use it, the pipeline will also work if you install [perl](https://www.perl.org),  [bedtools](http://bedtools.readthedocs.org/en/latest/), [samtools](http://www.htslib.org), vcfoverlay from [vcflib](https://github.com/ekg/vcflib) and Rscript from [R](https://www.r-project.org) and put them in your path (executables are assumed to be respectively called `perl`, `bedtools`, `samtools`, `vcflib` and `Rscript`). In this case, remove the `-with-docker` option from step 5 above.

The exact same pipeline can be run on your computer or on a HPC cluster, by adding a [nextflow configuration file](http://www.nextflow.io/docs/latest/config.html) to choose an appropriate [executor](http://www.nextflow.io/docs/latest/executor.html). For example to work on a cluster using [SGE scheduler](https://en.wikipedia.org/wiki/Oracle_Grid_Engine), simply add a file named `nextflow.config` in the current directory (or `~/.nextflow/config` to make global changes) containing:  
```java
process.executor = 'sge'
```

Other popular schedulers such as LSF, SLURM, PBS, TORQUE etc. are also compatible. See the nextflow documentation [here](http://www.nextflow.io/docs/latest/executor.html) for more details. Also have a look at the [other parameters for the executors](http://www.nextflow.io/docs/latest/config.html#scope-executor), in particular `queueSize` that defines the number of tasks the executor will handle in a parallel manner.  

The default number of tasks the executor will handle in a parallel is 100, which is certainly too high if you are executing it on your local machine. In this case a good idea is to set it to the number of computing cores your local machine has. Following is an example to create a config file with this information automatically (works on Linux and Mac OS X):
```bash
echo "executor.\$local.queueSize = "`getconf _NPROCESSORS_ONLN` > ~/.nextflow/config
```

Replace `>` by `>>` if you want to add the argument line to an existing nextflow config file.


`--bed`, `--bam_folder` and `--fasta_ref` are compulsary. The optional parameters with default values are:

| Parameter | Default value | Description |
|-----------|--------------:|-------------| 
| min_dp    |            50 | Minimum coverage in at least one sample to consider a site |
| min_ao | 5 | Minimum number of non-ref reads in at least one sample to consider a site |
| nsplit | 1 | Split the bed file in nsplit pieces and run in parallel |
| min_qval | 50 | qvalue threshold in [Phred scale](https://en.wikipedia.org/wiki/Phred_quality_score) to consider a variant |
| sb_type | SOR | Strand bias measure, either SOR or RVSB |
| sb_snv | 100 | Strand bias threshold for SNVs (100 =no filter) |
| sb_indel | 100 | Strand bias threshold for indels (100 = no filter)|
| map_qual | 20 | Min mapping quality (passed to samtools) |
| base_qual | 20 | Min base quality (passed to samtools) |
| max_DP | 30000 | Downsample coverage per sample (passed to samtools) |
| use_file_name |   | Put this argument to use the bam file names as sample names. By default the sample name is extracted from the bam file SM tag. |
| all_SNVs |   | Put this argument to output all SNVs, even when no variant is detected |
| no_plots |   | Put this argument to remove pdf plots of regressions from the output |
| no_indels |   | Put this argument to do not perform the variant calling on insertions and deletions |
| out_folder | --bam_folder | Output folder, by default equals to the input bam folder |
| out_vcf | all_variants.vcf | File name of final VFC. |

Simply add the parameters you want in the command line like `--min_dp 1000` for exmaple to change the min coverage.

[Recommended values](http://gatkforums.broadinstitute.org/discussion/5533/strandoddsratio-computation) for SOR strand bias are SOR < 4 for SNVs and < 10 for indels. For RVSB, a good starting point is to filter out variant with RVSB>0.85. There is no hard filter by default as this is easy to do afterward using [bcftools filter](http://samtools.github.io/bcftools/bcftools.html#filter) command.

A good practice is to keep (and publish) the `.nextflow.log` file create during the pipeline process, as it contains useful information for reproducibility (full command line, software versions etc.). You should also add the option `-with-trace` in the `nextflow run` command line that will create an additional `trace.csv` file containing even more information to keep for records. The option `-with-timeline` also creates a nice processes execution timeline file (a web page).

All in all, our recommended alias for the full command line to run needlestack is:
```bash
alias needlestack='nextflow run iarcbioinfo/needlestack -with-docker iarcbioinfo/needlestack -latest -with-trace --with-timeline'
```

[![Join the chat at https://gitter.im/iarcbioinfo/needlestack](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/iarcbioinfo/needlestack?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

[![Circle CI](https://circleci.com/gh/IARCbioinfo/needlestack/tree/master.svg?style=shield&circle-token=402d456a7c50af352bb4e1a52425ce0fe645f78f)](https://circleci.com/gh/IARCbioinfo/needlestack/tree/master)

[![Docker Hub](https://img.shields.io/badge/docker-ready-blue.svg)](https://hub.docker.com/r/iarcbioinfo/needlestack/)
