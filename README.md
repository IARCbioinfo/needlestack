# needlstack: multi-sample somatic variant caller

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
	```bash
	cd NGS_data_test/1000G_CEU_TP53/
	nextflow run mfoll/needlestack -with-docker mfoll/robust-regression-caller \
	         --bed TP53_all.bed --bam_folder BAM/ --fasta_ref 17.fasta.gz
	```
	
	You will find a [VCF file](https://samtools.github.io/hts-specs/) called `all_variants.vcf` in the `BAM/` folder once done.
	
	The first time it will take more time as the pipeline will be downloaded from github and the docker container from [dockerhub](https://hub.docker.com/r/mfoll/robust-regression-caller/).

	Creating an alias for the long command above can be useful. For example:
	```bash
	alias rrcaller='nextflow run mfoll/needlestack -with-docker mfoll/robust-regression-caller'
	```
	
	If you want to permanantly add this alias (and not just for your current session), add the above  line to your `~/.bashrc` file (assuming you are using bash).
	
	Will allow you to do this:
	```bash
	rrcaller --bed TP53_all.bed --bam_folder BAM/ --fasta_ref 17.fasta.gz
	```
6. Update the pipeline

	You can update the nextflow sofware and the pipeline itself simply using:
	```bash
	nextflow -self-update
	nextflow pull mfoll/needlestack
	```

	You can also automatically update the pipeline when you run it by adding the option `-latest` in the `nextflow run` command. Doing so you will always run the latest version from [Github](https://github.com/mfoll/needlestack).

	Official releases can be found [here](https://github.com/mfoll/needlestack/releases/). There is a corresponding official [docker container](https://hub.docker.com/r/mfoll/robust-regression-caller/) for each release and one can run a particular version using (for example for v0.1):
	```bash
	nextflow run mfoll/needlestack -r v0.1 -with-docker mfoll/robust-regression-caller:v0.1 \
	         --bed TP53_all.bed --bam_folder BAM/ --fasta_ref 17.fasta.gz
	```

## Detailed instructions

If you can't install [docker](https://www.docker.com) or don't want to use it, the pipeline will also work if you install [perl](https://www.perl.org),  [bedtools](http://bedtools.readthedocs.org/en/latest/), [samtools](http://www.htslib.org), vcfoverlay from [vcflib](https://github.com/ekg/vcflib) and Rscript from [R](https://www.r-project.org) and put them in your path (executables are assumed to be respectively called `perl`, `bedtools`, `samtools`, `vcflib` and `Rscript`). In this case, remove the `-with-docker` option from step 5 above.

The exact same pipeline can be run on your computer or on a HPC cluster, by adding a [nextflow configuration file](http://www.nextflow.io/docs/latest/config.html) to choose an appropriate [executor](http://www.nextflow.io/docs/latest/executor.html). For example to work on a cluster using [SGE scheduler](https://en.wikipedia.org/wiki/Oracle_Grid_Engine), simply add a file named `nextflow.config` in the current directory (or `~/.nextflow/config` to make global changes) containing:  
```java
process.executor = 'sge'
```

Other popular schedulers such as LSF, SLURM, PBS, TORQUE etc. are also compatible. See the nextflow documentation [here](http://www.nextflow.io/docs/latest/executor.html) for more details. Also have a look at the [other parameters for the executors](http://www.nextflow.io/docs/latest/config.html#scope-executor), in particular `queueSize` that defines the number of tasks the executor will handle in a parallel manner (default is 100 which is certainly too high if you are executing it on your local machine).

Following is an example of creating a global nextflow config file, with `queueSize` equals to your number of processors (works on linux platforms and Mac OS X), replace `>` by `>>` if you want to add the argument line to an existing nextflow config file :
```bash
echo "executor.\$local.queueSize = "`getconf _NPROCESSORS_ONLN` > ~/.nextflow/config
```

`--bed`, `--bam_folder` and `--fasta_ref` are compulsary. The optional parameters with default values are:

| Parameter | Default value | Description |
|-----------|--------------:|-------------| 
| min_dp    |            50 | Minimum coverage in at least one sample to consider a site |
| min_ao | 5 | Minimum number of non-ref reads in at least one sample to consider a site|
| nsplit | 1 | Split the bed file in nsplit pieces and run in parallel |
| min_qval | 50 | qvalue in Phred scale to consider a variant |
| sb_type | "SOR" | Strand bias measure, either "SOR" or "RVSB" |
| sb_snv | 100 | Strand bias threshold for SNVs (100 =no filter) |
| sb_indel | 100 | Strand bias threshold for indels (100 = no filter)|
| map_qual | 20 | Min mapping quality (passed to samtools) |
| base_qual | 20 | Min base quality (passed to samtools) |
| max_DP | 30000 | Downsample coverage per sample (passed to samtools) |
| use_file_name |   | Put this argument to use the bam file names as sample names and do not to use the sample name filed from the bam file (SM tag) |
| all_SNVs |   | Put this argument to output all SNVs sites, even when no variant is detected |
| no_plots |   | Put this argument to do not produce pdf plots of regressions |
| out_folder | --bam_folder | Output folder, by default equals to the input bam folder |

Simply add the parameters you want in the command line like `--min_dp 1000` for exmaple to change the min coverage.

[Recommended values](http://gatkforums.broadinstitute.org/discussion/5533/strandoddsratio-computation) for SOR strand bias are SOR < 4 for SNVs and < 10 for indels. For RVSB, a good starting point is to filter out variant with RVSB>0.85. There is no hard filter by default as this is easy to do afterward using [bcftools filter](http://samtools.github.io/bcftools/bcftools.html#filter) command.

A good practice is to keep (and publish) the `.nextflow.log` file create during the pipeline process, as it contains useful information for reproducibility (full command line, software versions etc.). You should also add the option `-with-trace` in the `nextflow run` command line that will create an additional `trace.csv` file containing even more information to keep for records.

[![Join the chat at https://gitter.im/mfoll/needlestack](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/mfoll/needlestack?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

[![Circle CI](https://circleci.com/gh/mfoll/needlestack/tree/master.svg?style=shield)](https://circleci.com/gh/mfoll/needlestack/tree/master) 

[![Docker Hub](https://img.shields.io/badge/docker-ready-blue.svg)](https://hub.docker.com/r/mfoll/robust-regression-caller/)
