#! /usr/bin/env nextflow

//vim: syntax=groovy -*- mode: groovy;-*-

// needlestack: a multi-sample somatic variant caller
// Copyright (C) 2017 IARC/WHO

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.


params.help = null
params.input_vcf = null
params.out_vcf = null
params.region = null
params.bed = null
params.out_annotated_vcf = null
params.ref_genome = null
params.min_dp = 30 // minimum median coverage to consider a site
params.min_ao = 3 // minimum number of non-ref reads in at least one sample to consider a site
params.nsplit = 1 // split the positions for calling in nsplit pieces and run in parallel
params.min_qval = 50 // qvalue in Phred scale to consider a variant
params.sb_type = "SOR" // strand bias measure to be used: "SOR" or "RVSB"
if(params.sb_type in ["SOR", "RVSB"] ) {
  params.sb_snv = 100 // strand bias threshold for snv
  params.sb_indel = 100 // strand bias threshold for indels
} else {
  params.sb_snv = 1000 // strand bias threshold for snv
  params.sb_indel = 1000 // strand bias threshold for indels
}
params.power_min_af = -1 // minimum allelic fraction for power computations
params.sigma_normal = 0.1 // sigma parameter for negative binomial modeling germline mutations
params.map_qual = 0 // min mapping quality (passed to samtools)
params.base_qual = 13 // min base quality (passed to samtools)
params.max_dp = 50000 // downsample coverage per sample (passed to samtools)
params.use_file_name = false //put these argument to use the bam file names as sample names and do not to use the sample name filed from the bam files (SM tag)
params.all_SNVs = false //  output all sites, even when no variant is detected
params.extra_robust_gl = false //  perform an extra robust regression basically for germline variants

params.tn_pairs = "FALSE" // by default R will get a false boolean value for tn_pairs option
assert (params.tn_pairs != true) : "please enter a file name when using --tn_pairs option"
if (params.tn_pairs != "FALSE") { try { assert file(params.tn_pairs).exists() : "\n ERROR : input tumor-normal pairs file not located in execution directory, exit" } catch (AssertionError e) { println e.getMessage() ; System.exit(0)} }
pairs_file = file(params.tn_pairs)

if (params.tn_pairs != "FALSE") {
  params.plots = "SOMATIC"  // produce pdf plots of regressions for somatic variants
}else {
  params.plots = "ALL" // produce pdf plots of regressions for all variants
}

params.do_alignments = false // do not produce alignment plots in the pdf
params.no_indels = false // do not skip indels
params.no_labels = false // label outliers
params.no_contours = false // add contours to the plots and plot min(AF)~DP

/* If --help in parameters, print software usage */

if (params.help) {
    log.info ''
    log.info '--------------------------------------------------------'
    log.info 'NEEDLESTACK v1.0b: A MULTI-SAMPLE SOMATIC VARIANT CALLER'
    log.info '--------------------------------------------------------'
    log.info 'Copyright (C) 2015 Matthieu Foll and Tiffany Delhomme'
    log.info 'This program comes with ABSOLUTELY NO WARRANTY; for details see LICENSE.txt'
    log.info 'This is free software, and you are welcome to redistribute it'
    log.info 'under certain conditions; see LICENSE.txt for details.'
    log.info '--------------------------------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info '    nextflow run iarcbioinfo/needlestack [-with-docker] --bed bedfile.bed --bam_folder BAM/ --ref reference.fasta [other options]'
    log.info ''
    log.info 'Mandatory arguments:'
    log.info '    --bam_folder     BAM_DIR                  BAM files directory.'
    log.info '    --ref      REF_IN_FASTA             Reference genome in fasta format.'
    log.info '    OR '
    log.info '    --input_vcf      VCF FILE                 VCF file (basically from GATK pipeline) to annotate.'
    log.info 'Options:'
    log.info '    --nsplit         INTEGER                  Split the region for calling in nsplit pieces and run in parallel.'
    log.info '    --min_dp         INTEGER                  Minimum median coverage (in addition, min_dp in at least 10 samples).'
    log.info '    --min_ao         INTEGER                  Minimum number of non-ref reads in at least one sample to consider a site.'
    log.info '    --min_qval       VALUE                    Qvalue in Phred scale to consider a variant.'
    log.info '    --sb_type        SOR or RVSB              Strand bias measure.'
    log.info '    --sb_snv         VALUE                    Strand bias threshold for SNVs.'
    log.info '    --sb_indel       VALUE                    Strand bias threshold for indels.'
    log.info '    --power_min_af   VALUE                    Minimum allelic fraction for power computations.'
    log.info '    --sigma_normal   VALUE                    Sigma parameter for negative binomial modeling germline mutations.'
    log.info '    --map_qual       VALUE                    Samtools minimum mapping quality.'
    log.info '    --base_qual      VALUE                    Samtools minimum base quality.'
    log.info '    --max_dp         INTEGER                  Samtools maximum coverage before downsampling.'
    log.info '    --use_file_name                           Sample names are taken from file names, otherwise extracted from the bam file SM tag.'
    log.info '    --all_SNVs                                Output all SNVs, even when no variant found.'
    log.info '    --extra_robust_gl                         Perform an extra robust regression, basically for germline variants'
    log.info '    --plots                                Output PDF regression plots.'
    log.info '    --do_alignments                           Add alignment plots.'
    log.info '    --no_labels                               Do not add labels to outliers in regression plots.'
    log.info '    --no_indels                               Do not call indels.'
    log.info '    --no_contours                             Do not add contours to plots and do not plot min(AF)~DP.'
    log.info '    --output_folder     OUTPUT FOLDER            Output directory, by default input bam folder.'
    log.info '    --bed            BED FILE                 A BED file for calling.'
    log.info '    --region         CHR:START-END            A region for calling.'
    log.info '    --tn_pairs     TEXT FILE                A tab-delimited file containing two columns (normal and tumor sample name) for each sample in line.'
    log.info '    --ref_genome                     Reference genome for alignments plot'
    log.info ''
    exit 1
}


/* Software information */
log.info ''
log.info '--------------------------------------------------------'
log.info 'NEEDLESTACK v1.0b: A MULTI-SAMPLE SOMATIC VARIANT CALLER'
log.info '--------------------------------------------------------'
log.info 'Copyright (C) 2015 Matthieu Foll and Tiffany Delhomme'
log.info 'This program comes with ABSOLUTELY NO WARRANTY; for details see LICENSE.txt'
log.info 'This is free software, and you are welcome to redistribute it'
log.info 'under certain conditions; see LICENSE.txt for details.'
log.info '--------------------------------------------------------'
log.info(params.tn_pairs == "FALSE" ? "Perform a tumor-normal somatic variant calling (--tn_pairs)   : no"  : "Perform a tumor-normal somatic variant calling (--tn_pairs)   : yes (file ${params.tn_pairs})" )
log.info "To consider a site for calling:"
log.info "     minimum median coverage (--min_dp)                         : ${params.min_dp}"
log.info "     minimum of alternative reads (--min_ao)                    : ${params.min_ao}"
log.info "Phred-scale qvalue threshold (--min_qval)                       : ${params.min_qval}"

if(params.plots == "ALL"){
	log.info "PDF regression plots (--plots)                               : ALL"
} else if (params.plots == "SOMATIC"){
	log.info "PDF regression plots (--plots)                               : SOMATIC"
} else {
	log.info "PDF regression plots (--plots)                               : NONE"
}

log.info(params.do_alignments == true ? "Alignment plots (--do_alignments)                               : yes"  : "Alignment plots (--do_alignments)                               : no" )
log.info(params.no_labels == true ? "Labeling outliers in regression plots (--no_labels)             : no"  : "Labeling outliers in regression plots (--no_labels)             : yes" )
log.info(params.no_contours == true ? "Add contours in plots and plot min(AF)~DP (--no_contours)       : no"  : "Add contours in plots and plot min(AF)~DP (--no_contours)       : yes" )

if(params.input_vcf) {

  params.output_folder = "annotated_vcf"
  params.chunk_size = 10000
  input_vcf = file(params.input_vcf)
  params.out_annotated_vcf = null
  out_annotated_vcf = params.out_annotated_vcf ? params.out_annotated_vcf : "annotated.vcf"
  assert params.extra_robust_gl in [true,false] : "do not assign a value to --extra_robust_gl"

  log.info "Number of vcf chunks for parallel computing (--nsplit)          : ${params.nsplit}"
  log.info "Size of read chunks by VariantAnnotation (--chunk_size)         : ${params.chunk_size}"
  log.info "Input vcf for annotation by needlestack (--input_vcf)           : ${params.input_vcf}"
  log.info "Output annotated file (--out_annotated_vcf)                     : ${out_annotated_vcf}"
  log.info(params.extra_robust_gl == true ? "Perform an extra-robust regression (--extra_robust_gl)          : yes" : "Perform an extra-robust regression (--extra_robust_gl)          : no" )
  log.info "output folder (--output_folder)                                    : ${params.output_folder}"
  log.info "\n"

  process split_vcf {

    input:
    file input_vcf

    output:
    file 'split*' into splitted_vcf mode flatten

    shell:
    '''
    zcat !{input_vcf} | sed '/^#CHROM/q' | grep -v "<redacted>" > header
    ((nb_total_lines= $((`zcat !{input_vcf} | wc -l`)) ))
    ((core_lines = $nb_total_lines - $((`cat header | wc -l`)) ))
    ((lines_per_file = ( $core_lines + !{params.nsplit} - 1) / !{params.nsplit}))
    ((start=( $((`cat header | wc -l`)) +1 ) ))

    # this works only with split (coreutils) version > 8.13 but is much faster
    zcat !{input_vcf} | tail -n+$start | split -l $lines_per_file -a 10 --filter='{ cat header; cat; } | bgzip > $FILE.gz' - split_

    ## slower but does not require any specific version of split
    #for i in `seq 1 !{params.nsplit}`;
    #    do
    #    if ((start < nb_total_lines)); then
    #        { cat header && zcat !{input_vcf} | tail -n+$start  | head -n$lines_per_file ; } | bgzip > split${i}.vcf.gz
    #        ((start=start+lines_per_file))
    #    fi
    #done
    '''
  }

  process annotate_vcf {

    if(!params.no_plots) {
          publishDir params.output_folder+'/PDF/', mode: 'move', pattern: '*.pdf'
    }

    input:
    file svcf from splitted_vcf

    output:
    file '*.pdf' optional true into PDF
    file '*.vcf' into annotated

    shell:
    '''
    tabix -p vcf !{svcf}
    Rscript !{baseDir}/bin/annotate_vcf.r --source_path=!{baseDir}/bin/ --input_vcf=!{svcf} --chunk_size=!{params.chunk_size} --do_plots=!{params.plots} --plot_labels=!{!params.no_labels} --add_contours=!{!params.no_contours} --min_coverage=!{params.min_dp} --min_reads=!{params.min_ao} --GQ_threshold=!{params.min_qval} --extra_rob=!{params.extra_robust_gl}
    '''
  }

  process merge_vcf {

    publishDir params.output_folder, mode: 'move'

    input:
    val out_annotated_vcf
    file all_vcf from annotated.toList()

    output:
    file "$out_annotated_vcf" into merged_vcf

    shell:
    '''
    # Extract the header from the first VCF
    grep '^#' !{all_vcf[0]} > !{out_annotated_vcf}

    # Add version numbers in the VCF header just after fileformat line
    echo '##NeedlestackCommand=!{workflow.commandLine}' > versions.txt
    echo '##NeedlestackRepository=!{workflow.repository}' >> versions.txt
    echo '##NeedlestackCommitId=!{workflow.commitId}' >> versions.txt
    echo '##NeedlestackRevision=!{workflow.revision}' >> versions.txt
    echo '##NeedlestackContainer=!{workflow.container}' >> versions.txt
    echo '##nextflow=v!{workflow.nextflow.version}' >> versions.txt
    echo '##Rscript='$(Rscript --version 2>&1) >> versions.txt
    sed -i '/##fileformat=.*/ r versions.txt' !{out_annotated_vcf}

    # this is only for the split_vcf process when using the split linux command that ensures files are in the right order
    grep -h -v '^#' split_*.vcf >> !{out_annotated_vcf}
    # this is for the slow version of the split_vcf process
    #for i in `seq 1 !{params.nsplit}`;
    #    do
    #        grep -v '^#' split${i}_annotated_needlestack.vcf >> !{out_annotated_vcf}
    #    done
    '''
  }

} else {

  params.output_folder = params.bam_folder // if not provided, outputs will be held on the input bam folder
  assert (params.ref != true) && (params.ref != null) : "please specify --ref option (--ref reference.fasta(.gz))"
  assert (params.bam_folder != true) && (params.bam_folder != null) : "please specify --bam_folder option (--bam_folder bamfolder)"

  fasta_ref = file( params.ref )
  fasta_ref_fai = file( params.ref+'.fai' )
  fasta_ref_gzi = file( params.ref+'.gzi' )

  /* Verify user inputs are correct */

  assert params.sb_type in ["SOR","RVSB","FS"] : "--sb_type must be SOR, RVSB or FS "
  assert params.all_SNVs in [true,false] : "do not assign a value to --all_SNVs"
  assert params.extra_robust_gl in [true,false] : "do not assign a value to --extra_robust_gl"
  assert params.plots in ["SOMATIC","ALL","NONE"] : "option not reconized for --plots (SOMATIC,ALL or NONE)"
  assert params.do_alignments in [true,false] : "do not assign a value to --no_alignments"
  assert params.no_indels in [true,false] : "do not assign a value to --no_indels"
  assert params.use_file_name in [true,false] : "do not assign a value to --use_file_name"
  if ( (params.plots == "SOMATIC") && (params.tn_pairs == "FALSE") ) {
      println "\n ERROR : --plots can not be set to SOMATIC since no tn_pairs was provided (--tn_pairs option), exit."; System.exit(0)
  }
  if ( (params.plots == "NONE") && (params.do_alignments == true) ) {
      println "\n ERROR : --do_alignments can not be true since --plots is set to NONE, exit."; System.exit(0)
  }
  if ( (params.do_alignments == true) && (params.ref_genome == null) ) {
      println "\n ERROR : --do_alignments is true, --ref_genome can not be null, exit."; System.exit(0)
  }
  if ( (params.do_alignments == false) && (params.ref_genome != null) ) {
    println "\n WARNING : value assign to --ref_genome although do_alignments is false."
  }
  if (params.bed) { try { assert file(params.bed).exists() : "\n WARNING : input bed file not located in execution directory" } catch (AssertionError e) { println e.getMessage() } }
  try { assert fasta_ref.exists() : "\n WARNING : fasta reference not located in execution directory. Make sure reference index is in the same folder as fasta reference" } catch (AssertionError e) { println e.getMessage() }
  if (fasta_ref.exists()) {assert fasta_ref_fai.exists() : "input fasta reference does not seem to have a .fai index (use samtools faidx)"}
  if (fasta_ref.exists() && params.ref.tokenize('.')[-1] == 'gz') {assert fasta_ref_gzi.exists() : "input gz fasta reference does not seem to have a .gzi index (use samtools faidx)"}
  try { assert file(params.bam_folder).exists() : "\n WARNING : input BAM folder not located in execution directory" } catch (AssertionError e) { println e.getMessage() }
  assert file(params.bam_folder).listFiles().findAll { it.name ==~ /.*bam/ }.size() > 0 : "BAM folder contains no BAM"
  if (file(params.bam_folder).exists()) {
      if (file(params.bam_folder).listFiles().findAll { it.name ==~ /.*bam/ }.size() < 10) { println "\n ERROR : BAM folder contains less than 10 BAM, exit."; System.exit(0) }
      else if (file(params.bam_folder).listFiles().findAll { it.name ==~ /.*bam/ }.size() < 20) { println "\n WARNING : BAM folder contains less than 20 BAM, method accuracy not warranted." }
      bamID = file(params.bam_folder).listFiles().findAll { it.name ==~ /.*bam/ }.collect { it.getName() }.collect { it.replace('.bam','') }
      baiID = file(params.bam_folder).listFiles().findAll { it.name ==~ /.*bam.bai/ }.collect { it.getName() }.collect { it.replace('.bam.bai','') }
      assert baiID.containsAll(bamID) : "check that every bam file has an index (.bam.bai)"
  }
  assert (params.min_dp >= 0) : "minimum coverage must be higher than or equal to 0 (--min_dp)"
  assert (params.max_dp > 1) : "maximum coverage before downsampling must be higher than 1 (--max_dp)"
  assert (params.min_ao >= 0) : "minimum alternative reads must be higher than or equal to 0 (--min_ao)"
  assert (params.nsplit > 0) : "number of regions to split must be higher than 0 (--nsplit)"
  assert (params.min_qval >= 0) : "minimum Phred-scale qvalue must be higher than or equal to 0 (--min_qval)"
  assert ( (params.power_min_af > 0 && params.power_min_af <= 1) || params.power_min_af == -1 ) : "minimum allelic fraction for power computations must be in [0,1] (--power_min_af)"
  assert (params.sigma_normal >= 0) : " sigma parameter for negative binomial must be positive (--sigma_normal)"
  if(params.sb_type in ["SOR", "RVSB"] ) {
  assert (params.sb_snv > 0 && params.sb_snv < 101) : "strand bias (SOR or RVSB) for SNVs must be in [0,100]"
  assert (params.sb_indel > 0 && params.sb_indel < 101) : "strand bias (SOR or RVSB) for indels must be in [0,100]"
  } else {
  assert (params.sb_snv > 0 && params.sb_snv < 1001) : "strand bias (FS) for SNVs must be in [0,1000]"
  assert (params.sb_indel > 0 && params.sb_indel < 1001) : "strand bias (FS) for indels must be in [0,1000]"
  }

  assert (params.map_qual >= 0) : "minimum mapping quality (samtools) must be higher than or equal to 0"
  assert (params.base_qual >= 0) : "minimum base quality (samtools) must be higher than or equal to 0"

  sample_names = params.use_file_name ? "FILE" : "BAM"
  out_vcf = params.out_vcf ? params.out_vcf : "all_variants.vcf"

  /* manage input positions to call (bed or region or whole-genome) */
  if (params.region) {
      input_region = 'region'
  } else if (params.bed) {
      input_region = 'bed'
      bed = file(params.bed)
  } else {
      input_region = 'whole_genome'
  }

  log.info "Input BAM folder (--bam_folder)                                 : ${params.bam_folder}"
  log.info "output folder (--output_folder)                                    : ${params.output_folder}"
  log.info "Reference in fasta format (--ref)                         : ${params.ref}"
  log.info "Intervals for calling (--bed)                                   : ${input_region}"
  log.info "Number of regions to split (--nsplit)                           : ${params.nsplit}"
  log.info "Strand bias measure (--sb_type)                                 : ${params.sb_type}"
  log.info "Strand bias threshold for SNVs (--sb_snv)                       : ${params.sb_snv}"
  log.info "Strand bias threshold for indels (--sb_indel)                   : ${params.sb_indel}"
  log.info "Minimum allelic fraction for power computations (--power_min_af): ${params.power_min_af}"
  log.info "Sigma parameter for germline (--sigma)                          : ${params.sigma_normal}"
  log.info "Samtools minimum mapping quality (--map_qual)                   : ${params.map_qual}"
  log.info "Samtools minimum base quality (--base_qual)                     : ${params.base_qual}"
  log.info "Samtools maximum coverage before downsampling (--max_dp)        : ${params.max_dp}"
  log.info "Sample names definition (--use_file_name)                       : ${sample_names}"
  log.info(params.all_SNVs == true ? "Output all SNVs (--all_SNVs)                                    : yes" : "Output all SNVs (--all_SNVs)                                    : no" )
  log.info(params.extra_robust_gl == true ? "Perform an extra-robust regression (--extra_robust_gl)          : yes" : "Perform an extra-robust regression (--extra_robust_gl)          : no" )
  log.info(params.no_indels == true ? "Skip indels (--no_indels)                                       : yes" : "Skip indels (--no_indels)                                       : no" )
  log.info "\n"

  bam = Channel.fromPath( params.bam_folder+'/*.bam' ).toList()
  bai = Channel.fromPath( params.bam_folder+'/*.bam.bai' ).toList()

  /* Building the bed file where calling would be done */
  process bed {
      output:
      file "temp.bed" into outbed

      shell:
      if (input_region == 'region')
      '''
      echo !{params.region} | sed -e 's/[:|-]/	/g' > temp.bed
      '''

      else if (input_region == 'bed')
      '''
      ln -s !{bed} temp.bed
      '''

      else if (input_region == 'whole_genome')
      '''
      cat !{fasta_ref_fai} | awk '{print $1"	"0"	"$2 }' > temp.bed
      '''
  }


  /* split bed file into nsplit regions */
  process split_bed {

      input:
      file bed from outbed

      output:
      file '*_regions' into split_bed mode flatten

      shell:
      '''
      grep -v '^track' !{bed} | sort -k1,1 -k2,2n | bedtools merge -i stdin | awk '{print $1" "$2" "$3}' | bed_cut.r !{params.nsplit}
      '''
  }


  // create mpileup file + parse mpileup file to send it to Rscript
  process mpileup2vcf {

	  if(params.plots) {
          publishDir params.output_folder+'/PDF/', mode: 'move', pattern: '*.pdf'
      }

      tag { region_tag }

      input:
      file split_bed
      file 'BAM/*' from bam
      file 'BAM/*' from bai
      file fasta_ref
      file fasta_ref_fai
      file fasta_ref_gzi
      file pairs_file

      output:
      file "${region_tag}.vcf" into vcf
      file '*.pdf' optional true into PDF


      shell:
      region_tag = split_bed.baseName
      if ( params.no_indels ) {
          indel_par = "true"
      } else {
          indel_par = "false"
      }

      '''
      for cur_bam in BAM/*.bam
      do
          if [ "!{sample_names}" == "FILE" ]; then
              # use bam file name as sample name
              bam_file_name=$(basename "${cur_bam%.*}")
              # remove whitespaces from name
              SM1="$(echo -e "${bam_file_name}" | tr -d '[[:space:]]')"
              SM2="$(echo -e "${bam_file_name}" | tr -d '[[:space:]]')"
          else
              # get bam file names
	      bam_file_name=$(basename "${cur_bam%.*}")
	      # remove whitespaces from name
	      SM1="$(echo -e "${bam_file_name}" | tr -d '[[:space:]]')"
	      # extract sample name from bam file read group info field
	      SM2=$(samtools view -H $cur_bam | grep "^@RG" | tail -n1 | sed "s/.*SM:\\([^	]*\\).*/\\1/" | tr -d '[:space:]')
          fi


          printf "$SM1	$SM2\\n" >> names.txt
      done

      if [ "!{params.tn_pairs}" != "FALSE" ]; then
          abs_pairs_file=$(readlink -f !{pairs_file})
      else
          abs_pairs_file="FALSE"
      fi
      set -o pipefail
      i=1

      { while read bed_line; do
          samtools mpileup --fasta-ref !{fasta_ref} --region $bed_line --ignore-RG --min-BQ !{params.base_qual} --min-MQ !{params.map_qual} --max-idepth 1000000 --max-depth !{params.max_dp} BAM/*.bam | sed 's/		/	*	*/g'
          i=$((i+1))
      done < !{split_bed}
      } | mpileup2readcounts 0 -5 !{indel_par} !{params.min_ao} | Rscript !{baseDir}/bin/needlestack.r --pairs_file=${abs_pairs_file} --source_path=!{baseDir}/bin/ --out_file=!{region_tag}.vcf --fasta_ref=!{fasta_ref} --bam_folder=BAM/ --ref_genome=!{params.ref_genome} --GQ_threshold=!{params.min_qval} --min_coverage=!{params.min_dp} --min_reads=!{params.min_ao} --SB_type=!{params.sb_type} --SB_threshold_SNV=!{params.sb_snv} --SB_threshold_indel=!{params.sb_indel} --output_all_SNVs=!{params.all_SNVs} --do_plots=!{params.plots} --do_alignments=!{params.do_alignments} --plot_labels=!{!params.no_labels} --add_contours=!{!params.no_contours} --extra_rob=!{params.extra_robust_gl} --afmin_power=!{params.power_min_af} --sigma=!{params.sigma_normal}
      '''
  }


  // merge all vcf files in one big file
  process collect_vcf_result {

      publishDir params.output_folder, mode: 'move'

      input:
      val out_vcf
      file all_vcf from vcf.toList()
      file fasta_ref_fai

      output:
      file "$out_vcf" into big_vcf

      when:
      !all_vcf.empty

      shell:
      '''
      mkdir VCF
      mv *.vcf VCF
      # Extract the header from the first VCF
      sed '/^#CHROM/q' VCF/!{all_vcf[0]} > header.txt

      # Add contigs in the VCF header
      cat !{fasta_ref_fai} | cut -f1,2 | sed -e 's/^/##contig=<ID=/' -e 's/[	 ][	 ]*/,length=/' -e 's/$/>/' > contigs.txt
      sed -i '/##reference=.*/ r contigs.txt' header.txt

      # Add version numbers in the VCF header
      echo '##command=!{workflow.commandLine}' > versions.txt
      echo '##repository=!{workflow.repository}' >> versions.txt
      echo '##commitId=!{workflow.commitId}' >> versions.txt
      echo '##revision=!{workflow.revision}' >> versions.txt
      echo '##container=!{workflow.container}' >> versions.txt
      echo '##nextflow=v!{workflow.nextflow.version}' >> versions.txt
      echo '##samtools='$(samtools --version | tr '\n' ' ') >> versions.txt
      echo '##bedtools='$(bedtools --version) >> versions.txt
      echo '##Rscript='$(Rscript --version 2>&1) >> versions.txt
      sed -i '/##source=.*/ r versions.txt' header.txt

      # Check if sort command allows sorting in natural order (chr1 chr2 chr10 instead of chr1 chr10 chr2)
      if [ `sort --help | grep -c 'version-sort' ` == 0 ]
      then
          sort_ops="-k1,1d"
      else
          sort_ops="-k1,1V"
      fi
      # Add all VCF contents and sort
      grep --no-filename -v '^#' VCF/*.vcf | LC_ALL=C sort -t '	' $sort_ops -k2,2n >> header.txt
      mv header.txt !{out_vcf}
      '''
  }
}
