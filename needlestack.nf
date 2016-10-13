#! /usr/bin/env nextflow

// needlestack: a multi-sample somatic variant caller
// Copyright (C) 2015 Matthieu Foll and Tiffany Delhomme

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

// requirement:
// - bedtools
// - samtools
// - Rscript (R)
// - bed_cut.r (in bin folder)
// - needlestack.r (in bin folder)
// - pileup2baseindel.pl (in bin folder) (+ perl)

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
params.map_qual = 0 // min mapping quality (passed to samtools)
params.base_qual = 13 // min base quality (passed to samtools)
params.max_DP = 50000 // downsample coverage per sample (passed to samtools)
params.use_file_name = false //put these argument to use the bam file names as sample names and do not to use the sample name filed from the bam files (SM tag)
params.all_SNVs = false //  output all sites, even when no variant is detected
params.extra_robust_gl = false //  perform an extra robust regression basically for germline variants
params.no_plots = false  // do not produce pdf plots of regressions
params.no_indels = false // do not skip indels
params.no_labels = false // label outliers
params.no_contours = false // add contours to the plots and plot min(AF)~DP
params.pairs_file = "FALSE" // by default R will get a false boolean value for pairs_file option
assert (params.pairs_file != true) : "please enter a file name when using --pairs_file option"
if (params.pairs_file != "FALSE") { try { assert file(params.pairs_file).exists() : "\n WARNING : input tumor-normal pairs file not located in execution directory" } catch (AssertionError e) { println e.getMessage() } }
pairs_file = file(params.pairs_file)

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
    log.info '    nextflow run iarcbioinfo/needlestack [-with-docker] --bed bedfile.bed --bam_folder BAM/ --fasta_ref reference.fasta [other options]'
    log.info ''
    log.info 'Mandatory arguments:'
    log.info '    --bam_folder     BAM_DIR                  BAM files directory.'
    log.info '    --fasta_ref      REF_IN_FASTA             Reference genome in fasta format.'
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
    log.info '    --map_qual       VALUE                    Samtools minimum mapping quality.'
    log.info '    --base_qual      VALUE                    Samtools minimum base quality.'
    log.info '    --max_DP         INTEGER                  Samtools maximum coverage before downsampling.'
    log.info '    --use_file_name                           Sample names are taken from file names, otherwise extracted from the bam file SM tag.'
    log.info '    --all_SNVs                                Output all SNVs, even when no variant found.'
    log.info '    --extra_robust_gl                         Perform an extra robust regression, basically for germline variants'
    log.info '    --no_plots                                Do not output PDF regression plots.'
    log.info '    --no_labels                               Do not add labels to outliers in regression plots.'
    log.info '    --no_indels                               Do not call indels.'
    log.info '    --no_contours                             Do not add contours to plots and do not plot min(AF)~DP.'
    log.info '    --out_folder     OUTPUT FOLDER            Output directory, by default input bam folder.'
    log.info '    --bed            BED FILE                 A BED file for calling.'
    log.info '    --region         CHR:START-END            A region for calling.'
    log.info '    --pairs_file     TEXT FILE                A tab-delimited file containing two columns (normal and tumor sample name) for each sample in line.'
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
log.info(params.pairs_file == "FALSE" ? "Perform a tumor-normal somatic variant calling (--pairs_file)   : no"  : "Perform a tumor-normal somatic variant calling (--pairs_file)   : yes (file ${params.pairs_file})" )
log.info "To consider a site for calling:"
log.info "     minimum median coverage (--min_dp)                         : ${params.min_dp}"
log.info "     minimum of alternative reads (--min_ao)                    : ${params.min_ao}"
log.info "Phred-scale qvalue threshold (--min_qval)                       : ${params.min_qval}"
log.info(params.no_plots == true ? "PDF regression plots (--no_plots)                               : no"  : "PDF regression plots (--no_plots)                               : yes" )
log.info(params.no_labels == true ? "Labeling outliers in regression plots (--no_labels)             : no"  : "Labeling outliers in regression plots (--no_labels)             : yes" )
log.info(params.no_contours == true ? "Add contours in plots and plot min(AF)~DP (--no_contours)       : no"  : "Add contours in plots and plot min(AF)~DP (--no_contours)       : yes" )

if(params.input_vcf) {

  params.out_folder = "annotated_vcf"
  params.chunk_size = 10000
  input_vcf = file(params.input_vcf)
  out_annotated_vcf = params.out_annotated_vcf ? params.out_annotated_vcf : "annotated.vcf"
  assert params.extra_robust_gl in [true,false] : "do not assign a value to --extra_robust_gl"

  log.info "Number of vcf chunks for parallel computing (--nsplit)          : ${params.nsplit}"
  log.info "Size of read chunks by VariantAnnotation (--chunk_size)         : ${params.chunk_size}"
  log.info "Input vcf for annotation by needlestack (--input_vcf)           : ${params.input_vcf}"
  log.info "Output annotated file (--out_annotated_vcf)                     : ${out_annotated_vcf}"
  log.info(params.extra_robust_gl == true ? "Perform an extra-robust regression (--extra_robust_gl)          : yes" : "Perform an extra-robust regression (--extra_robust_gl)          : no" )
  log.info "output folder (--out_folder)                                    : ${params.out_folder}"
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

    publishDir params.out_folder+'/PDF/', mode: 'move',  pattern: "*[ATCG-]*.pdf"

    input:
    file svcf from splitted_vcf

    output:
    file '*.pdf' into PDF
    file '*.vcf' into annotated

    shell:
    '''
    touch empty.pdf
    tabix -p vcf !{svcf}
    Rscript !{baseDir}/bin/annotate_vcf.r --source_path=!{baseDir}/bin/ --input_vcf=!{svcf} --chunk_size=!{params.chunk_size} --do_plots=!{!params.no_plots} --plot_labels=!{!params.no_labels} --add_contours=!{!params.no_contours} --min_coverage=!{params.min_dp} --min_reads=!{params.min_ao} --GQ_threshold=!{params.min_qval} --extra_rob=!{params.extra_robust_gl}
    '''
  }

  process merge_vcf {

    publishDir params.out_folder, mode: 'move'

    input:
    val out_annotated_vcf
    file all_vcf from annotated.toList()

    output:
    file "$out_annotated_vcf" into merged_vcf

    shell:
    '''
    # Extract the header from the first VCF
    grep '^#' !{all_vcf[0]} > !{out_annotated_vcf}
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

  params.out_folder = params.bam_folder // if not provided, outputs will be held on the input bam folder
  assert (params.fasta_ref != true) && (params.fasta_ref != null) : "please specify --fasta_ref option (--fasta_ref reference.fasta(.gz))"
  assert (params.bam_folder != true) && (params.bam_folder != null) : "please specify --bam_folder option (--bam_folder bamfolder)"

  fasta_ref = file( params.fasta_ref )
  fasta_ref_fai = file( params.fasta_ref+'.fai' )
  fasta_ref_gzi = file( params.fasta_ref+'.gzi' )

  /* Verify user inputs are correct */

  assert params.sb_type in ["SOR","RVSB","FS"] : "--sb_type must be SOR, RVSB or FS "
  assert params.all_SNVs in [true,false] : "do not assign a value to --all_SNVs"
  assert params.extra_robust_gl in [true,false] : "do not assign a value to --extra_robust_gl"
  assert params.no_plots in [true,false] : "do not assign a value to --no_plots"
  assert params.no_indels in [true,false] : "do not assign a value to --no_indels"
  assert params.use_file_name in [true,false] : "do not assign a value to --use_file_name"
  if (params.bed) { try { assert file(params.bed).exists() : "\n WARNING : input bed file not located in execution directory" } catch (AssertionError e) { println e.getMessage() } }
  try { assert fasta_ref.exists() : "\n WARNING : fasta reference not located in execution directory. Make sure reference index is in the same folder as fasta reference" } catch (AssertionError e) { println e.getMessage() }
  if (fasta_ref.exists()) {assert fasta_ref_fai.exists() : "input fasta reference does not seem to have a .fai index (use samtools faidx)"}
  if (fasta_ref.exists() && params.fasta_ref.tokenize('.')[-1] == 'gz') {assert fasta_ref_gzi.exists() : "input gz fasta reference does not seem to have a .gzi index (use samtools faidx)"}
  try { assert file(params.bam_folder).exists() : "\n WARNING : input BAM folder not located in execution directory" } catch (AssertionError e) { println e.getMessage() }
  assert file(params.bam_folder).listFiles().findAll { it.name ==~ /.*bam/ }.size() > 0 : "BAM folder contains no BAM"
  if (file(params.bam_folder).exists()) {
      if (file(params.bam_folder).listFiles().findAll { it.name ==~ /.*bam/ }.size() < 10) { println "\n ERROR : BAM folder contains less than 10 BAM, exit."; System.exit(0) }
      else if (file(params.bam_folder).listFiles().findAll { it.name ==~ /.*bam/ }.size() < 20) { println "\n WARNING : BAM folder contains less than 20 BAM, method accuracy not warranted." }
      bamID = file(params.bam_folder).listFiles().findAll { it.name ==~ /.*bam/ }.collect { it.getName() }.collect { it.replace('.bam','') }
      baiID = file(params.bam_folder).listFiles().findAll { it.name ==~ /.*bam.bai/ }.collect { it.getName() }.collect { it.replace('.bam.bai','') }
      assert baiID.containsAll(bamID) : "check that every bam file has an index (.bam.bai)"
  }
  assert (params.min_dp > 0) : "minimum coverage must be higher than 0 (--min_dp)"
  assert (params.max_DP > 1) : "maximum coverage before downsampling must be higher than 1 (--max_DP)"
  assert (params.min_ao >= 0) : "minimum alternative reads must be higher than or equal to 0 (--min_ao)"
  assert (params.nsplit > 0) : "number of regions to split must be higher than 0 (--nsplit)"
  assert (params.min_qval >= 0) : "minimum Phred-scale qvalue must be higher than or equal to 0 (--min_qval)"
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
  log.info "output folder (--out_folder)                                    : ${params.out_folder}"
  log.info "Reference in fasta format (--fasta_ref)                         : ${params.fasta_ref}"
  log.info "Intervals for calling (--bed)                                   : ${input_region}"
  log.info "Number of regions to split (--nsplit)                           : ${params.nsplit}"
  log.info "Strand bias measure (--sb_type)                                 : ${params.sb_type}"
  log.info "Strand bias threshold for SNVs (--sb_snv)                       : ${params.sb_snv}"
  log.info "Strand bias threshold for indels (--sb_indel)                   : ${params.sb_indel}"
  log.info "Samtools minimum mapping quality (--map_qual)                   : ${params.map_qual}"
  log.info "Samtools minimum base quality (--base_qual)                     : ${params.base_qual}"
  log.info "Samtools maximum coverage before downsampling (--max_DP)        : ${params.max_DP}"
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

  // create mpileup file + sed to have "*" when there is no coverage (otherwise pileup2baseindel.pl is unhappy)
  process samtools_mpileup {

      tag { region_tag }

      input:
      file split_bed
      file 'BAM/*' from bam
      file 'BAM/*' from bai
      file fasta_ref
      file fasta_ref_fai
      file fasta_ref_gzi

      output:
      set val(region_tag), file("${region_tag}.pileup") into pileup

      shell:
      region_tag = split_bed.baseName
      '''
      set -o pipefail
      while read bed_line; do
          samtools mpileup --fasta-ref !{fasta_ref} --region $bed_line --ignore-RG --min-BQ !{params.base_qual} --min-MQ !{params.map_qual} --max-idepth 1000000 --max-depth !{params.max_DP} BAM/*.bam | sed 's/		/	*	*/g' >> !{region_tag}.pileup
      done < !{split_bed}
      '''
  }

  // split mpileup file and convert to table
  process mpileup2table {

      tag { region_tag }

      input:
      set val(region_tag), file("${region_tag}.pileup") from pileup.filter { tag, file -> !file.isEmpty() }
      file 'BAM/*' from bam
      val sample_names

      output:
      set val(region_tag), file('TABLE/sample*.txt'), file('names.txt') into table

      shell:
      if ( params.no_indels ) {
          indel_par = "-no-indels"
      } else {
          indel_par = " "
      }
      '''
      nb_pos=$(wc -l < !{region_tag}.pileup)
      if [ $nb_pos -gt 0 ]; then
          # split and convert pileup file
          pileup2baseindel.pl -i !{region_tag}.pileup !{indel_par}
          mkdir TABLE
          mv sample*.txt TABLE
          i=1
          for cur_bam in BAM/*.bam
          do
              if [ "!{sample_names}" == "FILE" ]; then
                  # use bam file name as sample name
                  bam_file_name=$(basename "${cur_bam%.*}")
                  # remove whitespaces from name
                  SM="$(echo -e "${bam_file_name}" | tr -d '[[:space:]]')"
              else
                  # extract sample name from bam file read group info field
                  SM=$(samtools view -H $cur_bam | grep @RG | head -1 | sed "s/.*SM:\\([^	]*\\).*/\\1/" | tr -d '[:space:]')
              fi
              printf "sample$i	$SM\\n" >> names.txt
              i=$((i+1))
          done
      fi
      '''
  }

  // perform regression in R
  process R_regression {

      publishDir params.out_folder+'/PDF/', mode: 'move', pattern: "*[ATCG-]*.pdf"

      tag { region_tag }

      input:
      set val(region_tag), file('TABLE/*'), file('names.txt') from table
      file fasta_ref
      file fasta_ref_fai
      file fasta_ref_gzi
      file pairs_file

      output:
      file "${region_tag}.vcf" into vcf
      file '*.pdf' into PDF

      shell:
      '''
      # create a dummy empty pdf to avoid an error in the process when no variant is found
      touch empty.pdf
      needlestack.r --pairs_file=!{params.pairs_file} --source_path=!{baseDir}/bin/ --out_file=!{region_tag}.vcf --fasta_ref=!{fasta_ref} --GQ_threshold=!{params.min_qval} --min_coverage=!{params.min_dp} --min_reads=!{params.min_ao} --SB_type=!{params.sb_type} --SB_threshold_SNV=!{params.sb_snv} --SB_threshold_indel=!{params.sb_indel} --output_all_SNVs=!{params.all_SNVs} --do_plots=!{!params.no_plots} --plot_labels=!{!params.no_labels} --add_contours=!{!params.no_contours} --extra_rob=!{params.extra_robust_gl}
      '''
  }

  // merge all vcf files in one big file
  process collect_vcf_result {

      publishDir params.out_folder, mode: 'move'

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
      echo '##perl=v'$(perl -e 'print substr($^V, 1)') >> versions.txt
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
