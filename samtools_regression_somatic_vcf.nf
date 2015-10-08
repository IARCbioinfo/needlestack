#! /usr/bin/env nextflow
// run using for ex.: 
// samtools_regression_somatic_vcf.nf --bed my_bed_file.bed --nsplit 20 --fasta_ref /scratch/appli57_local_duplicates/reference/hg19_torrentserver.fasta --bam_folder BAM/

// requirement:
// - bedtools
// - samtools
// - Rscript (R)
// - bed_large_cut.r
// - pileup_nbrr_caller_vcf.r 
// - pileup2baseindel.pl (+ perl)
// - vcfoverlay from vcflib

params.min_dp = 50 // minimum coverage in at least one sample to consider a site
params.min_ao = 5 // minimum number of non-ref reads in at least one sample to consider a site
params.nsplit = 1 // split the bed file in nsplit pieces and run in parallel 
params.min_qval = 50 // qvalue in Phred scale to consider a variant 
// http://gatkforums.broadinstitute.org/discussion/5533/strandoddsratio-computation filter out SOR > 4 for SNVs and > 10 for indels
// filter out RVSB > 0.85 (maybe less stringent for SNVs)
params.sb_type = "SOR" // strand bias measure to be used: "SOR" or "RVSB"
params.sb_snv = 100 // strand bias threshold for snv 
params.sb_indel = 100 // strand bias threshold for indels
params.map_qual = 20 // min mapping quality (passed to samtools)
params.base_qual = 20 // min base quality (passed to samtools)
params.max_DP = 30000 // downsample coverage per sample (passed to samtools)
params.sample_names = "BAM" // put FILE to use the bam file names as sample names and BAM to use the sample name filed from the bam files
params.all_sites = "FALSE" //  output all sites, even when no variant is detected
params.do_plots = "TRUE" // produce pdf plots of regressions 
params.out_folder = params.bam_folder // if not provided, outputs will be held on the input bam folder

/* If --help in parameters, print software usage */

if (params.help) {
    log.info ''
    log.info '----------------------------'
    log.info 'ROBUST REGRESSION VARIANT CALLER'
    log.info '----------------------------'
    log.info 'somatic variant calling pipeline using multi-sampling from (ultra)deep next-generation sequencing'
    log.info '----------------------------'
    log.info ''
    log.info 'Usage: '
    log.info '    nextflow run mfoll/robust-regression-caller -with-docker mfoll/robust-regression-caller --bed your_bedfile.bed --bam_folder BAM/ --fasta_ref hg19.fasta [other options]'
    log.info ''
    log.info 'Mandatory arguments:'
    log.info '    --bed            BED_FILE                 Query intervals in bed format.'
    log.info '    --bam_folder     BAM_DIR                  Query bam files directory.'
    log.info '    --fasta_ref      REF_IN_FASTA             Reference genome in fasta.'
    log.info 'Options:'
    log.info '    --nsplit         INTEGER                  Number of splitted regions for parallel computations'
    log.info '    --min_dp         INTEGER                  Minimum coverage for a site to be considered as alternative.'
    log.info '    --min_ao         INTEGER                  Minimum of alternative reads to be considered as alternative.'
    log.info '    --min_qval       VALUE                    Phred-scale qvalue threshold to consider a sample as variant.'
    log.info '    --sb_type        SOR or RVSB              Strand bias measure.'
    log.info '    --sb_snv         VALUE                    Strand bias threshold for SNV.'
    log.info '    --sb_indel       VALUE                    Strand bias threshold for indel.'
    log.info '    --map_qual       VALUE                    Samtools minimum mapping quality.'
    log.info '    --base_qual      VALUE                    Samtools minimum base quality.'
    log.info '    --max_DP         INTEGER                  Samtools maximum coverage before downsampling.'
    log.info '    --sample_names   FILE or BAM              Sample names definition.'
    log.info '    --all_sites      BOOLEAN                  Output all sites, even when no variant found.'
    log.info '    --do_plots       BOOLEAN                  PDF regression plots in the output.'
    log.info '    --out_folder     OUTPUT FOLDER            Output directory, bu default input bam folder.'
    log.info ''
    exit 1
}

bed = file( params.bed )
fasta_ref = file( params.fasta_ref )
fasta_ref_fai = file( params.fasta_ref+'.fai' )
fasta_ref_gzi = file( params.fasta_ref+'.gzi' )


/* Verify user inputs are correct */

assert params.sb_type in ["SOR","RVSB"]
assert params.all_sites in ["TRUE","FALSE"]
assert params.do_plots in ["TRUE","FALSE"]
assert params.sample_names in ["FILE","BAM"]

/* Software information */

log.info "============================================"
log.info "ROBUST REGRESSION VARIANT CALLER"
log.info "somatic variant calling pipeline using multi-sampling from (ultra)deep next-generation sequencing"
log.info "============================================"
log.info "query bam folder                                                : ${params.bam_folder}"
log.info "reference in fasta format                                       : ${params.fasta_ref}"
log.info "intervals for calling                                           : ${params.bed}"
log.info "number of splitted regions (--nsplit)                           : ${params.nsplit}"
log.info "to consider a site as alternative                               : "
log.info "	minimum coverage (--min_dp)                             : ${params.min_dp}"
log.info "	minimum of alternative reads (--min_ao)                 : ${params.min_ao}"
log.info "to consider a sample as variant                                 : "
log.info "	Phred-scale qvalue threshold (--min_qval)               : ${params.min_qval}"
log.info "strand bias measure (--sb_type)                                 : ${params.sb_type}"
log.info "strand bias threshold for SNV (--sb_snv)                        : ${params.sb_snv}"
log.info "strand bias threshold for indel (--sb_indel)                    : ${params.sb_indel}"
log.info "samtools minimum mapping quality (--map_qual)                   : ${params.map_qual}"
log.info "samtools minimum base quality (--base_qual)                     : ${params.base_qual}"
log.info "samtools maximum coverage before downsampling (--max_DP)        : ${params.max_DP}"          
log.info "sample names definition (--sample_names)                        : ${params.sample_names}"
log.info "output all sites (--all_sites)                                  : ${params.all_sites}"
log.info "pdf regression plots (--do_plots)                               : ${params.do_plots}"
log.info "output folder (--out_folder)                                    : ${params.out_folder}"
log.info "\n"


bam = Channel.fromPath( params.bam_folder+'/*.bam' ).toList()   
bai = Channel.fromPath( params.bam_folder+'/*.bam.bai' ).toList()

/* split bed file into nsplit regions */
process split_bed {
     
     storeDir { params.out_folder+'/BED_REGIONS/' } 
     
     intput:
     file bed
        
	output:
	file '*_regions' into split_bed mode flatten 

	shell:
	'''
	grep -v '^track' !{bed} | sort -k1,1 -k2,2n | bedtools merge -i stdin | awk '{print $1" "$2" "$3}' | bed_large_cut.r !{params.nsplit}
	'''
}

// create mpileup file + sed to have "*" when there is no coverage (otherwise pileup2baseindel.pl is unhappy)
process samtools_mpileup {

     storeDir { params.out_folder+'/PILEUP/' }   
        
     tag { region_tag }   
        
     input:
     file split_bed
	file bam
	file bai  
	file fasta_ref
	file fasta_ref_fai
	file fasta_ref_gzi
     
     output:
     set val(region_tag), file("${region_tag}.pileup") into pileup
        
 	shell:
 	region_tag = split_bed.baseName
	'''
	while read bed_line; do
		samtools mpileup --fasta-ref !{fasta_ref} --region $bed_line --ignore-RG --min-BQ !{params.base_qual} --min-MQ !{params.map_qual} --max-idepth 1000000 --max-depth !{params.max_DP} !{bam} | sed 's/		/	*	*/g' >> !{region_tag}.pileup
	done < !{split_bed}
	'''
}

// split mpileup file and convert to table 
process mpileup2table {
     
     errorStrategy 'ignore'
     
     storeDir { params.out_folder+'/PILEUP/'+region_tag }   
        
     tag { region_tag }   
        
     input:
     set val(region_tag), file("${region_tag}.pileup") from pileup
     file bam
     
     output:
     set val(region_tag), file('sample*.txt'), file('names.txt') into table
        
 	shell:
 	'''
 	nb_pos=$(wc -l < !{region_tag}.pileup)
	if [ $nb_pos -gt 0 ]; then 
		# split and convert pileup file
		pileup2baseindel.pl -i !{region_tag}.pileup
		# rename the output (the converter call files sample{i}.txt)
		i=1
		for cur_bam in !{bam} 
		do 
			if [ "!{params.sample_names}" == "FILE" ]; then
				# use bam file name as sample name
				bam_file_name="${cur_bam%.*}"
				# remove whitespaces from name
				SM="$(echo -e "${bam_file_name}" | tr -d '[[:space:]]')"
			else
				# extract sample name from bam file read group info field
				SM=$(samtools view -H $cur_bam | grep @RG | head -1 | sed "s/.*SM:\\([^\\t]*\\).*/\\1/" | tr -d '[:space:]')
			fi
			printf "sample$i\\t$SM\\n" >> names.txt
			i=$((i+1)) 
		done
	fi
	'''
}

// perform regression in R
process R_regression {
       
     storeDir { params.out_folder+'/VCF/' }   
        
     tag { region_tag }   
        
     input:
     set val(region_tag), file(table_file), file('names.txt') from table
     file fasta_ref
     file fasta_ref_fai
	file fasta_ref_gzi
     
     output:
     file "${region_tag}.vcf" into vcf
     file '*.pdf' into PDF
        
 	shell:
 	'''
 	# create a dummy empty pdf to avoid an error in the process when no variant is found 
 	touch !{region_tag}_empty.pdf
	pileup_nbrr_caller_vcf.r --out_file=!{region_tag}.vcf --fasta_ref=!{fasta_ref} --GQ_threshold=!{params.min_qval} --min_coverage=!{params.min_dp} --min_reads=!{params.min_ao} --SB_type=!{params.sb_type} --SB_threshold_SNV=!{params.sb_snv} --SB_threshold_indel=!{params.sb_indel} --output_all_sites=!{params.all_sites} --do_plots=!{params.do_plots}
	'''
}
PDF.flatten().filter { it.size() == 0 }.subscribe { it.delete() }

// merge all vcf files in one big file 
process collect_vcf_result {

	storeDir { params.out_folder }

	input:
	file '*.vcf' from vcf.toList()
     file fasta_ref_fai
  
	output:
	file 'all_variants.vcf' into big_vcf

	shell:
	'''
	nb_vcf=$(find . -maxdepth 1 -name '*vcf' | wc -l)
	if [ $nb_vcf -gt 1 ]; then
		vcfoverlay *.vcf > all_variants.vcf
	else 
		cp .vcf all_variants.vcf
	fi
	# Add contigs in the VCF header
	cat !{fasta_ref_fai} | cut -f1,2 | sed -e 's/^/##contig=<ID=/' -e 's/[	 ][	 ]*/,length=/' -e 's/$/>/' > contigs.txt
	sed -i '/##reference=.*/ r contigs.txt' all_variants.vcf
	'''
}
