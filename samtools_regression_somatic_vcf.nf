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
params.sor_snv = 4 // strand bias SOR threshold for snv 
params.sor_indel = 10 // strand bias SOR threshold for indels
params.map_qual = 20 // min mapping quality (passed to samtools)
params.base_qual = 20 // min base quality (passed to samtools)
params.max_DP = 30000 // downsample coverage per sample (passed to samtools)
params.sample_names = "BAM" // put FILE to use the bam file names as sample names and BAM to use the sample name filed from the bam files
params.all_sites = "FALSE" //  output all sites, even when no variant is detected (but still affected by min_dp and min_ao)
params.do_plots = "TRUE" // produce pdf plots of regressions 

bed = file( params.bed )
fasta_ref = file( params.fasta_ref )
fasta_ref_fai = file( params.fasta_ref+'.fai' )
fasta_ref_gzi = file( params.fasta_ref+'.gzi' )

bam = Channel.fromPath( params.bam_folder+'/*.bam' ).toList()   
bai = Channel.fromPath( params.bam_folder+'/*.bam.bai' ).toList()

/* split bed file into nsplit regions */
process split_bed {
     
     storeDir { params.bam_folder+'/BED_REGIONS/' } 
     
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

     storeDir { params.bam_folder+'/PILEUP/' }   
        
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
     
     storeDir { params.bam_folder+'/PILEUP/'+region_tag }   
        
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
     
     errorStrategy 'ignore'
     
     storeDir { params.bam_folder+'/VCF/' }   
        
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
	pileup_nbrr_caller_vcf.r !{region_tag}.vcf !{fasta_ref} !{params.min_qval} !{params.min_dp} !{params.min_ao} !{params.sor_snv} !{params.sor_indel} !{params.all_sites} !{params.do_plots}
	'''
}

// merge all vcf files in one big file 
process collect_vcf_result {

	storeDir { params.bam_folder }

	input:
	file '*.vcf' from vcf.toList()
  
	output:
	file 'all_variants.vcf' into big_vcf

	'''
	nb_vcf=$(find . -maxdepth 1 -name '*vcf' | wc -l)
	if [ $nb_vcf -gt 1 ]; then
		vcfoverlay *.vcf > all_variants.vcf
	else 
		cp *.vcf all_variants.vcf
	fi
	'''
}
