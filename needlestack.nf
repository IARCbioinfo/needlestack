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

// run using for ex.: 
// needlestack.nf --bed my_bed_file.bed --nsplit 20 --fasta_ref reference.fasta --bam_folder BAM/

// requirement:
// - bedtools
// - samtools
// - Rscript (R)
// - bed_large_cut.r (in bin folder)
// - needlestack.r (in bin folder)
// - pileup2baseindel.pl (in bin folder) (+ perl)
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
params.use_file_name = false //put these argument to use the bam file names as sample names and do not to use the sample name filed from the bam files (SM tag)
params.all_SNVs = false //  output all sites, even when no variant is detected
params.no_plots = false  // do not produce pdf plots of regressions 
params.out_folder = params.bam_folder // if not provided, outputs will be held on the input bam folder
params.no_indels = false // do not skip indels

/* If --help in parameters, print software usage */

if (params.help) {
    log.info ''
    log.info '--------------------------------------------------'
    log.info 'NEEDLESTACK: A MULTI-SAMPLE SOMATIC VARIANT CALLER'
    log.info '--------------------------------------------------'
    log.info 'Copyright (C) 2015 Matthieu Foll and Tiffany Delhomme'
    log.info 'This program comes with ABSOLUTELY NO WARRANTY; for details see LICENSE.txt'
    log.info 'This is free software, and you are welcome to redistribute it'
    log.info 'under certain conditions; see LICENSE.txt for details.'
    log.info '--------------------------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info '    nextflow run iarcbioinfo/needlestack -with-docker iarcbioinfo/needlestack --bed bedfile.bed --bam_folder BAM/ --fasta_ref reference.fasta [other options]'
    log.info ''
    log.info 'Mandatory arguments:'
    log.info '    --bed            BED_FILE                 Restrict the calling to regions listed in the BED file'
    log.info '    --bam_folder     BAM_DIR                  BAM files directory.'
    log.info '    --fasta_ref      REF_IN_FASTA             Reference genome in fasta format.'
    log.info 'Options:'
    log.info '    --nsplit         INTEGER                  Split the bed file in nsplit pieces and run in parallel.'
    log.info '    --min_dp         INTEGER                  Minimum coverage in at least one sample to consider a site.'
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
    log.info '    --no_plots                                Do not output PDF regression plots.'
    log.info '    --no_indels                               Do not call indels.'
    log.info '    --out_folder     OUTPUT FOLDER            Output directory, by default input bam folder.'
    log.info ''
    exit 1
}

assert (params.fasta_ref != true) && (params.fasta_ref != null) : "please specify --fasta_ref option (--fasta_ref reference.fasta(.gz))"
assert (params.bed != true) && (params.bed != null) : "please specify --bed option (--bed bedfile.bed)"
assert (params.bam_folder != true) && (params.bam_folder != null) : "please specify --bam_folder option (--bam_folder bamfolder)"

bed = file( params.bed )
fasta_ref = file( params.fasta_ref )
fasta_ref_fai = file( params.fasta_ref+'.fai' )
fasta_ref_gzi = file( params.fasta_ref+'.gzi' )

/* Verify user inputs are correct */

assert params.sb_type in ["SOR","RVSB"] : "--sb_type must be SOR or RVSB "
assert params.all_SNVs in [true,false] : "do not assign a value to --all_SNVs"
assert params.no_plots in [true,false] : "do not assign a value to --no_plots"
assert params.no_indels in [true,false] : "do not assign a value to --no_indels"
assert params.use_file_name in [true,false] : "do not assign a value to --use_file_name"
try { assert fasta_ref.exists() : "\n WARNING : fasta reference not located in execution directory. Make sure reference index is in the same folder as fasta reference" } catch (AssertionError e) { println e.getMessage() }
if (fasta_ref.exists()) {assert fasta_ref_fai.exists() : "input fasta reference does not seem to have a .fai index (use samtools faidx)"}
try { assert bed.exists() : "\n WARNING : bed file not located in execution directory" } catch (AssertionError e) { println e.getMessage() }
try { assert file(params.bam_folder).exists() : "\n WARNING : input BAM folder not located in execution directory" } catch (AssertionError e) { println e.getMessage() }
assert (params.min_dp > 0) : "minimum coverage must be higher than 0 (--min_dp)"
assert (params.max_DP > 1) : "maximum coverage before downsampling must be higher than 1 (--max_DP)"
assert (params.min_ao > 0) : "minimum alternative reads must be higher than 0 (--min_ao)"
assert (params.nsplit > 0) : "number of regions to split must be higher than 0 (--nsplit)"
assert (params.min_qval > 0) : "minimum Phred-scale qvalue must be higher than 0 (--min_qval)"
assert (params.sb_snv > 0 && params.sb_snv < 101) : "strand bias for SNVs must be in [0,100]"
assert (params.sb_indel > 0 && params.sb_indel < 101) : "strand bias for indels must be in [0,100]"
assert (params.map_qual >= 0) : "minimum mapping quality (samtools) must be higher than or equals to 0"
assert (params.base_qual >= 0) : "minimum base quality (samtools) must be higher than or equals to 0"

if(params.use_file_name == true){
  sample_names = "FILE"
} else { sample_names = "BAM" }

/* Software information */

log.info ''
log.info '--------------------------------------------------'
log.info 'NEEDLESTACK: A MULTI-SAMPLE SOMATIC VARIANT CALLER'
log.info '--------------------------------------------------'
log.info 'Copyright (C) 2015 Matthieu Foll and Tiffany Delhomme'
log.info 'This program comes with ABSOLUTELY NO WARRANTY; for details see LICENSE.txt'
log.info 'This is free software, and you are welcome to redistribute it'
log.info 'under certain conditions; see LICENSE.txt for details.'
log.info '--------------------------------------------------'
log.info "Input BAM folder (--bam_folder)                                 : ${params.bam_folder}"
log.info "Reference in fasta format (--fasta_ref)                         : ${params.fasta_ref}"
log.info "Intervals for calling (--bed)                                   : ${params.bed}"
log.info "Number of regions to split (--nsplit)                           : ${params.nsplit}"
log.info "To consider a site for calling:"
log.info "     minimum coverage (--min_dp)                                : ${params.min_dp}"
log.info "     minimum of alternative reads (--min_ao)                    : ${params.min_ao}"
log.info "Phred-scale qvalue threshold (--min_qval)                       : ${params.min_qval}"
log.info "Strand bias measure (--sb_type)                                 : ${params.sb_type}"
log.info "Strand bias threshold for SNVs (--sb_snv)                       : ${params.sb_snv}"
log.info "Strand bias threshold for indels (--sb_indel)                   : ${params.sb_indel}"
log.info "Samtools minimum mapping quality (--map_qual)                   : ${params.map_qual}"
log.info "Samtools minimum base quality (--base_qual)                     : ${params.base_qual}"
log.info "Samtools maximum coverage before downsampling (--max_DP)        : ${params.max_DP}"          
log.info "Sample names definition (--use_file_name)                       : ${sample_names}"
log.info(params.all_SNVs == true ? "Output all SNVs (--all_SNVs)                                    : yes" : "Output all SNVs (--all_SNVs)                                    : no" ) 
log.info(params.no_plots == true ? "PDF regression plots (--no_plots)                               : no"  : "PDF regression plots (--no_plots)                               : yes" )
log.info(params.no_indels == true ? "Skip indels (--no_indels)                                       : yes" : "Skip indels (--no_indels)                                       : no" )  
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
     val sample_names
     
     output:
     set val(region_tag), file('sample*.txt'), file('names.txt') into table
        
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
		# rename the output (the converter call files sample{i}.txt)
		i=1
		for cur_bam in !{bam} 
		do 
			if [ "!{sample_names}" == "FILE" ]; then
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
	needlestack.r --out_file=!{region_tag}.vcf --fasta_ref=!{fasta_ref} --GQ_threshold=!{params.min_qval} --min_coverage=!{params.min_dp} --min_reads=!{params.min_ao} --SB_type=!{params.sb_type} --SB_threshold_SNV=!{params.sb_snv} --SB_threshold_indel=!{params.sb_indel} --output_all_SNVs=!{params.all_SNVs} --do_plots=!{!params.no_plots}
	'''
}
//PDF.flatten().filter { it.size() == 0 }.subscribe { it.delete() }

// merge all vcf files in one big file 
process collect_vcf_result {

	storeDir { params.out_folder }

	input:
	file '*.vcf' from vcf.toList()
     file fasta_ref_fai
  
	output:
	file 'all_variants.vcf' into big_vcf

	shell:
	empty_files = file(params.out_folder+'/VCF/*_empty.pdf')
	for( def file : empty_files ) {
    		file.delete()
	}
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
