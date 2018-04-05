#! /usr/bin/env Rscript

# needlestack: a multi-sample somatic variant caller
# Copyright (C) 2017 IARC/WHO
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


print(Sys.time())
############################## ARGUMENTS SECTION ##############################
args <- commandArgs(TRUE)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsL <- as.list(as.character(as.data.frame(do.call("rbind", parseArgs(args)))$V2))
names(argsL) <- as.data.frame(do.call("rbind", parseArgs(args)))$V1
args <- argsL;rm(argsL)

if("--help" %in% args | is.null(args$out_file) | is.null(args$fasta_ref) ) {
  cat("
      The R Script arguments_section.R
      Mandatory arguments:
      --writeHeader                  - write header if true
      --out_file=file_name           - name of output vcf
      --fasta_ref=path               - path of fasta ref
      --help                         - print this text
      Optionnal arguments:
      --source_path=path             - path to source files (glm_rob_nb.r, plot_rob_nb.r, plot_alignments.r)
      --samtools=path                - path of samtools, default=samtools
      --SB_type=SOR, RVSB or FS      - strand bias measure, default=SOR
      --SB_threshold_SNV=value       - strand bias threshold for SNV, default=100
      --SB_threshold_indel=value     - strand bias threshold for indel, default=100
      --min_coverage=value           - minimum coverage in at least one sample to consider a site, default=50
      --min_reads=value              - minimum number of non-ref reads in at least one sample to consider a site, default=5
      --min_af=value                 - minimum allelic fraction in at least one sample to consider a site, default=0
      --GQ_threshold=value           - phred scale qvalue threshold for variants, default=50
      --output_all_SNVs=boolean      - output all SNVs, even when no variant is detected, default=FALSE
      --do_plots=boolean             - output regression plots, default=TRUE
      --do_alignments=boolean        - output alignment plots, default=FALSE
      --ref_genome=string            - reference genome for alignments plot, examples : Hsapiens.UCSC.hg19, Hsapiens.UCSC.hg18, Mmusculus.UCSC.mm10...
      --extra_rob=boolean            - perform an extra-robust regression, default=FALSE
      --min_af_extra_rob=value       - minimum allelic fraction to exclude a sample at a position for extra-robust regression, default=0.2
      --min_prop_extra_rob=value     - minimum proportion of samples having an allelic fraction to be excluded from extra-robust regression, default=0.1
      --max_prop_extra_rob=value     - maximum proportion of samples having an allelic fraction to be excluded from extra-robust regression, default=0.1
      --pairs_file=file_name         - name of file containing the list of matched Tumor/Normal pairs for somatic variant calling
      --afmin_power=value            - minimum allelic fraction in mean for somatic mutations, default=0.01
      --sigma=value                  - sigma parameter for negative binomial modeling germline mutations, default=0.1
      WARNING : by default samtools has to be in your path
      Example:
      needlestack.r --out_file=test.vcf --fasta_ref=~/Documents/References/ \n\n")

  q(save="no")
}
if(is.null(args$samtools)) {args$samtools="samtools"}
if(is.null(args$SB_type)) {args$SB_type="SOR"}
if(is.null(args$SB_threshold_SNV)) {args$SB_threshold_SNV=100} else {args$SB_threshold_SNV=as.numeric(args$SB_threshold_SNV)}
if(is.null(args$SB_threshold_indel)) {args$SB_threshold_indel=100} else {args$SB_threshold_indel=as.numeric(args$SB_threshold_indel)}
if(is.null(args$min_coverage)) {args$min_coverage=30} else {args$min_coverage=as.numeric(args$min_coverage)}
if(is.null(args$min_reads)) {args$min_reads=3} else {args$min_reads=as.numeric(args$min_reads)}
if(is.null(args$min_af)) {args$min_af=0} else {args$min_af=as.numeric(args$min_af)}
if(is.null(args$GQ_threshold)) {args$GQ_threshold=50} else {args$GQ_threshold=as.numeric(args$GQ_threshold)}
if(is.null(args$output_all_SNVs)) {args$output_all_SNVs=FALSE} else {args$output_all_SNVs=as.logical(args$output_all_SNVs)}
if(is.null(args$pairs_file)) {args$pairs_file=FALSE}

if(is.null(args$do_plots)) {
  if(args$pairs_file!=FALSE){
    args$do_plots="SOMATIC"
  }else{
    args$do_plots="ALL"
  }
} else {
  args$do_plots=args$do_plots
}
if(is.null(args$do_alignments)) {args$do_alignments=FALSE} else {args$do_alignments=as.logical(args$do_alignments)}
if(is.null(args$plot_labels)) {args$plot_labels=FALSE} else {args$plot_labels=as.logical(args$plot_labels)}
if(is.null(args$add_contours)) {args$add_contours=FALSE} else {args$add_contours=as.logical(args$add_contours)}
if(is.null(args$extra_rob)) {args$extra_rob=FALSE} else {args$extra_rob=as.logical(args$extra_rob)}
if(is.null(args$min_af_extra_rob)) {args$min_af_extra_rob=0.2} else {args$min_af_extra_rob=as.numeric(args$min_af_extra_rob)}
if(is.null(args$min_prop_extra_rob)) {args$min_prop_extra_rob=0.1} else {args$min_prop_extra_rob=as.numeric(args$min_prop_extra_rob)}
if(is.null(args$max_prop_extra_rob)) {args$max_prop_extra_rob=0.5} else {args$max_prop_extra_rob=as.numeric(args$max_prop_extra_rob)}
if(is.null(args$afmin_power)) {args$afmin_power=-1} else {args$afmin_power=as.numeric(args$afmin_power)}
if(is.null(args$sigma)) {args$sigma=0.1} else {args$sigma=as.numeric(args$sigma)}

samtools=args$samtools
out_file=args$out_file
fasta_ref=args$fasta_ref
bam_folder=args$bam_folder
GQ_threshold=args$GQ_threshold
min_coverage=args$min_coverage
min_reads=args$min_reads
min_af=args$min_af
print(min_af)
# http://gatkforums.broadinstitute.org/discussion/5533/strandoddsratio-computation filter out SOR > 4 for SNVs and > 10 for indels
# filter out RVSB > 0.85 (maybe less stringent for SNVs)
SB_type=args$SB_type
SB_threshold_SNV=args$SB_threshold_SNV
SB_threshold_indel=args$SB_threshold_indel
output_all_SNVs=args$output_all_SNVs
do_plots=args$do_plots
do_alignments=args$do_alignments

ref_genome=args$ref_genome

plot_labels=args$plot_labels
add_contours=args$add_contours
extra_rob=args$extra_rob
min_af_extra_rob=args$min_af_extra_rob
min_prop_extra_rob=args$min_prop_extra_rob
max_prop_extra_rob=args$max_prop_extra_rob
pairs_file=args$pairs_file
sigma=args$sigma
afmin_power=args$afmin_power


if( do_plots==FALSE & do_alignments==TRUE){
  cat("Error : do_alignments can not be TRUE since do_plots=FALSE")
  q(save="no")
}

#source R files
getScriptPath <- function(){
  cmd.args <- commandArgs()
  m <- regexpr("(?<=^--file=).+", cmd.args, perl=TRUE)
  script.dir <- dirname(regmatches(cmd.args, m))
  if(length(script.dir) == 0) stop("can't determine script dir: please call the script with Rscript")
  if(length(script.dir) > 1) stop("can't determine script dir: more than one '--file' argument detected")
  return(script.dir)
}

if(is.null(args$source_path)) {
  source_path=paste0(getScriptPath(),"/")
} else {
  source_path=args$source_path
}

source(paste(source_path,'glm_rob_nb.r',sep=""))
source(paste(source_path,"plot_rob_nb.r",sep=""))
source(paste(source_path,"plot_alignments.r",sep=""))


############################################################

options("scipen"=100)

indiv_run=read.table("names.txt",stringsAsFactors=F,colClasses = "character")
#check if samples names exist in names.txt file
if(ncol(indiv_run)!=2){cat("Error : names.txt contains only one column : samples names can't be obtained from the bam files : use --use_file_name option"); q(save="no")}
indiv_run[,2]=make.unique(indiv_run[,2],sep="_")
nindiv=nrow(indiv_run)
id_samples=seq(1,nindiv) #samples IDs (integers)

isTNpairs = FALSE #activates Tumor-Normal pairs mode
if(pairs_file != FALSE) { #if user gives a pairs_file to needlestack
  isTNpairs = file.exists(pairs_file) #checks existence of tumour-normal pairs file => will be redundant once this is checked in the workflow
  if( isTNpairs ){
    if(afmin_power==-1) afmin_power = 0.01 #put default value
    pairsname = scan(pairs_file,nmax = 2,what = "character")
    TNpairs=read.table(pairs_file,h=T)
    names(TNpairs)[grep("TU",pairsname,ignore.case =T)] = "TUMOR" #set columns names (to avoid problems due to spelling variations or typos)
    names(TNpairs)[grep("NO",pairsname,ignore.case =T)] = "NORMAL"

    #Check that all sample names in the pairs_file file actually exist
    names = c(as.vector(subset(TNpairs$NORMAL,TNpairs$NORMAL!="NA")),as.vector(subset(TNpairs$TUMOR,TNpairs$TUMOR!="NA")))
    if(length(which((names %in% indiv_run[,2])==FALSE))!=0){cat("Error : sample name not present in names.txt"); q(save = "no")}

    onlyNindex = which( sapply(indiv_run[,2], function(x) x%in%TNpairs$NORMAL[is.na(TNpairs$TUMOR)] ) )
    onlyTindex = which( sapply(indiv_run[,2], function(x) x%in%TNpairs$TUMOR[is.na(TNpairs$NORMAL)] ) )
    TNpairs.complete = TNpairs[!(is.na(TNpairs$TUMOR)|is.na(TNpairs$NORMAL) ),] # all complete T-N pairs
    Tindex = sapply( 1:nrow(TNpairs.complete) , function(k) return(which( indiv_run[,2]==TNpairs.complete$TUMOR[k])) )
    Nindex = sapply( 1:nrow(TNpairs.complete) , function(k) return(which( indiv_run[,2]==TNpairs.complete$NORMAL[k])) )
  }else{
    cat("Error : pairs_file not located in execution directory ")
    q(save="no")
  }
}


depth =(1:nindiv)*11 - 10 +3
A_cols=(1:nindiv)*11 - 9 + 3
T_cols=(1:nindiv)*11 - 8 + 3
C_cols=(1:nindiv)*11 - 7 + 3
G_cols=(1:nindiv)*11 - 6 + 3
a_cols=(1:nindiv)*11 - 5 + 3
t_cols=(1:nindiv)*11 - 4 + 3
c_cols=(1:nindiv)*11 - 3 + 3
g_cols=(1:nindiv)*11 - 2 + 3
insertion_cols=(1:nindiv)*11 - 1 + 3
deletion_cols=(1:nindiv)*11 - 0 + 3


############################## FUNCTION SECTION ##############################
non_ref_bases=function(ref) {
  setdiff(c("A","T","C","G"),ref)
}

SOR=function(RO_forward,AO_forward,RO_reverse,AO_reverse) {
  X00=RO_forward
  X10=AO_forward
  X01=RO_reverse
  X11=AO_reverse
  refRatio=pmax(X00,X01)/pmin(X00,X01)
  altRatio=pmax(X10,X11)/pmin(X10,X11)
  R=(X00/X01)*(X11/X10)
  log((R+1/R)/(refRatio/altRatio))
}

common_annot=function() {
  Rp<<-as.numeric(linepos[eval(as.name(paste(linepos[3],"_cols",sep="")))])
  Rm<<-as.numeric(linepos[eval(as.name(paste(tolower(linepos[3]),"_cols",sep="")))])
  Cp<<-Vp+Rp
  Cm<<-Vm+Rm
  tmp_rvsbs=pmax(Vp * Cm, Vm * Cp) / (Vp * Cm + Vm * Cp)
  tmp_rvsbs[which(is.na(tmp_rvsbs))]=-1
  rvsbs<<-tmp_rvsbs
  all_rvsb<<-max(as.numeric(sum(Vp)) * sum(Cm), as.numeric(sum(Vm)) * sum(Cp)) / (as.numeric(sum(Vp)) * sum(Cm) + as.numeric(sum(Vm)) * sum(Cp))
  if (is.na(all_rvsb)) all_rvsb<<- (-1)
  tmp_sors<<-SOR(Rp,Vp,Rm,Vm)
  tmp_sors[which(is.na(tmp_sors))]=-1
  tmp_sors[which(is.infinite(tmp_sors))]=99
  sors<<-tmp_sors
  all_sor<<-SOR(sum(Rp),sum(Vp),sum(Rm),sum(Vm))
  if (is.na(all_sor)) all_sor<<- (-1)
  if (is.infinite(all_sor)) all_sor<<- 99
  FisherStrand<<- -10*log10(unlist(lapply(1:nindiv, function(indiv) fisher.test(matrix(c(Vp[indiv],Vm[indiv],Rp[indiv],Rm[indiv]), nrow=2))$p.value ))) #col1=alt(V), col2=ref(R)
  FisherStrand[which(FisherStrand>1000)] = 1000
  FisherStrand_all<<--10*log10(fisher.test(matrix(c(sum(Vp),sum(Vm),sum(Rp),sum(Rm)), nrow=2))$p.value)
  if(FisherStrand_all>1000) FisherStrand_all=1000
  all_RO<<-sum(Rp+Rm)
}

## functions to compute the Qvalue of a function with a given AF
toQvalueN <- function(x,rob_nb_res,sigma){ #change sigma value to change the departure from binomial distribution with parameter 0.5
  y = qnbinom(0.01,size = 1/sigma,mu = 0.5*x)
  unlist(-10*log10(p.adjust((dnbinom(c(rob_nb_res$ma_count,y),size=1/rob_nb_res$coef[[1]],mu=rob_nb_res$coef[[2]]*c(rob_nb_res$coverage,x)) +
                               pnbinom(c(rob_nb_res$ma_count,y),size=1/rob_nb_res$coef[[1]],mu=rob_nb_res$coef[[2]]*c(rob_nb_res$coverage,x),lower.tail = F)))[length(rob_nb_res$coverage)+1]))
}

toQvalueT <- function(x,rob_nb_res,afmin_power){ #change afmin_power to
  y = qbinom(0.01,x,afmin_power,lower.tail = T)
  unlist(-10*log10(p.adjust((dnbinom(c(rob_nb_res$ma_count,y),size=1/rob_nb_res$coef[[1]],mu=rob_nb_res$coef[[2]]*c(rob_nb_res$coverage,x)) +
                               pnbinom(c(rob_nb_res$ma_count,y),size=1/rob_nb_res$coef[[1]],mu=rob_nb_res$coef[[2]]*c(rob_nb_res$coverage,x),lower.tail = F)))[length(rob_nb_res$coverage)+1]))
}

###############################################################################################


############################## WRITEOUT SECTION ##############################
write_out=function(...) {
  cat(paste(...,sep=""),"\n",file=out_file,sep="",append=T)
}


if(file.exists(out_file)) file.remove(out_file)

write_out("##fileformat=VCFv4.1")
write_out("##fileDate=",format(Sys.Date(), "%Y%m%d"))
write_out("##source=needlestack v1.0b")
write_out("##reference=",fasta_ref)
write_out("##phasing=none")
write_out("##filter=\"QVAL > ",GQ_threshold," & ",SB_type,"_SNV < ",SB_threshold_SNV," & ",SB_type,"_INDEL < ",SB_threshold_indel," & min(AO) >= ",min_reads, " & min(AF) >= ",min_af," & min(DP) >= ",min_coverage,"\"")

write_out("##INFO=<ID=TYPE,Number=1,Type=String,Description=\"The type of allele, either snp, ins or del\">")
write_out("##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of samples with data\">")
write_out("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total read depth\">")
write_out("##INFO=<ID=RO,Number=1,Type=Integer,Description=\"Total reference allele observation count\">")
write_out("##INFO=<ID=AO,Number=1,Type=Integer,Description=\"Total alternate allele observation count\">")
write_out("##INFO=<ID=AF,Number=1,Type=Float,Description=\"Estimated allele frequency in the range [0,1]\">")
write_out("##INFO=<ID=SRF,Number=1,Type=Integer,Description=\"Total number of reference observations on the forward strand\">")
write_out("##INFO=<ID=SRR,Number=1,Type=Integer,Description=\"Total number of reference observations on the reverse strand\">")
write_out("##INFO=<ID=SAF,Number=1,Type=Integer,Description=\"Total number of alternate observations on the forward strand\">")
write_out("##INFO=<ID=SAR,Number=1,Type=Integer,Description=\"Total number of alternate observations on the reverse strand\">")
write_out("##INFO=<ID=SOR,Number=1,Type=Float,Description=\"Total Symmetric Odds Ratio of 2x2 contingency table to detect strand bias\">")
write_out("##INFO=<ID=RVSB,Number=1,Type=Float,Description=\"Total Relative Variant Strand Bias\">")
write_out("##INFO=<ID=FS,Number=1,Type=Float,Description=\"Total Fisher Exact Test p-value for detecting strand bias (Phred-scale)\">")
write_out("##INFO=<ID=ERR,Number=1,Type=Float,Description=\"Estimated error rate for the alternate allele\">")
write_out("##INFO=<ID=SIG,Number=1,Type=Float,Description=\"Estimated overdispersion for the alternate allele\">")
write_out("##INFO=<ID=CONT,Number=1,Type=String,Description=\"Context of the reference sequence\">")
write_out("##INFO=<ID=WARN,Number=1,Type=String,Description=\"Warning message when position is processed specifically\">")

write_out("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")
write_out("##FORMAT=<ID=QVAL,Number=1,Type=Float,Description=\"Phred-scaled qvalue for not being an error\">")
write_out("##FORMAT=<ID=QVAL_INV,Number=1,Type=Float,Description=\"Phred-scaled qvalue for not being an error in position where reference is switched\">")
write_out("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">")
write_out("##FORMAT=<ID=RO,Number=1,Type=Integer,Description=\"Reference allele observation count\">")
write_out("##FORMAT=<ID=AO,Number=1,Type=Integer,Description=\"Alternate allele observation count\">")
write_out("##FORMAT=<ID=AF,Number=1,Type=Float,Description=\"Allele fraction of the alternate allele with regard to reference\">")
write_out("##FORMAT=<ID=SB,Number=4,Type=Integer,Description=\"Per-sample component statistics to detect strand bias as SRF,SRR,SAF,SAR\">")
write_out("##FORMAT=<ID=SOR,Number=1,Type=Float,Description=\"Symmetric Odds Ratio of 2x2 contingency table to detect strand bias\">")
write_out("##FORMAT=<ID=RVSB,Number=1,Type=Float,Description=\"Relative Variant Strand Bias\">")
write_out("##FORMAT=<ID=FS,Number=1,Type=Float,Description=\"Fisher Exact Test p-value for detecting strand bias (Phred-scale)\">")
write_out("##FORMAT=<ID=FS,Number=1,Type=Float,Description=\"Fisher Exact Test p-value for detecting strand bias (Phred-scale)\">")
write_out("##FORMAT=<ID=QVAL_minAF,Number=1,Type=Float,Description=\"Phred-scaled qvalue for an allelic fraction of minAF\">")
write_out("##FORMAT=<ID=STATUS,Number=1,Type=String,Description=\"Somatic status (SOMATIC, GERMLINE_CONFIRMED, GERMLINE_UNCONFIRMED, GERMLINE_UNCONFIRMABLE, , or UNKNOWN) of variants\">")

write_out("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t",paste(indiv_run[,2],collapse = "\t"))


###############################################################################################

############################## REGRESSION SECTION ##############################

#Reading input from STDIN
f <- file("stdin")
open(f)
i=1
while(length(line <- readLines(f,n=1, warn = FALSE)) > 0) {
  if(i!=1){
    linepos=unlist(strsplit(line,"\t"))
    if (is.element(linepos[3],c("A","T","C","G"))) {
      # SNV
      for (alt in non_ref_bases(linepos[3])) { #for each base different from the reference
        #base on the reference
        ref=linepos[3]

        #count the number of alternative bases (nonreference)
        Vp=as.numeric(linepos[eval(as.name(paste(alt,"_cols",sep="")))])
        Vm=as.numeric(linepos[eval(as.name(paste(tolower(alt),"_cols",sep="")))])
        ma_count=Vp+Vm
        DP=as.numeric(linepos[eval(depth)])
        if( sum( (ma_count/DP) > 0.8 , na.rm = T) > 0.5*length(ma_count) ){  #here we need to reverse alt and ref for the regression (use ma_count of the ref)
          alt_inv=linepos[3]
          Vp_inv=as.numeric(linepos[eval(as.name(paste(alt_inv,"_cols",sep="")))])
          Vm_inv=as.numeric(linepos[eval(as.name(paste(tolower(alt_inv),"_cols",sep="")))])
          ma_count_inv=Vp_inv+Vm_inv
          ref_inv=TRUE
          reg_res=glmrob.nb(x=DP,y=ma_count_inv,min_coverage=min_coverage,min_reads=min_reads,min_af=min_af,extra_rob=extra_rob,min_af_extra_rob=min_af_extra_rob,min_prop_extra_rob=min_prop_extra_rob,max_prop_extra_rob=max_prop_extra_rob)
        } else {ref_inv=FALSE; reg_res=glmrob.nb(x=DP,y=ma_count,min_coverage=min_coverage,min_reads=min_reads,min_af=min_af,extra_rob=extra_rob,min_af_extra_rob=min_af_extra_rob,min_prop_extra_rob=min_prop_extra_rob,max_prop_extra_rob=max_prop_extra_rob) }

        # compute Qval for minAF
        qval_minAF= rep(0,nindiv)
        somatic_status = rep(".",nindiv)
        if(isTNpairs){#pairs file supplied
          qval_minAF[Nindex] = sapply(1:length(Nindex),function(ii) toQvalueN(DP[Nindex][ii],reg_res,sigma) )
          if(sum(reg_res$GQ[Nindex]>GQ_threshold)>0) qval_minAF[Tindex] = sapply(1:length(Tindex),function(ii) toQvalueN(DP[Tindex][ii],reg_res,sigma) ) #GERMLINE VARIANT, check if power to CONFIRM in tumors
          else qval_minAF[Tindex] = sapply(1:length(Tindex),function(ii) toQvalueT(DP[Tindex][ii],reg_res,afmin_power) ) #no germline variant, check if power to call SOMATIC in tumors
          if( length(onlyNindex)>0 ) qval_minAF[onlyNindex] = sapply(1:length(onlyNindex),function(ii) toQvalueN(DP[onlyNindex][ii],reg_res,sigma) )
          if( length(onlyTindex)>0 ) qval_minAF[onlyTindex] = sapply(1:length(onlyTindex),function(ii) toQvalueT(DP[onlyTindex][ii],reg_res,afmin_power) )
          #no matching normal -> UNKNOWN (impossible to call somatic status)
          somatic_status[onlyTindex][(reg_res$GQ[onlyTindex] > GQ_threshold) ] = "UNKNOWN"
          #no matching tumor -> GERMLINE_UNCONFIRMABLE (impossible to call somatic status)
          somatic_status[onlyNindex][(reg_res$GQ[onlyNindex] > GQ_threshold) ] = "GERMLINE_UNCONFIRMABLE"
          #no tumor variant, no normal variant despite good power in both -> ".", confirmed; already set by default
          #no tumor variant, no normal variant but low power in N and T-> ".", unknown; already set by default
          #no tumor variant, no normal variant but low power in N or T-> ".", unconfirmed ; already set by default
          #no tumor variant, normal variant with good power in T (with or without good power in N)-> "GERMLINE_UNCONFIRMED"
          somatic_status[Tindex][(reg_res$GQ[Tindex] < GQ_threshold)&(qval_minAF[Tindex]>GQ_threshold)&(reg_res$GQ[Nindex]>GQ_threshold) ] = "GERMLINE_UNCONFIRMED"
          somatic_status[Nindex][(reg_res$GQ[Tindex] < GQ_threshold)&(qval_minAF[Tindex]>GQ_threshold)&(reg_res$GQ[Nindex]>GQ_threshold) ] = "GERMLINE_UNCONFIRMED"
          #no tumor variant, normal variant with low power in T (with or without good power in N)-> "GERMLINE_UNCONFIRMABLE"
          somatic_status[Tindex][(reg_res$GQ[Tindex] < GQ_threshold)&(qval_minAF[Tindex]<GQ_threshold)&(reg_res$GQ[Nindex]>GQ_threshold) ] = "GERMLINE_UNCONFIRMABLE"
          somatic_status[Nindex][(reg_res$GQ[Tindex] < GQ_threshold)&(qval_minAF[Tindex]<GQ_threshold)&(reg_res$GQ[Nindex]>GQ_threshold) ] = "GERMLINE_UNCONFIRMABLE"

          #tumor variant, no normal variant but without enough power -> UNKNOWN for both Tumor and Normal
          somatic_status[Tindex][(reg_res$GQ[Tindex] > GQ_threshold)&(qval_minAF[Nindex]<GQ_threshold)&(reg_res$GQ[Nindex]<GQ_threshold) ] = "UNKNOWN"
          somatic_status[Nindex][(reg_res$GQ[Tindex] > GQ_threshold)&(qval_minAF[Nindex]<GQ_threshold)&(reg_res$GQ[Nindex]<GQ_threshold) ] = "UNKNOWN"
          #tumor variant, normal variant (with or without good power in T and N) -> GERMLINE_CONFIRMED
          somatic_status[Tindex][ (reg_res$GQ[Tindex] > GQ_threshold)&(reg_res$GQ[Nindex]>GQ_threshold) ] = "GERMLINE_CONFIRMED"
          somatic_status[Nindex][ (reg_res$GQ[Tindex] > GQ_threshold)&(reg_res$GQ[Nindex]>GQ_threshold) ] = "GERMLINE_CONFIRMED"
          #tumor variant, no normal variant despite good power -> SOMATIC
          somatic_status[Tindex][(reg_res$GQ[Tindex] > GQ_threshold)&(qval_minAF[Nindex]>GQ_threshold)&(reg_res$GQ[Nindex]<GQ_threshold) ] = "SOMATIC"

          #flag possible contamination
          wh.germ = grep("GERMLINE|UNKNOWN",somatic_status)
          if( length(wh.germ)>0 ) somatic_status[Tindex][somatic_status[Tindex] == "SOMATIC"] = paste("POSSIBLE_CONTAMINATION_FROM", paste(indiv_run[wh.germ,2],sep="_",collapse="_"),sep="_")
        }else{# no pairs file supplied
          if(afmin_power==-1 ){#no minimum frequency supplied -> use a negative binomial distribution to check the power
            qval_minAF = sapply(1:nindiv,function(ii) toQvalueN(DP[ii],reg_res,sigma) )
          }else{#minimum frequency supplied -> use a binomial distribution centered around afmin_power to check the power
            qval_minAF = sapply(1:nindiv,function(ii) toQvalueT(DP[ii],reg_res,afmin_power) )
          }
        }
        if (output_all_SNVs | (!is.na(reg_res$coef["slope"]) & sum(reg_res$GQ>=GQ_threshold,na.rm=TRUE)>0)) {
          all_AO=sum(ma_count)
          all_DP=sum(as.numeric(linepos[eval(depth)]))
          common_annot()
          all_RO=sum(Rp+Rm)
          if (SB_type=="SOR") sbs=sors else { if (SB_type=="RVSB") sbs=rvsbs else {sbs=FisherStrand} }
          if (output_all_SNVs | (sum(sbs<=SB_threshold_SNV & reg_res$GQ>=GQ_threshold,na.rm=TRUE)>0)) {
            #to get context around the base
            con=pipe(paste(samtools," faidx ",fasta_ref," ",linepos[1],":",as.numeric(linepos[2])-3,"-",as.numeric(linepos[2])-1," | tail -n1",sep=""))
            before=readLines(con)
            close(con)
            con=pipe(paste(samtools," faidx ",fasta_ref," ",linepos[1],":",as.numeric(linepos[2])+1,"-",as.numeric(linepos[2])+3," | tail -n1",sep=""))
            after=readLines(con)
            close(con)

            cat(linepos[1],"\t",linepos[2],"\t",".","\t",linepos[3],"\t",alt,"\t",max(reg_res$GQ),"\t",".",sep = "",file=out_file,append=T)
            # INFO field
            cat("\t","TYPE=snv;NS=",sum(as.numeric(linepos[eval(depth)])>0),";AF=",sum(reg_res$GQ>=GQ_threshold)/sum(as.numeric(linepos[eval(depth)])>0),";DP=",all_DP,";RO=",all_RO,";AO=",all_AO,";SRF=",sum(Rp),";SRR=",sum(Rm),";SAF=",sum(Vp),";SAR=",sum(Vm),";SOR=",all_sor,";RVSB=",all_rvsb,";FS=",FisherStrand_all,";ERR=",reg_res$coef["slope"],";SIG=",reg_res$coef["sigma"],";CONT=",paste(before,after,sep="x"),ifelse(reg_res$extra_rob & !(ref_inv),";WARN=EXTRA_ROBUST_GL",""),ifelse(reg_res$extra_rob & ref_inv,";WARN=EXTRA_ROBUST_GL/INV_REF",""),ifelse(ref_inv & !(reg_res$extra_rob),";WARN=INV_REF",""),sep="",file=out_file,append=T)
            # FORMAT field
            cat("\t",paste("GT:",ifelse(ref_inv,"QVAL_INV","QVAL"),":DP:RO:AO:AF:SB:SOR:RVSB:FS:QVAL_minAF:STATUS",sep=""),sep = "",file=out_file,append=T)

            # all samples
            if(ref_inv) { genotype=rep("1/1",l=nindiv) } else { genotype=rep("0/0",l=nindiv) }
            heterozygotes=which(reg_res$GQ>=GQ_threshold & sbs<=SB_threshold_SNV & reg_res$ma_count/reg_res$coverage < 0.75)
            genotype[heterozygotes]="0/1"
            homozygotes=which(reg_res$GQ>=GQ_threshold & sbs<=SB_threshold_SNV & reg_res$ma_count/reg_res$coverage >= 0.75)
            if(ref_inv) { genotype[homozygotes]="0/0" }  else { genotype[homozygotes]="1/1" }
            #no tumor variant but low power -> "./."
            genotype[(reg_res$GQ < GQ_threshold)&(qval_minAF<GQ_threshold) ]="./."
            for (cur_sample in 1:nindiv) {
              cat("\t",genotype[cur_sample],":",reg_res$GQ[cur_sample],":",DP[cur_sample],":",(Rp+Rm)[cur_sample],":",ma_count[cur_sample],":",(ma_count/DP)[cur_sample],":",Rp[cur_sample],",",Rm[cur_sample],",",Vp[cur_sample],",",Vm[cur_sample],":",sors[cur_sample],":",rvsbs[cur_sample],":",FisherStrand[cur_sample],":",qval_minAF[cur_sample],":",somatic_status[cur_sample],sep = "",file=out_file,append=T)
            }
            cat("\n",sep = "",file=out_file,append=T)

            if (do_plots=="ALL") { #output all variants
              if(do_alignments==TRUE){ #add alignment plots
                if(isTNpairs){
                  pdf(paste(linepos[1],"_",linepos[2],"_",linepos[2],"_",ref,"_",alt,ifelse(ref_inv,"_inv_ref",""),ifelse(reg_res$extra_rob,"_extra_robust",""),".pdf",sep=""),11,12)
                  par(mar=c(12,7,12,7))
                }else{
                  pdf(paste(linepos[1],"_",linepos[2],"_",linepos[2],"_",ref,"_",alt,ifelse(ref_inv,"_inv_ref",""),ifelse(reg_res$extra_rob,"_extra_robust",""),".pdf",sep=""),9,8)
                  par(mar=c(8,7,8,7))
                }
                #add regression plots to the pdf output
                plot_rob_nb(reg_res, 10^-(GQ_threshold/10), plot_title=bquote(paste(.(linepos[1]),":",.(linepos[2])," (",.(ref) %->% .(alt),")",.(ifelse(ref_inv," INV REF","")),.(ifelse(reg_res$extra_rob," EXTRA ROBUST","")),sep="")), sbs=sbs, SB_threshold=SB_threshold_SNV,plot_labels=plot_labels,add_contours=add_contours,names=indiv_run[,2])
                par(mar=c(1,1,1,1))
                #add alignment plots to the pdf output
                plotGviz(isTNpairs,sTrack,ref_genome,txdb,annotation,indiv_run,linepos,genotype,somatic_status,do_plots,Tindex,Nindex,onlyTindex,onlyNindex,bam_folder)
                dev.off()
              }else{
                #output only regression plots
                pdf(paste(linepos[1],"_",linepos[2],"_",linepos[2],"_",ref,"_",alt,ifelse(ref_inv,"_inv_ref",""),ifelse(reg_res$extra_rob,"_extra_robust",""),".pdf",sep=""),7,6)
                plot_rob_nb(reg_res, 10^-(GQ_threshold/10), plot_title=bquote(paste(.(linepos[1]),":",.(linepos[2])," (",.(ref) %->% .(alt),")",.(ifelse(ref_inv," INV REF","")),.(ifelse(reg_res$extra_rob," EXTRA ROBUST","")),sep="")), sbs=sbs, SB_threshold=SB_threshold_SNV,plot_labels=plot_labels,add_contours=add_contours,names=indiv_run[,2])
                dev.off()
              }
            }else if(do_plots=="SOMATIC"){ #output only somatic variants
              if("SOMATIC" %in% somatic_status){
                if(do_alignments==TRUE){ #add alignment plots
                  pdf(paste(linepos[1],"_",linepos[2],"_",linepos[2],"_",ref,"_",alt,ifelse(ref_inv,"_inv_ref",""),ifelse(reg_res$extra_rob,"_extra_robust",""),".pdf",sep=""),11,12)
                  par(mar=c(12,7,12,7))
                  #add regression plots to the pdf output
                  plot_rob_nb(reg_res, 10^-(GQ_threshold/10), plot_title=bquote(paste(.(linepos[1]),":",.(linepos[2])," (",.(ref) %->% .(alt),")",.(ifelse(ref_inv," INV REF","")),.(ifelse(reg_res$extra_rob," EXTRA ROBUST","")),sep="")), sbs=sbs, SB_threshold=SB_threshold_SNV,plot_labels=plot_labels,add_contours=add_contours,names=indiv_run[,2])
                  par(mar=c(1,1,1,1))
                  #add alignment plots to the pdf output
                  plotGviz(isTNpairs,sTrack,ref_genome,txdb,annotation,indiv_run,linepos,genotype,somatic_status,do_plots,Tindex,Nindex,onlyTindex,onlyNindex,bam_folder)
                  dev.off()
                }else{
                  #output only regression plots
                  pdf(paste(linepos[1],"_",linepos[2],"_",linepos[2],"_",ref,"_",alt,ifelse(ref_inv,"_inv_ref",""),ifelse(reg_res$extra_rob,"_extra_robust",""),".pdf",sep=""),7,6)
                  plot_rob_nb(reg_res, 10^-(GQ_threshold/10), plot_title=bquote(paste(.(linepos[1]),":",.(linepos[2])," (",.(ref) %->% .(alt),")",.(ifelse(ref_inv," INV REF","")),.(ifelse(reg_res$extra_rob," EXTRA ROBUST","")),sep="")), sbs=sbs, SB_threshold=SB_threshold_SNV,plot_labels=plot_labels,add_contours=add_contours,names=indiv_run[,2])
                  dev.off()
                }
              }
            }

          }
        }
      }

      #DEL
      ldel=linepos[deletion_cols]
      names(ldel)=id_samples
      all_del=ldel[!(ldel=="NA")]
      if (length(all_del)>0) {
        all_del=as.data.frame(strsplit(unlist(strsplit(paste(all_del),split = "|",fixed=T)),split = ":",fixed=T),stringsAsFactors=F,)
        uniq_del=unique(toupper(as.character(all_del[2,])))
        coverage_all_del=sum(as.numeric(all_del[1,]))
        for (cur_del in uniq_del) {
          Vp=rep(0,nindiv)
          names(Vp)=id_samples
          Vm=Vp
          all_cur_del=strsplit(grep(cur_del,ldel,ignore.case=T,value=T),split="|",fixed=T)
          all_cur_del_p=lapply(all_cur_del,function(x) {grep(paste(":",cur_del,"$",sep=""),x,value=T)})
          all_cur_del_m=lapply(all_cur_del,function(x) {grep(paste(":",tolower(cur_del),"$",sep=""),x,value=T)})
          ma_p_cur_del=unlist(lapply(all_cur_del_p,function(x) {as.numeric(gsub(paste("([0-9]+):",cur_del,"$",sep=""),"\\1",x,ignore.case = T))}))
          ma_m_cur_del=unlist(lapply(all_cur_del_m,function(x) {as.numeric(gsub(paste("([0-9]+):",cur_del,"$",sep=""),"\\1",x,ignore.case = T))}))
          Vp[names(ma_p_cur_del)]=ma_p_cur_del
          Vm[names(ma_m_cur_del)]=ma_m_cur_del
          ma_count=Vp+Vm
          DP=as.numeric(linepos[eval(depth)])
          if( sum( (ma_count/DP) > 0.8 , na.rm = T) > 0.5*length(ma_count) ){  #here we need to reverse alt and ref for the regression (use ma_count of the ref)
            ref=linepos[3]
            cur_del_inv=linepos[3]
            Vp_inv=as.numeric(linepos[eval(as.name(paste(cur_del_inv,"_cols",sep="")))])
            Vm_inv=as.numeric(linepos[eval(as.name(paste(tolower(cur_del_inv),"_cols",sep="")))])
            ma_count_inv=Vp_inv+Vm_inv
            ref_inv=TRUE
            reg_res=glmrob.nb(x=DP,y=ma_count_inv,min_coverage=min_coverage,min_reads=min_reads,min_af=min_af,extra_rob=extra_rob,min_af_extra_rob=min_af_extra_rob,min_prop_extra_rob=min_prop_extra_rob,max_prop_extra_rob=max_prop_extra_rob)
          } else { ref_inv=FALSE; reg_res=glmrob.nb(x=DP,y=ma_count,min_coverage=min_coverage,min_reads=min_reads,min_af=min_af,extra_rob=extra_rob,min_af_extra_rob=min_af_extra_rob,min_prop_extra_rob=min_prop_extra_rob,max_prop_extra_rob=max_prop_extra_rob) }
          # compute Qval20pc
          qval_minAF = rep(0,nindiv)
          somatic_status = rep(".",nindiv)
          if(isTNpairs){
            qval_minAF[Nindex] = sapply(1:length(Nindex),function(ii) toQvalueN(DP[Nindex][ii],reg_res,sigma) )
            if(sum(reg_res$GQ[Nindex]>GQ_threshold)>0) qval_minAF[Tindex] = sapply(1:length(Tindex),function(ii) toQvalueN(DP[Tindex][ii],reg_res,sigma) ) #GERMLINE VARIANT, check if power to CONFIRM in tumors
            else qval_minAF[Tindex] = sapply(1:length(Tindex),function(ii) toQvalueT(DP[Tindex][ii],reg_res,afmin_power) ) #no germline variant, check if power to call SOMATIC in tumors
            if( length(onlyNindex)>0 ) qval_minAF[onlyNindex] = sapply(1:length(onlyNindex),function(ii) toQvalueN(DP[onlyNindex][ii],reg_res,sigma) )
            if( length(onlyTindex)>0 ) qval_minAF[onlyTindex] = sapply(1:length(onlyTindex),function(ii) toQvalueT(DP[onlyTindex][ii],reg_res,afmin_power) )

            #no matching normal -> UNKNOWN (impossible to call somatic status)
            somatic_status[onlyTindex][(reg_res$GQ[onlyTindex] > GQ_threshold) ] = "UNKNOWN"
            #no matching tumor -> GERMLINE_UNCONFIRMABLE (impossible to call somatic status)
            somatic_status[onlyNindex][(reg_res$GQ[onlyNindex] > GQ_threshold) ] = "GERMLINE_UNCONFIRMABLE"
            #no tumor variant, no normal variant despite good power in both -> ".", confirmed; already set by default
            #no tumor variant, no normal variant but low power in N and T-> ".", unknown; already set by default
            #no tumor variant, no normal variant but low power in N or T-> ".", unconfirmed ; already set by default
            #no tumor variant, normal variant with good power in T (with or without good power in N)-> "GERMLINE_UNCONFIRMED"
            somatic_status[Tindex][(reg_res$GQ[Tindex] < GQ_threshold)&(qval_minAF[Tindex]>GQ_threshold)&(reg_res$GQ[Nindex]>GQ_threshold) ] = "GERMLINE_UNCONFIRMED"
            somatic_status[Nindex][(reg_res$GQ[Tindex] < GQ_threshold)&(qval_minAF[Tindex]>GQ_threshold)&(reg_res$GQ[Nindex]>GQ_threshold) ] = "GERMLINE_UNCONFIRMED"
            #no tumor variant, normal variant with low power in T (with or without good power in N)-> "GERMLINE_UNCONFIRMABLE"
            somatic_status[Tindex][(reg_res$GQ[Tindex] < GQ_threshold)&(qval_minAF[Tindex]<GQ_threshold)&(reg_res$GQ[Nindex]>GQ_threshold) ] = "GERMLINE_UNCONFIRMABLE"
            somatic_status[Nindex][(reg_res$GQ[Tindex] < GQ_threshold)&(qval_minAF[Tindex]<GQ_threshold)&(reg_res$GQ[Nindex]>GQ_threshold) ] = "GERMLINE_UNCONFIRMABLE"

            #tumor variant, no normal variant but without enough power -> UNKNOWN for both Tumor and Normal
            somatic_status[Tindex][(reg_res$GQ[Tindex] > GQ_threshold)&(qval_minAF[Nindex]<GQ_threshold)&(reg_res$GQ[Nindex]<GQ_threshold) ] = "UNKNOWN"
            somatic_status[Nindex][(reg_res$GQ[Tindex] > GQ_threshold)&(qval_minAF[Nindex]<GQ_threshold)&(reg_res$GQ[Nindex]<GQ_threshold) ] = "UNKNOWN"
            #tumor variant, normal variant (with or without good power in T and N) -> GERMLINE_CONFIRMED
            somatic_status[Tindex][ (reg_res$GQ[Tindex] > GQ_threshold)&(reg_res$GQ[Nindex]>GQ_threshold) ] = "GERMLINE_CONFIRMED"
            somatic_status[Nindex][ (reg_res$GQ[Tindex] > GQ_threshold)&(reg_res$GQ[Nindex]>GQ_threshold) ] = "GERMLINE_CONFIRMED"
            #tumor variant, no normal variant despite good power -> SOMATIC
            somatic_status[Tindex][(reg_res$GQ[Tindex] > GQ_threshold)&(qval_minAF[Nindex]>GQ_threshold)&(reg_res$GQ[Nindex]<GQ_threshold) ] = "SOMATIC"

            #flag possible contamination
            wh.germ = grep("GERMLINE|UNKNOWN",somatic_status)
            if( length(wh.germ)>0 ) somatic_status[Tindex][somatic_status[Tindex] == "SOMATIC"] = paste("POSSIBLE_CONTAMINATION_FROM", paste(indiv_run[wh.germ,2],sep="_",collapse="_"),sep="_")
          }else{# no pairs file supplied
            if(afmin_power==-1 ){#no minimum frequency supplied -> use a negative binomial distribution to check the power
              qval_minAF = sapply(1:nindiv,function(ii) toQvalueN(DP[ii],reg_res,sigma) )
            }else{#minimum frequency supplied -> use a binomial distribution centered around afmin_power to check the power
              qval_minAF = sapply(1:nindiv,function(ii) toQvalueT(DP[ii],reg_res,afmin_power) )
            }
          }

          if (!is.na(reg_res$coef["slope"]) & sum(reg_res$GQ>=GQ_threshold,na.rm=TRUE)>0) {
            all_AO=sum(ma_count)
            all_DP=sum(as.numeric(linepos[eval(depth)]))+sum(ma_count)
            common_annot()
            all_RO=sum(Rp+Rm)
            if (SB_type=="SOR") sbs=sors else { if (SB_type=="RVSB") sbs=rvsbs else {sbs=FisherStrand} }
            if (sum(sbs<=SB_threshold_indel & reg_res$GQ>=GQ_threshold,na.rm=TRUE)>0) {
              con=pipe(paste(samtools," faidx ",fasta_ref," ",linepos[1],":",as.numeric(linepos[2])+1-3,"-",as.numeric(linepos[2])+1-1," | tail -n1",sep=""))
              before=toupper(readLines(con))
              close(con)
              con=pipe(paste(samtools," faidx ",fasta_ref," ",linepos[1],":",as.numeric(linepos[2])+1+nchar(cur_del),"-",as.numeric(linepos[2])+1+3+nchar(cur_del)-1," | tail -n1",sep=""))
              after=toupper(readLines(con))
              close(con)
              prev_bp=substr(before,3,3)
              next_bp=substr(after,1,1)
              cat(linepos[1],"\t",as.numeric(linepos[2]),"\t",".","\t",paste(prev_bp,cur_del,sep=""),"\t",prev_bp,"\t",max(reg_res$GQ),"\t",".",sep = "",file=out_file,append=T)
              # INFO field
              cat("\t","TYPE=del;NS=",sum(as.numeric(linepos[eval(depth)])>0),";AF=",sum(reg_res$GQ>=GQ_threshold)/sum(as.numeric(linepos[eval(depth)])>0),";DP=",all_DP,";RO=",all_RO,";AO=",all_AO,";SRF=",sum(Rp),";SRR=",sum(Rm),";SAF=",sum(Vp),";SAR=",sum(Vm),";SOR=",all_sor,";RVSB=",all_rvsb,";FS=",FisherStrand_all,";ERR=",reg_res$coef["slope"],";SIG=",reg_res$coef["sigma"],";CONT=",paste(before,after,sep="x"),ifelse(reg_res$extra_rob & !(ref_inv),";WARN=EXTRA_ROBUST_GL",""),ifelse(reg_res$extra_rob & ref_inv,";WARN=EXTRA_ROBUST_GL/INV_REF",""),ifelse(ref_inv & !(reg_res$extra_rob),";WARN=INV_REF",""),sep="",file=out_file,append=T)
              # FORMAT field
              cat("\t",paste("GT:",ifelse(ref_inv,"QVAL_INV","QVAL"),":DP:RO:AO:AF:SB:SOR:RVSB:FS:QVAL_minAF:STATUS",sep=""),sep = "",file=out_file,append=T)
              # all samples
              if(ref_inv) { genotype=rep("1/1",l=nindiv) } else { genotype=rep("0/0",l=nindiv) }
              heterozygotes=which(reg_res$GQ>=GQ_threshold & sbs<=SB_threshold_SNV & reg_res$ma_count/reg_res$coverage < 0.75)
              genotype[heterozygotes]="0/1"
              homozygotes=which(reg_res$GQ>=GQ_threshold & sbs<=SB_threshold_SNV & reg_res$ma_count/reg_res$coverage >= 0.75)
              if(ref_inv) { genotype[homozygotes]="0/0" }  else { genotype[homozygotes]="1/1" }
              #no tumor variant but low power -> "./."
              genotype[(reg_res$GQ < GQ_threshold)&(qval_minAF<GQ_threshold) ]="./."

              for (cur_sample in 1:nindiv) {
                cat("\t",genotype[cur_sample],":",reg_res$GQ[cur_sample],":",DP[cur_sample],":",(Rp+Rm)[cur_sample],":",ma_count[cur_sample],":",(ma_count/DP)[cur_sample],":",Rp[cur_sample],",",Rm[cur_sample],",",Vp[cur_sample],",",Vm[cur_sample],":",sors[cur_sample],":",rvsbs[cur_sample],":",FisherStrand[cur_sample],":",qval_minAF[cur_sample],":",somatic_status[cur_sample],sep = "",file=out_file,append=T)
              }

              cat("\n",sep = "",file=out_file,append=T)

              if(!ref_inv & nchar(cur_del)>50) cur_del = paste("DEL",nchar(cur_del),sep="") #paste(substr(cur_del,1,5+match(cur_del,uniq_del)),substr(cur_del,nchar(cur_del)-(5+match(cur_del,uniq_del)),nchar(cur_del)),sep="...")
              if(ref_inv & nchar(ref)>50) ref = paste("DEL",nchar(ref),sep="") #paste(substr(ref,1,5+match(ref,uniq_del)),substr(ref,nchar(ref)-(5+match(ref,uniq_del)),nchar(ref)),sep="...")

              if (do_plots=="ALL") { #output all variants
                if(do_alignments==TRUE){ #add alignment plots

                  # deletions are shifted in samtools mpileup by 1bp, so put them at the right place by adding + to pos_ref[i,"loc"] everywhere in what follows

                  if(isTNpairs){
                    pdf(paste(linepos[1],"_",linepos[2],"_",as.numeric(linepos[2])+nchar(cur_del)-1,"_",paste(prev_bp,cur_del,sep=""),"_",prev_bp,ifelse(ref_inv,"_inv_ref",""),ifelse(reg_res$extra_rob,"_extra_robust",""),".pdf",sep=""),11,12)
                    par(mar=c(12,7,12,7))
                  }else{
                    pdf(paste(linepos[1],"_",linepos[2],"_",as.numeric(linepos[2])+nchar(cur_del)-1,"_",paste(prev_bp,cur_del,sep=""),"_",prev_bp,ifelse(ref_inv,"_inv_ref",""),ifelse(reg_res$extra_rob,"_extra_robust",""),".pdf",sep=""),9,8)
                    par(mar=c(8,7,8,7))
                  }
                  #add regression plots to the pdf output
                  plot_rob_nb(reg_res, 10^-(GQ_threshold/10), plot_title=bquote(paste(.(linepos[1]),":",.(as.numeric(linepos[2]))," (",.(paste(prev_bp,cur_del,sep="")) %->% .(prev_bp),")",.(ifelse(ref_inv," INV REF","")),.(ifelse(reg_res$extra_rob," EXTRA ROBUST","")),sep="")),sbs=sbs, SB_threshold=SB_threshold_indel,plot_labels=plot_labels,add_contours=add_contours,names=indiv_run[,2])
                  par(mar=c(1,1,1,1))
                  #add alignment plots to the pdf output
                  plotGviz(isTNpairs,sTrack,ref_genome,txdb,annotation,indiv_run,linepos,genotype,somatic_status,do_plots,Tindex,Nindex,onlyTindex,onlyNindex,bam_folder)
                  dev.off()
                }else{
                  #output only regression plots
                  pdf(paste(linepos[1],"_",linepos[2],"_",as.numeric(linepos[2])+nchar(cur_del)-1,"_",paste(prev_bp,cur_del,sep=""),"_",prev_bp,ifelse(ref_inv,"_inv_ref",""),ifelse(reg_res$extra_rob,"_extra_robust",""),".pdf",sep=""),7,6)
                  plot_rob_nb(reg_res, 10^-(GQ_threshold/10), plot_title=bquote(paste(.(linepos[1]),":",.(as.numeric(linepos[2]))," (",.(paste(prev_bp,cur_del,sep="")) %->% .(prev_bp),")",.(ifelse(ref_inv," INV REF","")),.(ifelse(reg_res$extra_rob," EXTRA ROBUST","")),sep="")),sbs=sbs, SB_threshold=SB_threshold_indel,plot_labels=plot_labels,add_contours=add_contours,names=indiv_run[,2])
                  dev.off()
                }
              }else if(do_plots=="SOMATIC"){ #output only somatic variants
                if("SOMATIC" %in% somatic_status){
                  if(do_alignments==TRUE){ #add alignment plots
                    pdf(paste(linepos[1],"_",linepos[2],"_",as.numeric(linepos[2])+nchar(cur_del)-1,"_",paste(prev_bp,cur_del,sep=""),"_",prev_bp,ifelse(ref_inv,"_inv_ref",""),ifelse(reg_res$extra_rob,"_extra_robust",""),".pdf",sep=""),11,12)
                    par(mar=c(12,7,12,7))
                    #add regression plots to the pdf output
                    plot_rob_nb(reg_res, 10^-(GQ_threshold/10), plot_title=bquote(paste(.(linepos[1]),":",.(as.numeric(linepos[2]))," (",.(paste(prev_bp,cur_del,sep="")) %->% .(prev_bp),")",.(ifelse(ref_inv," INV REF","")),.(ifelse(reg_res$extra_rob," EXTRA ROBUST","")),sep="")),sbs=sbs, SB_threshold=SB_threshold_indel,plot_labels=plot_labels,add_contours=add_contours,names=indiv_run[,2])
                    par(mar=c(1,1,1,1))
                    #add alignment plots to the pdf output
                    plotGviz(isTNpairs,sTrack,ref_genome,txdb,annotation,indiv_run,linepos,genotype,somatic_status,do_plots,Tindex,Nindex,onlyTindex,onlyNindex,bam_folder)
                    dev.off()
                  }else{
                    #output only regression plots
                    pdf(paste(linepos[1],"_",linepos[2],"_",as.numeric(linepos[2])+nchar(cur_del)-1,"_",paste(prev_bp,cur_del,sep=""),"_",prev_bp,ifelse(ref_inv,"_inv_ref",""),ifelse(reg_res$extra_rob,"_extra_robust",""),".pdf",sep=""),7,6)
                    plot_rob_nb(reg_res, 10^-(GQ_threshold/10), plot_title=bquote(paste(.(linepos[1]),":",.(as.numeric(linepos[2]))," (",.(paste(prev_bp,cur_del,sep="")) %->% .(prev_bp),")",.(ifelse(ref_inv," INV REF","")),.(ifelse(reg_res$extra_rob," EXTRA ROBUST","")),sep="")),sbs=sbs, SB_threshold=SB_threshold_indel,plot_labels=plot_labels,add_contours=add_contours,names=indiv_run[,2])
                    dev.off()
                  }
                }
              }


            }
          }
        }
      }


      # INS
      lins=linepos[insertion_cols]
      names(lins)=id_samples
      all_ins=lins[!(lins=="NA")]
      if (length(all_ins)>0) {
        all_ins=as.data.frame(strsplit(unlist(strsplit(paste(all_ins),split = "|",fixed=T)),split = ":",fixed=T),stringsAsFactors=F,)
        uniq_ins=unique(toupper(as.character(all_ins[2,])))
        coverage_all_ins=sum(as.numeric(all_ins[1,]))
        for (cur_ins in uniq_ins) {
          Vp=rep(0,nindiv)
          names(Vp)=id_samples
          Vm=Vp
          all_cur_ins=strsplit(grep(cur_ins,lins,ignore.case=T,value=T),split="|",fixed=T)
          all_cur_ins_p=lapply(all_cur_ins,function(x) {grep(paste(":",cur_ins,"$",sep=""),x,value=T)})
          all_cur_ins_m=lapply(all_cur_ins,function(x) {grep(paste(":",tolower(cur_ins),"$",sep=""),x,value=T)})
          ma_p_cur_ins=unlist(lapply(all_cur_ins_p,function(x) {as.numeric(gsub(paste("([0-9]+):",cur_ins,"$",sep=""),"\\1",x,ignore.case = T))}))
          ma_m_cur_ins=unlist(lapply(all_cur_ins_m,function(x) {as.numeric(gsub(paste("([0-9]+):",cur_ins,"$",sep=""),"\\1",x,ignore.case = T))}))
          Vp[names(ma_p_cur_ins)]=ma_p_cur_ins
          Vm[names(ma_m_cur_ins)]=ma_m_cur_ins
          ma_count=Vp+Vm
          DP=as.numeric(linepos[eval(depth)])
          if( sum( (ma_count/DP) > 0.8 , na.rm = T) > 0.5*length(ma_count) ){  #here we need to reverse alt and ref for the regression (use ma_count of the ref)
            ref=linepos[3]
            cur_ins_inv=linepos[3]
            Vp_inv=as.numeric(linepos[eval(as.name(paste(cur_ins_inv,"_cols",sep="")))])
            Vm_inv=as.numeric(linepos[eval(as.name(paste(tolower(cur_ins_inv),"_cols",sep="")))])
            ma_count_inv=Vp_inv+Vm_inv
            ref_inv=TRUE
            reg_res=glmrob.nb(x=DP,y=ma_count_inv,min_coverage=min_coverage,min_reads=min_reads,min_af=min_af,extra_rob=extra_rob,min_af_extra_rob=min_af_extra_rob,min_prop_extra_rob=min_prop_extra_rob,max_prop_extra_rob=max_prop_extra_rob)
          } else { ref_inv=FALSE; reg_res=glmrob.nb(x=DP,y=ma_count,min_coverage=min_coverage,min_reads=min_reads,min_af=min_af,extra_rob=extra_rob,min_af_extra_rob=min_af_extra_rob,min_prop_extra_rob=min_prop_extra_rob,max_prop_extra_rob=max_prop_extra_rob) }
          # compute Qval20pc
          qval_minAF = rep(0,nindiv)
          somatic_status = rep(".",nindiv)

          if(isTNpairs){
            qval_minAF[Nindex] = sapply(1:length(Nindex),function(ii) toQvalueN(DP[Nindex][ii],reg_res,sigma) )
            if(sum(reg_res$GQ[Nindex]>GQ_threshold)>0) qval_minAF[Tindex] = sapply(1:length(Tindex),function(ii) toQvalueN(DP[Tindex][ii],reg_res,sigma) ) #GERMLINE VARIANT, check if power to CONFIRM in tumors
            else qval_minAF[Tindex] = sapply(1:length(Tindex),function(ii) toQvalueT(DP[Tindex][ii],reg_res,afmin_power) ) #no germline variant, check if power to call SOMATIC in tumors
            if( length(onlyNindex)>0 ) qval_minAF[onlyNindex] = sapply(1:length(onlyNindex),function(ii) toQvalueN(DP[onlyNindex][ii],reg_res,sigma) )
            if( length(onlyTindex)>0 ) qval_minAF[onlyTindex] = sapply(1:length(onlyTindex),function(ii) toQvalueT(DP[onlyTindex][ii],reg_res,afmin_power) )

            #no matching normal -> UNKNOWN (impossible to call somatic status)
            somatic_status[onlyTindex][(reg_res$GQ[onlyTindex] > GQ_threshold) ] = "UNKNOWN"
            #no matching tumor -> GERMLINE_UNCONFIRMABLE (impossible to call somatic status)
            somatic_status[onlyNindex][(reg_res$GQ[onlyNindex] > GQ_threshold) ] = "GERMLINE_UNCONFIRMABLE"
            #no tumor variant, no normal variant despite good power in both -> ".", confirmed; already set by default
            #no tumor variant, no normal variant but low power in N and T-> ".", unknown; already set by default
            #no tumor variant, no normal variant but low power in N or T-> ".", unconfirmed ; already set by default
            #no tumor variant, normal variant with good power in T (with or without good power in N)-> "GERMLINE_UNCONFIRMED"
            somatic_status[Tindex][(reg_res$GQ[Tindex] < GQ_threshold)&(qval_minAF[Tindex]>GQ_threshold)&(reg_res$GQ[Nindex]>GQ_threshold) ] = "GERMLINE_UNCONFIRMED"
            somatic_status[Nindex][(reg_res$GQ[Tindex] < GQ_threshold)&(qval_minAF[Tindex]>GQ_threshold)&(reg_res$GQ[Nindex]>GQ_threshold) ] = "GERMLINE_UNCONFIRMED"
            #no tumor variant, normal variant with low power in T (with or without good power in N)-> "GERMLINE_UNCONFIRMABLE"
            somatic_status[Tindex][(reg_res$GQ[Tindex] < GQ_threshold)&(qval_minAF[Tindex]<GQ_threshold)&(reg_res$GQ[Nindex]>GQ_threshold) ] = "GERMLINE_UNCONFIRMABLE"
            somatic_status[Nindex][(reg_res$GQ[Tindex] < GQ_threshold)&(qval_minAF[Tindex]<GQ_threshold)&(reg_res$GQ[Nindex]>GQ_threshold) ] = "GERMLINE_UNCONFIRMABLE"

            #tumor variant, no normal variant but without enough power -> UNKNOWN for both Tumor and Normal
            somatic_status[Tindex][(reg_res$GQ[Tindex] > GQ_threshold)&(qval_minAF[Nindex]<GQ_threshold)&(reg_res$GQ[Nindex]<GQ_threshold) ] = "UNKNOWN"
            somatic_status[Nindex][(reg_res$GQ[Tindex] > GQ_threshold)&(qval_minAF[Nindex]<GQ_threshold)&(reg_res$GQ[Nindex]<GQ_threshold) ] = "UNKNOWN"
            #tumor variant, normal variant (with or without good power in T and N) -> GERMLINE_CONFIRMED
            somatic_status[Tindex][ (reg_res$GQ[Tindex] > GQ_threshold)&(reg_res$GQ[Nindex]>GQ_threshold) ] = "GERMLINE_CONFIRMED"
            somatic_status[Nindex][ (reg_res$GQ[Tindex] > GQ_threshold)&(reg_res$GQ[Nindex]>GQ_threshold) ] = "GERMLINE_CONFIRMED"
            #tumor variant, no normal variant despite good power -> SOMATIC
            somatic_status[Tindex][(reg_res$GQ[Tindex] > GQ_threshold)&(qval_minAF[Nindex]>GQ_threshold)&(reg_res$GQ[Nindex]<GQ_threshold) ] = "SOMATIC"

            #flag possible contamination
            wh.germ = grep("GERMLINE|UNKNOWN",somatic_status)
            if( length(wh.germ)>0 ) somatic_status[Tindex][somatic_status[Tindex] == "SOMATIC"] = paste("POSSIBLE_CONTAMINATION_FROM", paste(indiv_run[wh.germ,2],sep="_",collapse="_"),sep="_")
          }else{# no pairs file supplied
            if(afmin_power==-1 ){#no minimum frequency supplied -> use a negative binomial distribution to check the power
              qval_minAF = sapply(1:nindiv,function(ii) toQvalueN(DP[ii],reg_res,sigma) )
            }else{#minimum frequency supplied -> use a binomial distribution centered around afmin_power to check the power
              qval_minAF = sapply(1:nindiv,function(ii) toQvalueT(DP[ii],reg_res,afmin_power) )
            }
          }

          if (!is.na(reg_res$coef["slope"]) & sum(reg_res$GQ>=GQ_threshold,na.rm=TRUE)>0) {
            all_AO=sum(ma_count)
            all_DP=sum(as.numeric(linepos[eval(depth)]))+sum(ma_count)
            common_annot()
            all_RO=sum(Rp+Rm)
            if (SB_type=="SOR") sbs=sors else { if (SB_type=="RVSB") sbs=rvsbs else {sbs=FisherStrand} }
            if (sum(sbs<=SB_threshold_indel & reg_res$GQ>=GQ_threshold,na.rm=TRUE)>0) {
              con=pipe(paste(samtools," faidx ",fasta_ref," ",linepos[1],":",as.numeric(linepos[2])-2,"-",as.numeric(linepos[2])," | tail -n1",sep=""))
              before=toupper(readLines(con))
              close(con)
              con=pipe(paste(samtools," faidx ",fasta_ref," ",linepos[1],":",as.numeric(linepos[2])+1,"-",as.numeric(linepos[2])+3," | tail -n1",sep=""))
              after=toupper(readLines(con))
              close(con)
              prev_bp=substr(before,3,3)
              cat(linepos[1],"\t",linepos[2],"\t",".","\t",prev_bp,"\t",paste(prev_bp,cur_ins,sep=""),"\t",max(reg_res$GQ),"\t",".",sep = "",file=out_file,append=T)
              # INFO field
              cat("\t","TYPE=ins;NS=",sum(as.numeric(linepos[eval(depth)])>0),";AF=",sum(reg_res$GQ>=GQ_threshold)/sum(as.numeric(linepos[eval(depth)])>0),";DP=",all_DP,";RO=",all_RO,";AO=",all_AO,";SRF=",sum(Rp),";SRR=",sum(Rm),";SAF=",sum(Vp),";SAR=",sum(Vm),";SOR=",all_sor,";RVSB=",all_rvsb,";FS=",FisherStrand_all,";ERR=",reg_res$coef["slope"],";SIG=",reg_res$coef["sigma"],";CONT=",paste(before,after,sep="x"),ifelse(reg_res$extra_rob & !(ref_inv),";WARN=EXTRA_ROBUST_GL",""),ifelse(reg_res$extra_rob & ref_inv,";WARN=EXTRA_ROBUST_GL/INV_REF",""),ifelse(ref_inv & !(reg_res$extra_rob),";WARN=INV_REF",""),sep="",file=out_file,append=T)
              # FORMAT field
              cat("\t",paste("GT:",ifelse(ref_inv,"QVAL_INV","QVAL"),":DP:RO:AO:AF:SB:SOR:RVSB:FS:QVAL_minAF:STATUS",sep=""),sep = "",file=out_file,append=T)
              # all samples
              if(ref_inv) { genotype=rep("1/1",l=nindiv) } else { genotype=rep("0/0",l=nindiv) }
              heterozygotes=which(reg_res$GQ>=GQ_threshold & sbs<=SB_threshold_SNV & reg_res$ma_count/reg_res$coverage < 0.75)
              genotype[heterozygotes]="0/1"
              homozygotes=which(reg_res$GQ>=GQ_threshold & sbs<=SB_threshold_SNV & reg_res$ma_count/reg_res$coverage >= 0.75)
              if(ref_inv) { genotype[homozygotes]="0/0" }  else { genotype[homozygotes]="1/1" }
              #no tumor variant but low power -> "./."
              genotype[(reg_res$GQ < GQ_threshold)&(qval_minAF<GQ_threshold) ]="./."

              for (cur_sample in 1:nindiv) {
                cat("\t",genotype[cur_sample],":",reg_res$GQ[cur_sample],":",DP[cur_sample],":",(Rp+Rm)[cur_sample],":",ma_count[cur_sample],":",(ma_count/DP)[cur_sample],":",Rp[cur_sample],",",Rm[cur_sample],",",Vp[cur_sample],",",Vm[cur_sample],":",sors[cur_sample],":",rvsbs[cur_sample],":",FisherStrand[cur_sample],":",qval_minAF[cur_sample],":",somatic_status[cur_sample],sep = "",file=out_file,append=T)
              }

              cat("\n",sep = "",file=out_file,append=T)

              if(!ref_inv & nchar(cur_ins)>50) cur_ins = paste("INS",nchar(cur_ins),sep="") #paste(substr(cur_ins,1,5+match(cur_ins,uniq_ins)),substr(cur_ins,nchar(cur_ins)-(5+match(cur_ins,uniq_ins)),nchar(cur_ins)),sep="...")
              if(ref_inv & nchar(ref)>50) ref = paste("INS",nchar(ref),sep="") #paste(substr(ref,1,5+match(ref,uniq_ins)),substr(ref,nchar(ref)-(5+match(ref,uniq_ins)),nchar(ref)),sep="...")

              if (do_plots=="ALL") { #output all variants

                if(do_alignments==TRUE){ #add alignment plots

                  if(isTNpairs){
                    pdf(paste(linepos[1],"_",linepos[2],"_",as.numeric(linepos[2]),"_",prev_bp,"_",paste(prev_bp,cur_ins,sep=""),ifelse(ref_inv,"_inv_ref",""),ifelse(reg_res$extra_rob,"_extra_robust",""),".pdf",sep=""),11,12)
                    par(mar=c(12,7,12,7))
                  }else{
                    pdf(paste(linepos[1],"_",linepos[2],"_",as.numeric(linepos[2]),"_",prev_bp,"_",paste(prev_bp,cur_ins,sep=""),ifelse(ref_inv,"_inv_ref",""),ifelse(reg_res$extra_rob,"_extra_robust",""),".pdf",sep=""),9,8)
                    par(mar=c(8,7,8,7))
                  }
                  #add regression plots to the pdf output
                  plot_rob_nb(reg_res, 10^-(GQ_threshold/10), plot_title=bquote(paste(.(linepos[1]),":",.(linepos[2])," (",.(prev_bp) %->% .(paste(prev_bp,cur_ins,sep="")),")",.(ifelse(ref_inv," INV REF","")),.(ifelse(reg_res$extra_rob," EXTRA ROBUST","")),sep="")),sbs=sbs, SB_threshold=SB_threshold_indel,plot_labels=plot_labels,add_contours=add_contours,names=indiv_run[,2])
                  par(mar=c(1,1,1,1))
                  #add alignment plots to the pdf output
                  plotGviz(isTNpairs,sTrack,ref_genome,txdb,annotation,indiv_run,linepos,genotype,somatic_status,do_plots,Tindex,Nindex,onlyTindex,onlyNindex,bam_folder)
                  dev.off()
                }else{
                  #output only regression plots
                  pdf(paste(linepos[1],"_",linepos[2],"_",as.numeric(linepos[2]),"_",prev_bp,"_",paste(prev_bp,cur_ins,sep=""),ifelse(ref_inv,"_inv_ref",""),ifelse(reg_res$extra_rob,"_extra_robust",""),".pdf",sep=""),7,6)
                  plot_rob_nb(reg_res, 10^-(GQ_threshold/10), plot_title=bquote(paste(.(linepos[1]),":",.(linepos[2])," (",.(prev_bp) %->% .(paste(prev_bp,cur_ins,sep="")),")",.(ifelse(ref_inv," INV REF","")),.(ifelse(reg_res$extra_rob," EXTRA ROBUST","")),sep="")),sbs=sbs, SB_threshold=SB_threshold_indel,plot_labels=plot_labels,add_contours=add_contours,names=indiv_run[,2])
                  dev.off()
                }
              }else if(do_plots=="SOMATIC"){ #output only somatic variants
                if("SOMATIC" %in% somatic_status){
                  if(do_alignments==TRUE){ #add alignment plots
                    pdf(paste(linepos[1],"_",linepos[2],"_",as.numeric(linepos[2]),"_",prev_bp,"_",paste(prev_bp,cur_ins,sep=""),ifelse(ref_inv,"_inv_ref",""),ifelse(reg_res$extra_rob,"_extra_robust",""),".pdf",sep=""),11,12)
                    par(mar=c(12,7,12,7))
                    #add regression plots to the pdf output
                    plot_rob_nb(reg_res, 10^-(GQ_threshold/10), plot_title=bquote(paste(.(linepos[1]),":",.(linepos[2])," (",.(prev_bp) %->% .(paste(prev_bp,cur_ins,sep="")),")",.(ifelse(ref_inv," INV REF","")),.(ifelse(reg_res$extra_rob," EXTRA ROBUST","")),sep="")),sbs=sbs, SB_threshold=SB_threshold_indel,plot_labels=plot_labels,add_contours=add_contours,names=indiv_run[,2])
                    par(mar=c(1,1,1,1))
                    #add alignment plots to the pdf output
                    plotGviz(isTNpairs,sTrack,ref_genome,txdb,annotation,indiv_run,linepos,genotype,somatic_status,do_plots,Tindex,Nindex,onlyTindex,onlyNindex,bam_folder)
                    dev.off()
                  }else{
                    #output only regression plots
                    pdf(paste(linepos[1],"_",linepos[2],"_",as.numeric(linepos[2]),"_",prev_bp,"_",paste(prev_bp,cur_ins,sep=""),ifelse(ref_inv,"_inv_ref",""),ifelse(reg_res$extra_rob,"_extra_robust",""),".pdf",sep=""),7,6)
                    plot_rob_nb(reg_res, 10^-(GQ_threshold/10), plot_title=bquote(paste(.(linepos[1]),":",.(linepos[2])," (",.(prev_bp) %->% .(paste(prev_bp,cur_ins,sep="")),")",.(ifelse(ref_inv," INV REF","")),.(ifelse(reg_res$extra_rob," EXTRA ROBUST","")),sep="")),sbs=sbs, SB_threshold=SB_threshold_indel,plot_labels=plot_labels,add_contours=add_contours,names=indiv_run[,2])
                    dev.off()
                  }
                }
              }

            }
          }
        }
      }


    }
  }
  i=i+1
}

print(Sys.time())
###############################################################################################
