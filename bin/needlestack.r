#! /usr/bin/env Rscript

# needlestack: a multi-sample somatic variant caller
# Copyright (C) 2015 Matthieu Foll and Tiffany Delhomme
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

############################## ARGUMENTS SECTION ##############################
args <- commandArgs(TRUE)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsL <- as.list(as.character(as.data.frame(do.call("rbind", parseArgs(args)))$V2))
names(argsL) <- as.data.frame(do.call("rbind", parseArgs(args)))$V1
args <- argsL;rm(argsL)

if("--help" %in% args | is.null(args$out_file) | is.null(args$fasta_ref) | is.null(args$source_path) ) {
  cat("
      The R Script arguments_section.R

      Mandatory arguments:
      --out_file=file_name           - name of output vcf
      --source_path=path             - path to source files (glm_rob_nb.r, plot_rob_nb.r)
      --fasta_ref=path               - path of fasta ref
      --help                         - print this text

      Optionnal arguments:
      --samtools=path                - path of samtools, default=samtools
      --SB_type=SOR, RVSB or FISH    - strand bias measure, default=SOR
      --SB_threshold_SNV=value       - strand bias threshold for SNV, default=100
      --SB_threshold_indel=value     - strand bias threshold for indel, default=100
      --min_coverage=value           - minimum coverage in at least one sample to consider a site, default=50
      --min_reads=value              - minimum number of non-ref reads in at least one sample to consider a site, default=5
      --GQ_threshold=value           - phred scale qvalue threshold for variants, default=50
      --output_all_SNVs=boolean     - output all SNVs, even when no variant is detected, default=FALSE
      --do_plots=boolean              - output regression plots, default=TRUE

      WARNING : by default samtools has to be in your path

      Example:
      pileup_nbrr_caller_vcf.r --out_file=test.vcf --fasta_ref=~/Documents/References/ \n\n")

  q(save="no")
}

if(is.null(args$samtools)) {args$samtools="samtools"}
if(is.null(args$SB_type)) {args$SB_type="SOR"}
if(is.null(args$SB_threshold_SNV)) {args$SB_threshold_SNV=100} else {args$SB_threshold_SNV=as.numeric(args$SB_threshold_SNV)}
if(is.null(args$SB_threshold_indel)) {args$SB_threshold_indel=100} else {args$SB_threshold_indel=as.numeric(args$SB_threshold_indel)}
if(is.null(args$min_coverage)) {args$min_coverage=50} else {args$min_coverage=as.numeric(args$min_coverage)}
if(is.null(args$min_reads)) {args$min_reads=5} else {args$min_reads=as.numeric(args$min_reads)}
if(is.null(args$GQ_threshold)) {args$GQ_threshold=50} else {args$GQ_threshold=as.numeric(args$GQ_threshold)}
if(is.null(args$output_all_SNVs)) {args$output_all_SNVs=FALSE} else {args$output_all_SNVs=as.logical(args$output_all_SNVs)}
if(is.null(args$do_plots)) {args$do_plots=TRUE} else {args$do_plots=as.logical(args$do_plots)}
if(is.null(args$plot_labels)) {args$plot_labels=FALSE} else {args$plot_labels=as.logical(args$plot_labels)}
if(is.null(args$add_contours)) {args$add_contours=FALSE} else {args$add_contours=as.logical(args$add_contours)}

samtools=args$samtools
out_file=args$out_file
fasta_ref=args$fasta_ref
GQ_threshold=args$GQ_threshold
min_coverage=args$min_coverage
min_reads=args$min_reads
# http://gatkforums.broadinstitute.org/discussion/5533/strandoddsratio-computation filter out SOR > 4 for SNVs and > 10 for indels
# filter out RVSB > 0.85 (maybe less stringent for SNVs)
SB_type=args$SB_type
SB_threshold_SNV=args$SB_threshold_SNV
SB_threshold_indel=args$SB_threshold_indel
output_all_SNVs=args$output_all_SNVs
do_plots=args$do_plots
plot_labels=args$plot_labels
add_contours=args$add_contours

source(paste(args$source_path,"glm_rob_nb.r",sep=""))
source(paste(args$source_path,"plot_rob_nb.r",sep=""))

############################################################

options("scipen"=100)

indiv_run=read.table("names.txt",stringsAsFactors=F,colClasses = "character")
indiv_run[,2]=make.unique(indiv_run[,2],sep="_")

pileups_files=paste("TABLE/",indiv_run[,1],".txt",sep="")
nindiv=nrow(indiv_run)

npos=length(readLines(pileups_files[1]))-1
atcg_matrix=matrix(nrow=npos,ncol=8*nindiv)
coverage_matrix=matrix(nrow=npos,ncol=nindiv)

ins=as.data.frame(setNames(replicate(nindiv,rep(NA,npos), simplify = F), indiv_run[,1]),optional=T)
del=ins
for (k in 1:nindiv) {
  cur_data=read.table(pileups_files[k],header = T,stringsAsFactors = F,sep="\t",colClasses = c("character","numeric","character","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","character","character"))
  if (k==1) {
  	pos_ref=cur_data[,1:3]
  	pos_ref[,"ref"]=toupper(pos_ref[,"ref"])
  }
  atcg_matrix[,((k-1)*8+1):(k*8)]=as.matrix(cur_data[,4:11])
  ins[,indiv_run[k,1]]=cur_data[,12]
  del[,indiv_run[k,1]]=cur_data[,13]
  coverage_matrix[,k]=rowSums(matrix(atcg_matrix[,((k-1)*8+1):(k*8)],nrow=npos))
}

A_cols=(1:nindiv)*8 - 7
T_cols=(1:nindiv)*8 - 6
C_cols=(1:nindiv)*8 - 5
G_cols=(1:nindiv)*8 - 4
a_cols=(1:nindiv)*8 - 3
t_cols=(1:nindiv)*8 - 2
c_cols=(1:nindiv)*8 - 1
g_cols=(1:nindiv)*8 - 0

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
  Rp<<-atcg_matrix[i,eval(as.name(paste(pos_ref[i,"ref"],"_cols",sep="")))]
  Rm<<-atcg_matrix[i,eval(as.name(paste(tolower(pos_ref[i,"ref"]),"_cols",sep="")))]
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
  FisherStrand<<-unlist(lapply(1:nindiv, function(indiv) fisher.test(matrix(c(Vp[indiv],Vm[indiv],Rp[indiv],Rm[indiv]), nrow=2))$p.value )) #col1=alt(V), col2=ref(R)
  FisherStrand_all<<--10*log10(fisher.test(matrix(c(sum(Vp),sum(Vm),sum(Rp),sum(Rm)), nrow=2))$p.value)
  all_RO<<-sum(Rp+Rm)
}

if(file.exists(out_file)) file.remove(out_file)
write_out=function(...) {
  cat(paste(...,sep=""),"\n",file=out_file,sep="",append=T)
}
write_out("##fileformat=VCFv4.1")
write_out("##fileDate=",format(Sys.Date(), "%Y%m%d"))
write_out("##source=needlestack v1.0b")
write_out("##reference=",fasta_ref)
write_out("##phasing=none")
write_out("##filter=\"QVAL > ",GQ_threshold," & ",SB_type,"_SNV < ",SB_threshold_SNV," & ",SB_type,"_INDEL < ",SB_threshold_indel," & min(AO) >= ",min_reads," & min(DP) >= ",min_coverage,"\"")

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
write_out("##INFO=<ID=FISH,Number=1,Type=Float,Description=\"Total Fisher Exact Test p-value for detecting strand bias (Phred-scale)\">")
write_out("##INFO=<ID=ERR,Number=1,Type=Float,Description=\"Estimated error rate for the alternate allele\">")
write_out("##INFO=<ID=SIG,Number=1,Type=Float,Description=\"Estimated overdispersion for the alternate allele\">")
write_out("##INFO=<ID=CONT,Number=1,Type=String,Description=\"Context of the reference sequence\">")

write_out("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")
write_out("##FORMAT=<ID=QVAL,Number=1,Type=Float,Description=\"Phred-scaled qvalue for not being an error\">")
write_out("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">")
write_out("##FORMAT=<ID=RO,Number=1,Type=Integer,Description=\"Reference allele observation count\">")
write_out("##FORMAT=<ID=AO,Number=1,Type=Integer,Description=\"Alternate allele observation count\">")
write_out("##FORMAT=<ID=AF,Number=1,Type=Float,Description=\"Allele fraction of the alternate allele with regard to reference\">")
write_out("##FORMAT=<ID=SB,Number=4,Type=Integer,Description=\"Per-sample component statistics to detect strand bias as SRF,SRR,SAF,SAR\">")
write_out("##FORMAT=<ID=SOR,Number=1,Type=Float,Description=\"Symmetric Odds Ratio of 2x2 contingency table to detect strand bias\">")
write_out("##FORMAT=<ID=RVSB,Number=1,Type=Float,Description=\"Relative Variant Strand Bias\">")
write_out("##FORMAT=<ID=FISH,Number=1,Type=Float,Description=\"Fisher Exact Test p-value for detecting strand bias (Phred-scale)\">")

write_out("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t",paste(indiv_run[,2],collapse = "\t"))

for (i in 1:npos) {
  if (is.element(pos_ref[i,"ref"],c("A","T","C","G"))) {
    # SNV
    for (alt in non_ref_bases(pos_ref[i,"ref"])) {
      Vp=atcg_matrix[i,eval(as.name(paste(alt,"_cols",sep="")))]
      Vm=atcg_matrix[i,eval(as.name(paste(tolower(alt),"_cols",sep="")))]
      ma_count=Vp+Vm
      DP=coverage_matrix[i,]
      reg_res=glmrob.nb(x=DP,y=ma_count,min_coverage=min_coverage,min_reads=min_reads)
      if (output_all_SNVs | (!is.na(reg_res$coef["slope"]) & sum(reg_res$GQ>=GQ_threshold,na.rm=TRUE)>0)) {
        all_AO=sum(ma_count)
        all_DP=sum(coverage_matrix[i,])
        common_annot()
        all_RO=sum(Rp+Rm)
        if (SB_type=="SOR") sbs=sors else { if (SB_type=="RVSB") sbs=rvsbs else {sbs=FisherStrand} }
        if (output_all_SNVs | (sum(sbs<=SB_threshold_SNV & reg_res$GQ>=GQ_threshold,na.rm=TRUE)>0)) {
          con=pipe(paste(samtools," faidx ",fasta_ref," ",pos_ref[i,"chr"],":",pos_ref[i,"loc"]-3,"-",pos_ref[i,"loc"]-1," | tail -n1",sep=""))
  		    before=readLines(con)
          close(con)
          con=pipe(paste(samtools," faidx ",fasta_ref," ",pos_ref[i,"chr"],":",pos_ref[i,"loc"]+1,"-",pos_ref[i,"loc"]+3," | tail -n1",sep=""))
          after=readLines(con)
  		    close(con)
  		    cat(pos_ref[i,"chr"],"\t",pos_ref[i,"loc"],"\t",".","\t",pos_ref[i,"ref"],"\t",alt,"\t",max(reg_res$GQ),"\t",".",sep = "",file=out_file,append=T)
          # INFO field
  		    cat("\t","TYPE=snv;NS=",sum(coverage_matrix[i,]>0),";AF=",sum(reg_res$GQ>=GQ_threshold)/sum(coverage_matrix[i,]>0),";DP=",all_DP,";RO=",all_RO,";AO=",all_AO,";SRF=",sum(Rp),";SRR=",sum(Rm),";SAF=",sum(Vp),";SAR=",sum(Vm),";SOR=",all_sor,";RVSB=",all_rvsb,";FISH=",FisherStrand_all,";ERR=",reg_res$coef["slope"],";SIG=",reg_res$coef["sigma"],";CONT=",paste(before,after,sep="x"),sep="",file=out_file,append=T)
  		    # FORMAT field
  		    cat("\t","GT:QVAL:DP:RO:AO:AF:SB:SOR:RVSB:FISH",sep = "",file=out_file,append=T)
          # all samples
  		    genotype=rep("0/0",l=nindiv)
  		    heterozygotes=which(reg_res$GQ>=GQ_threshold & sbs<=SB_threshold_SNV & reg_res$ma_count/reg_res$coverage < 0.75)
  		    genotype[heterozygotes]="0/1"
          homozygotes=which(reg_res$GQ>=GQ_threshold & sbs<=SB_threshold_SNV & reg_res$ma_count/reg_res$coverage >= 0.75)
  		    genotype[homozygotes]="1/1"
          for (cur_sample in 1:nindiv) {
            cat("\t",genotype[cur_sample],":",reg_res$GQ[cur_sample],":",DP[cur_sample],":",(Rp+Rm)[cur_sample],":",ma_count[cur_sample],":",(ma_count/DP)[cur_sample],":",Rp[cur_sample],",",Rm[cur_sample],",",Vp[cur_sample],",",Vm[cur_sample],":",sors[cur_sample],":",rvsbs[cur_sample],":",FisherStrand[cur_sample],sep = "",file=out_file,append=T)
          }
  		    cat("\n",sep = "",file=out_file,append=T)
          if (do_plots) {
            pdf(paste(pos_ref[i,"chr"],"_",pos_ref[i,"loc"],"_",pos_ref[i,"loc"],"_",pos_ref[i,"ref"],"_",alt,".pdf",sep=""),7,6)
            plot_rob_nb(reg_res, 10^-(GQ_threshold/10), plot_title=bquote(paste(.(pos_ref[i,"chr"]),":",.(pos_ref[i,"loc"])," (",.(pos_ref[i,"ref"]) %->% .(alt),")",sep="")), sbs=sbs, SB_threshold=SB_threshold_SNV,plot_labels=plot_labels,add_contours=add_contours,names=indiv_run[,2])
            dev.off()
          }
        }
      }
    }

    # DEL
    all_del=del[i,!is.na(del[i,])]
    if (length(all_del)>0) {
      all_del=as.data.frame(strsplit(unlist(strsplit(paste(all_del),split = "|",fixed=T)),split = ":",fixed=T),stringsAsFactors=F,)
      uniq_del=unique(toupper(as.character(all_del[2,])))
      coverage_all_del=sum(as.numeric(all_del[1,]))
      for (cur_del in uniq_del) {
        Vp=rep(0,nindiv)
        names(Vp)=indiv_run[,1]
        Vm=Vp
        all_cur_del=strsplit(grep(cur_del,del[i,],ignore.case=T,value=T),split="|",fixed=T)
        all_cur_del_p=lapply(all_cur_del,function(x) {grep(paste(":",cur_del,"$",sep=""),x,value=T)})
        all_cur_del_m=lapply(all_cur_del,function(x) {grep(paste(":",tolower(cur_del),"$",sep=""),x,value=T)})
        ma_p_cur_del=unlist(lapply(all_cur_del_p,function(x) {as.numeric(gsub(paste("([0-9]+):",cur_del,"$",sep=""),"\\1",x,ignore.case = T))}))
        ma_m_cur_del=unlist(lapply(all_cur_del_m,function(x) {as.numeric(gsub(paste("([0-9]+):",cur_del,"$",sep=""),"\\1",x,ignore.case = T))}))
        Vp[names(ma_p_cur_del)]=ma_p_cur_del
        Vm[names(ma_m_cur_del)]=ma_m_cur_del
        ma_count=Vp+Vm
        DP=coverage_matrix[i,]+ma_count
        reg_res=glmrob.nb(x=DP,y=ma_count,min_coverage=min_coverage,min_reads=min_reads)
        if (!is.na(reg_res$coef["slope"]) & sum(reg_res$GQ>=GQ_threshold,na.rm=TRUE)>0) {
          all_AO=sum(ma_count)
          all_DP=sum(coverage_matrix[i,])+sum(ma_count)
          common_annot()
          all_RO=sum(Rp+Rm)
          if (SB_type=="SOR") sbs=sors else { if (SB_type=="RVSB") sbs=rvsbs else {sbs=FisherStrand} }
          if (sum(sbs<=SB_threshold_indel & reg_res$GQ>=GQ_threshold,na.rm=TRUE)>0) {
            con=pipe(paste(samtools," faidx ",fasta_ref," ",pos_ref[i,"chr"],":",pos_ref[i,"loc"]+1-3,"-",pos_ref[i,"loc"]+1-1," | tail -n1",sep=""))
            before=readLines(con)
            close(con)
            con=pipe(paste(samtools," faidx ",fasta_ref," ",pos_ref[i,"chr"],":",pos_ref[i,"loc"]+1+nchar(cur_del),"-",pos_ref[i,"loc"]+1+3+nchar(cur_del)-1," | tail -n1",sep=""))
            after=readLines(con)
            close(con)
            prev_bp=substr(before,3,3)
            next_bp=substr(after,1,1)
            cat(pos_ref[i,"chr"],"\t",pos_ref[i,"loc"],"\t",".","\t",paste(prev_bp,cur_del,sep=""),"\t",prev_bp,"\t",max(reg_res$GQ),"\t",".",sep = "",file=out_file,append=T)
            # INFO field
            cat("\t","TYPE=del;NS=",sum(coverage_matrix[i,]>0),";AF=",sum(reg_res$GQ>=GQ_threshold)/sum(coverage_matrix[i,]>0),";DP=",all_DP,";RO=",all_RO,";AO=",all_AO,";SRF=",sum(Rp),";SRR=",sum(Rm),";SAF=",sum(Vp),";SAR=",sum(Vm),";SOR=",all_sor,";RVSB=",all_rvsb,";FISH=",FisherStrand_all,";ERR=",reg_res$coef["slope"],";SIG=",reg_res$coef["sigma"],";CONT=",paste(before,after,sep="x"),sep="",file=out_file,append=T)
            # FORMAT field
            cat("\t","GT:QVAL:DP:RO:AO:AF:SB:SOR:RVSB:FISH",sep = "",file=out_file,append=T)
            # all samples
            genotype=rep("0/0",l=nindiv)
            heterozygotes=which(reg_res$GQ>=GQ_threshold & sbs<=SB_threshold_indel & reg_res$ma_count/reg_res$coverage < 0.75)
    		    genotype[heterozygotes]="0/1"
            homozygotes=which(reg_res$GQ>=GQ_threshold & sbs<=SB_threshold_indel & reg_res$ma_count/reg_res$coverage >= 0.75)
    		    genotype[homozygotes]="1/1"
            for (cur_sample in 1:nindiv) {
              cat("\t",genotype[cur_sample],":",reg_res$GQ[cur_sample],":",DP[cur_sample],":",(Rp+Rm)[cur_sample],":",ma_count[cur_sample],":",(ma_count/DP)[cur_sample],":",Rp[cur_sample],",",Rm[cur_sample],",",Vp[cur_sample],",",Vm[cur_sample],":",sors[cur_sample],":",rvsbs[cur_sample],":",FisherStrand[cur_sample],sep = "",file=out_file,append=T)
            }
            cat("\n",sep = "",file=out_file,append=T)
            if (do_plots) {
              # deletions are shifted in samtools mpileup by 1bp, so put them at the right place by adding + to pos_ref[i,"loc"] everywhere in what follows
              pdf(paste(pos_ref[i,"chr"],"_",pos_ref[i,"loc"]+1,"_",pos_ref[i,"loc"]+1+nchar(cur_del)-1,"_",cur_del,"_","-",".pdf",sep=""),7,6)
              plot_rob_nb(reg_res, 10^-(GQ_threshold/10), plot_title=bquote(paste(.(pos_ref[i,"chr"]),":",.(pos_ref[i,"loc"]+1)," (",.(cur_del) %->% .("-"),")",sep="")),sbs=sbs, SB_threshold=SB_threshold_indel,plot_labels=plot_labels,add_contours=add_contours,names=indiv_run[,2])
              dev.off()
            }
          }
        }
     }
   }

    # INS
    all_ins=ins[i,!is.na(ins[i,])]
    if (length(all_ins)>0) {
      all_ins=as.data.frame(strsplit(unlist(strsplit(paste(all_ins),split = "|",fixed=T)),split = ":",fixed=T),stringsAsFactors=F,)
      uniq_ins=unique(toupper(as.character(all_ins[2,])))
      coverage_all_ins=sum(as.numeric(all_ins[1,]))
      for (cur_ins in uniq_ins) {
        Vp=rep(0,nindiv)
        names(Vp)=indiv_run[,1]
        Vm=Vp
        all_cur_ins=strsplit(grep(cur_ins,ins[i,],ignore.case=T,value=T),split="|",fixed=T)
        all_cur_ins_p=lapply(all_cur_ins,function(x) {grep(paste(":",cur_ins,"$",sep=""),x,value=T)})
        all_cur_ins_m=lapply(all_cur_ins,function(x) {grep(paste(":",tolower(cur_ins),"$",sep=""),x,value=T)})
        ma_p_cur_ins=unlist(lapply(all_cur_ins_p,function(x) {as.numeric(gsub(paste("([0-9]+):",cur_ins,"$",sep=""),"\\1",x,ignore.case = T))}))
        ma_m_cur_ins=unlist(lapply(all_cur_ins_m,function(x) {as.numeric(gsub(paste("([0-9]+):",cur_ins,"$",sep=""),"\\1",x,ignore.case = T))}))
        Vp[names(ma_p_cur_ins)]=ma_p_cur_ins
        Vm[names(ma_m_cur_ins)]=ma_m_cur_ins
        ma_count=Vp+Vm
        DP=coverage_matrix[i,]+ma_count[]
        reg_res=glmrob.nb(x=DP,y=ma_count,min_coverage=min_coverage,min_reads=min_reads)
        if (!is.na(reg_res$coef["slope"]) & sum(reg_res$GQ>=GQ_threshold,na.rm=TRUE)>0) {
          all_AO=sum(ma_count)
          all_DP=sum(coverage_matrix[i,])+sum(ma_count)
          common_annot()
          all_RO=sum(Rp+Rm)
          if (SB_type=="SOR") sbs=sors else { if (SB_type=="RVSB") sbs=rvsbs else {sbs=FisherStrand} }
          if (sum(sbs<=SB_threshold_indel & reg_res$GQ>=GQ_threshold,na.rm=TRUE)>0) {
            con=pipe(paste(samtools," faidx ",fasta_ref," ",pos_ref[i,"chr"],":",pos_ref[i,"loc"]-2,"-",pos_ref[i,"loc"]," | tail -n1",sep=""))
            before=readLines(con)
            close(con)
            con=pipe(paste(samtools," faidx ",fasta_ref," ",pos_ref[i,"chr"],":",pos_ref[i,"loc"]+1,"-",pos_ref[i,"loc"]+3," | tail -n1",sep=""))
            after=readLines(con)
            close(con)
            prev_bp=substr(before,3,3)
            cat(pos_ref[i,"chr"],"\t",pos_ref[i,"loc"],"\t",".","\t",prev_bp,"\t",paste(prev_bp,cur_ins,sep=""),"\t",max(reg_res$GQ),"\t",".",sep = "",file=out_file,append=T)
            # INFO field
            cat("\t","TYPE=ins;NS=",sum(coverage_matrix[i,]>0),";AF=",sum(reg_res$GQ>=GQ_threshold)/sum(coverage_matrix[i,]>0),";DP=",all_DP,";RO=",all_RO,";AO=",all_AO,";SRF=",sum(Rp),";SRR=",sum(Rm),";SAF=",sum(Vp),";SAR=",sum(Vm),";SOR=",all_sor,";RVSB=",all_rvsb,";FISH=",FisherStrand_all,";ERR=",reg_res$coef["slope"],";SIG=",reg_res$coef["sigma"],";CONT=",paste(before,after,sep="x"),sep="",file=out_file,append=T)
            # FORMAT field
            cat("\t","GT:QVAL:DP:RO:AO:AF:SB:SOR:RVSB:FISH",sep = "",file=out_file,append=T)
            # all samples
            genotype=rep("0/0",l=nindiv)
            heterozygotes=which(reg_res$GQ>=GQ_threshold & sbs<=SB_threshold_indel & reg_res$ma_count/reg_res$coverage < 0.75)
    		    genotype[heterozygotes]="0/1"
            homozygotes=which(reg_res$GQ>=GQ_threshold & sbs<=SB_threshold_indel & reg_res$ma_count/reg_res$coverage >= 0.75)
    		    genotype[homozygotes]="1/1"
            for (cur_sample in 1:nindiv) {
              cat("\t",genotype[cur_sample],":",reg_res$GQ[cur_sample],":",DP[cur_sample],":",(Rp+Rm)[cur_sample],":",ma_count[cur_sample],":",(ma_count/DP)[cur_sample],":",Rp[cur_sample],",",Rm[cur_sample],",",Vp[cur_sample],",",Vm[cur_sample],":",sors[cur_sample],":",rvsbs[cur_sample],":",FisherStrand[cur_sample],sep = "",file=out_file,append=T)
            }
            cat("\n",sep = "",file=out_file,append=T)
            if (do_plots) {
              pdf(paste(pos_ref[i,"chr"],"_",pos_ref[i,"loc"],"_",pos_ref[i,"loc"],"_","-","_",cur_ins,".pdf",sep=""),7,6)
              plot_rob_nb(reg_res, 10^-(GQ_threshold/10), plot_title=bquote(paste(.(pos_ref[i,"chr"]),":",.(pos_ref[i,"loc"])," (",.("-") %->% .(cur_ins),")",sep="")),sbs=sbs, SB_threshold=SB_threshold_indel,plot_labels=plot_labels,add_contours=add_contours,names=indiv_run[,2])
              dev.off()
            }
          }
        }
      }
    }
  }
}
