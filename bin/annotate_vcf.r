#! /usr/bin/env Rscript

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

args <- commandArgs(TRUE)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsL <- as.list(as.character(as.data.frame(do.call("rbind", parseArgs(args)))$V2))
names(argsL) <- as.data.frame(do.call("rbind", parseArgs(args)))$V1
args <- argsL;rm(argsL)

if(is.null(args$input_vcf))            {stop("no input VCF file")} else {input_vcf = args$input_vcf}
if(is.null(args$out_vcf))              {args$out_vcf = out_vcf = paste(gsub(".vcf.bgz","",input_vcf),"_annotated_needlestack.vcf",sep="")} else {out_vcf=args$out_vcf}
if(nchar(args$out_vcf)==0)             {out_vcf = paste(gsub(".vcf.bgz","",input_vcf),"_annotated_needlestack.vcf",sep="")}
if(is.null(args$chunk_size))           {chunk_size = 1000} else {chunk_size = as.numeric(args$chunk_size)}
if(is.null(args$do_plots))             {do_plots = TRUE} else {do_plots = as.logical(args$do_plots)}
if(is.null(args$plot_labels))          {plot_labels = TRUE} else {plot_labels = as.logical(args$plot_labels)}
if(is.null(args$add_contours))         {add_contours = TRUE} else {add_contours = as.logical(args$add_contours)}
if(is.null(args$min_coverage))         {min_coverage = 50} else {min_coverage = as.numeric(args$min_coverage)}
if(is.null(args$min_reads))            {min_reads = 5} else {min_reads = as.numeric(args$min_reads)}
if(is.null(args$GQ_threshold))         {GQ_threshold=50} else {GQ_threshold = as.numeric(args$GQ_threshold)}
if(is.null(args$SB_threshold))         {SB_threshold=100} else {SB_threshold = as.numeric(args$SB_threshold)}
if(is.null(args$extra_rob))            {extra_rob=FALSE} else {extra_rob=as.logical(args$extra_rob)}

source(paste(args$source_path,"glm_rob_nb.r",sep=""))
source(paste(args$source_path,"plot_rob_nb.r",sep=""))
library(VariantAnnotation)

#initiate the first chunk
vcf <- open(VcfFile(input_vcf,  yieldSize=chunk_size))
vcf_chunk = readVcf(vcf, "hg19")

#and continue
while(dim(vcf_chunk)[1] != 0) {
  # coverage (matrix of integers)
  DP_matrix = geno(vcf_chunk,"DP")
  # AO counts (matrix of lists of integers)
  AD_matrix = geno(vcf_chunk,"AD")

  #compute regressions and qvals,err,sig
  reg_list = lapply(1:dim(vcf_chunk)[1], function(var_line) { #for each line of the chunk return a list of reg for each AD
    print(var_line)
    print(start(ranges(rowRanges(vcf_chunk,"seqnames"))[var_line]))
    # replace NAs and integer(0) by correct number of 0 ADs
    AD_matrix[var_line, which(is.na(AD_matrix[var_line,]))] = lapply(AD_matrix[var_line, which(is.na(AD_matrix[var_line,]))], function(x) x=as.vector(rep(0,max(lengths(AD_matrix[var_line,])))))
    AD_matrix[var_line,] = lapply(AD_matrix[var_line,], function(x) {if(length(x)==0) { x=as.vector(rep(0, ifelse(max(lengths(AD_matrix[var_line,]),na.rm = T) >0, max(lengths(AD_matrix[var_line,]),na.rm = T), 2) )) } else {x=x} } )
    lapply(2:max(lengths(AD_matrix[var_line,])), function(AD_index) { #for each alternative
      DP=DP_matrix[var_line,]
      AO=unlist(lapply(AD_matrix[var_line,],"[[",AD_index)) #AD_matrix[var_line,] is a list of AD for each sample, here return list of ADs(i) for alt i
      reg_res=glmrob.nb(x=DP,y=AO,min_coverage=min_coverage,min_reads=min_reads,extra_rob=extra_rob)
      if (do_plots) {
        chr=as.character(seqnames(rowRanges(vcf_chunk,"seqnames"))[var_line])
        loc=start(ranges(rowRanges(vcf_chunk,"seqnames"))[var_line])
        ref=as.character(ref(vcf_chunk)[[var_line]])
        alt=paste(as.character(alt(vcf_chunk)[[var_line]]),collapse = ",")
        sbs=rep(NA,dim(vcf_chunk)[2])
        pdf(paste(chr,"_",loc,"_",loc,"_",ref,"_",alt,"_",AD_index-1,ifelse(reg_res$extra_rob,"_extra_robust",""),".pdf",sep=""),7,6)
        plot_rob_nb(reg_res, 10^-(GQ_threshold/10), plot_title=bquote(paste(.(chr),":",.(loc)," (",.(ref) %->% .(alt),"[",.(AD_index-1),"]",")",.(ifelse(reg_res$extra_rob," EXTRA ROBUST","")),sep="")), sbs=sbs, SB_threshold=SB_threshold,plot_labels=T,add_contours=T,names=samples(header(vcf_chunk)))
        dev.off()
      }
      reg_res
    })
  })
  qvals = lapply(reg_list, function(regs) {
    lapply(regs, function(reg) (unlist(reg["GQ"])+0)) #here add +0 to avoid sprintf wrinting -0
  })
  err = lapply(reg_list, function(regs) {
    lapply(regs, function(reg) unlist(reg$coef["slope"]))
  })
  sig = lapply(reg_list, function(regs) {
    lapply(regs, function(reg) unlist(reg$coef["sigma"]))
  })
  extra_robust_gl = unlist(lapply(reg_list, function(regs) {
    lapply(regs, function(reg) unlist(reg$extra_rob))
  }))

  #annotate the header of the chunk
  info(header(vcf_chunk))["ERR",]=list("A","Float","Error rate estimated by needlestack")
  info(header(vcf_chunk))["SIG",]=list("A","Float","Dispertion parameter estimated by needlestack")
  info(header(vcf_chunk))["WARN",]=list("A","Character","Warning message when position is processed specifically by needlestack")
  geno(header(vcf_chunk))["QVAL",]=list("A","Float","Phred q-values computed by needlestack")

  #annotate the chunk with computed values
  info(vcf_chunk)$ERR = NumericList(err)
  info(vcf_chunk)$SIG = NumericList(sig)
  geno(vcf_chunk)$QVAL = matrix(data = unlist(lapply(qvals, function(q) as.list(data.frame(t(mapply(c,q))))),recursive = FALSE),
                                nrow = dim(vcf_chunk)[1],
                                byrow = TRUE)
  info(vcf_chunk)$WARN = rep(NA, length(extra_robust_gl))
  info(vcf_chunk)$WARN[which(extra_robust_gl==TRUE)]="EXTRA_ROBUST_GL"

  #write out the annotated VCF
  con = file(out_vcf, open = "a")
  writeVcf(vcf_chunk, con)
  vcf_chunk = readVcf(vcf, "hg19")
  close(con)
}
