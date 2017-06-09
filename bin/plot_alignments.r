#! /usr/bin/env Rscript

# needlestack: a multi-sample somatic variant caller
# Copyright (C) 2017 Matthieu Foll, Tiffany Delhomme, Nicolas Alcala and Aurelie Gabriel
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

#loading libraries
if(do_alignments==TRUE){
  
  library("Gviz")
  #load libraries for the annotations
  library(paste0("TxDb.",ref_genome,".knownGene"),character.only=TRUE)
  assign("txdb",get(paste0("TxDb.",ref_genome,".knownGene")))
  
  annotation=paste0("org.",substr(ref_genome,1,2),".eg.db")
  library(annotation,character.only=TRUE)
  assign("annotation",get(annotation))
  
  options(ucscChromosomeNames=FALSE) #in case some chromosome names are not recognized as UCSC chromosome names
  
  
  #read ref genome
  library(Biostrings)
  ref=readDNAStringSet(fasta_ref)
  sTrack<-SequenceTrack(ref,cex=0.6)
}

#get UCSC chromosome correspondances
get_UCSC_correpondance<- function(chr,ref_genome){
  convert_table=read.table(paste(source_path,paste0(unlist(strsplit(ref_genome,"[.]"))[3],"_chromosomeNames2UCSC.txt"),sep=""),head=FALSE,fill = TRUE)
  if(length(which(convert_table$V1==chr))==0){
    return("")
  }else{
    return(as.character(unique(convert_table[which(convert_table$V1==chr),2])))
  }
}

#create alignments tracks
get_htTrack <- function(isTNpairs,i,Tindex,Nindex,onlyTindex,onlyNindex,indiv_run,bam_folder,plot_grtracks,UCSC,sTrack,grtrack,chr,pos,paired=FALSE){
  #get bam alignments tracks (AlignmentsTrack())
  if(isTNpairs){
    if(i %in% Tindex){
      pair_index=which(Tindex==i)
      alTrack_T <- AlignmentsTrack(paste0(bam_folder,indiv_run[,1][Tindex[pair_index]],".bam"), isPaired = paired,stacking = "squish",alpha=0.95,chromosome=chr,cex.mismatch=0.5,name="Tumor",cex.title=1.5)
      alTrack_N <- AlignmentsTrack(paste0(bam_folder,indiv_run[,1][Nindex[pair_index]],".bam"), isPaired = paired,stacking = "squish",alpha=0.95,chromosome=chr,cex.mismatch=0.5,name="Normal",cex.title=1.5)
      unique=FALSE #two alignments on the same plot
    }else if(i %in% Nindex) {
      pair_index=which(Nindex==i)
      alTrack_T <- AlignmentsTrack(paste0(bam_folder,indiv_run[,1][Tindex[pair_index]],".bam"), isPaired = paired,stacking = "squish",alpha=0.95,chromosome=chr,cex.mismatch=0.5,name="Tumor",cex.title=1.5)
      alTrack_N <- AlignmentsTrack(paste0(bam_folder,indiv_run[,1][Nindex[pair_index]],".bam"), isPaired = paired,stacking = "squish",alpha=0.95,chromosome=chr,cex.mismatch=0.5,name="Normal",cex.title=1.5)
      unique=FALSE
    }else if((i %in% onlyTindex)){
      alTrack_unique <- AlignmentsTrack(paste0(bam_folder,indiv_run[,1][i],".bam"), isPaired = paired,stacking = "squish",alpha=0.95,chromosome=chr,cex.mismatch=0.5,name="Normal",cex.title=1.5)
      unique=TRUE #only one alignment 
    }else if((i %in% onlyNindex)){
      alTrack_unique <- AlignmentsTrack(paste0(bam_folder,indiv_run[,1][i],".bam"), isPaired = paired,stacking = "squish",alpha=0.95,chromosome=chr,cex.mismatch=0.5,name="Tumor",cex.title=1.5)
      unique=TRUE
    }
    
    if(plot_grtracks){
      if(unique){
        ht <- HighlightTrack(trackList = list(alTrack_unique, sTrack, grtrack), start = c(pos), width =0,chromosome = chr) #highlight the position of the variant
        s=c(0.05,0.1,0.72,0.05,0.08) #tracks proportions
      }else{
        ht <- HighlightTrack(trackList = list(alTrack_T,alTrack_N, sTrack, grtrack), start = c(pos), width =0,chromosome = chr)
        s=c(0.05,0.1,0.36,0.36,0.05,0.08)
      }
    }else if(plot_grtracks==FALSE){
      if(unique){
        ht <- HighlightTrack(trackList = list(alTrack_unique, sTrack), start = c(pos), width =0,chromosome = chr)
        s=c(0.05,0.1,0.8,0.05)
      }else{
        ht <- HighlightTrack(trackList = list(alTrack_T,alTrack_N, sTrack), start = c(pos), width =0,chromosome = chr)
        s=c(0.05,0.1,0.4,0.4,0.05)
      }
    }
  }else{
  alTrack <- AlignmentsTrack(paste0(bam_folder,indiv_run[,1][i],".bam"), isPaired = paired,stacking = "squish",alpha=0.95,chromosome=chr,cex.mismatch=0.5,name="Alignment",cex.title=1.5)
  if(plot_grtracks){
      ht <- HighlightTrack(trackList = list(alTrack, sTrack, grtrack), start = c(pos), width =0,chromosome = chr)
      s=c(0.05,0.1,0.72,0.05,0.08)
    }else{
      ht <- HighlightTrack(trackList = list(alTrack, sTrack), start = c(pos), width =0,chromosome = chr)
      s=c(0.05,0.1,0.8,0.05)
    }
  }
  
  return(list(ht,s))
}

#remove repetition : to avoid having two identical figures for each tumor/normal pair
check_pair_repetition <- function(s,Tindex,Nindex){
  stock_pairs_index=vector()
  i=1
  while(i <= length(s)){
    if(s[i] %in% Tindex){
      pair_index=which(Tindex==s[i])
      if ( (pair_index %in% stock_pairs_index) ){
        s=s[-i]
        i=i-1
      }else{
        stock_pairs_index = append(stock_pairs_index,pair_index)
      }
    }else if(s[i] %in% Nindex) {
      pair_index=which(Nindex==s[i])
      if ( (pair_index %in% stock_pairs_index) ){
        s=s[-i]
        i=i-1
      }else{
        stock_pairs_index = append(stock_pairs_index,pair_index)
      }
    }
    i=i+1
  }
  return(s)
}

#plot alignments
plotGviz <- function(isTNpairs,sTrack,ref_genome,txdb,annotation,UCSC,indiv_run,linepos,genotype,somatic_status,do_plots,Tindex,Nindex,onlyTindex,onlyNindex,bam_folder,w=50,w_zoomout=1000,paired=FALSE,nb_toplot=5){
  chr=linepos[1]
  pos=as.numeric(linepos[2])
  #select samples with the variant
  if(do_plots=="SOMATIC"){ #select only the sample with a somatic variant
    samples_with_var=which(somatic_status=="SOMATIC")
    samples_with_var = check_pair_repetition(samples_with_var,Tindex,Nindex) #remove repetitions
    samples_without_var=vector()
  }else{
    samples_with_var=which((genotype=="0/1")|(genotype=="1/1"))
    if(isTNpairs){samples_with_var = check_pair_repetition(samples_with_var,Tindex,Nindex)} #remove repetitions
    samples_without_var=which((genotype=="./.")|(genotype=="0/0"))
    if(isTNpairs){samples_without_var = check_pair_repetition(samples_without_var,Tindex,Nindex)} #remove repetitions
  }
  
  chr_annotation=get_UCSC_correpondance(chr,ref_genome)
  if(chr_annotation==""){
    only_alignment=TRUE
  }else{
    only_alignment=FALSE
  }
  
  if(only_alignment==FALSE){
    #Define tracks common to all samples (reference sequence, chromosome representation, genome annotation)
    gtrack <- GenomeAxisTrack() #genomic axis
    sTrack@sequence@ranges@NAMES <- sapply( 1:length(sTrack@sequence@ranges@NAMES), function(i) unlist(strsplit(sTrack@sequence@ranges@NAMES[i], " "))[1] )
    sTrack@chromosome <-chr
    ideoTrack <- IdeogramTrack(genome = unlist(strsplit(ref_genome,".",fixed=TRUE))[3], chromosome = chr_annotation) #chromosome representation
    levels(ideoTrack@bandTable$chrom)[which(levels(ideoTrack@bandTable$chrom)==chr_annotation)] <- chr
    ideoTrack@chromosome<-chr
    
    grtrack <- GeneRegionTrack(txdb,chromosome = chr_annotation,start = pos-w_zoomout, end = pos+w_zoomout,exonAnnotation = "exon",shape = "arrow",showTitle=FALSE,alpha=0.95,cex=0.7)
    if( summary(is.na(grtrack@range@elementMetadata$gene))[2]>0 ){ #if no genes in the annotation it is not possible to use the collapseTranscripts option
      displayPars(grtrack)$collapseTranscripts <- "longest"
    }
    grtrack@range@seqinfo@seqnames<-chr
    levels(grtrack@range@seqnames)<-chr
    grtrack@chromosome<-chr
    
    #complete gene annotation 
    if(length(unique(gene(grtrack)))>=6){
      sampling=sample(unique(grtrack@range@elementMetadata$gene),5)
      grtrack@range=grtrack@range[grtrack@range@elementMetadata$gene %in% sampling,]
      displayPars(grtrack)$cex <- 0.6
    }
    displayPars(grtrack) <- list(background.title = "white")
    #zoom out on the gene
    grtrack_zoomout=grtrack
    displayPars(grtrack_zoomout)$transcriptAnnotation <- "symbol"
    displayPars(grtrack_zoomout)$showExonId <- FALSE
    displayPars(grtrack_zoomout)$shape <- c("smallArrow","box")
    if(length(gene(grtrack_zoomout))!=0){
      if( length( which( unique(gene(grtrack_zoomout)) %in% keys(annotation,keytype="ENTREZID") == TRUE))==length(unique(gene(grtrack_zoomout))) ){
        symbols <- unlist(mapIds(annotation, gene(grtrack_zoomout), "SYMBOL", "ENTREZID", multiVals = "first"))
        symbol(grtrack_zoomout) <- symbols[gene(grtrack_zoomout)]
      }
      plot_grtracks=TRUE
    }else{ #genome annotation can not be added if non UCSC genome
      plot_grtracks=FALSE
    }
    if(length(unique(gene(grtrack_zoomout)))>=6){
      grtrack_zoomout@range=grtrack_zoomout@range[grtrack_zoomout@range@elementMetadata$gene %in% sampling,]
    }
    ht_zoomout <- HighlightTrack(trackList = list(grtrack_zoomout,gtrack), start = c(pos), width =0,chromosome = chr) #highlight the position of the variant on the genomic axis and the genome annotation
    
  }else{
    plot_grtracks=FALSE
    sTrack@sequence@ranges@NAMES <- sapply( 1:length(sTrack@sequence@ranges@NAMES), function(i) unlist(strsplit(sTrack@sequence@ranges@NAMES[i], " "))[1] )
    sTrack@chromosome <-chr
  }
  
  #plot alignments for the samples with the variant
  for(i in samples_with_var){
    
    if(only_alignment=="FALSE"){
      res=get_htTrack(isTNpairs,i,Tindex,Nindex,onlyTindex,onlyNindex,indiv_run,bam_folder,plot_grtracks,UCSC,sTrack,grtrack,chr,pos,paired) #create alignments tracks
      if(plot_grtracks){
        tryCatch(
          {
            grid.newpage()
            pushViewport(viewport(x=0,y=1, height=0.85, width=1, just=c("left","top"))) 
            plotTracks(c(ideoTrack,gtrack,res[1]),sizes=unlist(res[2]),from = pos-w, to = pos+w,add = TRUE, add53=TRUE,min.height=4, main=paste0(indiv_run[i,2]),title.width=0.7,littleTicks = TRUE,cex.main=1.5)
            popViewport(1)
            pushViewport(viewport(x=0,y=0, height=0.15, width=1, just=c("left","bottom")))
            plotTracks(list(ht_zoomout),chromosome = chr, add = TRUE)
            popViewport(0)
          },
          error=function(e) print(paste0("Warning: error for this AlignmentsPlot : chr -> ",chr," pos -> ",pos," sample -> ",indiv_run[i,2]))
        )
      }else{
        tryCatch(
          {
            plotTracks(c(ideoTrack,gtrack,res[1]),sizes=unlist(res[2]),from = pos-w, to = pos+w, add53=TRUE,min.height=4, main=paste0(indiv_run[i,2]),title.width=0.7,littleTicks = TRUE,cex.main=1.5)
          },
          error=function(e) print(paste0("Warning: error for this AlignmentsPlot : chr -> ",chr," pos -> ",pos," sample -> ",indiv_run[i,2]))
        )
      }
      
    }else {
      res=get_htTrack(isTNpairs,i,Tindex,Nindex,onlyTindex,onlyNindex,indiv_run,bam_folder,plot_grtracks,UCSC,sTrack,grtrack,chr,pos,paired) #create alignments tracks
      tryCatch(
        {
          plotTracks(c(gtrack,res[1]),sizes=unlist(res[2])[-1],from = pos-w, to = pos+w, add53=TRUE,min.height=4, main=paste0(indiv_run[i,2]),title.width=0.7,littleTicks = TRUE,cex.main=1.5)
        },
        error=function(e) print(paste0("Warning: error for this AlignmentsPlot : chr -> ",chr," pos -> ",pos," sample -> ",indiv_run[i,2]))
      )
    }
  }
  
  #plot alignments for the samples without the variant
  if(length(samples_without_var)!=0){
    #sampling of the samples without the variant
    if(length(samples_without_var)<nb_toplot){
      samples_without_var_toplot=samples_without_var
    }else{
      samples_without_var_toplot=sample(samples_without_var,nb_toplot)
    }
    
    for(i in samples_without_var_toplot){
      if(plot_grtracks){
        res=get_htTrack(isTNpairs,i,Tindex,Nindex,onlyTindex,onlyNindex,indiv_run,bam_folder,plot_grtracks,UCSC,sTrack,grtrack,chr,pos,paired) #create alignments tracks
        if(plot_grtracks){
          tryCatch(
            {
              grid.newpage()
              pushViewport(viewport(x=0,y=1, height=0.85, width=1, just=c("left","top"))) 
              plotTracks(c(ideoTrack,gtrack,res[1]),sizes=unlist(res[2]),from = pos-w, to = pos+w,add = TRUE, add53=TRUE,min.height=4, main=paste0(indiv_run[i,2],"*"),title.width=0.7,littleTicks = TRUE,cex.main=1.5)
              popViewport(1)
              pushViewport(viewport(x=0,y=0, height=0.15, width=1, just=c("left","bottom")))
              plotTracks(list(ht_zoomout),chromosome = chr, add = TRUE)
              popViewport(0)
            },
            error=function(e) print(paste0("Warning: error for this AlignmentsPlot : chr -> ",chr," pos -> ",pos," sample* -> ",indiv_run[i,2]))
          )
        }else{
          tryCatch(
            {
              plotTracks(c(ideoTrack,gtrack,res[1]),sizes=unlist(res[2]),from = pos-w, to = pos+w,add53=TRUE,min.height=4, main=paste0(indiv_run[i,2],"*"),title.width=0.7,littleTicks = TRUE,cex.main=1.5)
            },
            error=function(e) print(paste0("Warning: error for this AlignmentsPlot : chr -> ",chr," pos -> ",pos," sample* -> ",indiv_run[i,2]))
          )
        }
      }else if(only_alignment){
        res=get_htTrack(isTNpairs,i,Tindex,Nindex,onlyTindex,onlyNindex,indiv_run,bam_folder,plot_grtracks,UCSC,sTrack,grtrack,chr,pos,paired) #create alignments tracks
        tryCatch(
          {
            plotTracks(c(gtrack,res[1]),sizes=unlist(res[2])[-1],from = pos-w, to = pos+w, add53=TRUE,min.height=4, main=paste0(indiv_run[i,2],"*"),title.width=0.7,littleTicks = TRUE,cex.main=1.5)
          },
          error=function(e) print(paste0("Warning: error for this AlignmentsPlot : chr -> ",chr," pos -> ",pos," sample* -> ",indiv_run[i,2]))
        )
      }
    }
  }
  
}


