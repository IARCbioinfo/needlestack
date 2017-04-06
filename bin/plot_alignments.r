if(do_alignments==TRUE){
  
  library("Gviz")
  ref_string=paste0("BSgenome.",ref_genome)
  library(ref_string,character.only=TRUE)
  
  assign("g",get(ref_string))
  
  if(ref_genome!="Hsapiens.1000genomes.hs37d5"){
    
    #define SequenceTrack (reference genome)
    library(paste0("TxDb.",ref_genome,".knownGene"),character.only=TRUE)
    assign("txdb",get(paste0("TxDb.",ref_genome,".knownGene")))
    sTrack <- SequenceTrack(g,cex = 0.6)
    
    annotation=paste0("org.",substr(ref_genome,1,2),".eg.db")
    library(annotation,character.only=TRUE)
    assign("annotation",get(annotation))
    UCSC=TRUE #UCSC reference genome
    
  }else if(ref_genome=="Hsapiens.1000genomes.hs37d5"){
    UCSC=FALSE #non UCSC reference genome
    options(ucscChromosomeNames=FALSE)
    
    #define SequenceTrack (reference genome)
    sTrack <- SequenceTrack(g,cex = 0.6)
    
  }
}

get_htTrack <- function(isTNpairs,i,Tindex,Nindex,onlyTindex,onlyNindex,indiv_run,bam_folder,plot_grtracks,UCSC,sTrack,grtrack,chr,pos,paired=TRUE){
  if(isTNpairs){
    if(i %in% Tindex){
      pair_index=which(Tindex==i)
      alTrack_T <- AlignmentsTrack(paste0(bam_folder,indiv_run[,1][Tindex[pair_index]],".bam"), isPaired = paired,stacking = "squish",alpha=0.95,chromosome=chr,cex.mismatch=0.5,name="Tumor",cex.title=1.5)
      alTrack_N <- AlignmentsTrack(paste0(bam_folder,indiv_run[,1][Nindex[pair_index]],".bam"), isPaired = paired,stacking = "squish",alpha=0.95,chromosome=chr,cex.mismatch=0.5,name="Normal",cex.title=1.5)
      unique=FALSE
    }else if(i %in% Nindex) {
      pair_index=which(Nindex==i)
      alTrack_T <- AlignmentsTrack(paste0(bam_folder,indiv_run[,1][Tindex[pair_index]],".bam"), isPaired = paired,stacking = "squish",alpha=0.95,chromosome=chr,cex.mismatch=0.5,name="Tumor",cex.title=1.5)
      alTrack_N <- AlignmentsTrack(paste0(bam_folder,indiv_run[,1][Nindex[pair_index]],".bam"), isPaired = paired,stacking = "squish",alpha=0.95,chromosome=chr,cex.mismatch=0.5,name="Normal",cex.title=1.5)
      unique=FALSE
    }else if((i %in% onlyTindex)){
      alTrack_unique <- AlignmentsTrack(paste0(bam_folder,indiv_run[,1][i],".bam"), isPaired = paired,stacking = "squish",alpha=0.95,chromosome=chr,cex.mismatch=0.5,name="Normal",cex.title=1.5)
      unique=TRUE
    }else if((i %in% onlyNindex)){
      alTrack_unique <- AlignmentsTrack(paste0(bam_folder,indiv_run[,1][i],".bam"), isPaired = paired,stacking = "squish",alpha=0.95,chromosome=chr,cex.mismatch=0.5,name="Tumor",cex.title=1.5)
      unique=TRUE
    }
    
    
    if(UCSC & plot_grtracks){
      if(unique){
        ht <- HighlightTrack(trackList = list(alTrack_unique, sTrack, grtrack), start = c(pos), width =0,chromosome = chr)
        s=c(0.05,0.1,0.72,0.05,0.08)
      }else{
        ht <- HighlightTrack(trackList = list(alTrack_T,alTrack_N, sTrack, grtrack), start = c(pos), width =0,chromosome = chr)
        s=c(0.05,0.1,0.36,0.36,0.05,0.08)
      }
    }else if(UCSC==FALSE | plot_grtracks==FALSE){
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
    if(UCSC & plot_grtracks){
      ht <- HighlightTrack(trackList = list(alTrack, sTrack, grtrack), start = c(pos), width =0,chromosome = chr)
      s=c(0.05,0.1,0.72,0.05,0.08)
    }else if(UCSC==FALSE | plot_grtracks==FALSE){
      ht <- HighlightTrack(trackList = list(alTrack, sTrack), start = c(pos), width =0,chromosome = chr)
      s=c(0.05,0.1,0.8,0.05)
    }
  }
  
  return(list(ht,s))
}

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

plotGviz <- function(isTNpairs,sTrack,ref_genome,txdb,annotation,UCSC,indiv_run,linepos,genotype,somatic_status,do_plots,Tindex,Nindex,onlyTindex,onlyNindex,bam_folder,w=50,w_zoomout=1000,paired=TRUE,nb_toplot=5){
  chr=linepos[1]
  pos=as.numeric(linepos[2])
  #select samples with the variant
  if(do_plots=="SOMATIC"){
    samples_with_var=which(somatic_status=="SOMATIC")
    samples_with_var = check_pair_repetition(samples_with_var,Tindex,Nindex)
    samples_without_var=vector()
  }else{
    samples_with_var=which((genotype=="0/1")|(genotype=="1/1"))
    if(isTNpairs){samples_with_var = check_pair_repetition(samples_with_var,Tindex,Nindex)}
    samples_without_var=which((genotype=="./.")|(genotype=="0/0"))
    if(isTNpairs){samples_without_var = check_pair_repetition(samples_without_var,Tindex,Nindex)}
  }
  
  
  gtrack <- GenomeAxisTrack()
  if(UCSC){
    sTrack@chromosome <- chr
    ideoTrack <- IdeogramTrack(genome = unlist(strsplit(ref_genome,".",fixed=TRUE))[3], chromosome = chr)
    grtrack <- GeneRegionTrack(txdb,chromosome = chr,start = pos-w, end = pos-w,exonAnnotation = "exon",collapseTranscripts = "longest",shape = "arrow",showTitle=FALSE,alpha=0.95)
    displayPars(grtrack) <- list(background.title = "white")
    grtrack_zoomout <- GeneRegionTrack(txdb,chromosome = chr,start = pos-w_zoomout, end = pos+w_zoomout,transcriptAnnotation = "symbol",collapseTranscripts = "longest",alpha=0.95,showTitle=FALSE)
    if(length(gene(grtrack_zoomout))!=0){
      symbols <- unlist(mapIds(annotation, gene(grtrack_zoomout), "SYMBOL", "ENTREZID", multiVals = "first"))
      symbol(grtrack_zoomout) <- symbols[gene(grtrack_zoomout)]
      plot_grtracks=TRUE
    }else{
      plot_grtracks=FALSE
    }
    ht_zoomout <- HighlightTrack(trackList = list(grtrack_zoomout,gtrack), start = c(pos), width =0,chromosome = chr)
  }else{
    sTrack@chromosome <- chr
    ideoTrack <- IdeogramTrack(genome = "hg19", chromosome = paste0("chr",chr))
    levels(ideoTrack@bandTable$chrom) <- sub("^chr", "", levels(ideoTrack@bandTable$chrom), ignore.case=T)
    ideoTrack@chromosome<-chr
    plot_grtracks=FALSE
  }
  
  
  for(i in samples_with_var){
    
    if(UCSC & plot_grtracks){
      res=get_htTrack(isTNpairs,i,Tindex,Nindex,onlyTindex,onlyNindex,indiv_run,bam_folder,plot_grtracks,UCSC,sTrack,grtrack,chr,pos,paired)
      grid.newpage()
      pushViewport(viewport(x=0,y=1, height=0.85, width=1, just=c("left","top"))) 
      plotTracks(c(ideoTrack,gtrack,res[1]),sizes=unlist(res[2]),from = pos-w, to = pos+w,add = TRUE, add53=TRUE,min.height=4, main=paste0(indiv_run[i,2]),title.width=0.7,littleTicks = TRUE,cex.main=1.5)
      popViewport(1)
      pushViewport(viewport(x=0,y=0, height=0.15, width=1, just=c("left","bottom")))
      plotTracks(list(ht_zoomout),chromosome = chr, add = TRUE)
      popViewport(0)
    }else if(UCSC==FALSE | plot_grtracks==FALSE){
      res=get_htTrack(isTNpairs,i,Tindex,Nindex,onlyTindex,onlyNindex,indiv_run,bam_folder,plot_grtracks,UCSC,sTrack,grtrack,chr,pos,paired)
      plotTracks(c(ideoTrack,gtrack,res[1]),sizes=unlist(res[2]),from = pos-w, to = pos+w, add53=TRUE,min.height=4, main=paste0(indiv_run[i,2]),title.width=0.7,littleTicks = TRUE,cex.main=1.5)
      
    }
  }
  
  if(length(samples_without_var)!=0){
    if(length(samples_without_var)<nb_toplot){
      samples_without_var_toplot=samples_without_var
    }else{
      samples_without_var_toplot=sample(samples_without_var,nb_toplot)
    }
    
    for(i in samples_without_var_toplot){
      if(UCSC & plot_grtracks){
        res=get_htTrack(isTNpairs,i,Tindex,Nindex,onlyTindex,onlyNindex,indiv_run,bam_folder,plot_grtracks,UCSC,sTrack,grtrack,chr,pos,paired)
        grid.newpage()
        pushViewport(viewport(x=0,y=1, height=0.85, width=1, just=c("left","top"))) 
        plotTracks(c(ideoTrack,gtrack,res[1]),sizes=unlist(res[2]),from = pos-w, to = pos+w,add = TRUE, add53=TRUE,min.height=4, main=paste0(indiv_run[i,2],"*"),title.width=0.7,littleTicks = TRUE,cex.main=1.5)
        popViewport(1)
        pushViewport(viewport(x=0,y=0, height=0.15, width=1, just=c("left","bottom")))
        plotTracks(list(ht_zoomout),chromosome = chr, add = TRUE)
        popViewport(0)
      }else if(UCSC==FALSE | plot_grtracks==FALSE){
        res=get_htTrack(isTNpairs,i,Tindex,Nindex,onlyTindex,onlyNindex,indiv_run,bam_folder,plot_grtracks,UCSC,sTrack,grtrack,chr,pos,paired)
        plotTracks(c(ideoTrack,gtrack,res[1]),sizes=unlist(res[2]),from = pos-w, to = pos+w, add53=TRUE,min.height=4, main=paste0(indiv_run[i,2],"*"),title.width=0.7,littleTicks = TRUE,cex.main=1.5)
      }
    }
  }
  
}

