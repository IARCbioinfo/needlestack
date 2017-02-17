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

plot_rob_nb <- function(rob_nb_res,qthreshold=0.01,plot_title=NULL,sbs,SB_threshold=Inf,names=NULL,plot_labels=FALSE,add_contours=FALSE){
  n=sum(rob_nb_res$qvalue>qthreshold)
  m=sum(rob_nb_res$qvalue<=qthreshold)
  
  cut_max_qvals=100
  cols=rep("black",length(rob_nb_res$coverage))
  palette=rev(rainbow(cut_max_qvals+1,start=0, end=4/6))
  col_indices = round(-10*log10(rob_nb_res$qvalues))+1
  col_indices[which(col_indices>cut_max_qvals)]=cut_max_qvals+1
  cols=palette[col_indices]
  outliers_color=cols
  outliers_color[which(sbs>SB_threshold)]="black"
  
  temp_title = bquote(e==.(format(rob_nb_res$coef[[2]],digits = 2)) ~ "," ~ sigma==.(format(rob_nb_res$coef[[1]],digits = 2)))
  plot(rob_nb_res$coverage, rob_nb_res$ma_count,
       pch=21,bg=cols,col=outliers_color,xlab="Coverage (DP)",ylab="Number of ALT reads (AO)",
       main=plot_title, xlim=c(0,max(rob_nb_res$coverage,na.rm = T)),ylim=c(0,max(rob_nb_res$ma_count, na.rm = T)))
  mtext(temp_title)
  #### labeling outliers
  if(!is.null(names) & plot_labels & length(names[which(rob_nb_res$qvalues<=qthreshold)]) > 0) {
    text(rob_nb_res$coverage[which(rob_nb_res$qvalues<=qthreshold)], rob_nb_res$ma_count[which(rob_nb_res$qvalues<=qthreshold)],
         labels=names[which(rob_nb_res$qvalues<=qthreshold)], cex= 0.6, pos=1)
  }
  
  #### plot the color palette
  plot_palette <- function(topright=FALSE) {
    xmin <- par("usr")[1]
    xmax <- par("usr")[2]
    ymin <- par("usr")[3]
    ymax <- par("usr")[4]
    xright=xmin+(xmax-xmin)*(1-0.9)
    xleft=xmin+(xmax-xmin)*(1-0.94)
    if(topright){xright=xmin+(xmax-xmin)*0.9;xleft=xmin+(xmax-xmin)*0.94}
    ybottom=ymin+(ymax-ymin)*0.72
    ytop=ymin+(ymax-ymin)*0.94
    
    rasterImage(as.raster(matrix(rev(palette), ncol=1)),xright ,ybottom ,xleft,ytop )
    rect(xright ,ybottom ,xleft,ytop )
    text(x=(xright+xleft)/2, y = ytop+(ytop-ybottom)*0.1, labels = "QVAL", cex=0.8)
    keep_labels=seq(0,cut_max_qvals,by=20)
    keep_labels_pos=seq(ybottom,ytop,l=length(keep_labels))
    tick_width=-(xleft-xright)/5
    for (i in 1:length(keep_labels)) {
      if (topright) {
        lines(c(xright,xright+tick_width),c(keep_labels_pos[i],keep_labels_pos[i]))
      } else {
        lines(c(xleft,xleft-tick_width),c(keep_labels_pos[i],keep_labels_pos[i]))
      }
    }
    if(topright) {
      text(x=xright+1.5*tick_width, y = keep_labels_pos, labels = keep_labels,adj = c(1,0.5), cex = 0.8)
    } else {
      text(x=(xright-(xright-xleft))*0.3, y = keep_labels_pos, labels = keep_labels,adj = c(1,0.5), cex = 0.8)
    }
  }
  plot_palette()
  ####
  
  if (!is.na(rob_nb_res$coef["slope"])) {
    ################### ADD CONTOURS ##################
    max_nb_grid_pts = 2500
    max_qvalue=100
    qlevels = c(10,30,50,70,100)
    # following function returns a qvalue for a new point by adding it in the set of observations, recomputing all qvalues and taking its corresponding value.
    toQvalue <- function(x,y){
      if(y>x) return(0)
      unlist(-10*log10(p.adjust((dnbinom(c(rob_nb_res$ma_count,y),size=1/rob_nb_res$coef[[1]],mu=rob_nb_res$coef[[2]]*c(rob_nb_res$coverage,x)) +
                                   pnbinom(c(rob_nb_res$ma_count,y),size=1/rob_nb_res$coef[[1]],mu=rob_nb_res$coef[[2]]*c(rob_nb_res$coverage,x),lower.tail = F)))
                       [length(rob_nb_res$coverage)+1]))
    }
    #here we compute the dimension of the qvalue grid (ylength*xlength), with min(ylength)=5 (this avoids a too "flat" grid)
    #if needed to sampling (too large grid if dimension=max(AO)*max(DP)), we verify two equations: equality of ratios ylength/xlength before and after sampling and ylength*xlength=max_nb_grid_pts
    #### compute zoom y limit
    if(add_contours){
      nb_pts_zoom_computation=1000
      if(max(rob_nb_res$coverage,na.rm=T) > nb_pts_zoom_computation){
        maxDP_AO = unique(sort(c(round(max(rob_nb_res$coverage,na.rm=T)*rbeta(nb_pts_zoom_computation,1,100)),runif(100,1,max(rob_nb_res$coverage,na.rm=T)))))
      } else {
        maxDP_AO = seq(1,max(rob_nb_res$coverage,na.rm=T),by=1)
      }
      maxDP_qvals = unlist(lapply(maxDP_AO,function(AO){ toQvalue(x=max(rob_nb_res$coverage,na.rm=T),y=AO) }))
      af_min_lim = log10(maxDP_AO[which(maxDP_qvals>=qlevels[1])[1]]/max(rob_nb_res$coverage,na.rm=T))
      if(is.na(af_min_lim)) af_min_lim=0 #if maxDP_qvals does not contain qlevels[1] 
      ylim_zoom = maxDP_AO[which(maxDP_qvals>=max_qvalue)[1]]
      ylim_zoom_cor=ifelse(is.na(ylim_zoom),max(rob_nb_res$coverage,na.rm=T),ylim_zoom)
      if(!is.na(ylim_zoom_cor)){ #ylim_zoom is na iff we found at least one qvalue >= max_qvalue (if error rate closed to 1, only qvalues closed to 0)
        #### compute dim of the qvalue grid
        if(ylim_zoom_cor*max(rob_nb_res$coverage,na.rm=T) <= max_nb_grid_pts){
          xgrid = seq(0,max(rob_nb_res$coverage,na.rm=T), by=1)
          ygrid = seq(0,ylim_zoom_cor,by=1) #use by=1 to have integer, if not dnbinom not happy
        } else {
          if(ylim_zoom_cor<=50){
            ygrid=seq(0,ylim_zoom_cor,by=1)
            xlength = round(max_nb_grid_pts/length(ygrid))
            xgrid = round(seq(0,max(rob_nb_res$coverage,na.rm=T),length=xlength))
          } else {
            ylength = 50
            xlength = round(max_nb_grid_pts/ylength)
            xgrid = round(seq(0,max(rob_nb_res$coverage,na.rm=T),length=xlength))
            ygrid = round(seq(0,ylim_zoom_cor,length=ylength))
          }
        }
        #here we initiate the grid with each case containing list=(DP,AO) from xgrid, ygrid
        matgrid = array(as.list(as.data.frame(t(expand.grid(xgrid,ygrid)))),dim=c(length(xgrid),length(ygrid)))
        #then we fill in the grid with qvalues for each pair of AO,DP taken from ygrid,xgrid vectors (we use initiated values to identify corresponding AO,DP). Finally we plot the contours.
        matgrid=matrix(sapply(matgrid,function(case) toQvalue(unlist(case)[1],unlist(case)[2])), length(xgrid),length(ygrid))
        #### plot the contour "by hands"
        for(qvalue in qlevels) {
          dat = na.omit(cbind(xgrid,unlist(lapply(xgrid,function(DP,ygrid,xgrid){
            if(sum(matgrid[match(DP,xgrid),]>=qvalue)==0) {
              qval = NA } else {
                qval = min(matgrid[match(DP,xgrid),which(matgrid[match(DP,xgrid),]>=qvalue)]) 
              } #NA iff no case>=qvalue in matgrid at DP
            if (is.na(qval)) {
              if( DP == 0 ) { AO = 0 } else { AO = NA }
            } else {
              AO = min(ygrid[which(matgrid[match(DP,xgrid),]==qval)])
            }
            AO },ygrid,xgrid))))
          lines(dat[,1],dat[,2],col=rev(rainbow(length(qlevels),start=0, end=4/6))[match(qvalue,qlevels)],lwd=1.3,lty=3)
        }
      }
    }
    #### plot confidence interval + error rate
    xi=max(rob_nb_res$coverage,na.rm=T)
    yi1=qnbinom(p=0.99, size=1/rob_nb_res$coef[[1]], mu=rob_nb_res$coef[[2]]*xi)
    yi2=qnbinom(p=0.01, size=1/rob_nb_res$coef[[1]], mu=rob_nb_res$coef[[2]]*xi)
    if(max(rob_nb_res$coverage,na.rm=T)>1000) { DP_for_IC = round(seq(0,max(rob_nb_res$coverage,na.rm=T),length=1000)) } else { DP_for_IC = seq(0,max(rob_nb_res$coverage,na.rm=T),by=1) }
    lines(DP_for_IC,qnbinom(p=0.99, size=1/rob_nb_res$coef[[1]], mu=rob_nb_res$coef[[2]]*DP_for_IC),col="black",lty=3,lwd=2)
    lines(DP_for_IC,qnbinom(p=0.01, size=1/rob_nb_res$coef[[1]], mu=rob_nb_res$coef[[2]]*DP_for_IC),col="black",lty=3,lwd=2)
    abline(a=0, b=rob_nb_res$coef[[2]], col="black")
    #### plot zoom on max qvalue if add_contours, otherwise on 2*IC
    if(!add_contours) ylim_zoom_cor = 2*yi1
    plot(rob_nb_res$coverage, rob_nb_res$ma_count,
         pch=21,bg=cols,col=outliers_color,xlab="Coverage (DP)",ylab="Number of ALT reads (AO)",
         main=plot_title, ylim=c(0,ylim_zoom_cor), xlim=c(0,max(rob_nb_res$coverage,na.rm=T)))
    #### plot the contour "by hands"
    if(add_contours){
      if(!is.na(ylim_zoom_cor)){
        for(qvalue in qlevels) {
          dat = na.omit(cbind(xgrid,unlist(lapply(xgrid,function(DP,ygrid,xgrid){
            if(sum(matgrid[match(DP,xgrid),]>=qvalue)==0) {
              qval = NA } else {
                qval = min(matgrid[match(DP,xgrid),which(matgrid[match(DP,xgrid),]>=qvalue)]) 
              } #NA iff no case>=qvalue in matgrid at DP
            if (is.na(qval)) {
              if( DP == 0 ) { AO = 0 } else { AO = NA }
            } else {
              AO = min(ygrid[which(matgrid[match(DP,xgrid),]==qval)])
            }
            AO },ygrid,xgrid))))
          lines(dat[,1],dat[,2],col=rev(rainbow(length(qlevels),start=0, end=4/6))[match(qvalue,qlevels)],lwd=1.3,lty=3)        }
      }
    }
    #contour(xgrid, ygrid, matgrid, levels=qlevels , col = rev(rainbow(length(qlevels),start=0, end=4/6)), add=T, lwd = 1.3, labcex = 0.8, lty=3)
    mtext(paste("zoom on maximum q-value =",max_qvalue))
    #### labeling outliers and plot confidence interval + error rate
    if(!is.null(names) & plot_labels & length(names[which(rob_nb_res$qvalues<=qthreshold)]) > 0) {
      text(rob_nb_res$coverage[which(rob_nb_res$qvalues<=qthreshold)], rob_nb_res$ma_count[which(rob_nb_res$qvalues<=qthreshold)],
           labels=names[which(rob_nb_res$qvalues<=qthreshold)], cex= 0.6, pos=1)
    }
    lines(DP_for_IC,qnbinom(p=0.99, size=1/rob_nb_res$coef[[1]], mu=rob_nb_res$coef[[2]]*DP_for_IC),col="black",lty=3,lwd=2)
    lines(DP_for_IC,qnbinom(p=0.01, size=1/rob_nb_res$coef[[1]], mu=rob_nb_res$coef[[2]]*DP_for_IC),col="black",lty=3,lwd=2)
    abline(a=0, b=rob_nb_res$coef[[2]], col="black")
    plot_palette()
    
    plot(rob_nb_res$GQ,log10(rob_nb_res$ma_count/rob_nb_res$coverage),pch=21,bg=cols,col=outliers_color,ylab=bquote("log"[10] ~ "[Allelic Fraction (AF)]"),xlab="QVAL",main="Allelic fraction effect")
    abline(v=-10*log10(qthreshold),col="red",lwd=2)
    plot_palette()
    ylim_zoom_af = ifelse(add_contours, ylim_zoom_cor/max(rob_nb_res$coverage,na.rm=T), (2*yi1)/xi)
    plot(rob_nb_res$GQ,log10(rob_nb_res$ma_count/rob_nb_res$coverage),pch=21,bg=cols,col=outliers_color,ylab=bquote("log"[10] ~ "[Allelic Fraction (AF)]"),xlab="QVAL",main="Allelic fraction effect",
         ylim=c(min(log10(rob_nb_res$ma_count/rob_nb_res$coverage)[is.finite(log10(rob_nb_res$ma_count/rob_nb_res$coverage))]),log10(ylim_zoom_af)), xlim=c(0,100))
    if(add_contours) { mtext(paste("zoom on maximum q-value =",max_qvalue)) } else { mtext("zoom on 99% confidence interval") }
    abline(v=-10*log10(qthreshold),col="red",lwd=2)
    plot_palette()
    if(add_contours){
      if(!is.na(ylim_zoom_cor)){
        #### plot min(af) ~ DP
        plot(1,type='n', ylim=c(1.1*af_min_lim,0), xlim=range(xgrid), xlab="DP", ylab=bquote("log"[10] ~ "[min(AF)]"))
        for(qvalue in qlevels) {
          lines(xgrid, unlist(lapply(xgrid,function(DP,ygrid,xgrid){
            if(sum(matgrid[match(DP,xgrid),]>=qvalue)==0) {
              qval = NA } else {
              qval = min(matgrid[match(DP,xgrid),which(matgrid[match(DP,xgrid),]>=qvalue)]) 
            } #NA iff no case>=qvalue in matgrid at DP
            if (is.na(qval)) {
              af = 0
            } else {
              af = min(ygrid[which(matgrid[match(DP,xgrid),]==qval)])/DP
            }
            if(DP==0 || af>1) { af=1 } #af>1 if min(...)>DP
            log10(af) },ygrid,xgrid)),col=rev(rainbow(length(qlevels),start=0, end=4/6))[match(qvalue,qlevels)])
        }
        plot_palette(topright = TRUE)
        #hist(rob_nb_res$pvalues,main="p-values distribution",ylab="Density",xlab="p-value",col="grey",freq=T,br=20,xlim=c(0,1))
        #hist(rob_nb_res$qvalues,main="q-values distribution",breaks=20,xlab="q-value",col="grey",freq=T,xlim=c(0,1))
      }
    }
  }
}
