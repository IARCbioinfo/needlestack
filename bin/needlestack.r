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

#library(qvalue)
#library(discreteMTP)
glmrob.nb <- function(y,x,bounding.func='T/T',c.tukey.beta=5,c.tukey.sig=3,c.by.beta=4,weights.on.x='none',
                      minsig=1e-3,maxsig=10,minmu=1e-10,maxmu=1e5,maxit=30,maxit.sig=50,sig.prec=1e-8,tol=1e-6,
                      n_ai.sig.tukey=100,n_xout=10^4,min_coverage=1,min_reads=1,size_min=10,...){

  if (max(x)<min_coverage | max(y)<min_reads) return(res=list("coverage"=x, "ma_count"=y, "coef"=c(sigma=NA,slope=NA), "pvalues"=rep(1,l=length(y)), "qvalues"=rep(1,l=length(y)),"GQ"=rep(0,l=length(y))))
    
  ### Written by William H. Aeberhard, February 2014
  ## Disclaimer: Users of these routines are cautioned that, while due care has been taken and they are
  ## believed accurate, they have not been rigorously tested and their use and results are
  ## solely the responsibilities of the user.
  ### Modified by Matthieu Foll and Tiffany Delhomme, 2015 (follm@iarc.fr)
  #-------------------------------------------------------------------
  # General set up
  #-------------------------------------------------------------------
  yy <- y
  xx <- x
  y <- y[which(x>0)]
  x <- x[which(x>0)]
  if(length(x)<10) return(res=list("coverage"=x, "ma_count"=y, "coef"=c(sigma=NA,slope=NA), "pvalues"=rep(1,l=length(y)), "qvalues"=rep(1,l=length(y)),"GQ"=rep(0,l=length(y))))
  n <- length(y)
  X <- matrix(log(x))
  if (dim(X)[1]!=n){stop('length(y) does not match dim(X)[1].')}
  onevec <- rep(1,n)
  if (identical(X[,1],onevec)){X <- X[,-1]}
  res <- list()
  #-------------------------------------------------------------------
  # initial estimates: MLEs for beta and sigma
  #-------------------------------------------------------------------
  invlink <- function(eta){exp(eta)}
  derivlink <- function(mu){1/mu}
  varfunc <- function(mu,sig){mu+sig*mu^2}
  # initial mu computed through median outliers method
  af=y/exp(X)
  af[which(y==0)] <- 0
  inislope <- mean(af[which(af<=quantile(af,0.75)+1.5*diff(quantile(af,c(0.25,0.75))))])
  if(inislope<=1/(10*mean(exp(X)))){inislope<-1/(10*mean(exp(X)))}
  mu <- inislope*exp(X)
  mu[which(mu>maxmu)] <- maxmu
  mu[which(mu<minmu)] <- minmu
  eta <- log(mu)
  # sig MLE based on initial mu, with starting value = moment based
  sig <- sum((y/mu-1)^2)/length(y)
  if (sig>maxsig){sig <- maxsig}
  if (sig<minsig){sig <- minsig}
  #-------------------------------------------------------------------
  # Robust estimations
  #-------------------------------------------------------------------
  derivinvlink <- function(etai){exp(etai)}
  if (weights.on.x=='none'){
    weights.x <- onevec
  } else if (weights.on.x=='hard'){
    require(MASS) # for cov.rob
    Xrc <- cov.rob(X,quantile.used=floor(0.8*n))
    D2 <- mahalanobis(X,center=Xrc$center,cov=Xrc$cov) # copied from robustbase:::wts_RobDist
    qchi2 <- qchisq(p=0.95,df=dim(X)[2])
    weights.x <- ifelse(D2<=qchi2,1,0)
  } else {stop('Only "hard" and "none" are implemented for weights.on.x.')}
  derivinvlink <- function(eta){exp(eta)}
  psi.sig.ML <- function(r,mu,sig){
    digamma(r*sqrt(mu*(sig*mu+1))+mu+1/sig)-sig*r*sqrt(mu/(sig*mu+1))-digamma(1/sig)-log(sig*mu+1)
  }
  if (bounding.func=='T/T'){
    ### estimations
    tukeypsi <- function(r,c.tukey){
      ifelse(abs(r)>c.tukey,0,((r/c.tukey)^2-1)^2*r)
    }
    E.tukeypsi.1 <- function(mui,sig,c.tukey){
      sqrtVmui <- sqrt(varfunc(mui,sig))
      j1 <- max(c(ceiling(mui-c.tukey*sqrtVmui),0))
      j2 <- floor(mui+c.tukey*sqrtVmui)
      if (j1>j2){0}
      else {
        if (j2-j1+1 > 2*n_ai.sig.tukey) {
          j12 <- round(seq(j1,j2,l=n_ai.sig.tukey))
          yj12 <- (((j12-mui)/(c.tukey*sqrtVmui))^2-1)^2*(j12-mui)*dnbinom(j12,mu=mui,size=1/sig)
          integrate(splinefun(x = j12, y=yj12,method = "natural"), lower = min(j12), upper = max(j12)+1)[[1]]/sqrtVmui
        } else {
          j12 <- j1:j2
          sum((((j12-mui)/(c.tukey*sqrtVmui))^2-1)^2*(j12-mui)*dnbinom(j12,mu=mui,size=1/sig))/sqrtVmui
        }
      }
    }
    E.tukeypsi.2 <- function(mui,sig,c.tukey){
      sqrtVmui <- sqrt(varfunc(mui,sig))
      j1 <- max(c(ceiling(mui-c.tukey*sqrtVmui),0))
      j2 <- floor(mui+c.tukey*sqrtVmui)
      if (j1>j2){0}
      else {
        if (j2-j1+1 > 2*n_ai.sig.tukey) {
          j12 <- round(seq(j1,j2,l=n_ai.sig.tukey))
          yj12 <- (((j12-mui)/(c.tukey*sqrtVmui))^2-1)^2*(j12-mui)^2*dnbinom(j12,mu=mui,size=1/sig)
          integrate(splinefun(x = j12, y=yj12,method = "natural"), lower = min(j12), upper = max(j12)+1)[[1]]/sqrtVmui
        } else {
          j12 <- j1:j2
          sum((((j12-mui)/(c.tukey*sqrtVmui))^2-1)^2*(j12-mui)^2*dnbinom(j12,mu=mui,size=1/sig))/sqrtVmui
        }
      }
    }
    ai.sig.tukey <- function(mui,sig,c.tukey){  
      psi.sig.ML.mod <- function(j,mui,invsig){
        digamma(j+invsig)-digamma(invsig)-log(mui/invsig+1)-(j-mui)/(mui+invsig)
      }
      sqrtVmui <- sqrt(mui*(sig*mui+1))
      invsig <- 1/sig
      j1 <- max(c(ceiling(mui-c.tukey*sqrtVmui),0))
      j2 <- floor(mui+c.tukey*sqrtVmui)
      if (j1>j2){0}
      else {
        if (j2-j1+1 > 2*n_ai.sig.tukey) {
          j12 <- round(seq(j1,j2,l=n_ai.sig.tukey))
          yj12 <- (((j12-mui)/(c.tukey*sqrtVmui))^2-1)^2*psi.sig.ML.mod(j=j12,mui=mui,invsig=invsig)*dnbinom(x=j12,mu=mui,size=invsig)
          integrate(splinefun(x = j12, y=yj12,method = "natural"), lower = min(j12), upper = max(j12)+1)[[1]] #xout=j1:j2
        } else {
          j12 <- j1:j2
          sum((((j12-mui)/(c.tukey*sqrtVmui))^2-1)^2*psi.sig.ML.mod(j=j12,mui=mui,invsig=invsig)*dnbinom(x=j12,mu=mui,size=invsig))
        }
      }
    }
    sig.rob.tukey <- function(sig,y,mu,c.tukey){
      r <- (y-mu)/sqrt(varfunc(mu,sig))
      wi <- tukeypsi(r=r,c.tukey=c.tukey)/r
      sum(wi*psi.sig.ML(r=r,mu=mu,sig=sig)-sapply(X=mu,FUN=ai.sig.tukey,sig=sig,c.tukey=c.tukey))
    }
    sig0 <- sig+1+tol
    beta11 <- log(inislope)
    #beta11 <- log(median(y[which(y<=1.5*quantile(y,0.75))]/x[which(y<=1.5*quantile(y,0.75))])) #slope estimated with median
    beta00 <- beta11+tol+1
    it <- 0
    while((abs(sig-sig0)>tol | abs(max(beta11-beta00))>tol) & it<maxit){
      sig0 <- sig
      beta00 <- beta11
      # estimate sigma given mu
      sig=tryCatch( { 
        uniroot(f=sig.rob.tukey,interval=c(minsig,maxsig),tol=sig.prec,maxiter=maxit.sig,mu=mu,y=y,c.tukey=c.tukey.sig)$root
      },error=function(cond){ return(NA) } )
      if (is.na(sig)){sig <- minsig}
      else {
        if (sig>maxsig){sig <- maxsig}
        if (sig<minsig){sig <- minsig}
      }
      # estimate mu given sigma
      beta1 <- log(inislope)
      beta0 <- beta1+tol+1
      it.mu <- 0
      while(abs(max(beta1-beta0))>tol & it.mu<maxit){
        beta0 <- beta1
        sap <- sapply(X=mu,FUN=E.tukeypsi.2,sig=sig,c.tukey=c.tukey.beta)
        bi <- sap*varfunc(mu,sig)^(-3/2)*derivinvlink(eta)^2
        ei <- (tukeypsi(r=(y-mu)/sqrt(varfunc(mu,sig)),c.tukey=c.tukey.beta)-sapply(X=mu,FUN=E.tukeypsi.1,sig=sig,c.tukey=c.tukey.beta))/
          sap*varfunc(mu,sig)*derivlink(mu)
        zi <- eta + ei
        wls <- lm(zi~offset(X),weights=bi) ## added offset
        beta1 <- c(coef(wls),1) # coef(wls)
        eta <- fitted(wls)
        mu <- invlink(eta)
        mu[which(mu>maxmu)] <- maxmu
        mu[which(mu==0)] <- minmu
        eta <- log(mu)
        it.mu <- it.mu+1
      }
      beta11 <- beta1
      it <- it+1
    }
    #build the result
    x <- xx
    y <- yy
    res$coverage <- x
    res$ma_count <- y
    res$coef <- c(sigma=sig,slope=exp(beta1[[1]]))
    #res$pvalues_bad <- 1-(pnbinom(y,size=1/res$coef[[1]],mu=res$coef[[2]]*x))
    #res$pvalues_bad[which((y-x*res$coef["slope"])<0)] <- 1
    #res$qvalues <- qvalue(res$pvalues_bad)$qvalues
    res$pvalues <- dnbinom(y,size=1/res$coef[[1]],mu=res$coef[[2]]*x) + pnbinom(y,size=1/res$coef[[1]],mu=res$coef[[2]]*x,lower.tail = F)
    res$qvalues=p.adjust(res$pvalues,method="BH")
    res$GQ=-log10(res$qvalues)*10
#     pCDFlist=list()
#     for (i in 1:length(x)) {
#       all_y=round(seq(min(min(y),qnbinom(10^-3,mu=x[i]*res$coef[[2]],size=1/res$coef[[1]])),max(max(y),qnbinom(-10^-9,mu=x[i]*res$coef[[2]],size=1/res$coef[[1]],log.p=T)),l=1000))
#       pCDFlist[[i]]=dnbinom(all_y,mu=x[i]*res$coef[[2]],size=1/res$coef[[1]])+pnbinom(all_y,mu=x[i]*res$coef[[2]],size=1/res$coef[[1]],lower.tail=F)
#       pCDFlist[[i]]=unique(rev(pCDFlist[[i]]))
#     }
#     res$qvalues_DBH=p.discrete.adjust(res$pvalues,pCDFlist,method="DBH")
#     res$qvalues_DBL=p.discrete.adjust(res$pvalues,pCDFlist,method="DBL")
  } else {stop('Available bounding.func is "T/T"')}
  return(res)
}

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
       main=plot_title, xlim=c(0,max(rob_nb_res$coverage)),ylim=c(0,max(rob_nb_res$ma_count)))
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
      unlist(-10*log10(p.adjust((dnbinom(c(rob_nb_res$ma_count,y),size=1/rob_nb_res$coef[[1]],mu=rob_nb_res$coef[[2]]*c(rob_nb_res$coverage,x)) + 
                                   pnbinom(c(rob_nb_res$ma_count,y),size=1/rob_nb_res$coef[[1]],mu=rob_nb_res$coef[[2]]*c(rob_nb_res$coverage,x),lower.tail = F)))
                       [length(rob_nb_res$coverage)+1]))
    }
    #here we compute the dimension of the qvalue grid (ylength*xlength), with min(ylength)=5 (this avoids a too "flat" grid)
    #if needed to sampling (too large grid if dimension=max(AO)*max(DP)), we verify two equations: equality of ratios ylength/xlength before and after sampling and ylength*xlength=max_nb_grid_pts
    #### compute zoom y limit
    if(add_contours){
      nb_pts_zoom_computation=1000
      if(max(rob_nb_res$coverage) > nb_pts_zoom_computation){
        maxDP_AO = unique(sort(c(round(max(rob_nb_res$coverage)*rbeta(nb_pts_zoom_computation,1,100)),runif(100,1,max(rob_nb_res$coverage)))))      
      } else {
        maxDP_AO = seq(1,max(rob_nb_res$coverage),by=1)
      }
      maxDP_qvals = unlist(lapply(maxDP_AO,function(AO){ toQvalue(x=max(rob_nb_res$coverage),y=AO) }))
      af_min_lim = log10(maxDP_AO[which(maxDP_qvals>=qlevels[1])[1]]/max(rob_nb_res$coverage))
      ylim_zoom = maxDP_AO[which(maxDP_qvals>=max_qvalue)[1]]
      ylim_zoom_cor=ifelse(is.na(ylim_zoom),max(rob_nb_res$ma_count),ylim_zoom)
      if(!is.na(ylim_zoom)){ #ylim_zoom is na iff we found at least one qvalue >= max_qvalue (if error rate closed to 1, only qvalues closed to 0)
        #### compute dim of the qvalue grid
        if(ylim_zoom*max(rob_nb_res$coverage) <= max_nb_grid_pts){
          xgrid = seq(0,max(rob_nb_res$coverage), by=1) 
          ygrid = seq(0,ylim_zoom,by=1) #use by=1 to have integer, if not dnbinom not happy
        } else {
          if(ylim_zoom<=50){
            ygrid=seq(0,ylim_zoom,by=1)
            xlength = round(max_nb_grid_pts/length(ygrid))
            xgrid = round(seq(0,max(rob_nb_res$coverage),length=xlength))
          } else {
            ylength = 50
            xlength = round(max_nb_grid_pts/ylength)
            xgrid = round(seq(0,max(rob_nb_res$coverage),length=xlength))
            ygrid = round(seq(0,ylim_zoom,length=ylength))    
          }
        }
        #here we initiate the grid with each case containing list=(DP,AO) from xgrid, ygrid
        matgrid = array(as.list(as.data.frame(t(expand.grid(xgrid,ygrid)))),dim=c(length(xgrid),length(ygrid)))
        #then we fill in the grid with qvalues for each pair of AO,DP taken from ygrid,xgrid vectors (we use initiated values to identify corresponding AO,DP). Finally we plot the contours.
        matgrid=matrix(sapply(matgrid,function(case) toQvalue(unlist(case)[1],unlist(case)[2])), length(xgrid),length(ygrid))
        #### plot the contour "by hands"
        for(qvalue in qlevels) {
          lines(xgrid, unlist(lapply(xgrid,function(DP,ygrid,xgrid){
            qval=min(matgrid[match(DP,xgrid),which(matgrid[match(DP,xgrid),]>=qvalue)])
            AO=min(ygrid[which(matgrid[match(DP,xgrid),]==qval)]) 
            AO },ygrid,xgrid)),col=rev(rainbow(length(qlevels),start=0, end=4/6))[match(qvalue,qlevels)],lwd=1.3,lty=3)
        }
      }
    }    
    #### plot confidence interval + error rate
    xi=max(rob_nb_res$coverage)
    yi1=qnbinom(p=0.99, size=1/rob_nb_res$coef[[1]], mu=rob_nb_res$coef[[2]]*xi)
    yi2=qnbinom(p=0.01, size=1/rob_nb_res$coef[[1]], mu=rob_nb_res$coef[[2]]*xi)
    abline(a=0, b=yi1/xi, lwd=2, lty=3, col="blue")
    abline(a=0, b=yi2/xi, lwd=2, lty=3, col="blue")
    abline(a=0, b=rob_nb_res$coef[[2]], col="blue")
    #### plot zoom on max qvalue if add_contours, otherwise on 2*IC
    if(!add_contours) ylim_zoom_cor = 2*yi1
    plot(rob_nb_res$coverage, rob_nb_res$ma_count,
         pch=21,bg=cols,col=outliers_color,xlab="Coverage (DP)",ylab="Number of ALT reads (AO)",
         main=plot_title, ylim=c(0,ylim_zoom_cor), xlim=c(0,max(rob_nb_res$coverage)))
    #### plot the contour "by hands"
    if(add_contours){
      if(!is.na(ylim_zoom)){
        for(qvalue in qlevels) {
          lines(xgrid, unlist(lapply(xgrid,function(DP,ygrid,xgrid){
            qval=min(matgrid[match(DP,xgrid),which(matgrid[match(DP,xgrid),]>=qvalue)])
            AO=min(ygrid[which(matgrid[match(DP,xgrid),]==qval)]) 
            AO },ygrid,xgrid)),col=rev(rainbow(length(qlevels),start=0, end=4/6))[match(qvalue,qlevels)],lwd=1.3,lty=3)
        }
      }
    }
    #contour(xgrid, ygrid, matgrid, levels=qlevels , col = rev(rainbow(length(qlevels),start=0, end=4/6)), add=T, lwd = 1.3, labcex = 0.8, lty=3)
    mtext(paste("zoom on maximum q-value =",max_qvalue))
    #### labeling outliers and plot confidence interval + error rate
    if(!is.null(names) & plot_labels & length(names[which(rob_nb_res$qvalues<=qthreshold)]) > 0) {
      text(rob_nb_res$coverage[which(rob_nb_res$qvalues<=qthreshold)], rob_nb_res$ma_count[which(rob_nb_res$qvalues<=qthreshold)],
         labels=names[which(rob_nb_res$qvalues<=qthreshold)], cex= 0.6, pos=1)
    }
    abline(a=0, b=yi1/xi, lwd=2, lty=3, col="blue")
    abline(a=0, b=yi2/xi, lwd=2, lty=3, col="blue")
    abline(a=0, b=rob_nb_res$coef[[2]], col="blue")
    plot_palette()  
    
    logqvals=log10(rob_nb_res$qvalues)
    logqvals[logqvals<=-9]=-9
    plot(logqvals,rob_nb_res$ma_count/rob_nb_res$coverage,pch=21,bg=cols,col=outliers_color,ylab="Allelic fraction (AF)",xlab=bquote("log"[10] ~ "(q-value)"),main="Allelic fraction effect")
    abline(v=log10(qthreshold),col="red",lwd=2)
    plot_palette(topright = TRUE)
    ylim_zoom_af = ifelse(add_contours, ylim_zoom_cor/max(rob_nb_res$coverage), (2*yi1)/xi)
    plot(logqvals,rob_nb_res$ma_count/rob_nb_res$coverage,pch=21,bg=cols,col=outliers_color,ylab="Allelic fraction (AF)",xlab=bquote("log"[10] ~ "(q-value)"),main="Allelic fraction effect",
         ylim=c(0,ylim_zoom_af))
    if(add_contours) { mtext(paste("zoom on maximum q-value =",max_qvalue)) } else { mtext("zoom on 99% confidence interval") }
    abline(v=log10(qthreshold),col="red",lwd=2)
    plot_palette(topright = TRUE)
    if(add_contours){
      if(!is.na(ylim_zoom)){
        #### plot min(af) ~ DP
        plot(1,type='n', ylim=c(1.1*af_min_lim,0), xlim=range(xgrid), xlab="DP", ylab=bquote("log"[10] ~ "[min(AF)]"))
        for(qvalue in qlevels) {
          lines(xgrid, unlist(lapply(xgrid,function(DP,ygrid,xgrid){
            qval=min(matgrid[match(DP,xgrid),which(matgrid[match(DP,xgrid),]>=qvalue)])
            af=min(ygrid[which(matgrid[match(DP,xgrid),]==qval)]) / DP
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

############################## ARGUMENTS SECTION ##############################
args <- commandArgs(TRUE)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsL <- as.list(as.character(as.data.frame(do.call("rbind", parseArgs(args)))$V2))
names(argsL) <- as.data.frame(do.call("rbind", parseArgs(args)))$V1
args <- argsL;rm(argsL)

if("--help" %in% args | is.null(args$out_file) | is.null(args$fasta_ref)) {
  cat("
      The R Script arguments_section.R
      
      Mandatory arguments:
      --out_file=file_name           - name of output vcf
      --fasta_ref=path               - path of fasta ref
      --help                         - print this text
      
      Optionnal arguments:
      --samtools=path                - path of samtools, default=samtools
      --SB_type=SOR or RVSB      - strand bias measure, default=SOR
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

############################################################

indiv_run=read.table("names.txt",stringsAsFactors=F,colClasses = "character")
indiv_run[,2]=make.unique(indiv_run[,2],sep="_")

pileups_files=paste(indiv_run[,1],".txt",sep="")
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
  #FisherStrand<<-fisher.test(matrix(c(Vp,Vm,Rp,Rm),nrow=2))$p.value 
  all_RO<<-sum(Rp+Rm)
}

if(file.exists(out_file)) file.remove(out_file) 
write_out=function(...) {
  cat(paste(...,sep=""),"\n",file=out_file,sep="",append=T)
}
write_out("##fileformat=VCFv4.1")
write_out("##fileDate=",format(Sys.Date(), "%Y%m%d"))
write_out("##source=needlestack v0.3")
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
        if (SB_type=="SOR") sbs=sors else sbs=rvsbs
        if (output_all_SNVs | (sum(sbs<=SB_threshold_SNV & reg_res$GQ>=GQ_threshold,na.rm=TRUE)>0)) {
          con=pipe(paste(samtools," faidx ",fasta_ref," ",pos_ref[i,"chr"],":",pos_ref[i,"loc"]-3,"-",pos_ref[i,"loc"]-1," | tail -n1",sep=""))
  		    before=readLines(con)
          close(con)
          con=pipe(paste(samtools," faidx ",fasta_ref," ",pos_ref[i,"chr"],":",pos_ref[i,"loc"]+1,"-",pos_ref[i,"loc"]+3," | tail -n1",sep=""))
          after=readLines(con)
  		    close(con)
  		    cat(pos_ref[i,"chr"],"\t",pos_ref[i,"loc"],"\t",".","\t",pos_ref[i,"ref"],"\t",alt,"\t",max(reg_res$GQ),"\t",".",sep = "",file=out_file,append=T)
          # INFO field
  		    cat("\t","TYPE=snv;NS=",sum(coverage_matrix[i,]>0),";AF=",sum(reg_res$GQ>=GQ_threshold)/sum(coverage_matrix[i,]>0),";DP=",all_DP,";RO=",all_RO,";AO=",all_AO,";SRF=",sum(Rp),";SRR=",sum(Rm),";SAF=",sum(Vp),";SAR=",sum(Vm),";SOR=",all_sor,";RVSB=",all_rvsb,";ERR=",reg_res$coef["slope"],";SIG=",reg_res$coef["sigma"],";CONT=",paste(before,after,sep="x"),sep="",file=out_file,append=T)
  		    # FORMAT field
  		    cat("\t","GT:QVAL:DP:RO:AO:AF:SB:SOR:RVSB",sep = "",file=out_file,append=T)
          # all samples
  		    genotype=rep("0/0",l=nindiv)
  		    variants=which(reg_res$GQ>=GQ_threshold & sbs<=SB_threshold_SNV)
  		    genotype[variants]="0/1"
          for (cur_sample in 1:nindiv) {
            cat("\t",genotype[cur_sample],":",reg_res$GQ[cur_sample],":",DP[cur_sample],":",(Rp+Rm)[cur_sample],":",ma_count[cur_sample],":",(ma_count/DP)[cur_sample],":",Rp[cur_sample],",",Rm[cur_sample],",",Vp[cur_sample],",",Vm[cur_sample],":",sors[cur_sample],":",rvsbs[cur_sample],sep = "",file=out_file,append=T)
          }
  		    cat("\n",sep = "",file=out_file,append=T)
          if (do_plots) {
            pdf(paste(pos_ref[i,"chr"],"_",pos_ref[i,"loc"],"_",pos_ref[i,"loc"],"_",pos_ref[i,"ref"],"_",alt,".pdf",sep=""),7,6)
            plot_rob_nb(reg_res, 10^-(GQ_threshold/10), plot_title=bquote(paste(.(pos_ref[i,"loc"])," (",.(pos_ref[i,"ref"]) %->% .(alt),")",sep="")), sbs=sbs, SB_threshold=SB_threshold_SNV,plot_labels=plot_labels,add_contours=add_contours,names=indiv_run[,2])
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
          if (SB_type=="SOR") sbs=sors else sbs=rvsbs
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
            cat("\t","TYPE=del;NS=",sum(coverage_matrix[i,]>0),";AF=",sum(reg_res$GQ>=GQ_threshold)/sum(coverage_matrix[i,]>0),";DP=",all_DP,";RO=",all_RO,";AO=",all_AO,";SRF=",sum(Rp),";SRR=",sum(Rm),";SAF=",sum(Vp),";SAR=",sum(Vm),";SOR=",all_sor,";RVSB=",all_rvsb,";ERR=",reg_res$coef["slope"],";SIG=",reg_res$coef["sigma"],";CONT=",paste(before,after,sep="x"),sep="",file=out_file,append=T)
            # FORMAT field
            cat("\t","GT:QVAL:DP:RO:AO:AF:SB:SOR:RVSB",sep = "",file=out_file,append=T)
            # all samples
            genotype=rep("0/0",l=nindiv)
            variants=which(reg_res$GQ>=GQ_threshold & sbs<=SB_threshold_indel)
            genotype[variants]="0/1"
            for (cur_sample in 1:nindiv) {
              cat("\t",genotype[cur_sample],":",reg_res$GQ[cur_sample],":",DP[cur_sample],":",(Rp+Rm)[cur_sample],":",ma_count[cur_sample],":",(ma_count/DP)[cur_sample],":",Rp[cur_sample],",",Rm[cur_sample],",",Vp[cur_sample],",",Vm[cur_sample],":",sors[cur_sample],":",rvsbs[cur_sample],sep = "",file=out_file,append=T)
            }
            cat("\n",sep = "",file=out_file,append=T)
            if (do_plots) {
              # deletions are shifted in samtools mpileup by 1bp, so put them at the right place by adding + to pos_ref[i,"loc"] everywhere in what follows
              pdf(paste(pos_ref[i,"chr"],"_",pos_ref[i,"loc"]+1,"_",pos_ref[i,"loc"]+1+nchar(cur_del)-1,"_",cur_del,"_","-",".pdf",sep=""),7,6)
              plot_rob_nb(reg_res, 10^-(GQ_threshold/10), plot_title=bquote(paste(.(pos_ref[i,"loc"]+1)," (",.(cur_del) %->% .("-"),")",sep="")),sbs=sbs, SB_threshold=SB_threshold_indel,plot_labels=plot_labels,add_contours=add_contours,names=indiv_run[,2])
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
          if (SB_type=="SOR") sbs=sors else sbs=rvsbs
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
            cat("\t","TYPE=ins;NS=",sum(coverage_matrix[i,]>0),";AF=",sum(reg_res$GQ>=GQ_threshold)/sum(coverage_matrix[i,]>0),";DP=",all_DP,";RO=",all_RO,";AO=",all_AO,";SRF=",sum(Rp),";SRR=",sum(Rm),";SAF=",sum(Vp),";SAR=",sum(Vm),";SOR=",all_sor,";RVSB=",all_rvsb,";ERR=",reg_res$coef["slope"],";SIG=",reg_res$coef["sigma"],";CONT=",paste(before,after,sep="x"),sep="",file=out_file,append=T)
            # FORMAT field
            cat("\t","GT:QVAL:DP:RO:AO:AF:SB:SOR:RVSB",sep = "",file=out_file,append=T)
            # all samples
            genotype=rep("0/0",l=nindiv)
            variants=which(reg_res$GQ>=GQ_threshold & sbs<=SB_threshold_indel)
            genotype[variants]="0/1"
            for (cur_sample in 1:nindiv) {
              cat("\t",genotype[cur_sample],":",reg_res$GQ[cur_sample],":",DP[cur_sample],":",(Rp+Rm)[cur_sample],":",ma_count[cur_sample],":",(ma_count/DP)[cur_sample],":",Rp[cur_sample],",",Rm[cur_sample],",",Vp[cur_sample],",",Vm[cur_sample],":",sors[cur_sample],":",rvsbs[cur_sample],sep = "",file=out_file,append=T)
            }
            cat("\n",sep = "",file=out_file,append=T)
            if (do_plots) {
              pdf(paste(pos_ref[i,"chr"],"_",pos_ref[i,"loc"],"_",pos_ref[i,"loc"],"_","-","_",cur_ins,".pdf",sep=""),7,6)
              plot_rob_nb(reg_res, 10^-(GQ_threshold/10), plot_title=bquote(paste(.(pos_ref[i,"loc"])," (",.("-") %->% .(cur_ins),")",sep="")),sbs=sbs, SB_threshold=SB_threshold_indel,plot_labels=plot_labels,add_contours=add_contours,names=indiv_run[,2])
              dev.off()
            }
          }
        }
      }
    }
  }
}
