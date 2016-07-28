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

glmrob.nb <- function(y,x,bounding.func='T/T',c.tukey.beta=5,c.tukey.sig=3,c.by.beta=4,weights.on.x='none',
                      minsig=1e-3,maxsig=10,minmu=1e-10,maxmu=1e5,maxit=30,maxit.sig=50,sig.prec=1e-8,tol=1e-6,
                      n_ai.sig.tukey=100,n_xout=10^4,min_coverage=1,min_reads=1,size_min=10,snps=NULL,...){

  if (max(x,na.rm = T)<min_coverage | max(y, na.rm = T)<min_reads | length(x[which(x>0)])<size_min ) {
    if(is.null(snps)) return(res=list("coverage"=x, "ma_count"=y, "coef"=c(sigma=NA,slope=NA), "pvalues"=rep(1,l=length(y)), "qvalues"=rep(1,l=length(y)),"GQ"=rep(0,l=length(y))))
    if(!is.null(snps) & ( max(c(x,snps$DP_snp),na.rm=T)<min_coverage | max(c(y,snps$ma_count_snp),na.rm=T)<min_reads ) ) {
      coverage = ma_count = rep(0,length(y)+length(snps$snp_pos))
      ma_count[setdiff(1:(length(y)+length(snps$snp_pos)),snps$snp_pos)]=y
      coverage[setdiff(1:(length(y)+length(snps$snp_pos)),snps$snp_pos)]=x
      ma_count[snps$snp_pos]=snps$ma_count_snp
      coverage[snps$snp_pos]=snps$DP_snp
      return(res=list("coverage"=coverage, "ma_count"=ma_count, "coef"=c(sigma=NA,slope=NA), "pvalues"=rep(1,l=length(ma_count)), "qvalues"=rep(1,l=length(ma_count)),"GQ"=rep(0,l=length(ma_count))))
    }
  }
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
    if(!is.null(snps)) {
      coverage = ma_count = rep(0,length(y)+length(snps$snp_pos))
      ma_count[setdiff(1:(length(y)+length(snps$snp_pos)),snps$snp_pos)]=y
      coverage[setdiff(1:(length(y)+length(snps$snp_pos)),snps$snp_pos)]=x
      ma_count[snps$snp_pos]=snps$ma_count_snp
      coverage[snps$snp_pos]=snps$DP_snp
      y=ma_count; x=coverage
    }
    res$coverage <- x
    res$ma_count <- y
    res$coef <- c(sigma=sig,slope=exp(beta1[[1]]))
    res$pvalues <- dnbinom(y,size=1/res$coef[[1]],mu=res$coef[[2]]*x) + pnbinom(y,size=1/res$coef[[1]],mu=res$coef[[2]]*x,lower.tail = F)
    res$qvalues=p.adjust(res$pvalues,method="BH")
    res$GQ=-log10(res$qvalues)*10
    res$GQ[res$GQ>1000]=1000 #here also manage qvalues=Inf
  } else {stop('Available bounding.func is "T/T"')}
  return(res)
}
