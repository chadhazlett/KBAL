#' Build Gaussian kernel matrix.
#' @description  Centers and rescales X then computes Guassian
#' kernel matrix. Entry {i,j} correspond to k(X_i,X_j) where k is the (Gaussian) kernel.
#' @param X A numeric matrix of data, typically a design matrix.
#' @param sigma The kernel ``bandwidth''. If NULL, defaults to ncol(X)
#' @examples
#' K=buildgauss(rnorm(10))
#' @export
buildgauss = function(X,sigma=NULL){
	X = as.matrix(X)
	if (is.numeric(X) == FALSE) {stop("X must be numeric")}
	if (sum(is.na(X)) > 0) {stop("X contains missing data")}
	n = nrow(X)
	d = ncol(X)
	X = scale(X, center = TRUE, scale = TRUE)
	if (is.null(sigma)) {	sigma = d}
	K=exp(-1 * as.matrix(dist(X)^2)/sigma)
  return(K)
}

#' Get the bound on the bias due to incomplete balance
#' @description XXX
#' @param D xxx .
#' @param w xxx .
#' @param V xxx .
#' @param a xxx .
#' @param hilbertnorm xxx .
#' @examples
#' biasbound=(D=D, w=w, V=svd.out$u, a=svd.out$d, hilbertnorm=1)
#' @export
biasbound=function(D,w,V,a, hilbertnorm=1){
  w1=w[D==1]/sum(D==1)
  w0=w[D==0]/sum(D==0)

  #optionally look only at remaining imbalance assuming first numdims
  #components were perfectly balanced.
  #V=V[,(numdimsbalanced+1):length(w)]
  #a=a[(numdimsbalanced+1):length(w)]

  V1=V[D==1, , drop=FALSE]
  V0=V[D==0, , drop=FALSE]
  eigenimbal=t(w1)%*%V1 -t(w0)%*%V0

  #eigenimbal*(a^.5)%*%t(V)
  effectiveimbal=(eigenimbal*(a^.5)) #%*%t(V)
  biasbound=sqrt(hilbertnorm)*sqrt(effectiveimbal%*%t(effectiveimbal))
  return(biasbound)
}

#' Function to multiply a square matrix, X, with a diagonal matrix, diag(d)
#' @description  Fast multiplication by a diagonal matrix
#' @param X A numeric matrix of data, typically a design matrix.
#' @param d a vector, such that you want to know X%*%diag(d)
#' @examples
#' Xd=multdiag(X,d)
multdiag <- function(X,d){
  R=matrix(NA,nrow=dim(X)[1],ncol=dim(X)[2])
  for (i in 1:dim(X)[2]){
    R[,i]=X[,i]*d[i]
  }
  return(R)
}


#' Kernel balancing function.
#' @description Chooses weights on control units that produces equal means on a kernel matrix, K, rather than the original design matrix, X.
#' @param X The original covariate data, as a numeric matrix.
#' @param D The treatment assignment variable taking values of 1 for treated units and 0 for control units.
#' @param K Optional user-provided kernel matrix. Typically this is not user-specified, but rather is computed internally by a call to \code{buildgauss}.
#' @param whiten Optional pre-whitening of the data prior to construction of K. If used, rotates the data by \code{solve(chol(var(X)))}, then centers and rescales.
#' @param trimratio Optional \code{logical}
#' @param numdims Optional user-specified number of projectionss of \code{K} to balance upon.
#' @param minnumdims Optional user-specified choice for the minimum number of projections of \code{K}. Defualts to 2.
#' @param maxnumdims Optional user-specified choice for the maximum number of projectsion ffo \code{K} to attempt balance on. Defaults to the number of control units.
#' @param sigma Optional user-specificied paramater for the Gaussian kernel. If blank, defaults to \code{nrow(X)}.
#' @param method "ebal" or "el". Whether balance should be obtained on each projection of \code{K} using entropy balancing ("ebal", default) or empirical likelihood ("el")
#' @return \item{w}{The weights, taking value of 1 for treated unit and the estimated weight for each control.}
#' \item{L1_orig}{The L1 imbalance metric on the original data}
#' \item{L1_kbal}{The L1 imbalance metric on the weighted data}
#' \item{dist.record}{The record of L1 imbalances at each number of principap components of K balanced upon}
#' \item{K}{The kernel matrix used}
#' \item{pX_D0}{The estimated density measure for the control, as measured at each X-coordinate observed (for both treated and control units)}
#' \item{pX_D1}{The estimated density measure for the treated, as measured at each X-coordinate observed (for both treated and control units)}
#' \item{sigma}{The choice of kernel bandwidth used}
#' \item{biasbound}{The maximum bias that can be due to incomplete balancing, or not balancing at all. For a function with RKHS norm of gamma, the maximal bias possible due to incomplete balance is sqrt(gamma)*biasbound }
#' @examples #Run Lalonde example as in paper:
#' data(lalonde)
#' lalonde$nodegr=as.numeric(lalonde$educ<=11)
#' xvars=c("age","black","educ","hisp","married","re74","re75","nodegr","u74","u75")
#' attach(lalonde)
#'
#' #Raw diff-in-means: way off, -$15205
#' mean(re78[nsw==1])-mean(re78[nsw==0])
#'
#' #OLS with covariates:
#' summary(lm(re78~nsw+., data=lalonde[,xvars]))
#'
#' #Kbal at defaults: $1806
#' kbalout=kbal(D=nsw,X=lalonde[,xvars])
#' summary(lm(re78~nsw,w=kbalout$w))
#' plot(x=seq(2:41),kbalout$dist.record[2:41],
#' ylab="L1 imbalance", xlab="Num. dims of K balanced")
#'
#' #Kbal with mean balance ensured first, at defaults: $2455
#' kbalout_mean=kbal_mean(D=nsw,X=lalonde[,xvars])
#' summary(lm(re78~nsw,w=kbalout_mean$w))
#' plot(x=seq(2:41),kbalout$dist.record[2:41], pch=16,
#'     ylab="L1 imbalance", xlab="Num. dims of K balanced")
#' points(kbalout_mean$dist.record[2:41], col=2, pch=16)
#' legend("topright", col=c(1,2), pch=16, legend=c("full kbal","kbal after mean"))
#' @export

kbal=function(X,D, K=NULL, whiten=FALSE, trimratio=NULL, numdims=NULL,
          maxnumdims=NULL, minnumdims=NULL, sigma=NULL, method="ebal", linkernel=FALSE){
	N=dim(X)[1]
  P=dim(X)[2]
	X=as.matrix(X)

	if (method=="el") library(glmc)

	if (is.null(maxnumdims)) maxnumdims=sum(D==0)
	if (maxnumdims>sum(D==0)) maxnumdims=sum(D==0)
	if (is.null(minnumdims)){minnumdims=1}

	if (linkernel==TRUE){
	  X = scale(X, center = FALSE, scale = TRUE)
	  maxnumdims=ncol(X)
	  }

	#Option to pre-whiten X, as if using Mahalanobis distance in the kernel
	if (whiten){ X=X%*%solve(chol(var(X)))}


  if (is.null(sigma)){
		sigma=2*dim(X)[2]
	}

	if (linkernel==TRUE){
	 	K=X%*%t(X)
	}

	if (linkernel==FALSE){
	  X=scale(X, center=TRUE, scale=TRUE)
	  K=buildgauss(X,sigma=sigma)
	}

	#For readability, construct these first:
	K_c=K[,D==0]
	K_t=K[,D==1, drop=FALSE]
	K_t_bar=as.matrix(apply(K_t,1,mean))
	N_t=sum(D==1)
	N_c=sum(D==0)

	#Disabling this here - 15 Sept. I think it's elsewhere.
	#Pseudo-Density of treated, p(X|D=1), taken at all points
  #pX_D1=K_t%*%matrix(1/(sqrt(sigma*pi)*sum(D==1)),sum(D==1),1)
	#Pseudo-Density of controls, p(X|D=0), taken at all points
  #pX_D0=K_c%*%matrix(1/(sqrt(sigma*pi)*sum(D==0)),sum(D==0),1)
  #pX_D1=pX_D1/sum(pX_D1)
  #pX_D0=pX_D0/sum(pX_D0)

	badbals=NULL
	treatdrop=NULL

	if (!is.null(trimratio)) {
		badbals=which((pX_D1/pX_D0)>trimratio)
		treatdrop= intersect(badbals,which(D==1))
	}

	if (length(treatdrop)>0) {
		X=as.matrix(X[-treatdrop,])
		D=as.matrix(D[-treatdrop])
		K=K[-treatdrop,-treatdrop]
	}

  #Optional improving of rank by adding to diagonal. Set to 0 for now.
	lambda=0
  Klambda=K+diag(lambda,N)

  #SVD version: 11 Oct 2017
  #svd.out=svd(scale(Klambda,center = TRUE, scale=TRUE)) #Scaling messes it up.
  svd.out=svd(Klambda)

  # Could just use eigenvectors directly using svd.out$u and
  # balance on these. However, performance better when projecting K,
  # which is really just rescaling the eigenvector by their
  # eigenvalues. I think this helps ebal in case of minor imperfections.
  #Kpc=multdiag(svd.out$u, svd.out$d)
  Kpc=svd.out$u
  cum.var.pct=cumsum(svd.out$d)[1:N]/N

  #Eigen? 11 Oct 2017
  #eig.out=eigen(Klambda, symmetric = TRUE)
  #Kpc=eig.out$vectors
  #cum.var.pct=cumsum(eig.out$values)[1:N]/N

  #Prcomp -- used in drafts
  #prcomp.out=prcomp(Klambda, retx = TRUE)
	#Kpc=prcomp.out$x
  #cum.var.pct=cumsum(prcomp.out$sdev^2)[1:N]/sum(prcomp.out$sdev^2)

  if (is.null(numdims)){
    thisnumdims=minnumdims
    dist.record=NULL
    #rep(NA, N_c+1)
    keepgoing=TRUE
    wayover=FALSE
    mindistsofar=998
    dist.now=998

    while (keepgoing==TRUE){
      #keepgoing=(dist.now!=999) & thisnumdims<=maxnumdims & wayover==FALSE
      get.dist.out=get.dist(numdims=thisnumdims, D=D,
                            X=X, Kpc=Kpc, K=K, K_t=K_t, K_c=K_c,
                            method=method, treatdrop=treatdrop, linkernel=linkernel,
                            svd.out=svd.out)
      dist.now=get.dist.out$dist
      #if(thisnumdims==minnumdims){paste("Starting with ",P, "mean balance dims, plus...")}
      print(paste("Trying",thisnumdims,"dims of K; distance at", round(dist.now, 5)))
      dist.record=c(dist.record,dist.now)
      thisnumdims=thisnumdims+1

      if (dist.now<mindistsofar){mindistsofar=dist.now}

      wayover=(dist.now/mindistsofar)>1.25
      keepgoing=(dist.now!=999) & thisnumdims<=maxnumdims & wayover==FALSE
    }

    dimseq=seq(minnumdims,maxnumdims,1)
    numdims=dimseq[which(dist.record==min(dist.record,na.rm=TRUE))]
  }   #end for null numdims

	#get pctvarK
	pctvarK=cum.var.pct[numdims]

  #Recover optimal answer:
	#most of the goodies will be in here:
  best.out=get.dist(X=X, numdims = numdims, D=D, Kpc = Kpc, K=K, K_t=K_t,
                    K_c=K_c, method=method, treatdrop=treatdrop, linkernel=linkernel, svd.out)

  if(!is.null(numdims)){
    dist.record=best.out$dist
  }

  L1_orig=sum(abs(best.out$pX_D1-best.out$pX_D0)) #removing 0.5
	L1_kbal=best.out$L1  #.5*sum(abs(best.out$pX_D1-best.out$pX_D0w))
	L2_kbal=sqrt(sum((best.out$pX_D1-best.out$pX_D0w)^2)) #removing 0.5

	biasbound_orig = 	biasbound(D = D, w = rep(1,N), V = svd.out$u, a = svd.out$d)
	biasbound_kbal = best.out$biasbound

	r=list()
	r$K=K
	r$w=best.out$w
	#r$treatdrop=treatdrop
	r$L1_orig=L1_orig
	r$L1_kbal=L1_kbal
	r$L2_kbal=L2_kbal
	r$pX_D1=best.out$pX_D1
	r$pX_D0=best.out$pX_D0
	r$pX_D0w=best.out$pX_D0w
	r$pctvarK=pctvarK
	r$numdims=numdims
	r$sigma=sigma
	r$min90=min(which(cumsum(sort(best.out$w[D==0],decreasing=TRUE))/sum(D==0)>=.90))
  r$dist.record=dist.record[1:max(which(!is.na(dist.record)))]
  r$biasbound_orig=biasbound_orig
  r$biasbound_kbal=biasbound_kbal
	return(r)
}


#' Find weights and compute L1 distance.
#' @description  Get's the weights at the desired settings and computes
#' the objective function, L1.
#' @export
get.dist= function(numdims, D, Kpc, K, K_t, K_c, method, treatdrop, linkernel, X, svd.out, ...){
  R=list()
  K2=Kpc[,1:numdims, drop=FALSE]
  N=nrow(K2)
  if (method=="ebal"){bal.out.pc=try(ebal::ebalance(Treatment=as.vector(D),X=K2, print.level=-1),silent=TRUE)}

  if (method=="el"){
    yfake=rnorm(sum(D==0))
    #get mean row of K2 among the treated, which will be our target.
    meanK2tx= apply(K2[D==1, , drop=FALSE],2,mean)
    K2_0=K2[D==0,]
    z=t(t(K2_0)-meanK2tx)
    bal.out.pc=try(glmc::glmc(yfake~+1, Amat=z))
  }

  if (method=="el_custom"){
    meanK2tx= apply(K2[D==1, , drop=FALSE],2,mean)
    K2_0=K2[D==0, , drop=F]
    z=t(t(K2_0)-meanK2tx)
    bal.out.pc=try(el_custom(z=z, eps=1/N))
  }

  #if (class(bal.out.pc)=="try-error"){
  if ("try-error"%in%class(bal.out.pc)){
    dist=999
    #w=rep(1,N)
    R$dist=dist
  }
  if (class(bal.out.pc)[1]!="try-error"){
    w=rep(1,N)
    w[D==0]=bal.out.pc$w
    #rescale to mean=1 among the controls
    w[D==0]=w[D==0]/mean(w[D==0])
    w[treatdrop]=0

    if (linkernel==FALSE){
      pX_D1=K_t%*%matrix(1,sum(D==1),1)/sum(D==1)
      pX_D0=K_c%*%matrix(1,sum(D==0),1)/sum(D==0)
      pX_D0w=K_c%*%w[D==0]/sum(w[D==0])

      #Commenting out this rescaling, instead using a rescaling
      #above by the number of treated and control, and sum of weights
      #so that these are all like averages. This will look more
      #like the biasbound scaling. 20 Oct 2017
      pX_D1=pX_D1/sum(pX_D1)
      pX_D0=pX_D0/sum(pX_D0)
      pX_D0w=pX_D0w/sum(pX_D0w)
    }

    if (linkernel==TRUE){
      pX_D1=colMeans(X[D==1, , drop=FALSE])
      pX_D0=colMeans(X[D==0, , drop=FALSE])
      pX_D0w=w[D==0]%*%X[D==0,]/sum(D==0)
      L1=sum(abs(pX_D1-pX_D0w))
    }

    if (linkernel==FALSE){
      L1 = sum(abs(pX_D1-pX_D0w)) #removed the 0.5 factor -- Oct 20 2017
      biasbound = biasbound(D = D, w=w, V=svd.out$v, a = svd.out$d, hilbertnorm = 1)
    }

    R$dist= biasbound  ##experimenting with using biasbound instead of L1
    R$L1=L1
    R$biasbound=biasbound
    R$w=w
    R$pX_D1=pX_D1
    R$pX_D0=pX_D0
    R$pX_D0w=pX_D0w
  }
  return(R)
} ## end of get.dist

#--------------------------------------------

#' Experimental KBAL that insists on mean balance first.
#' @description This is an experimental version of KBAL that begins by achieving
#' exact mean balance on X, before adding components from K as per the usual KBAL.
#' The multivariate distance metrics are still based on balance in K, not the mean balance.
#' The mean balance must be made exact first, then it attempts to add components of K until optimized.
#' @param X The original covariate data, as a numeric matrix.
#' @param D The treatment assignment variable taking values of 1 for treatet units and 0 for control units.
#' @param K Optional user-provided kernel matrix. Typically this is not user-specified, but rather is computed internally by a call to \code{buildgauss}.
#' @param whiten Optional pre-whitening of the data prior to construction of K. If used, rotates the data by \code{solve(chol(var(X)))}, then centers and rescales.
#' @param trimratio Optional \code{logical}
#' @param numdims Optional user-specified number of projectionss of \code{K} to balance upon.
#' @param minnumdims Optional user-specified choice for the minimum number of projections of \code{K}. Defualts to 2.
#' @param maxnumdims Optional user-specified choice for the maximum number of projectsion ffo \code{K} to attempt balance on. Defaults to the number of control units.
#' @param sigma Optional user-specificied paramater for the Gaussian kernel. If blank, defaults to \code{nrow(X)}.
#' @param method "ebal" or "el". Whether balance should be obtained on each projection of \code{K} using entropy balancing ("ebal", default) or empirical likelihood ("el")
#' @return \item{w}{The weights, taking value of 1 for treated unit and the estimated weight for each control.}
#' \item{L1_orig}{The L1 imbalance metric on the original data}
#' \item{L1_kbal}{The L1 imbalance metric on the weighted data}
#' \item{dist.record}{The record of L1 imbalances at each number of principap components of K balanced upon}
#' \item{K}{The kernel matrix used}
#' \item{pX_D0}{The estimated density measure for the control, as measured at each X-coordinate observed (for both treated and control units)}
#' \item{pX_D1}{The estimated density measure for the treated, as measured at each X-coordinate observed (for both treated and control units)}
#' \item{sigma}{The choice of kernel bandwidth used}
#' @export
# kbal_meanfirst=function(X,D, K=NULL, whiten=FALSE, trimratio=NULL,numdims=NULL,maxnumdims=NULL,minnumdims=NULL,sigma=NULL, method="ebal"){
#   N=dim(X)[1]
#   P=dim(X)[2]
#   X=as.matrix(X)
#   X=scale(X, center=TRUE, scale=TRUE)
#
#   #Option to pre-whiten X, as if using Mahalanobis distance in the kernel
#   if (whiten){ X=X%*%solve(chol(var(X)))}
#
#   if (is.null(sigma)){
#     sigma=2*dim(X)[2]
#   }
#
#   K=buildgauss(X,sigma=sigma)
#
#   #For readability, construct these first:
#   K_c=K[,D==0]
#   K_t=K[,D==1]
#   K_t_bar=as.matrix(apply(K_t,1,mean))
#   N_t=sum(D==1)
#   N_c=sum(D==0)
#
#   #Pseudo-Density of treated, p(X|D=1), taken at all points
#   pX_D1=K_t%*%matrix(1/(sqrt(sigma*pi)*sum(D==1)),sum(D==1),1)
#   #Pseudo-Density of controls, p(X|D=0), taken at all points
#   pX_D0=K_c%*%matrix(1/(sqrt(sigma*pi)*sum(D==0)),sum(D==0),1)
#
#   pX_D1=pX_D1/sum(pX_D1)
#   pX_D0=pX_D0/sum(pX_D0)
#
#   badbals=NULL
#   treatdrop=NULL
#
#   if (!is.null(trimratio)) {
#     badbals=which((pX_D1/pX_D0)>trimratio)
#     treatdrop= intersect(badbals,which(D==1))
#   }
#
#   if (length(treatdrop)>0) {
#     X=as.matrix(X[-treatdrop,])
#     D=as.matrix(D[-treatdrop])
#     K=K[-treatdrop,-treatdrop]
#   }
#
#   #Return to these options later...
#   lambda=0
#   Kreg=K+diag(lambda,N)
#   prcomp.out=prcomp(Kreg, retx = TRUE)
#   Kpc=prcomp.out$x
#
#   #Kpc=K2%*%prcomp.out$rotation
#   cum.var.pct=cumsum(prcomp.out$sdev^2)[1:N]/sum(prcomp.out$sdev^2)
#
#   #Function to return distance for given numdims balanced upon. And everything else.
#   get.dist= function(numdims, X, Kpc, K, ...){
#     R=list()
#     XKpc=cbind(X,Kpc)
#     trydims=P+numdims
#     K2=XKpc[,1:trydims]
#     if (method=="ebal"){bal.out.pc=try(ebal::ebalance(Treatment=as.vector(D),X=K2, print.level=-1),silent=TRUE)}
#
#     if (method=="el"){
#       yfake=rnorm(sum(D==0))
#       #get mean row of K2 among the treated, which will be our target.
#       meanK2tx=apply(K2[D==1,],2,mean)
#       K2_0=K2[D==0,]
#       z=t(t(K2_0)-meanK2tx)
#       bal.out.pc=try(glmc::glmc(yfake~+1, Amat=z))
#     }
#
#     #if (class(bal.out.pc)=="try-error"){
#     if ("try-error"%in%class(bal.out.pc)){
#       dist=999
#       #w=rep(1,N)
#       R$dist=dist
#     }
#     if (class(bal.out.pc)[1]!="try-error"){
#       w=rep(1,N)
#       w[D==0]=bal.out.pc$w
#       #rescale to mean=1 among the controls
#       w[D==0]=w[D==0]/mean(w[D==0])
#       w[treatdrop]=0
#
#       pX_D1=K_t%*%matrix(1,sum(D==1),1)
#       pX_D0=K_c%*%matrix(1,sum(D==0),1)
#       pX_D0w=K_c%*%w[D==0]
#
#       pX_D1=pX_D1/sum(pX_D1)
#       pX_D0=pX_D0/sum(pX_D0)
#       pX_D0w=pX_D0w/sum(pX_D0w)
#
#       L1=.5*sum(abs(pX_D1-pX_D0w))
#       dist=L1
#       R$dist=dist
#       R$w=w
#       R$pX_D1=pX_D1
#       R$pX_D0=pX_D0
#       R$pX_D0w=pX_D0w
#     }
#     return(R)
#   } ## End of get.dist function
#
#
#   if (is.null(numdims)){
#     if (is.null(minnumdims)){minnumdims=0}
#     thisnumdims=minnumdims
#     dist.record=NULL
#     #rep(NA, N_c+1)
#     if (is.null(maxnumdims)){maxnumdims=N_c}
#     keepgoing=TRUE
#     wayover=FALSE
#     mindistsofar=1
#     dist.now=1
#
#     while (keepgoing==TRUE){
#       #keepgoing=(dist.now!=999) & thisnumdims<=maxnumdims & wayover==FALSE
#       get.dist.out=get.dist(numdims=thisnumdims,
#                             X=X, Kpc=Kpc, K=K)
#       dist.now=get.dist.out$dist
#       if(thisnumdims==minnumdims){paste("Starting with ",P, "mean balance dims, plus...")}
#       print(paste(thisnumdims,"dims of K; L1 dist. at", round(dist.now, 5)))
#       dist.record=c(dist.record,dist.now)
#       thisnumdims=thisnumdims+1
#
#       if (dist.now<mindistsofar){mindistsofar=dist.now}
#       wayover=(dist.now/mindistsofar)>2
#       keepgoing=(dist.now!=999) & thisnumdims<=maxnumdims & wayover==FALSE
#     }
#
#      dimseq=seq(minnumdims,maxnumdims,1)
#      numdims=dimseq[which(dist.record==min(dist.record,na.rm=TRUE))]
#   }   #end for null numdims
#
#   #get pctvarK
#
#   pctvarK=cum.var.pct[numdims]
#
#   #Recover optimal answer:
#   #most of the goodies will be in here:
#   best.out=get.dist(numdims=numdims, X=X, Kpc=Kpc, K=K)
#
#   L1_orig=.5*sum(abs(best.out$pX_D1-best.out$pX_D0))
#   L1_kbal=.5*sum(abs(best.out$pX_D1-best.out$pX_D0w))
#
#   r=list()
#   r$K=K
#   #r$dist=best.out$dist
#   r$w=best.out$w
#   #r$treatdrop=treatdrop
#   r$L1_orig=L1_orig
#   r$L1_kbal=L1_kbal
#   r$pX_D1=best.out$pX_D1
#   r$pX_D0=best.out$pX_D0
#   r$pX_D0w=best.out$pX_D0w
#   r$pctvarK=pctvarK
#   r$numdims=numdims
#   r$sigma=sigma
#   r$min90=min(which(cumsum(sort(best.out$w[D==0],decreasing=TRUE))/sum(D==0)>=.90))
#   r$dist.record=dist.record[1:max(which(!is.na(dist.record)))]
#   return(r)
# }
