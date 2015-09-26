require(ebal, glmc)
#' Build Gaussian kernel matrix.
#' @description  Centers and rescales X then computes Guassian
#' kernel matrix. Entry {i,j} correspond to k(X_i,X_j) where k is the (Gaussian) kernel.
#' @param X A numeric matrix of data, typically a design matrix.
#' @param sigma The kernel ``bandwidth''. If NULL, defaults to ncol(X)
#' @examples
#' K=buildgauss(rnorm(10))
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

#' Kernel balancing function.
#' @description Chooses weights on control units that produces equal means on a kernel matrix, K, rather than the original design matrix, X.
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
kbal=function(X,D, K=NULL, whiten=FALSE, trimratio=NULL,numdims=NULL,maxnumdims=NULL,minnumdims=NULL,sigma=NULL, method="ebal"){
	N=dim(X)[1]
  P=dim(X)[2]
	X=as.matrix(X)

	#Option to pre-whiten X, as if using Mahalanobis distance in the kernel
	if (whiten){ X=X%*%solve(chol(var(X)))}

	X=scale(X, center=TRUE, scale=TRUE)

  if (is.null(sigma)){
		sigma=2*dim(X)[2]
	}

	K=buildgauss(X,sigma=sigma)

	#For readability, construct these first:
	K_c=K[,D==0]
	K_t=K[,D==1]
	K_t_bar=as.matrix(apply(K_t,1,mean))
	N_t=sum(D==1)
	N_c=sum(D==0)

	#Pseudo-Density of treated, p(X|D=1), taken at all points
  pX_D1=K_t%*%matrix(1/(sqrt(sigma*pi)*sum(D==1)),sum(D==1),1)
	#Pseudo-Density of controls, p(X|D=0), taken at all points
  pX_D0=K_c%*%matrix(1/(sqrt(sigma*pi)*sum(D==0)),sum(D==0),1)

  pX_D1=pX_D1/sum(pX_D1)
  pX_D0=pX_D0/sum(pX_D0)

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

  #Return to these options later...
	lambda=0
  K2=K+diag(lambda,N)
  prcomp.out=prcomp(K2, retx = TRUE)
	Kpc=prcomp.out$x
  #Kpc=K2%*%prcomp.out$rotation
  cum.var.pct=cumsum(prcomp.out$sdev^2)[1:N]/sum(prcomp.out$sdev^2)

  #Function to return distance for given numdims balanced upon. And everything else.
  get.dist= function(numdims, Kpc, K, ...){
    R=list()
    K2=Kpc[,1:numdims]
    if (method=="ebal"){bal.out.pc=try(ebalance(Treatment=as.vector(D),X=K2, print.level=-1),silent=TRUE)}

    if (method=="el"){
      yfake=rnorm(sum(D==0))
      #get mean row of K2 among the treated, which will be our target.
      meanK2tx=apply(K2[D==1,],2,mean)
      K2_0=K2[D==0,]
      z=t(t(K2_0)-meanK2tx)
      bal.out.pc=try(glmc(yfake~+1, Amat=z))
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

	    pX_D1=K_t%*%matrix(1,sum(D==1),1)
	    pX_D0=K_c%*%matrix(1,sum(D==0),1)
	    pX_D0w=K_c%*%w[D==0]

      pX_D1=pX_D1/sum(pX_D1)
	    pX_D0=pX_D0/sum(pX_D0)
	    pX_D0w=pX_D0w/sum(pX_D0w)

	    L1=.5*sum(abs(pX_D1-pX_D0w))
      dist=L1
      R$dist=dist
      R$w=w
      R$pX_D1=pX_D1
      R$pX_D0=pX_D0
      R$pX_D0w=pX_D0w
	  }
	  return(R)
	}

	if (is.null(numdims)){
    if (is.null(minnumdims)){minnumdims=2}
    thisnumdims=2
		dist.record=rep(NA, N_c)
		if (is.null(maxnumdims)){maxnumdims=N_c}
    keepgoing=TRUE
		wayover=FALSE
		mindistsofar=1
		dist.now=1

    while (keepgoing==TRUE){
			#keepgoing=(dist.now!=999) & thisnumdims<=maxnumdims & wayover==FALSE
			get.dist.out=get.dist(numdims = thisnumdims,Kpc = Kpc, K=K)
			dist.now=get.dist.out$dist
      print(dist.now)
			dist.record[thisnumdims]=dist.now
			thisnumdims=thisnumdims+1

			if (dist.now<mindistsofar){mindistsofar=dist.now}
			wayover=(dist.now/mindistsofar)>2
			keepgoing=(dist.now!=999) & thisnumdims<=maxnumdims & wayover==FALSE
    }
		numdims=which(dist.record==min(dist.record,na.rm=TRUE))
	}   #end for null numdims

	#get pctvarK

	pctvarK=cum.var.pct[numdims]

  #Recover optimal answer:
	#most of the goodies will be in here:
  best.out=get.dist(numdims = numdims,Kpc = Kpc, K=K)

	L1_orig=.5*sum(abs(best.out$pX_D1-best.out$pX_D0))
	L1_kbal=.5*sum(abs(best.out$pX_D1-best.out$pX_D0w))

	r=list()
	r$K=K
  r$dist=best.out$dist
	r$w=best.out$w
	#r$treatdrop=treatdrop
	r$L1_orig=L1_orig
	r$L1_kbal=L1_kbal
	r$pX_D1=best.out$pX_D1
	r$pX_D0=best.out$pX_D0
	r$pX_D0w=best.out$pX_D0w
	r$pctvarK=pctvarK
	r$numdims=numdims
	r$sigma=sigma
	r$min90=min(which(cumsum(sort(best.out$w[D==0],decreasing=TRUE))/sum(D==0)>=.90))
  r$dist.record=dist.record[1:max(which(!is.na(dist.record)))]
	return(r)
}

#--------------------------------------------

