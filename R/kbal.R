require(ebal)

#' Build Gaussian kernel matrix.
#' Centers and rescales X then computes Guassian
#' kernel matrix. The $j^{th}$ row and $i^th$
#' @examples
#' K=buildgauss(X,sigma=ncol(X))
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

#' Main kernel balancing function.
#' @X love Do you love cats? Defaults to TRUE.
#' @keywords cats
#' @export
#' @examples
#' cat_function()
kbal=function(X,D, K=NULL, whiten=TRUE, trimratio=NULL,numdims=NULL,maxnumdims=NULL,minnumdims=NULL,sigma=NULL){
	N=dim(X)[1]
  P=dim(X)[2]
	X=as.matrix(X)

	#Rather than Mahalanobis distance in Gaussian kernel, give option to pre-whiten X
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

  cum.var.pct=cumsum(prcomp.out$sdev^2)[1:N]/sum(prcomp.out$sdev^2)
  #plot(cum.var.pct)
  #abline(v=min(which(cum.var.pct>.99)))

  #eig=eigen(K)
	#prcomp.out=list()
	#prcomp.out$sdev=sqrt(eig$values)
	#prcomp.out$rotation=eig$vectors



  #Function to return distance for given numdims balanced upon. And everything else.
  get.dist= function(numdims, Kpc, K, ...){
    K2=Kpc[,1:numdims]

	  ebal.out.pc=try(ebalance(Treatment=as.vector(D),X=K2, print.level=-1),silent=TRUE)

    if (class(ebal.out.pc)=="try-error") {dist=999}
	  if (class(ebal.out.pc)!="try-error"){
	    w=rep(1,N)
	    w[D==0]=ebal.out.pc$w
      #rescale to mean=1 among the controls
      w[D==0]=w[D==0]/mean(w[D==0])
	    w[treatdrop]=0


	    pX_D1=K_t%*%matrix(1,sum(D==1),1)
	    pX_D0=K_c%*%matrix(1,sum(D==0),1)
	    pX_D0w=K_c%*%w[D==0]

      pX_D1=pX_D1/sum(pX_D1)
	    pX_D0=pX_D0/sum(pX_D0)
	    pX_D0w=pX_D0w/sum(pX_D0w)

	    #Try a new distance relating to how far pscores are from flat.
      dist= mean((pX_D1/pX_D0w-1)^2)
      #L1=.5*sum(abs(density_treat-density_weighted_control))
      #dist=L1
	  }

    R=list()
    R$dist=dist
    R$w=w
    R$pX_D1=pX_D1
    R$pX_D0=pX_D0
    R$pX_D0w=pX_D0w
	  return(R)
	}

	if (is.null(numdims)){
    if (is.null(minnumdims)){minnumdims=2}

    thisnumdims=2
		dist.record=rep(NA, N)
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

		numdims=which(dist.record==min(dist.record,na.rm=T))
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
  r$dist.record=dist.record
	return(r)
}

#--------------------------------------------

