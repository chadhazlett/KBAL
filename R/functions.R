### Helper functions

### Build kernel, possibly non-square
#' Build the Kernel Matrix
#' 
#' @description Builds the kernel matrix using Rcpp.
#' @param allx A data matrix containing all observations where rows are units and columns are covariates.
#' @param useasbases Binary vector argument to specify what bases to use when constructing the kernel matrix and finding weights. While the number of observations is under 2000, the default maximum is to use all observations. Due to the computation burden, when the number of observations is over 2000, the default is to use sampled units.
#' @param b Scaling factor in the calculation of gaussian kernel distance equivalent to the entire denominator \eqn{2\sigma^2} of the exponent. Default is twice the number of covariates or columns in \code{allx}.
#' @return \item{K}{The Kernel matrix}
#' @examples
#' #load and clean data a bit
#' data(lalonde)
#' lalonde$nodegr=as.numeric(lalonde$educ<=11)
#' xvars=c("age","black","educ","hisp","married","re74","re75","nodegr","u74","u75")
#' 
#' #note that lalonde$nsw is the treatment vector, so the observered is 1-lalodne$nsw
#' #running makeK with the sampled units as the bases given the large size of the data
#' #and with b as twice the number of covariates
#' K = makeK(allx = lalonde[,xvars], 
#' useasbases = as.numeric(1-lalonde$nsw), 
#' b = 2*ncol(lalonde[,xvars]))
#' @useDynLib KBAL
#' @importFrom Rcpp sourceCpp
#' @export
makeK = function(allx, useasbases=NULL, b=NULL){
  N=nrow(allx)
  
  # If no "useasbasis" given, assume all observations are to be used.
  #default b is set to 2ncol to match kbal for now
  if (is.null(b)){ b=2*ncol(allx) }
  bases = allx[useasbases==1, ]
  
  Xmeans <- colMeans(bases)
  Xsds <- apply(bases,2,sd)
  bases <- scale(bases, center = Xmeans, scale = Xsds)
  allx <- scale(allx, center = Xmeans, scale = Xsds)
  
  K = new_gauss_kern(newx = allx, oldx = bases, b = b) 
  return(K)
  #K=KRLS2::newKernel(X = bases , newData = allx , b = b) #if use this then do not rescale above
}

### Get bias bound, which will be the distance
### function to drive the optimization
#similar to biasboud in kbal
#' Get the bound on the bias due to incomplete balance
#' @description Calculate the upper bound on the bias induced by approximate balance  
#' @param observed a numeric vector of length equal to the total number of units where sampled units take a value of 1 and population units take a value of 0. 
#' @param target a numeric vector of length equal to the total number of units where population units take a value of 1 and sample units take a value of 0.
#' @param svd.out the object output from \code{\link{svd}} performed on the kernel matrix.
#' @param w numeric vector of length equal to the total number of units containing the weight for each corresponding unit. Note that these weights should sum to the total number of units not to one. They are divided by the number of units internally within the function. 
#' @param hilbertnorm numeric value of the hilbertnorm.
#' @examples
#' #load and clean data a bit
#' data(lalonde)
#' lalonde$nodegr=as.numeric(lalonde$educ<=11)
#' xvars=c("age","black","educ","hisp","married","re74","re75","nodegr","u74","u75")
#' 
#' #need a kernel matrix to run SVD on and pass in so get that first with makeK
#' #running makeK with the sampled units as the bases
#' K = makeK(allx = lalonde[,xvars], 
#' useasbases = as.numeric(1-lalonde$nsw), 
#' b = 2*ncol(lalonde[,xvars]))
#' 
#' #svd on this kernel
#' svd_pass = svd(K)
#' #let's use the original weights of 1/number of sampled units, and 1/numer of target units
#' #this is the default if we pass in w as all 1's
#' biasbound(observed=(1-lalonde$nsw),
#'  target=lalonde$nsw, 
#'  svd.out = svd_pass, 
#'  w = rep(1,nrow(lalonde)), hilbertnorm=1)
#' @export
biasbound=function(observed, target, svd.out, w, hilbertnorm=1){
  wtarget=w[target==1]/sum(target==1)
  wobserved=w[observed==1]/sum(observed==1)

  V=svd.out$u
  eigenvals=svd.out$d

  V1=V[target==1, , drop=FALSE]
  V0=V[observed==1, , drop=FALSE]
  eigenimbal=as.vector(t(wtarget)%*%V1 - t(wobserved)%*%V0)

  effectiveimbal=(eigenimbal*(eigenvals^.5)) #%*%t(V)
  biasbound=sqrt(hilbertnorm)*sqrt(t(effectiveimbal)%*%(effectiveimbal))
  return(biasbound)
}

# Simple diffference in mean and in weighted means
# (Not actually used right now but can be convenient)

#' Difference in means and Difference in weighted means
#'
#' Calcuates the simple difference in means or weighted difference in means between the treated/sample population and the control/target popultion.
#' @param X matrix of data where rows are observations and columns are covariates.
#' @param w numeric vector of weights for each observation.
#' @param target numeric vector of length equal to the total number of units where population units take a value of 1 and sample units take a value of 0. 
#'
#' @param X matrix of combined treated/sample population and control/target population data. Rows are observations, columns are covariates.
#' @param w numeric vector of weights for every obervation
#' @param target vector taking values of 0 or 1's indicating which observations (whuch rows  of X and w) are in the treated/sample population and which are in the control/target population.
#' @return \code{dim} the simple, unweighted difference in means
#' @return \code{dimw} the weighted difference in means
#' @examples 
#' #let's say we want to get the unweighted DIM and the weighted DIM using weights from the kpop
#' #function with the lalonde data:
#' #load and clean data a bit
#' data(lalonde)
#' lalonde$nodegr=as.numeric(lalonde$educ<=11)
#' xvars=c("age","black","educ","hisp","married","re74","re75","nodegr","u74","u75")
#' 
#' #get the kpop weights
#' kpopout= kpop(allx=lalonde[,xvars],
#'                useasbases=NULL, b=NULL,
#'                sampled=NULL, sampledinpop=FALSE,
#'                treatment=lalonde$nsw,
#'                ebal.tol=1e-6, numdims=NULL, 
#'                minnumdims=NULL, maxnumdims=NULL, 
#'                incrementby=1,
#'                printprogress =TRUE)
#'  #now use dimw to get the DIMs
#'  dimw(X = lalonde[,xvars], w = kpopout$w, target = lalonde$nsw)
#' @export
dimw = function(X,w,target){
  w1=w[target==1]/sum(w[target==1])
  w0=w[target!=1]/sum(w[target!=1])

  X1=as.matrix(X[target==1, , drop=FALSE])
  X0=as.matrix(X[target!=1, , drop=FALSE])

  R=list()
  R$dim=colMeans(X1)-colMeans(X0)
  R$dimw=t(as.vector(w1))%*%X1-t(w0)%*%X0
  return(R)
}

# Function to get the moments that solve the constraints
# we have setup.
# Currently just uses ebal, but we could to others
# or give options.
# Will need some work to better carry error messages from ebal.

#' Find weights using entropy balancing.
#' @description Using entropy balancing to find and return the weights such that XXXX
#' 
#' @param target a numeric vector of length equal to the total number of units where population units take a value of 1 and sample units take a value of 0. 
#' @param observed a numeric vector of length equal to the total number of units where sampled units take a value of 1 and population units take a value of 0. 
#' @param allrows Matrix whose columns contain the left singular vectors of the kernel matrix 
#' @param ebal.tol tolerance level used by \code{ebalance}
#' @return \item{w}{numeric vector of weights}
#' @examples
#' #load and clean data a bit
#' data(lalonde)
#' lalonde$nodegr=as.numeric(lalonde$educ<=11)
#' xvars=c("age","black","educ","hisp","married","re74","re75","nodegr","u74","u75")
#' 
#' #need a kernel matrix to run SVD on then find weights with so get that first with makeK
#' #running makeK with the sampled units as the bases
#' K = makeK(allx = lalonde[,xvars], 
#' useasbases = as.numeric(1-lalonde$nsw), 
#' b = 2*ncol(lalonde[,xvars]))
#' 
#' #svd on this kernel and get matrix with left singular values
#' Kpc = svd(K)$u
#' #usually we are getting weights using different number of columns of this matrix, the finding
#' # the bias and looking for the minimum. For now let's just use the first 10
#' Kpc2=Kpc[,1:10, drop=FALSE]
#' getw.out=getw(target=lalonde$nsw, observed=1-lalonde$nsw, allrows=Kpc2)
#' @export
getw = function(target, observed, allrows, ebal.tol=1e-6){

  # To trick ebal into using a control group that corresponds to the
  # observed and a treated that corresponds to the "target" group,
  # (1) anybody who is "observed" but also considered part of the target
  # group has to get a duplicate row in the data
  # (2) construct a variables target_ebal that = 1 when target=1
  # but =0 in the appended data, i.e. for the observed who are
  # entering a second time.
  Xappended = rbind(allrows,  allrows[observed==1 & target==1, , drop=FALSE] )
  target_ebal = c(target, rep(0, sum(observed==1 & target==1)))

    bal.out.pc=try(ebal::ebalance(Treatment=target_ebal,X=Xappended,
        constraint.tolerance=ebal.tol, print.level=-1),
        silent=TRUE)
  N=nrow(allrows)
  
  if ("try-error"%in%class(bal.out.pc)){
      if(ncol(allrows) <= 2) {
          stop("ebalance convergence failed within first two dimensions")
      }
    w=rep(1,N)
    #R$dist="ebalerror"
  }

  if (class(bal.out.pc)[1]!="try-error"){
    w=rep(1,N)
    w[observed==1]=bal.out.pc$w
    #rescale to mean=1 among the donors
    w[observed==1]=w[observed==1]/mean(w[observed==1])
    #biasbound.out = biasbound(D = D, w=w, V=svd.out$v, a = svd.out$d, hilbertnorm = 1)
    #R$dist= biasbound.out  ##experimenting with using biasbound instead of L1
    #R$biasbound = biasbound.out
  }
  return(w)
} # end of getw.


# The main event: Actual kpop function!
#' Kernel Balancing
#' 
#' @description The kernel balancing function. XXXXX
#' 
#' @param allx a data matrix containing all observations where rows are units and columns are covariates.
#' @param useasbases optional binary vector argument to specify what bases to use when constructing the kernel matrix and finding weights. While the number of observations is under 2000, the default maximum is to use all observations. Due to the computation burden, when the number of observations is over 2000, the default is to use sampled units.}
#' @param b scaling factor in the calculation of gaussian kernel distance equivalent to the entire denominator \eqn{2\sigma^2} of the exponent. 
#' @param sampled a numeric vector of length equal to the total number of units where sampled units take a value of 1 and population units take a value of 0.
#' @param sampledinpop a logical to be used in compination with input \code{sampled} that when \code{TRUE} indicates that sampled units should be included in the target popultion XXX.
#' @param treatment an alternative input to \code{sampled} and \code{sampledinpop} that is a numeric vector of length equal to the total number of units where population units take a value of 1 and sample units take a value of 0. Note that this input corresponds to \code{sampledinpop} as \code{FALSE} only XXX.
#' @param ebal.tol tolerance level used by \code{ebalance}
#' @param numdims optional numeric argument to specify the number of dimensions of the kernel matrix to find balance on rather than searching for the number of dimensions which minimize the bias.
#' @param minnumdims optional numeric argument to specify the minimum number of dimensions of the SVD of the kernel matrix to find balance on in the search for the number of dimesions which minimize the bias. Default minimum is 1.
#' @param maxnumdims optional numeric argument to specify the maximum number of dimensions of the SVD of the kernel matrix to find balance on in the search for the number of dimesions which minimize the bias. While the number of observations is under 2000, the default maximum is the total number of observations. Due to the computation burden, when the number of observations is over 2000, the default is the number of sampled units. 
#' @param incrementby optional argument to specify the number of dimesions to increase by from \code{minnumdims} to \code{maxnumdims} in each iteration of the search for the number of dimensions which minimizes the bias. Default is 1.
#' @param printprogress optional logical argument to print current number of dimensions and bias. 
#' 
#' @return \item{dist.record}{a numeric matrix recording the bias bound corresponding to balance on increasing dimesions of the SVD of the kernel matrix starting from \code{minnumdims} increasing by \code{incrementby} to \code{maxnumdims} or until the bias grows to be 1.25 times the minimal bias found.}
#'  \item{biasbound.orig}{the bias bound found when all sampled units have a weight of one over the number of sampled units and all target units have a weight of one over the number of target units.}
#'  \item{numdims}{the optimal number of dimensions of the SVD of the kernel matrix which minimizes the bias bound.}
#'  \item{w}{the weights found using entropy balancing on \code{numdims} dimensions of the SVD of the kernel matrix}
#'  \item{biasbound.opt}{the minimal bias bound found using \code{numdims} as the number of dimestions of the SVD of the kernel matrix. When \code{numdims} is user-specified, the bias bound using this number of dimensions of the kernel matrix.}
#' \item{K}{the kernel matrix.}
#' @examples 
#' #Run Lalonde example as in paper:
#' data(lalonde)
#' lalonde$nodegr=as.numeric(lalonde$educ<=11)
#' xvars=c("age","black","educ","hisp","married","re74","re75","nodegr","u74","u75")
#' #kpop at defaults: 
#' kpopout= kpop(allx=lalonde[,xvars],
#'                useasbases=NULL, b=NULL,
#'                sampled=NULL, sampledinpop=FALSE,
#'                treatment=lalonde$nsw,
#'                ebal.tol=1e-6, numdims=NULL, 
#'                minnumdims=NULL, maxnumdims=NULL, 
#'                incrementby=1,
#'                printprogress =TRUE)
#'  summary(lm(re78~nsw,w=kpopout$w, data = lalonde))
#' @export
kpop = function(allx, useasbases=NULL, b=NULL, 
                sampled=NULL, sampledinpop=NULL,
                treatment=NULL,
                ebal.tol=1e-6, numdims=NULL, 
                minnumdims=NULL, maxnumdims=NULL, 
                incrementby=1, 
                printprogress = TRUE){

    
  #need to throw error if try to pass both sample and target
  if(!is.null(sampled) & !is.null(treatment)) {
       stop("\"sampled\" and \"treatment\" arguments can not be specified simultaneously")
  }

  #throw warning for using default sampledinpop=TRUE
  if(is.null(sampledinpop)) {
      warning("using default parameter \"sampledinpop\" = TRUE", immediate. = TRUE)
      sampledinpop = TRUE
  }

  N=nrow(allx)
  user_numdims = NULL
  if(!is.null(sampled) & sampledinpop==FALSE) {
      observed = sampled
      target = 1-sampled
  } else if(!is.null(sampled)) {
      observed = sampled
      target = rep(1,N)
  } else if(is.null(treatment)) { #error for passing in both as null
      stop("either the \"sampled\" or \"treatment\" argument is required")
  } else{
      observed = 1-treatment
      target = treatment
      #adding warning that sampledinpop==FALSE with treatment passed in
      #NB: this may change later if we add ATE option
      if(sampledinpop == TRUE) { 
          warning("\"treatment\" input requires \"sampledinpop\" to be FALSE. Proceeding with \"sampledinpop\" as FALSE.", immediate. = TRUE) 
          }
  }
  
  #error catch for passing in all 1's for either sampled or treated
  if(sum(observed != 1) == 0) {
      stop("\"sampled\" contains only sampled units (all entries are 1)")
  } else if(sum(target != 1) == 0 ){
      stop("\"target\" contains only target units (all entries are 1)")
  }
  
  #error catch for covariates with no variance:
  if(0 %in% apply(allx, 2, sd)) {
       stop("One or more column in \"allx\" have zero variance")
  }
  
  #error catch for NAs in data
  if(sum(is.na(allx) != 0)) {
      stop("\"allx\" contains missing values")
  }
  

  # If we don't specify which observations to use as bases, 
  # use all as default unless K is very large, then use sample set. 
  if (is.null(useasbases) & N <= 2000){
      useasbases = rep(1,N) 
  } else if(is.null(useasbases)) {
      warning("Dimensions of K greater than 2000, using sampled as default bases", 
              immediate. = TRUE)
      useasbases = as.numeric(observed==1)
  }

  # If we don't specify which observations to use as bases, 
  # use all as default unless K is very large, then use sample set. 
  if (is.null(useasbases) & N <= 2000){
      useasbases = rep(1,N) 
  } else if(is.null(useasbases)) {
      warning("Dimensions of K greater than 2000, using sampled as default bases",
              immediate. = TRUE)
      useasbases = as.numeric(observed==1)
  }

  if (is.null(minnumdims)){minnumdims=1}
  
  if(!is.null(minnumdims)) {
      if(minnumdims == 0) {
          stop("Minimum number of dimensions specified must be greater than zero")
      }
  }

  if(!is.null(numdims)) {
      if(numdims == 0) {
          stop("Number of dimensions specified must be greater than zero")
      }
  }
  
  
  # The most dims you can use is the number of bases
  if (!is.null(maxnumdims) && maxnumdims>sum(useasbases)){
    warning("Cannot allow dimensions of K to be greater than the number of bases. Reducing maxnumdims.", immediate. = TRUE)
    maxnumdims=sum(useasbases)
    }
  if (is.null(maxnumdims)){maxnumdims=sum(useasbases)}

  #adding default b within the kpop function rather than in makeK
  #changing default to be 2*ncol to match kbal
  if (is.null(b)){ b = 2*ncol(allx) }
  
  if(printprogress == TRUE) {print(paste0("Building Kernel matrix"))}
  K = makeK(allx = allx, useasbases = useasbases, b=b)
  if(printprogress == TRUE) {print(paste0("Running SVD on Kernel matrix"))}
  svd.out=svd(K)
  Kpc=svd.out$u

  # Get biasbound with no improvement in balance:
  biasbound_orig=biasbound(w = rep(1,N), observed=observed, target = target, 
                            svd.out = svd.out, hilbertnorm = 1)
  paste0("Without balancing, biasbound (norm=1) is ", round(biasbound_orig,3))

  # If numdims given, just get the weights in one shot:
  if (!is.null(numdims)){
    Kpc2=Kpc[,1:numdims, drop=FALSE]
    getw.out=getw(target=target, observed=observed, allrows=Kpc2)
    # XXX This would be place to add check for non-convergence of ebal.
    w=getw.out
    biasboundnow=biasbound( w = w, observed=observed,  target = target, 
                            svd.out = svd.out, hilbertnorm = 1)
    print(paste0("With ",numdims," dimensions, biasbound (norm=1) of ", 
                 round(biasboundnow,3)))
    user_numdims = TRUE
    }

  # If numdims not given, we search to minimize biasbound:
  if (is.null(numdims)){
    thisnumdims=minnumdims
    dist.record=NULL
    #rep(NA, N_c+1)
    keepgoing=TRUE
    wayover=FALSE
    mindistsofar=998

    while (keepgoing==TRUE){
      Kpc_try=Kpc[,1:thisnumdims, drop=FALSE]
      getw.out=getw(target = target, observed=observed, allrows = Kpc_try)
      w=getw.out
      # Need to work on case where ebal fails and flagging this in result.
      # For now just returns all even weights.

      biasboundnow=biasbound( w = w, observed=observed, target = target, svd.out = svd.out, hilbertnorm = 1)
      if(printprogress == TRUE) {
          print(paste0("With ",thisnumdims," dimensions, biasbound (norm=1) of ", 
                       round(biasboundnow,3)))
      }
     
      dist.record=c(dist.record,biasboundnow)

      dist.now=biasboundnow # To make more generic, distance could be any measure.
      dist.orig=biasbound_orig

      thisnumdims=thisnumdims+incrementby

      if (dist.now<mindistsofar){mindistsofar=dist.now}

      wayover=(dist.now/mindistsofar)>1.25

      keepgoing= (thisnumdims<=maxnumdims) & (wayover==FALSE) #& (dist.now<dist.orig)
      #& dist.now!="error"
      # (dist.now>mindistsofar)  # XXX this was in there, but needed?
      # XXX In above, need to work on "keepgoing" for case where ebal
      # is failing.
    } # End of while loop for "keepgoing"

    dimseq=seq(minnumdims,maxnumdims,incrementby)
    numdims=dimseq[which(dist.record==min(dist.record,na.rm=TRUE))]

    #if nothing improved balance, there will be "many" minima.
    # throw warning, and choose the fewest numdims.
    if (length(numdims)>1){
      warning("Lack of improvement in balance; choosing fewest dimensions to balance on among those with the same (lack of) improvement. But beware that balance is likely poor.", 
              immediate. = TRUE)
      numdims=min(numdims)
    }

    # Finally, we didn't save weights each time, so go back and re-run
    # at optimal  number of dimensions
    if(is.null(user_numdims)){
        paste0("Re-running at optimal choice of numdims, ", numdims) 
    } else {paste0("Running at user-specified choice of numdims, ", numdims)}
    Kpc2=Kpc[,1:numdims, drop=FALSE]
    getw.out=getw(target=target, observed=observed, allrows=Kpc2)
    w=getw.out
    biasbound_opt=biasbound( w = w, observed=observed, target = target, svd.out = svd.out, hilbertnorm = 1)
}   #end for "if null numdims"

  dist_pass = rbind(dimseq[1:length(dist.record)], dist.record)
  rownames(dist_pass) <- c("Dims", "BiasBound")
  
  R=list()
  R$dist.record= dist_pass
  R$biasbound.orig=dist.orig
  R$w=w
  R$numdims=numdims
  R$biasbound.opt=biasbound_opt
  R$K = K

  return(R)
} # end kpop main function



