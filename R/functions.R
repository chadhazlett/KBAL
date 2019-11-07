### Helper functions

### Build kernel, possibly non-square
makeK = function(allx, useasbases=NULL, b=NULL){
  N=nrow(allx)
  
  # If no "useasbasis" given, assume all observations are to be used. 
  if (is.null(b)){ b=ncol(allx) }
  bases = allx[useasbases==1, ]
  K=KRLS2::newKernel(X = bases , newData = allx , b = b)
}

### Get bias bound, which will be the distance
### function to drive the optimization
#similar to biasboud in kbal
#' Get the bound on the bias due to incomplete balance
#' @description XXX
#' @param observed xxx .
#' @param target xxx .
#' @param svd.out xxx .
#' @param w xxx .
#' @param hilbertnorm xxx .
#' @examples
#' biasbound=(observed=xxx, target=xxx, svd.out = xxx, w = xxx, hilbertnorm=1)
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
#' 
#' @param X Matrix of combined treated/sample population and control/target population data. Rows are observations, columns are covariates.
#' @param w vector of weights
#' @param target vector taking values of 0 or 1's indicating which observations (whuch rows  of X and w) are in the treated/sample population and which are in the control/target population.
#' @return \code{dim} the simple, unweighted difference in means
#' @return \code{dimw} the weighted difference in means
#' @example XXX 
dimw = function(X,w,target){
  w1=w[target==1]/sum(w[target==1])
  w0=w[target!=1]/sum(w[target!=1])
  
  X1=X[target==1, , drop=FALSE]
  X0=X[target!=1, , drop=FALSE]
  
  R=list()
  R$dim=colMeans(X1)-colMeans(X0)
  R$dimw=t(as.vector(w1))%*%X1-t(w0)%*%X0
}

# Function to get the moments that solve the constraints
# we have setup. 
# Currently just uses ebal, but we could to others
# or give options.
# Will need some work to better carry error messages from ebal.
#' Find weights and compute L1 distance.
#' @description  Get's the weights at the desired settings and computes
#' the objective function, L1.
#' @export
getw = function(target, observed, allrows, ebal.tol=1e-6,...){
  
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
  R=list()
  if ("try-error"%in%class(bal.out.pc)){
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
  R$w=w
  return(R)
} # end of getw.


# Oct 16 2019 Update on reasoning of names; two forms of input to kpop and 
#a third language change to the whole code
# 1. pass in "sampled" and a logical "sampledinpop" that when true indicates 
#that the sampled pop should be re-incorperated into the target population
# 2. pass in "treatment" only which specifies the target population;
#zeros in this matrix are inferred to be the "sampled" or observed popultion
# 3. throughout the code itself use "target" to refer to the target populationtion (or the treated)
#use "observed" to refer to the sampled population (or the controls)
#note that this just means that "sampled"  and "treatment" only appear as an inputs to the kpop function
#and "observed"  and "target" are used everywhere else internally in the code (biasbound, getw etc)

# 10 Oct 2019 -- re-reasoning about how to manage sample vs. target
# Please note this may not match what is implemented below!!!
# 1. Only specify "sampled" and a logical flag "sampleinpop" that indicates
# your intended target population should re-include units who happened to also be sampled.  Throw an error is sampleinpop is not specified (i.e. no default).
# There will no longer be a "target" vector at all.



# OLD stuff (pre 10 Oct 2019)
# Added "sampled", in addition to "target".
# If you tell it only who was sampled, it assumes the targets are 
# the non-sampled. 
# If you don't tell it who was sampled but do tell it who the targets are,
# it assumes the sampled are the non-targets.
# If you specify both, you can use sampled to tell it whose data you have
# and use target to specify who should inform the target group, without restriction. 
# But the usual use case for specifying both would be to say everybody
# is in the target, and only the sampled are sampled (which also the default 
# if you say nothing about the target.)


# The main event: Actual kpop function!
kpop = function(allx, useasbases=NULL, b=NULL, 
                sampled=NULL, sampledinpop=NULL,
                treatment=NULL,
                ebal.tol=1e-6, numdims=NULL, 
                minnumdims=NULL, maxnumdims=NULL, 
                incrementby=1){
    
  #need to throw error if try to pass both sample and target
  if(!is.null(sampled) & !is.null(treatment)) {
       stop("\"sampled\" and \"treatment\" arguments can not be specified simultaneously")
  }

  #throw warning for using default sampledinpop=TRUE
  if(is.null(sampledinpop)) {
      warning("using default parameter \"sampledinpop\" = TRUE")
      sampledinpop = TRUE
      
  }
    
  N=nrow(allx)
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
  }
  
  #error catch for passing in all 1's for either sampled or treated
  if(sum(observed != 1) == 0) {
      stop("\"sampled\" contains only sampled units (all entries are 1)")
  } else if(sum(target != 1) == 0 ){
      stop("\"target\" contains only target units (all entries are 1)")
  }
  
  
  # If we don't specify which observations to use as bases, 
  # use just the "sample" set, i.e. the non-targets. 
  if (is.null(useasbases)){
    useasbases=as.numeric(observed==1)
    # useasbases=rep(1,N)  #or if you want all obs as bases
  }
 
  if (is.null(minnumdims)){minnumdims=1}
  
  # The most dims you can use is the number of bases
  if (!is.null(maxnumdims) && maxnumdims>sum(useasbases)){
    warning("Cannot allow more dimensions of K than the number of bases. Reducing maxnumdims.")
    maxnumdims=sum(useasbases)
    }
  if (is.null(maxnumdims)){maxnumdims=sum(useasbases)}

  K = makeK(allx = allx, useasbases = useasbases, b=b)
  svd.out=svd(K)
  Kpc=svd.out$u
  
  # Get biasbound with no improvement in balance:
  biasbound_orig=biasbound( w = rep(1,N), observed=observed, target = target, svd.out = svd.out, hilbertnorm = 1)
  paste0("Without balancing, biasbound (norm=1) is ",round(biasbound_orig,3))
  
  # If numdims given, just get the weights in one shot:
  if (!is.null(numdims)){
    Kpc2=Kpc[,1:numdims, drop=FALSE]
    getw.out=getw(target=target, observed=observed, allrows=Kpc2)
    # XXX This would be place to add check for non-convergence of ebal.
    w=getw.out$w
    biasboundnow=biasbound( w = w, observed=observed,  target = target, svd.out = svd.out, hilbertnorm = 1)
    print(paste0("With ",numdims," dimensions, biasbound (norm=1) of ", round(biasboundnow,3)))
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
      w=getw.out$w
      # Need to work on case where ebal fails and flagging this in result. 
      # For now just returns all even weights. 
      
      biasboundnow=biasbound( w = w, observed=observed, target = target, svd.out = svd.out, hilbertnorm = 1)
      print(paste0("With ",thisnumdims," dimensions, biasbound (norm=1) of ", round(biasboundnow,3)))
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
      warning("Lack of improvement in balance; choosing fewest dimensions to balance on among those with the same (lack of) improvement. But beware that balance is likely poor.")
      numdims=min(numdims)
    }
    
    # Finally, we didn't save weights each time, so go back and re-run
    # at optimal  number of dimensions
    paste0("Re-running at optimal choice of numdims, ", numdims)
    Kpc2=Kpc[,1:numdims, drop=FALSE]
    getw.out=getw(target=target, observed=observed, allrows=Kpc2)
    w=getw.out$w
    biasbound_opt=biasbound( w = w, observed=observed, target = target, svd.out = svd.out, hilbertnorm = 1)
}   #end for "if null numdims"
  
  
  R=list()
  R$w=w
  R$dist.record=dist.record
  R$biasbound.orig=dist.orig
  R$dist.sequence=dimseq[1:length(dist.record)]
  R$numdims=numdims
  R$biasbound.opt=biasbound_opt

  return(R)
} # end kpop main function

