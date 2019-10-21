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
biasbound=function(sampled, target, svd.out, w, hilbertnorm=1){
  wtarget=w[target==1]/sum(target==1)
  wsample=w[sampled==1]/sum(sampled==1)
  
  V=svd.out$u
  eigenvals=svd.out$d
  
  V1=V[target==1, , drop=FALSE]
  V0=V[sampled==1, , drop=FALSE]
  eigenimbal=as.vector(t(wtarget)%*%V1 - t(wsample)%*%V0)
  
  effectiveimbal=(eigenimbal*(eigenvals^.5)) #%*%t(V)
  biasbound=sqrt(hilbertnorm)*sqrt(t(effectiveimbal)%*%(effectiveimbal))
  return(biasbound)
}

# Simple diffference in mean and in weighted means
# (Not actually used right now but can be convenient)
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
getw = function(target, sampled, allrows, ebal.tol=1e-6,...){
  
  # To trick ebal into using a control group that corresponds to the 
  # sampled and a treated that corresponds to the "target" group, 
  # (1) anybody who is "sampled" but also considered part of the target
  # group has to get a duplicate row in the data 
  # (2) construct a variables target_ebal that = 1 when target=1
  # but =0 in the appended data, i.e. for the sampled who are 
  # entering a second time. 
  Xappended = rbind(allrows,  allrows[sampled==1 & target==1, , drop=FALSE] )
  target_ebal = c(target, rep(0, sum(sampled==1 & target==1)))

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
    w[sampled==1]=bal.out.pc$w
    #rescale to mean=1 among the donors
    w[sampled==1]=w[sampled==1]/mean(w[sampled==1])
    #biasbound.out = biasbound(D = D, w=w, V=svd.out$v, a = svd.out$d, hilbertnorm = 1)
    #R$dist= biasbound.out  ##experimenting with using biasbound instead of L1
    #R$biasbound = biasbound.out
  }
  R$w=w
  return(R)
} # end of getw.


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
kpop = function(allx, useasbases=NULL, b=NULL, target=NULL, sampled=NULL, 
                ebal.tol=1e-6, numdims=NULL, 
                minnumdims=NULL, maxnumdims=NULL, 
                incrementby=1){
 
  N=nrow(allx)
  if (is.null(target)){target=rep(1,N)}
  
  if(is.null(sample)){sample=(target!=1)}
  
  # If we don't specify which observations to use as bases, 
  # use just the "sample" set, i.e. the non-targets. 
  if (is.null(useasbases)){
    useasbases=as.numeric(sampled==1)
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
  biasbound_orig=biasbound( w = rep(1,N), sampled=sampled, target = target, svd.out = svd.out, hilbertnorm = 1)
  paste0("Without balancing, biasbound (norm=1) is ",round(biasbound_orig,3))
  
  # If numdims given, just get the weights in one shot:
  if (!is.null(numdims)){
    Kpc2=Kpc[,1:numdims, drop=FALSE]
    getw.out=getw(target=target, sampled=sampled, allrows=Kpc2)
    # XXX This would be place to add check for non-convergence of ebal.
    w=getw.out$w
    biasboundnow=biasbound( w = w, sampled=sampled,  target = target, svd.out = svd.out, hilbertnorm = 1)
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
      getw.out=getw(target = target, sampled=sampled, allrows = Kpc_try)
      w=getw.out$w
      # Need to work on case where ebal fails and flagging this in result. 
      # For now just returns all even weights. 
      
      biasboundnow=biasbound( w = w, sampled=sampled, target = target, svd.out = svd.out, hilbertnorm = 1)
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
    getw.out=getw(target=target, sampled=sampled, allrows=Kpc2)
    w=getw.out$w
    biasbound_opt=biasbound( w = w, sampled=sampled, target = target, svd.out = svd.out, hilbertnorm = 1)
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


