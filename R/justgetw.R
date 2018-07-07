###======================
# Testing
###=======================
library(KBAL)
data(lalonde)
lalonde$nodegr=as.numeric(lalonde$educ<=11)
xvars=c("age","black","educ","hisp","married","re74","re75","nodegr","u74","u75")
kbalout=kbal(D=nsw,X=lalonde[,xvars])

# Prepare inputs:
D=lalonde$nsw
K=kbalout$K

#Choices:
numdims=kbalout$numdims  #so we can compare
ebal.tol=1e-4

w_check=getweights(K=K,D=D, numdims=numdims)

identical(kbalout$w, w_check)

###================================
### Simple function to get weights
###===============================
getweights = function(K,D, numdims, ebal.tol=1e-4){
  N=nrow(K)
  svd.out=svd(K)  #this is a bit slow.
  Kpc=svd.out$u
  K2=Kpc[,1:numdims, drop=FALSE]
  bal.out.pc=ebal::ebalance(Treatment=as.vector(D),
                          constraint.tolerance=ebal.tol,
                          X=K2,
                          print.level=-1)
  w=rep(1,N)
  w[D==0]=bal.out.pc$w
  #rescale to mean=1 among the controls
  w[D==0]=w[D==0]/mean(w[D==0])
  return(w)
} # end function
