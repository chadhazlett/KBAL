### Helper functions

### Build kernel, possibly non-square
#' Build the Gaussian Kernel Matrix
#'
#' @description Builds the Gaussian kernel matrix using Rcpp.
#' @param allx A data matrix containing all observations where rows are units and columns are covariates.
#' @param useasbases Vector argument containing one's and zero's with length equal to the number of obervations (rows in \code{allx}) to specify which bases to use when constructing the kernel matrix and finding weights. If not specified, the default is to use all observations.
#' @param b Scaling factor in the calculation of gaussian kernel distance equivalent to the entire denominator \eqn{2\sigma^2} of the exponent. Default is twice the number of covariates or columns in \code{allx}.
#' @param linkernel Indicates that user wants linear kernel, \eqn{K=XX'}, which in practice employs \eqn{K=X} and achieves (approximate) mean balance on \eqn{X}.  
#' @return \item{K}{The kernel matrix}
#' @examples
#' #load and clean data a bit
#' \donttest{
#' data(lalonde)
#' lalonde$nodegr=as.numeric(lalonde$educ<=11)
#' xvars=c("age","black","educ","hisp","married","re74","re75","nodegr","u74","u75")
#'
#' #note that lalonde$nsw is the treatment vector, so the observered is 1-lalodne$nsw
#' #running makeK with the sampled/observed units as the bases given the large size of the data
#' #and with b as twice the number of covariates
#' K = makeK(allx = lalonde[,xvars],
#' useasbases = 1-lalonde$nsw,
#' b = 2*ncol(lalonde[,xvars]))}
#' @useDynLib KBAL
#' @importFrom stats sd 
#' @importFrom Rcpp sourceCpp
#' @export
makeK = function(allx, useasbases=NULL, b=NULL, linkernel = FALSE){
  N=nrow(allx)
  # If no "useasbasis" given, assume all observations are to be used.
  if(is.null(useasbases)) {useasbases = rep(1, N)}
  
  #default b is set to 2ncol to match kbal for now
  if (is.null(b)){ b=2*ncol(allx) }
  
  bases = allx[useasbases==1, ]
  
  Xmeans.bases <- colMeans(bases)
  Xsds.bases <- apply(bases,2,sd)
  bases <- scale(bases, center = Xmeans.bases, scale = Xsds.bases)
  allx <- scale(allx, center = Xmeans.bases, scale = Xsds.bases)
  
  if(linkernel == TRUE) {
      K = allx
  } else {
      K = new_gauss_kern(newx = allx, oldx = bases, b = b)
  }
  return(K)
}

### Get bias bound, which will be the distance
### function to drive the optimization
#similar to biasboud in kbal
#' Worst-Case Bias Bound due to Incomplete Balance
#' @description Calculate the upper bound on the bias induced by approximate balance with a given \code{hilbertnorm}. Approximate balance is conducted in \code{kbal()} and uses only the first \code{numdims} dimensions of the singular value decomposition of the kernel matrix to generate weights \code{w} which produce mean balance between control or sampled units and treated or population units. The following function calculates the worse-case bias induced by this approximate balancing with weights \code{w} and a given \code{hilbertnorm.}
#' @param observed a numeric vector of length equal to the total number of units where sampled units take a value of 1 and population units take a value of 0.
#' @param target a numeric vector of length equal to the total number of units where population units take a value of 1 and sample units take a value of 0.
#' @param svd.out the list object output from \code{svd()} performed on the kernel matrix.
#' @param w numeric vector containing the weight for every corresponding unit. Note that these weights should sum to the total number of units, not to one. They are divided by the number of control or sample and treated or population units internally.
#' @param hilbertnorm numeric value of the hilbertnorm. Default value is 1.
#' @examples
#' \donttest{
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
#'  w = rep(1,nrow(lalonde)), hilbertnorm=1)}
#' @export
biasbound=function(observed, target, svd.out, w, w.pop = NULL, sampledinpop = FALSE,
                   hilbertnorm=1){
    N = nrow(svd.out$u)
    #error check for pop weights
    if(is.null(w.pop)) {
        w.pop = rep(1,N)
    } else {
        if(sum(sign(w.pop)) != length(observed)) {
            stop("\"w.pop\" must be non-negative")
        }
        if(length(w.pop) != length(observed)) {
            stop("\"w.pop\" must have the same length as the total number of units")
        }
        if(sum(target) == length(target) && !(sum(w.pop[observed]) == sum(observed)
                                             & sd(w.pop[observed]) == 0)) {
            stop("\"w.pop\" must the value 1 for all sampled/treated units")
        }
        #check population weights sum to num of treated/population units
        if(sum(w.pop[target ==1 ]) != sum(target)) {
            #allow user to pass in weights that sum to one and transform them here
            if(sum(w.pop[target ==1]) == 1) {
                w.pop[target==1] = w.pop[target ==1]/mean(w.pop[target==1])
            } else { #in this case they don't sum to N_t or 1 so ng
                stop("\"population.w\" must sum to either 1 or the number of treated/population units")
            }
        }
    }
        
    if(!sampledinpop) {
        #weights sum to N_0, N_1 so normalizing so sum to 1
        wtarget=w.pop[target==1]/sum(target==1) 
        wobserved=w[observed==1]/sum(observed==1)
    } else { #when sampledinpop true
        #3.2.20: idk why I had this here, seems unnecessary and also never used bc we set to be 1 if null above??
        if(is.null(w.pop)) {error("\"w.pop\" required when \"sampledinpop\"= TRUE. ")}
        N = length(observed)
        wtarget= w.pop/N #over N right because w.pop and length N
        #drop the first N entries of w because all 1's, only last N_observed are meaningful
        wobserved=w[observed ==1]/sum(observed==1)
   }
   
    U=svd.out$u
    #ISSUE for frankenstein
    eigenvals=svd.out$d #singular values (A)
    
    U1=U[target==1, , drop=FALSE]
    U0=U[observed==1, , drop=FALSE]
    eigenimbal=as.vector(t(wtarget)%*%U1 - t(wobserved)%*%U0)
    effectiveimbal=(eigenimbal*(eigenvals^.5))
    biasbound=sqrt(hilbertnorm)*sqrt(t(effectiveimbal)%*%(effectiveimbal))
    return(biasbound)
}

# Simple diffference in mean and in weighted means
# (Not actually used right now but can be convenient)

#' Difference in Means and Difference in Weighted Means
#'
#' Calcuates the simple difference in means or weighted difference in means between the control or sample population and the treated or target popultion.
#' @param X matrix of data where rows are observations and columns are covariates.
#' @param w numeric vector of weights for each observation.
#' @param target numeric vector of length equal to the total number of units where population units take a value of 1 and sample units take a value of 0.
#' @return \item{dim}{the simple, unweighted difference in means.}
#' \item{dimw}{the weighted difference in means.}
#' @examples
#' \donttest{
#' #let's say we want to get the unweighted DIM and the weighted DIM using weights from the kbal
#' #function with the lalonde data:
#' #load and clean data a bit
#' data(lalonde)
#' lalonde$nodegr=as.numeric(lalonde$educ<=11)
#' xvars=c("age","black","educ","hisp","married","re74","re75","nodegr","u74","u75")
#'
#' #get the kbal weights
#' kbalout= kbal(allx=lalonde[,xvars],
#'                sampledinpop=FALSE,
#'                treatment=lalonde$nsw)
#'  #now use dimw to get the DIMs
#'  dimw(X = lalonde[,xvars], w = kbalout$w, target = lalonde$nsw)}
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

#' Find Weights using Entropy Balancing.
#' @description Uses entropy balancing to find and return the weights that produce mean balance on \eqn{\phi(X_i)}, the expaned features of \eqn{X_i} using a given kernel \eqn{\phi(.)}, for the control or sample group and treated group or target population.
#'
#' @param target a numeric vector of length equal to the total number of units where population units take a value of 1 and sample units take a value of 0.
#' @param observed a numeric vector of length equal to the total number of units where sampled units take a value of 1 and population units take a value of 0.
#' @param svd.U matrix whose columns contain the left singular vectors of the kernel matrix.
#' @param ebal.tol tolerance level used by custom entropy balancing function \code{ebalance_custom}.
#' @return \item{w}{numeric vector of weights.}
#' \item{earlyfail}{boolean indicating if ebalance failed to converge within the first two dimensions of \code{svd.U}.}
#' @examples
#' \donttest{
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
#' U = svd(K)$u
#' #usually we are getting weights using different number of columns of this matrix, the finding
#' # the bias and looking for the minimum. For now let's just use the first 10
#' U2=U[,1:10, drop=FALSE]
#' getw.out=getw(target=lalonde$nsw, observed=1-lalonde$nsw, svd.U=U2)}
#' @export
getw = function(target, observed, svd.U, ebal.tol=1e-6){

  # To trick ebal into using a control group that corresponds to the
  # observed and a treated that corresponds to the "target" group,
  # (1) anybody who is "observed" but also considered part of the target
  # group has to get a duplicate row in the data
  # (2) construct a variables target_ebal that = 1 when target=1
  # but =0 in the appended data, i.e. for the observed who are
  # entering a second time.
  Xappended = rbind(svd.U,  svd.U[observed==1 & target==1, , drop=FALSE] )
  target_ebal = c(target, rep(0, sum(observed==1 & target==1)))

    bal.out.pc=try(ebalance_custom(Treatment = target_ebal,
                                   X = Xappended,
                                   base.weight = NULL,
                                   norm.constant  = NULL,
                                   coefs = NULL ,
                                   max.iterations = 200,
                                   constraint.tolerance = ebal.tol,
                                   print.level=0),
                   silent=TRUE)
  N=nrow(svd.U)
  earlyfail = FALSE
  if ("try-error"%in%class(bal.out.pc)){
      if(ncol(svd.U) <= 2) {
          warning("Kbal was unable to successfully find weights because ebalance convergence failed within first two dimensions of K. Returning equal weights.")
          earlyfail = TRUE
      }
    w=rep(1,N)
    
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
    out <- list(w = w, 
                earlyfail = earlyfail)
  return(out)
} # end of getw.

#' L1 Distance
#' @description Calculates the L1 distance between the treated or population units and the kernel balanced control or sampled units.
#' @param target a numeric vector of length equal to the total number of units where population units take a value of 1 and sample units take a value of 0.
#' @param observed a numeric vector of length equal to the total number of units where sampled units take a value of 1 and population units take a value of 0.
#' @param K the kernel matrix
#' @param svd.out the list object output from performing \code{svd()} on the kernel matrix.
#' @param w a numeric vector of weights for every obervation. If unspecified, these are found using \code{numdims} dimensions of the SVD of the kernel matrix \code{svd.out$u} with custom entropy balancing function \code{ebalance_custom()}. Note that these weights should sum to the total number of units, where treated or population units have a weight of 1 and control or sample units have appropriate weights dervied from kernel balancing with mean 1 which is consistent with the ouput of \code{getw()}.
#' @param numdims a numeric input specifying the number of columns of the singular value decomposition of the kernel matrix to use when finding weights in the case that \code{w} is not specified.
#' @param ebal.tol an optional numeric input speccifying the tolerance level used by custom entropy balancing function \code{ebalance_custom()} in the case that \code{w} is not specified. When not specified, the default is 1e-6/
#' @return \item{w}{numeric vector of weights used}
#' \item{L1}{a numeric giving the L1 distance, the absolute difference between \code{pX_D1} and \code{pX_D0w}}
#' \item{pX_D1}{a numeric vector of length equal to the total number of observations where the nth entry is the sum of the kernel distances from the nth unit to every treated or population unit.}
#' \item{pX_D0}{a numeric vector of length equal to the total number of observations where the nth entry is the sum of the kernel distances from the nth unit to every control or sampled unit.}
#' \item{pX_D0w}{a numeric vector of length equal to the total number of observations where the nth entry is the weighted sum of the kernel distances from the nth unit to every control or sampled unit. The weights are given by entropy balancing and produce mean balance on \eqn{\phi(X)}, the expaned features of \eqn{X} using a given kernel \eqn{\phi(.)}, for the control or sample group and treated group or target population.}
#' @examples
#' \donttest{
#' #loading and cleaning lalonde data
#' data(lalonde)
#' lalonde$nodegr=as.numeric(lalonde$educ<=11)
#' xvars=c("age","black","educ","hisp","married","re74","re75","nodegr","u74","u75")
#'
#' #need to first build gaussian kernel matrix
#' K_pass <- makeK(allx = lalonde[,xvars])
#' #also need the SVD of this matrix
#' svd.U_pass <- svd(K_pass)
#'
#' #running without passing weights in directly, using numdims=33
#' l1_lalonde <- getdist(target = lalonde$nsw,
#'                       observed = 1-lalonde$nsw,
#'                       K = K_pass,
#'                       svd.out = svd.U_pass,
#'                       numdims = 33)
#'
#'  #alternatively, we can get the weights ourselves and pass them in directly
#'  w_opt <- getw(target= lalonde$nsw,
#'                observed = 1-lalonde$nsw,
#'                svd.U = svd.U_pass$u[,1:33, drop=FALSE],
#'                ebal.tol=1e-6)$w
#'  l1_lalonde2 <- getdist(target = lalonde$nsw,
#'                   observed = 1-lalonde$nsw,
#'                   K = K_pass,
#'                   svd.out = svd.U_pass,
#'                   w = w_opt)}
#' @export
getdist <- function(target, observed, K,
                    w=NULL, numdims = NULL, w.pop = NULL, ebal.tol=NULL) {

        R=list()
        N=nrow(K)
        K_c=K[observed==1, ,drop = FALSE]
        K_t=K[target==1, ,drop=FALSE]
        if(is.null(w.pop)) {
            w.pop = rep(1,N)
        } else {
            if(sum(sign(w.pop)) != length(observed)) {
                stop("\"w.pop\" must be non-negative")
            }
            if(length(w.pop) != length(observed)) {
                stop("\"w.pop\" must have the same length as the total number of units")
            }
            if(!(sum(w.pop[observed==1]) == sum(observed) & sd(w.pop[observed==1]) == 0)) {
                stop("\"w.pop\" must the value 1 for all sampled/treated units")
            }
            #check population weights sum to num of treated/population units
            if(sum(w.pop[target ==1 ]) != sum(target)) {
                #allow user to pass in weights that sum to one and transform them here
                if(sum(w.pop[target ==1]) == 1) {
                    w.pop[target==1] = w.pop[target ==1]/mean(w.pop[target==1])
                } else { #in this case they don't sum to N_t or 1 so ng
                    stop("\"population.w\" must sum to either 1 or the number of treated/population units")
                }
            }
            
        }
        #if user does not provide weights, go get them
        # XXXX INCORPERATE w.pop XXX
        if(is.null(w)) {
            if(is.null(ebal.tol)) {ebal.tol = 1e-6}
            if(is.null(numdims)) {stop("If weights w input is not specified, numdims must be in order to calculate these weights internally")}
            U_w.pop <- w.pop*U
            w = getw(target = target, observed=observed,
                     svd.U = U_w.pop[,1:numdims, drop=FALSE],
                     ebal.tol=ebal.tol)$w

            #if ebal fails we get weights of 1 for everyone
            if (sum(w ==1) == length(w)){
                stop("ebalance failed to converge for this choice of numdims dimensions of the SVD of the kernel matrix")
            }
        }
        #just "average row Kt"
        pX_D1=(matrix(1,1, sum(target==1))/sum(target==1))%*% K_t 
        #average row Kc
        pX_D0=(matrix(1,1,sum(observed==1))/sum(observed==1))%*% K_c 
        #weighted average Kc with ebal weights
        pX_D0w=(w[observed==1]/sum(w[observed==1])) %*% K_c
        #weighted average Kt ONLT DIFF from pX_D1 if have pop weights not all equal to 1
        pX_D1wpop = (w.pop[target==1]/sum(w.pop[target==1])) %*% K_t 

            # A rescaling to ensure sum is 1 as would be an integral.
        pX_D1=pX_D1/sum(pX_D1)
        pX_D0=pX_D0/sum(pX_D0)
        pX_D0w=pX_D0w/sum(pX_D0w)
        pX_D1wpop =pX_D1wpop/sum(pX_D1wpop)
        L1 = sum(abs(pX_D1wpop-pX_D0w)) #removed the 0.5 factor -- Oct 20 2017


        R$L1=L1
        R$w=w
        R$pX_D1wpop=pX_D1wpop
        R$pX_D0w=pX_D0w
        R$pX_D1=pX_D1
        R$pX_D0=pX_D0
        R$pX_D0w=pX_D0w

        return(R)
} ## end of getdist


# The main event: Actual kbal function!
#' Kernel Balancing
#'
#' @description Kernel balancing (KBAL) is non-parametric weighting tool to make two groups have a similar distribution of covariates, not just in terms of means or marginal distributions but (i) on general smooth functions of the covariates, including (ii) on a smoothing estimator of the joint distribution of the covariates. It was originally designed (Hazlett, 2017) to make control and treated groups look alike, as desired when estimating causal effects under conditional ignorabiity. This package also facilitates use of this approach for more general distribution-alignment tasks, such as making a sampled group have a similar distribution of covariates as a target population, as in survey reweighting. The examples below provide an introduction to both settings.
#' 
#' To proceed in the causal effect setting, KBAL assumes that the expectation of the non-treatment potential outcome conditional on the covariates falls in a large, flexible space of functions associated with a kernel. It then constructs linear bases for this function space and achieves approximate balance on these bases. The approximation is one that minimizes the worst-case bias that could persist due to remaining imbalances. 
#' 
#' The \code{kbal} function implements kernel balancing using a gaussian kernel to expand the features of \eqn{X_i} to infinite dimensions.  It finds approximate mean balance for the control or sample group and treated group or target population in this expanded feature space by using the first \code{numdims} dimensions of the singular value decomposition of the gaussian kernel matrix. It employs entropy balancing to find the weights for each unit which produce this approximate balance. When \code{numdims} is not user-specified, it searches through increasing dimensions of the SVD of the kernel matrix to find the number of dimensions which produce weights that minimizes the worst-case bias bound with a given \code{hilbertnorm}. It then results these optimal weights, along with the minimized bias, the kernel matrix, a record of the number of dimensions used and the corresponding bais, as well as an original bias using naive group size weights for comparison.
#'
#' @references Hazlett, C. (2017), "Kernel Balancing: A flexible non-parametric weighting procedure for estimating causal effects." Forthcoming in Statistica Sinice. https://doi.org/10.5705/ss.202017.0555
#' 
#' @param allx a data matrix containing all observations where rows are units and columns are covariates.
#' @param useasbases optional vector of 0/1 or FALSE/TRUE to specify what observations are to be used in forming bases (columns of the kernel matrix) balanced upon.  If the number of observations is under 2000, the default is to use all observations. When the number of observations is over 2000, the default is to use the sampled (control) units only.
#' @param b scaling factor in the calculation of gaussian kernel distance equivalent to the entire denominator \eqn{2\sigma^2} of the exponent.
#' @param K optional user-supplied kernel matrix. If given, this kernel matrix is used rather than computing one.
#' @param sampled a numeric vector of length equal to the total number of units where sampled units take a value of 1 and population units take a value of 0.
#' @param sampledinpop a logical to be used in combination with input \code{sampled} that when \code{TRUE} indicates that sampled units should also be included in the target population.
#' @param treatment an alternative input to \code{sampled} and \code{sampledinpop} that is a numeric vector of length equal to the total number of units. Current version supports the ATT estimand. Accordingly, the treated units are the target population, and the control are equivalent to the sampled. Weights play the role of making the control groups (sampled) look like the target population (treated).  \code{sampledinpop} is forced to be \code{FALSE}.
#' @param linkernel if true, uses the linear kernel which is technicaly \eqn{K=XX'}. In practice this simply achieves mean balance on the original X. For speed purposes, the code effectively employs \eqn{K=X} instead, but this is equivalent to \eqn{K=XX'} for our purposes because they have the same left-singular vectors. It is thus nearly equivalent to entropy balancinng on means. The difference is that it employs SVD on X then seeks balanc on the left singular vectors, using the bias bound to determine how many dimensions to balance. Thus in cases where full balance may be infeasible, it automatically resorts to approximate balance.
#' @param ebal.tol tolerance level used by custom entropy balancing function \code{ebalance_custom()}.
#' @param numdims optional numeric argument to specify the number of dimensions of the kernel matrix to find balance on rather than searching for the number of dimensions which minimize the bias.
#' @param minnumdims optional numeric argument to specify the minimum number of dimensions of the SVD of the kernel matrix to find balance on in the search for the number of dimesions which minimize the bias. Default minimum is 1.
#' @param maxnumdims optional numeric argument to specify the maximum number of dimensions of the SVD of the kernel matrix to find balance on in the search for the number of dimesions which minimize the bias. While the number of observations is under 2000, the default maximum is the total number of observations. Due to the computation burden, when the number of observations is over 2000, the default is the number of sampled units.
#' @param incrementby optional argument to specify the number of dimesions to increase by from \code{minnumdims} to \code{maxnumdims} in each iteration of the search for the number of dimensions which minimizes the bias. Default is 1.
#' @param printprogress optional logical argument to print current number of dimensions and bias.
#'
#' @return \item{dist.record}{a numeric matrix recording the bias bound corresponding to balance on increasing dimesions of the SVD of the kernel matrix starting from \code{minnumdims} increasing by \code{incrementby} to \code{maxnumdims} or until the bias grows to be 1.25 times the minimal bias found.}
#'  \item{biasbound.orig}{the bias bound found when all sampled units have a weight of one over the number of sampled units and all target units have a weight of one over the number of target units.}
#'  \item{numdims}{the optimal number of dimensions of the SVD of the kernel matrix which minimizes the bias bound.}
#'  \item{w}{the weights found using entropy balancing on \code{numdims} dimensions of the SVD of the kernel matrix.}
#'  \item{earlyfail}{boolean to indicate whether ebalance failed to converge within the first two dimesions of the SVD of K} 
#'  \item{biasbound.opt}{the minimal bias bound found using \code{numdims} as the number of dimestions of the SVD of the kernel matrix. When \code{numdims} is user-specified, the bias bound using this number of dimensions of the kernel matrix.}
#' \item{K}{the kernel matrix.}
#' @examples
#' #----------------------------------------------------------------
#' # Example 1: Reweight a control group to a treated to esimate ATT. 
#' # Benchmark using Lalonde et al.
#' #----------------------------------------------------------------
#' data(lalonde)
#' lalonde$nodegr=as.numeric(lalonde$educ<=11)
#' xvars=c("age","black","educ","hisp","married","re74","re75","nodegr","u74","u75")
#'  \donttest{
#' # Rerun Lalonde example with settings as in the KBAL paper:
#' kbalout.full= kbal(allx=lalonde[,xvars], b=length(xvars),
#'                useasbases=rep(1,nrow(lalonde)),
#'                treatment=lalonde$nsw)
#' summary(lm(re78~nsw,w=kbalout.full$w, data = lalonde))  
#'  }
#' #----------------------------------------------------------------
#' # Example 1B: Reweight a control group to a treated to esimate ATT. 
#' # Benchmark using Lalonde et al. -- but just mean balancing now 
#' # via "linkernel".
#' #----------------------------------------------------------------
#'  
#' # Rerun Lalonde example with settings as in the KBAL paper:
#' kbalout.lin= kbal(allx=lalonde[,xvars], b=length(xvars),
#'                useasbases=rep(1,nrow(lalonde)),
#'                treatment=lalonde$nsw, linkernel=TRUE)
#' 
#' # Check balance with and without these weights:
#' dimw(X=lalonde[,xvars], w=kbalout.lin$w, target=lalonde$nsw)
#' 
#' summary(lm(re78~nsw,w=kbalout.lin$w, data = lalonde))  
#'  
#' #----------------------------------------------------------------
#' # Example 2: Reweight a sample to a target population.
#' #----------------------------------------------------------------
#' # Suppose a population consists of four groups in equal shares: 
#' # white republican, non-white republican, white non-republicans, 
#' # and non-white non-republicans. A given policy happens to be supported 
#' # by all white republicans, and nobody else. Thus the mean level of 
#' # support in the population should be 25%. 
#' #
#' # Further, the sample is surveyed in such a way that was careful 
#' # to quota on party and race, obtaining 50% republican and 50% white.
#' # However, among republicans three-quarters are white and among non-republicans,
#' # three quarters are non-white. This biases the average level of support
#' # despite having a sample that matches the population on its marginal distributions. #'
#' # We'd like to reweight the sample so it resembles the population not 
#' # just on the margins, but in the joint distribution of characteristics.
#' 
#' pop <- data.frame(
#' republican =  c(rep(0,400), rep(1,400)),
#' white = c(rep(1,200), rep(0,200), rep(1,200), rep(0,200)),
#' support = c(rep(1,200), rep(0,600)))
#'   
#' mean(pop$support)  # Target value
#'  
#' # Survey sample: correct margins/means, but wrong joint distribution
#' samp <- data.frame( republican = c(rep(1, 40), rep(0,40)),
#'    white = c(rep(1,30), rep(0,10), rep(1,10), rep(0,30)),
#'    support = c(rep(1,30), rep(0,50)))
#'   
#' mean(samp$support)  # Appears that support is 37.5% instead of 25%.
#'  
#' # Mean Balancing -----------------------------------------
#' # Sample is already mean-balanced to the population on each 
#' # characteristic. However for illustrative purposes, use ebal() 
#' dat <- rbind(pop,samp)
#' 
#' # Indicate which units are sampled (1) and which are population units(0)
#' sampled <- c(rep(0,800), rep(1,80))
#'  
#' # Run ebal (treatment = population units = 1-sampled)
#' ebal_out <- ebalance_custom(Treatment = 1-sampled, 
#'                             X=dat[,1:2],
#'                             constraint.tolerance=1e-6, 
#'                             print.level=-1)
#'                             
#'                             
#'  
#' # We can see everything gets even weights, since already mean balanced.
#' length(unique(ebal_out$w))
#' 
#' # And we end up with the same estimate we started with
#' weighted.mean(samp[,3], w = ebal_out$w)
#'  
#' # We see that, because the margins are correct, all weights are equal
#' unique(cbind(samp, e_bal_weight = ebal_out$w))
#' 
#' # Kernel balancing for weighting to a population (i.e. kpop) -------
#' kbalout = kbal(allx=dat[,1:2],
#'                 useasbases=rep(1,nrow(dat)), 
#'                 sampled = sampled, 
#'                 b = 1,
#'                 sampledinpop = FALSE)
#'                 
#' # The weights now vary:
#' plot(kbalout$w[sampled ==1], pch=16)
#' 
#' # And produce correct estimate:
#' weighted.mean(samp$support, w = kbalout$w[sampled==1])    
#'  
#' # KBAL correctly downweights white republicans and non-white non-republicans
#' # and upweights the non-white republicans and white non-republicans
#' unique(cbind(samp[,-3], k_bal_weight = kbalout$w[sampled==1]))
#' @export
kbal = function(allx, useasbases=NULL, b=NULL, 
                sampled=NULL, sampledinpop=NULL,
                treatment=NULL,
                population.w = NULL,
                K=NULL, K.svd = NULL,
                linkernel = FALSE, constraint = NULL,
                numdims=NULL,
                minnumdims=NULL, maxnumdims=NULL,
                fullSVD = FALSE,
                incrementby=1,
                ebal.tol=1e-6,
                printprogress = TRUE) {

    N=nrow(allx)
    
#####start of big error catch series to check if data is passed in correctly and
    #default setting/data set up
  #1. checking sampled and sampledinpop
  if(!is.null(sampled)) { 
      if(!(all(sampled %in% c(0,1)))) { #note this is still ok for logicals
          stop("\"sampled\" contains non-binary elements")
      }
      if(sd(sampled) == 0) { 
          stop("\"sampled\" has zero variance")
      }
      if(length(sampled) != N) {
          stop("Dimensions of \"sampled\" do not match data \"allx\"")
      }
      #now check sampledinpop
      if(is.null(sampledinpop)) { 
          #if pass in sampled and dn specify sampledinpop set default and give warning
          warning("using default parameter \"sampledinpop\" = TRUE", immediate. = TRUE)
          sampledinpop = TRUE
      } else if(!(sampledinpop %in% c(0,1))) { #if pass in sampledinpop check its binary
          stop("\"sampledinpop\" is not binary" )
      }
    #2. now checking treatment if dn pass in sampled 
  } else if(!is.null(treatment)) { 
      if(!(all(treatment %in% c(0,1)))) { 
          stop("\"treated\" contains non-binary elements")
      }
      if(sd(treatment) == 0) {
          stop("\"treated\" has zero variance")
      }
      if(length(treatment) != N) {
          stop("Dimensions of \"treatment\" do not match data \"allx\"")
      }
  } else { #only get here if both sampled and treated are null
      stop("Either \"sampled\" or \"treatment\" must be specified")
  } #end of sampled/treated if else check
    
    #3. checking if user tried to pass in both
    if(!is.null(sampled) & !is.null(treatment)) {
        stop("\"sampled\" and \"treatment\" arguments can not be specified simultaneously")
    }
    #4. For now we will only support ATT for "treatment" case.  This means sampledinpop is FALSE
    if(!is.null(treatment) & (is.null(sampledinpop) || sampledinpop == TRUE)) {
        sampledinpop=FALSE
        warning("Targeting ATT, which implies sampledinpop=FALSE.", immediate. = TRUE)
    }
    #5. checking for covariates with no variance:
    if(0 %in% apply(allx, 2, sd)) {
        stop("One or more column in \"allx\" has zero variance")
    }
    #6. error catch for NAs in data
    if(sum(is.na(allx) != 0)) {
        stop("\"allx\" contains missing values")
    }
    
#####Setting up data: build observed and target from inputs 
    if(!is.null(sampled) & sampledinpop==FALSE) {
        observed = sampled
        target = 1-sampled
    } else if(!is.null(sampled)) { 
        observed = sampled
        target = rep(1,N)
    } else{ #note that we checked above for both sampled/treated null, so this must be
        #treated passed in case
        observed = 1-treatment
        target = treatment
    }
    
    #7. checking useasbases if passed in
    if(!is.null(useasbases)) {
        if(length(useasbases) != N) {
            stop("Dimensions of \"useasbases\" do not match data \"allx\"")
        }
        if(!(all(useasbases %in% c(0,1)))) {
            stop("\"useasbases\" contains non-binary elements")
        }
        #check do not pass in K
        if(!is.null(K) | !is.null(K.svd)) {
            warning("\"useasbases\" argument only used in the construction of the kernel matrix \"K\" and should not be specified when \"K\" or \"K.svd\" is already user-supplied. Using all columns.", immediate. = TRUE)
        }
    }
    
    #Setting defaults - useasbases: If we don't specify which observations to use as bases,
    # use all as default unless K is very large, then use sample set.
    if (is.null(useasbases) & N <= 2000) {
            useasbases = rep(1,N)
    } else if(is.null(useasbases)) {
        if(is.null(K) & is.null(K.svd) ) {
            warning("Dimensions of K greater than 2000, using sampled as default bases",
                    immediate. = TRUE)
        }
          useasbases = as.numeric(observed==1)
    }
 
    #Population  weights
    #default for population weights is to have them all equal
    #population weights need to sum to N because the averaging among the treated occurs
    #within ebal (so we don't want to double average)
    if(is.null(population.w)) {
        w.pop = rep(1,N)
    } else {
        if(sum(target) != 1 && length(population.w) ==1) {
            stop("\"population.w\" has length one. Please ensure it is a vector.")
        }
        #check do not have any negative weights
        if(length(population.w) != sum(target)) {
            stop("\"population.w\" must have the same length as the number of population/treated units")
        }
        if(sum(sign(population.w)) != sum(target)) {
            stop("\"population.w\" must be non-negative")
        }
        #check population weights sum to num of treated/population units
        if(sum(population.w) != sum(target)) {
            #allow user to pass in weights that sum to one and transform them here
            if(sum(population.w) == 1) {
                population.w = population.w/mean(population.w)
            } else { #in this case they don't sum to N_t or 1 so ng
                stop("\"population.w\" must sum to either 1 or the number of treated/population units")
            }
        }
        #adding uniform weights for control/sampled units
        w.pop <- rep(1, N)
        w.pop[target==1] = population.w
        
    }
    
    #8. now checking maxnumdims: the most dims you can use is the number of bases
    if (!is.null(maxnumdims) && maxnumdims>sum(useasbases)) {
        warning("Cannot allow dimensions of K to be greater than the number of bases. Reducing \"maxnumdims\".", immediate. = TRUE)
        maxnumdims=sum(useasbases)
    }#make sure don't send in a neg
    if (!is.null(maxnumdims) && maxnumdims<=0) {
        stop("\"maxnumdims\" must be greater than zero")
    } #when linkernel ==TRUE ensure maxnumdims not greater than cols of X
    if(linkernel == TRUE && !is.null(maxnumdims) && maxnumdims > ncol(allx)) {
        warning("When using a linear kernel, cannot allow dimensions of K to be greater than the number of columns in \"allx\". Reducing to the number of coumns in \"allx\".", 
                immediate. = TRUE )
        maxnumdims = ncol(allx)
    }
    
    #Setting defaults: minnumdims, maxnumdims
    if (is.null(minnumdims)){minnumdims=1}
    if (is.null(maxnumdims)){
        if(linkernel == FALSE) {
            maxnumdims= min(200, sum(useasbases)) #reducing to use RSpectra to max 200
        } else { maxnumdims = min(ncol(allx), 200) } 
    }
    #setting defaults - b: dding default b within the kbal function rather than in makeK
    #changing default to be 2*ncol to match kbal
    if (is.null(b)){ b = 2*ncol(allx) }
    
    #9. now checking numdims if passed in
    if(!is.null(numdims) && numdims>maxnumdims) { #check not over max
        warning("\"numdims\" cannot exceed number of bases. Reducing to maximum allowed.",
                immediate. = TRUE)
        numdims=sum(useasbases)
    } else if(!is.null(numdims) && numdims <= 0) { #check not leq zero
        stop("Specified \"numdims\" must be greater than zero")
    }
    
    #10. incrementby not null and geq 1
    if(is.null(incrementby) || incrementby < 1){
        warning(" \"incrementby\" must be greater than or equal to 1. Setting \"incrementby\" to be 1.", immediate. = TRUE)
        incrementby = 1
    }
#####end of big error catch series and data setup

########### BUILDING K ################
#Setting Up K: Make the kernel or take in user K or user K.svd and check those work with dims
   #first check didn't pass both
  if(!is.null(K) & !is.null(K.svd)){
      stop("\"K\" and \"K.svd\" can not be specified simultaneously")
#CASE 1: if user pases K, always conduct SVD upto maxnumdims or numdims
  } else if(!is.null(K)) { 
      #error checks:
      #check maxnumdims
      if(maxnumdims > ncol(K) ) {
          warning("\"maxnumdims\" cannot exceed number of columns of \"K\". Reducing to maximum allowed.", immediate. = TRUE)
          maxnumdims = ncol(K)
      }
      #check numdims
      if(!is.null(numdims) && numdims > ncol(K) ) {
          warning("\"numdims\" cannot exceed number of columns of \"K\". Reducing to maximum allowed.", immediate. = TRUE)
          numdims = ncol(K)
      }
      #check minnumdims 
      if(minnumdims > maxnumdims) {
          warning("\"minnumdims\" cannot exceed \"maxnumdims\". Reducing \"minnumdims\" to 1.", immediate. = TRUE)
          minnumdims = 1
      }
      #warning that linkernel meaningless
      if(!is.null(linkernel) && linkernel) {
          warning("\"linkernel\" argument only used in the construction of the kernel matrix \"K\" and is not used when \"K\" or \"K.svd\" is already user-supplied.", immediate. = TRUE)
      }
      if(b != 2*ncol(allx)) {
          warning("\"b\" argument only used in the construction of the kernel matrix \"K\" and is not used when \"K\" or \"K.svd\" is already user-supplied. Using all columns.", immediate. = TRUE)
      }
      #provided we pass all those checks get svd with RSpectra
      if(printprogress == TRUE) {cat("Using user-supplied K \n")}
      #if user does not ask for fullsvd, and does not give numdims, get svd upto maxnumdims
      if(!fullSVD & is.null(numdims)) { 
          if(printprogress) {cat("Running SVD on kernel matrix upto maxnumdims", maxnumdims, "dimensions \n")}
          #for a symmetric K just do eigs_sym as is:
          if(nrow(K) == ncol(K)) {
              rspec.out= RSpectra::eigs_sym(K, maxnumdims)
              #this is really SAS' but since we use svd throughout pass out in svd lang
              svd.out = list(u = rspec.out$vectors, d = rspec.out$values,
                             v= rspec.out$vectors)
          } else { #use svds
              svd.out= RSpectra::svds(K, maxnumdims)
          }
          U=svd.out$u
      #if user does not want full svd but did pass numdims, get svd upto numdims
      } else if(!fullSVD) { 
          if(printprogress) {cat("Running SVD on kernel matrix upto numdims", numdims, "dimensions \n")}
          if(nrow(K) == ncol(K)) {
              rspec.out= RSpectra::eigs_sym(K, numdims)
              #this is really SAS' but since we use svd throughout pass out in svd lang
              svd.out = list(u = rspec.out$vectors, d = rspec.out$values, 
                             v= rspec.out$vectors)
          } else { #use svds
              svd.out = RSpectra::svds(K, numdims)
          }
          U=svd.out$u
    #if user asked for full SVD, get it
      } else { 
          if(printprogress) {cat("Running full SVD on kernel matrix \n")}
          svd.out = svd(K)
          U = svd.out$u
      }
#CASE 2: if user pases in K.svd never conduct SVD   
  } else if(!is.null(K.svd)) {
      #NB: we require K.svd to have $u and $d just as a real svd woulda
      #error catches
      #check maxnumdims
      if(maxnumdims > ncol(K.svd$u) ) {
          warning("\"maxnumdims\" cannot exceed number of columns of \"K.svd\". Reducing to maximum allowed.", immediate. = TRUE)
          maxnumdims = ncol(K)
      }
      #check numdims
      if(!is.null(numdims) && numdims > ncol(K.svd$u) ) {
          warning("\"numdims\" cannot exceed number of columns of \"K.svd\". Reducing to maximum allowed.", immediate. = TRUE)
          numdims = ncol(K)
      }
      #check minnumdims 
      if(minnumdims > maxnumdims) {
          warning("\"minnumdims\" cannot exceed \"maxnumdims\". Reducing \"minnumdims\" to 1.", immediate. = TRUE)
          minnumdims = 1
      }
      if(!is.null(linkernel) && linkernel) { #only if linkernel = TRUE
          warning("\"linkernel\" argument only used in the construction of the kernel matrix \"K\" and is not used when \"K\" or \"K.svd\" is already user-supplied.", immediate. = TRUE)
      }
      if(b != 2*ncol(allx)) {
          warning("\"b\" argument only used in the construction of the kernel matrix \"K\" and is not used when \"K\" or \"K.svd\" is already user-supplied. Using all columns.", immediate. = TRUE)
      }
      svd.out = K.svd
      U = K.svd$u
      K = U
#CASE 3: if user does not specify either K or K.svd, build K and get svd of K
  } else { 
      if(printprogress == TRUE) {cat("Building kernel matrix \n")}
      if(linkernel == FALSE) {
          K = makeK(allx = allx, useasbases = useasbases, b=b)
      } else {
          K = makeK(allx = allx,
                    useasbases = useasbases, 
                    linkernel = TRUE)
      }
     #if user does not ask for full svd, and does not pass in numdims, get svd upto maxnumdim
      if(!fullSVD & is.null(numdims)) {
          if(printprogress) {cat("Running SVD on kernel matrix upto maxnumdims", maxnumdims, "dimensions \n")}
          #for a symmetric K just do eigs_sym as is:
          if(nrow(K) == ncol(K)) {
              rspec.out= RSpectra::eigs_sym(K, maxnumdims)
              #this is really SAS' but since we use svd throughout pass out in svd lang
              svd.out = list(u = rspec.out$vectors, d = rspec.out$values,
                             v= rspec.out$vectors)
          } else { #use trucated svd
              svd.out= RSpectra::svds(K, maxnumdims)
          }
          U=svd.out$u
      #if user does not ask for full svd, but does pass in numdims, get svd upto numdims  
      } else if(!fullSVD) {
          if(printprogress) {cat("Running SVD on kernel matrix upto numdims", numdims, "dimensions \n")}
          if(nrow(K) == ncol(K)) {
              rspec.out= RSpectra::eigs_sym(K, numdims)
              #this is really SAS' but since we use svd throughout pass out in svd lang
              svd.out = list(u = rspec.out$vectors, d = rspec.out$values,
                             v= rspec.out$vectors)
          } else { 
              svd.out= RSpectra::svds(K, numdims)
          }
          U=svd.out$u
      }
      #if user askes for full svd, go get it
      else {
          if(printprogress) {cat("Running full SVD on kernel matrix \n")}
          svd.out = svd(K)
          U = svd.out$u
      }
  }
    
####### Adding Constraint to minimization: paste constraint vector to front of U
    if(!is.null(constraint)) {
        #check dims of constraint
        #this will cause problems if pass in one constraint as a vector, it needs to be a 1 column matrix
        if(nrow(constraint) != N) { 
            error("\"constraint\" must have the same number of rows as \"allx\"")
        } 
        minnumdims <- ncol(constraint) + 1
        U = cbind(constraint, U) 
        #if numdims given move it up to accomodate the constraint
        if(!is.null(numdims)) {
            numdims = numdims + ncol(constraint)
        }
    }

  # BASELINE
  # Get biasbound with no improvement in balance:
  #recall: w.pop is flat weights for sampled, user specified weights for population
  biasbound_orig=biasbound(w = rep(1,N),
                           w.pop = w.pop, 
                           observed=observed, target = target,
                           svd.out = svd.out, hilbertnorm = 1)
#XXXXX ISSUE not well defined when pass in K.svd not K
  getdist.orig = getdist(target=target, observed = observed,
                         w = rep(1,N), w.pop = w.pop, K=K)
  L1_orig = getdist.orig$L1

  if(printprogress==TRUE) {
      cat("Without balancing, biasbound (norm=1) is ", round(biasbound_orig,3), " and the L1 discrepancy is ", round(L1_orig,3), " \n")
  }

  # If numdims given, just get the weights in one shot:
  if (!is.null(numdims)){
    U2=U[,1:numdims, drop=FALSE]
    
    U2.w.pop <- w.pop*U2
    getw.out=getw(target=target, observed=observed, svd.U=U2.w.pop)
    biasboundnow=biasbound( w = getw.out$w,
                            observed=observed,  target = target,
                            svd.out = svd.out, 
                            w.pop = w.pop,
                            sampledinpop = sampledinpop, 
                            hilbertnorm = 1)
    if(printprogress == TRUE & is.null(constraint)) {
        cat("With user-specified ", numdims," dimension(s), biasbound (norm=1) of ",
                 round(biasboundnow,3), " \n")
    } else if(printprogress) {
        cat("With user-specified",numdims - ncol(constraint)," dimensions of K, biasbound (norm=1) of ",
            round(biasboundnow,3), " \n")
    }
    
    #stuff to set so we can skip the entire if statement below and just printout
    
    dist.record = biasboundnow
    if(is.null(constraint)) {
        dist_pass = rbind(numdims, dist.record)
    } else {
        dist_pass = rbind(numdims - ncol(constraint), dist.record)
    }
    biasbound_opt = biasboundnow
    dist.orig= biasbound_orig
    L1_optim = getdist(target=target, observed = observed,
                       w = getw.out$w, w.pop = w.pop, K=K)$L1
    
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
      U_try=U[,1:thisnumdims, drop=FALSE]
      U_try.w.pop <- w.pop*U_try
      getw.out=getw(target = target, observed=observed, svd.U = U_try.w.pop)
      # Need to work on case where ebal fails and flagging this in result.
      # For now just returns all even weights.
      
      biasboundnow=biasbound(w = getw.out$w,
                             observed=observed,  target = target,
                             svd.out = svd.out, 
                             w.pop = w.pop,
                             sampledinpop = sampledinpop, 
                             hilbertnorm = 1)
      if(printprogress == TRUE & is.null(constraint)) {
          cat("With ",thisnumdims," dimensions, biasbound (norm=1) of ",
                       round(biasboundnow,3), " \n")
      } else if(printprogress == TRUE) {
          cat("With ",thisnumdims - ncol(constraint)," dimensions of K, biasbound (norm=1) of ",
               round(biasboundnow,3), " \n")
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
      # Need to work on "keepgoing" for case where ebal fails.
      
      #if ebal fails to converge within first 2 dims we stop searching (1/22/19)
      if(getw.out$earlyfail == TRUE) { keepgoing = FALSE}
    } # End of while loop for "keepgoing"

    dimseq=seq(minnumdims,maxnumdims,incrementby)
    numdims=dimseq[which(dist.record==min(dist.record,na.rm=TRUE))]
    

    # If nothing improved balance, there will be multiple minima.
    # Throw warning, and choose the fewest numdims.
    if (length(numdims)>1){
      warning("Lack of improvement in balance; choosing fewest dimensions to balance on among those with the same (lack of) improvement. But beware that balance is likely poor.",
              immediate. = TRUE)
      numdims=min(numdims)
    }

    # Finally, we didn't save weights each time, so go back and re-run
    # at optimal  number of dimensions
    
    if(printprogress == TRUE) {
        cat("Re-running at optimal choice of numdims, ", numdims, "\n")
    }
    U_final.w.pop <- w.pop*U[,1:numdims, drop = FALSE]
    getw.out = getw(target= target, observed=observed, svd.U=U_final.w.pop)
    biasbound_opt= biasbound(w = getw.out$w, observed=observed, target = target, 
                             svd.out = svd.out, 
                             w.pop = w.pop,
                             sampledinpop = sampledinpop,
                             hilbertnorm = 1)
    #ISSUE IF K.svd passed in XXX
    L1_optim = getdist(target=target, observed = observed,
                       w = getw.out$w, w.pop = w.pop, K=K)$L1
    dist_pass = rbind(dimseq[1:length(dist.record)], dist.record)
  }
  
  rownames(dist_pass) <- c("Dims", "BiasBound")
  
  #for now a crude warning if pass K.svd in for L1 distance
  if(!is.null(K.svd)) {
      warning("Please note that the L1 distance is calculated on the svd of the kernel matrix passed in as \"K.svd\".", immediate. = TRUE)
  }
  R=list()
  R$w= getw.out$w
  R$biasbound.opt=biasbound_opt
  R$biasbound.orig=dist.orig
  R$dist.record= dist_pass
  R$numdims=numdims
  R$L1.orig = L1_orig
  R$L1.opt = L1_optim
  R$K = K
  R$earlyfail = getw.out$earlyfail
  R$linkernel = linkernel
  R$svdK = svd.out

  return(R)
} # end kbal main function



