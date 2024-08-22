### Helper functions

### Build kernel, possibly non-square
#' Build the Gaussian Kernel Matrix
#'
#' @description Builds the Gaussian kernel matrix using Rcpp.
#' @param allx a data matrix containing all observations where rows are units and columns are covariates.
#' @param useasbases a binary vector with length equal to the number of observations (rows in \code{allx}) to specify which bases to use when constructing the kernel matrix (columns of \eqn{K}). If not specified, the default is to use all observations.
#' @param b Scaling factor in the calculation of Gaussian kernel distance equivalent to the entire denominator \eqn{2\sigma^2} of the exponent. Default is twice the number of covariates or columns in \code{allx}.
#' @param linkernel a logical value indicating whether to use a linear kernel, \eqn{K=XX'}, which in practice employs \eqn{K=X}.  
#' @param scale a logical value indicating whether to standardize \code{allx} (demeaned with sd=1) before constructing the kernel matrix 
#' @return \item{K}{The kernel matrix}
#' @examples
#' #load and clean data a bit
#' \donttest{
#' data(lalonde)
#' xvars <- c("age","black","educ","hisp","married","re74","re75","nodegr","u74","u75")
#'
#' #note that lalonde$nsw is the treatment vector, so the observed is 1-lalonde$nsw
#' #running makeK with the sampled/control units as the bases given 
#' #the large size of the data
#' K <- makeK(allx = lalonde[,xvars], useasbases = 1-lalonde$nsw) 
#' }
#' @useDynLib kbal
#' @importFrom stats sd 
#' @importFrom Rcpp sourceCpp 
#' @importFrom RcppParallel RcppParallelLibs
#' @importFrom dplyr filter '%>%' group_by summarise n
#' @export
makeK = function(allx, useasbases=NULL, b=NULL, linkernel = FALSE, scale = TRUE){

  # Error handling
  if (!is.matrix(allx)) stop("`allx` must be a matrix.")
  if (!is.null(useasbases) && (!is.numeric(useasbases) || length(useasbases) != nrow(allx) || any(!useasbases %in% c(0, 1)))) {
    stop("`useasbases` must be a binary vector with the same length as the number of rows in `allx`.")
  }
  if (!is.null(b) && (!is.numeric(b) || length(b) != 1)) stop("`b` must be a single numeric value.")
  if (!is.logical(linkernel)) stop("`linkernel` must be a logical value.")
  if (!is.logical(scale)) stop("`scale` must be a logical value.")
  if (!is.null(useasbases) && (sum(useasbases) == 0)) stop("`useasbases` must have at least one element set to 1.")
  
  N=nrow(allx)
  # If no "useasbasis" given, assume all observations are to be used.
  if(is.null(useasbases)) {useasbases = rep(1, N)}

  single.base = ( sum(useasbases)==1 )
  if(single.base) {
    base.is.one = ( which(useasbases==1)==1 )
    if(base.is.one){
      useasbases[2] = 1
      }else{
      useasbases[1] = 1
      }}
  
  #default b is set to 2ncol to match kbal for now
  if (is.null(b)){ b=2*ncol(allx) }
  
  if(scale) {
      allx = scale(allx)
  } 
  bases = allx[useasbases==1, ]
  
  #removed scaling based on bases and just rescaled all of allx then subsetted
  #Xmeans.bases <- colMeans(bases)
  #Xsds.bases <- apply(bases,2,sd)
  #bases <- scale(bases, center = Xmeans.bases, scale = Xsds.bases)
  #allx <- scale(allx, center = Xmeans.bases, scale = Xsds.bases)
  #bases <- scale(bases)
  #allx <- scale(allx)
  
  if(linkernel == TRUE) {
      K = allx
  } else {
      if(sum(useasbases) == N) {
          #symmetric K, build faster using
          K = kernel_parallel(X = allx, b = b)
      } else {
          #updated to build only triangle and mirror (4x faster)
          K = kernel_parallel_2(X = allx, Y = bases, b = b)
          #K = kernel_parallel_old(X = allx, Y = bases, b = b)

          if(single.base) {
            if(base.is.one){
              K = as.matrix(K[,1])
            }else{
              K = as.matrix(K[,-1])
            }}
          
      }
            #old non parallel
          #new_gauss_kern(newx = allx, oldx = bases, b = b)
  }
  return(K)
}

### Get bias bound, which will be the distance
### function to drive the optimization
#' Worst-Case Bias Bound due to Incomplete Balance
#' @description Calculate the upper bound on the bias induced by approximate balance with a given \code{hilbertnorm}. Approximate balance is conducted in \code{kbal()} and uses only the first \code{numdims} dimensions of the singular value decomposition of the kernel matrix to generate weights \code{w} which produce mean balance between control or sampled units and treated or population units. The following function calculates the worse-case bias induced by this approximate balancing with weights \code{w} and a given \code{hilbertnorm}.
#' @param observed a numeric vector of length equal to the total number of units where sampled/control units take a value of 1 and population/treated units take a value of 0.
#' @param target a numeric vector of length equal to the total number of units where population/treated units take a value of 1 and sample/control units take a value of 0.
#' @param svd.out the list object output from \code{svd()} performed on the kernel matrix. Requires a list object with left singular vectors in \code{svd.out$u} and singular values in \code{svd.out$d}
#' @param w numeric vector containing the weight for every corresponding unit. Note that these weights should sum to the total number of units, not to one. They are divided by the number of control or sample and treated or population units internally.
#' @param w.pop an optional vector input to specify population weights. Must be of length equal to the total number of units (rows in \code{svd.out}) with all sampled units receiving a weight of 1. The sum of the weights for population units must be either 1 or the number of population units.
#' @param hilbertnorm numeric value of the hilbertnorm.
#' @return \item{biasbound}{value of worst-case bias bound due to incomplete balance with inputted weights}
#' @examples
#' \donttest{
#' #load and clean data a bit
#' data(lalonde)
#' xvars=c("age","black","educ","hisp","married","re74","re75","nodegr","u74","u75")
#'
#' #need a kernel matrix to run SVD on and pass in so get that first with makeK
#' #running makeK with the sampled units as the bases
#' K = makeK(allx = lalonde[,xvars], useasbases = 1-lalonde$nsw)
#'
#' #svd on this kernel
#' svd_pass = svd(K)  
#' #let's use the original weights of 1/number of sampled units, and 1/number of target units
#' #this is the default if we pass in w as all 1's
#' biasbound(observed=(1-lalonde$nsw),
#'           target=lalonde$nsw,
#'           svd.out = svd_pass,
#'           w = rep(1,nrow(lalonde)), hilbertnorm=1)
#'  }
#' @export
biasbound=function(observed, target, svd.out, w, w.pop = NULL,
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
        #sampledinpop == TRUE, check that w.pop = 1 for all sampled/control units
        if(sum(target) == length(target) && !(sum(w.pop[observed]) == sum(observed)
                                             & sd(w.pop[observed]) == 0)) {
            stop("\"w.pop\" must the value 1 for all sampled/control units")
        }
        #check population weights sum to num of treated/population units
        if(round(sum(w.pop[target ==1 ])) != sum(target)) {
            #allow user to pass in weights that sum to one and transform them here
            if(round(sum(w.pop[target ==1])) == 1) {
                w.pop[target==1] = w.pop[target ==1]/mean(w.pop[target==1])
            } else { #in this case they don't sum to N_t or 1 so ng
                stop("\"w.pop\" must sum to either 1 or the number of treated/population units")
            }
        }
    }
    #sampledinpop = FALSE   
    if(sum(target) != length(target)) {
        #weights sum to N_0, N_1 so normalizing so sum to 1
        wtarget=w.pop[target==1]/sum(target==1) 
        wobserved=w[observed==1]/sum(observed==1)
    } else { #when sampledinpop true
        wtarget= w.pop/N #over N because sum(target ==1) = 1
        wobserved=w[observed ==1]/sum(observed==1)
   }
   
    U=svd.out$u
    eigenvals=svd.out$d #singular values (A)
    if(sum(sign(eigenvals) == -1) != 0) {
        stop("Encountered negative eigenvalues. Cannot compute biasbound.")
    }
    
    U1=U[target==1, , drop=FALSE]
    U0=U[observed==1, , drop=FALSE]
    eigenimbal=as.vector(t(wtarget)%*%U1 - t(wobserved)%*%U0)
    effectiveimbal=(eigenimbal*(eigenvals^.5))
    biasbound=sqrt(hilbertnorm)*sqrt(t(effectiveimbal)%*%(effectiveimbal))
    return(biasbound)
}

# Simple difference in mean and in weighted means
# (Not actually used right now but can be convenient)
#' Difference in Means and Difference in Weighted Means
#'
#' Calculates the simple difference in means or weighted difference in means between the control or sample population and the treated or target population.
#' @param X matrix of data where rows are observations and columns are covariates.
#' @param w numeric vector of weights for each observation.
#' @param target numeric vector of length equal to the total number of units where population/treated units take a value of 1 and sample/control units take a value of 0.
#' @return \item{dim}{the simple, unweighted difference in means.}
#' \item{dimw}{the weighted difference in means.}
#' @examples
#' \donttest{
#' #let's say we want to get the unweighted DIM and the weighted DIM using weights from the kbal
#' #function with the lalonde data:
#' #load and clean data a bit
#' data(lalonde)
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
#' @description Uses entropy balancing to find and return the weights that produce mean balance on \eqn{\phi(X_i)}, the expanded features of \eqn{X_i} using a given kernel \eqn{\phi(.)}, for the control or sample group and treated group or target population.
#'
#' @param target a numeric vector of length equal to the total number of units where population/treated units take a value of 1 and sample/control units take a value of 0.
#' @param observed a numeric vector of length equal to the total number of units where sampled/control units take a value of 1 and population/treated units take a value of 0.
#' @param svd.U a matrix of left singular vectors from performing \code{svd()} on the kernel matrix.
#' @param ebal.tol tolerance level used by custom entropy balancing function \code{ebalance_custom}
#' @param ebal.maxit maximum number of iterations in optimization search used by \code{ebalance_custom}
#' @return \item{w}{numeric vector of weights.}
#' \item{converged}{boolean indicating if \code{ebalance_custom} converged}
#' \item{ebal_error}{returns error message if \code{ebalance_custom} encounters an error}
#' @examples
#' \donttest{
#' #load and clean data a bit
#' data(lalonde)
#' xvars=c("age","black","educ","hisp","married","re74","re75","nodegr","u74","u75")
#'
#' #need a kernel matrix to run SVD on then find weights with so get that first with makeK
#' #running makeK with the sampled units as the bases
#' K = makeK(allx = lalonde[,xvars], useasbases = 1-lalonde$nsw)
#'
#' #svd on this kernel and get matrix with left singular values
#' U = svd(K)$u
#' #usually search over all columns of U to find which produces weights with 
#' #the minimum bias bound; in this ex, search only first 10 dims
#' U2=U[,1:10]
#' getw.out=getw(target=lalonde$nsw, 
#'               observed=1-lalonde$nsw, 
#'               svd.U=U2)
#'  }
#' @export
getw = function(target, observed, svd.U, ebal.tol=1e-6, ebal.maxit = 500){
    
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
                                   max.iterations = ebal.maxit,
                                   constraint.tolerance = ebal.tol,
                                   print.level=0),
                   silent=TRUE)
  N=nrow(svd.U)
  converged = FALSE
  #earlyfail = FALSE
  error = NULL
  if ("try-error"%in%class(bal.out.pc)[1]){
          warning("\'ebalace_custom()\' encountered an error. Returning equal weights. See \"ebal_error\" for details. ", 
                  immediate. = T)
      error = bal.out.pc[1]
      
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
    converged = bal.out.pc$converged
  }
 
    out <- list(w = w, 
                converged=converged, 
                ebal_error = error)
  return(out)
} # end of getw.

#' L1 Distance
#' @description Calculates the L1 distance between the treated or population units and the kernel balanced control or sampled units.
#' @param target a numeric vector of length equal to the total number of units where population/treated units take a value of 1 and sample/control units take a value of 0.
#' @param observed a numeric vector of length equal to the total number of units where sampled/control units take a value of 1 and population/treated units take a value of 0.
#' @param K the kernel matrix
#' @param w.pop an optional vector input to specify population weights. Must be of length equal to the total number of units (rows in \code{svd.U}) with all sampled units receiving a weight of 1. The sum of the weights for population units must be either 1 or the number of population units.
#' @param w a optional numeric vector of weights for every observation. Note that these weights should sum to the total number of units, where treated or population units have a weight of 1 and control or sample units have appropriate weights derived from kernel balancing with mean 1, is consistent with the ouput of \code{getw()}. If unspecified, these weights are found internally using \code{numdims} dimensions of the SVD of the kernel matrix \code{svd.U} with \code{ebalance_custom()}. 
#' @param numdims an optional numeric input specifying the number of columns of the singular value decomposition of the kernel matrix to use when finding weights when \code{w} is not specified.
#' @param ebal.tol an optional numeric input specifying the tolerance level used by custom entropy balancing function \code{ebalance_custom()} in the case that \code{w} is not specified. 
#' @param ebal.maxit maximum number of iterations in optimization search used by \code{ebalance_custom} when \code{w} is not specified. 
#' @param svd.U an optional matrix of left singular vectors from performing \code{svd()} on the kernel matrix in the case that \code{w} is unspecified. If unspecified when \code{w} also not specified, internally computes the svd of \code{K}.
#' @return \item{L1}{a numeric giving the L1 distance, the absolute difference between \code{pX_D1} and \code{pX_D0w}}
#' \item{w}{numeric vector of weights used}
#' \item{pX_D1}{a numeric vector of length equal to the total number of observations where the nth entry is the sum of the kernel distances from the nth unit to every treated or population unit. If population units are specified, this sum is weighted by \code{w.pop} accordingly.}
#' \item{pX_D0}{a numeric vector of length equal to the total number of observations where the nth entry is the sum of the kernel distances from the nth unit to every control or sampled unit.}
#' \item{pX_D0w}{a numeric vector of length equal to the total number of observations where the nth entry is the weighted sum of the kernel distances from the nth unit to every control or sampled unit. The weights are given by entropy balancing and produce mean balance on \eqn{\phi(X)}, the expanded features of \eqn{X} using a given kernel \eqn{\phi(.)}, for the control or sample group and treated group or target population.}
#' @examples
#' \donttest{
#' #loading and cleaning lalonde data
#' data(lalonde)
#' xvars=c("age","black","educ","hisp","married","re74","re75","nodegr","u74","u75")
#'
#' #need to first build gaussian kernel matrix
#' K_pass <- makeK(allx = lalonde[,xvars])
#' #also need the SVD of this matrix
#' svd_pass <- svd(K_pass)
#'
#' #running without passing weights in directly, using numdims=33
#' l1_lalonde <- getdist(target = lalonde$nsw,
#'                       observed = 1-lalonde$nsw,
#'                       K = K_pass,
#'                       svd.U = svd_pass$u,
#'                       numdims = 33)
#'
#'  #alternatively, we can get the weights ourselves and pass them in directly
#'  #using the first 33 dims of svd_pass$u to match the above
#' w_opt <- getw(target= lalonde$nsw,
#'               observed = 1-lalonde$nsw,
#'               svd.U = svd_pass$u[,1:33])$w
#' l1_lalonde2 <- getdist(target = lalonde$nsw,
#'                  observed = 1-lalonde$nsw,
#'                  K = K_pass,
#'                  w = w_opt)
#' }
#' @export
getdist <- function(target, observed, K, w.pop = NULL, 
                    w=NULL, numdims = NULL, 
                    ebal.tol= 1e-6, 
                    ebal.maxit = 500,
                    svd.U = NULL) {

        
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
            if(!(sum(w.pop[observed==1]) == sum(observed) & 
                 (sum(observed) ==1 | sd(w.pop[observed==1]) == 0) ) ) {
                stop("\"w.pop\" must the value 1 for all sampled units")
            }
            #check population weights sum to num of treated/population units
            if(round(sum(w.pop[target ==1 ])) != sum(target)) {
                #allow user to pass in weights that sum to one and transform them here
                if(round(sum(w.pop[target ==1])) == 1) {
                    w.pop[target==1] = w.pop[target ==1]/mean(w.pop[target==1])
                } else { #in this case they don't sum to N_t or 1 so ng
                    stop("\"population.w\" must sum to either 1 or the number of treated/population units")
                }
            }
            
        }
        #if user does not provide weights, go get them
        if(is.null(w)) {
            if(is.null(numdims)) {stop("If weights w input is not specified, numdims must be in order to calculate these weights internally")}
            if(is.null(svd.U)) {
                svd.U = svd(K)$u
            }
            U_w.pop <- w.pop*svd.U
            w = suppressWarnings(getw(target = target, observed=observed,
                     svd.U = U_w.pop[,1:numdims, drop=FALSE],
                     ebal.tol=ebal.tol, ebal.maxit = ebal.maxit)$w)

            #if ebal fails we get weights of 1 for everyone
            if(sum(w ==1) == length(w)){
                warning("ebalance failed for this choice of numdims dimensions of the SVD of the kernel matrix returning equal weights (1) for all units", immediate. = T)
            }
        }
        #just "average row Kt"
        #pX_D1=(matrix(1,1, sum(target==1))/sum(target==1))%*% K_t 
        #average row Kc
        pX_D0=(matrix(1,1,sum(observed==1))/sum(observed==1))%*% K_c 
        #weighted average Kc with ebal weights
        pX_D0w=(w[observed==1]/sum(w[observed==1])) %*% K_c
        #weighted average Kt ONLY DIFF from pX_D1 if have pop weights not all equal to 1
        pX_D1wpop = (w.pop[target==1]/sum(w.pop[target==1])) %*% K_t 

        # A rescaling to ensure sum is 1 as would be an integral.
        #pX_D1=pX_D1/sum(pX_D1)
        pX_D0=pX_D0/sum(pX_D0)
        pX_D0w=pX_D0w/sum(pX_D0w)
        pX_D1wpop =pX_D1wpop/sum(pX_D1wpop)
        L1 = sum(abs(pX_D1wpop-pX_D0w)) #removed the 0.5 factor -- Oct 20 2017

        R=list()
        R$L1=L1
        R$w=w
        R$pX_D1=pX_D1wpop
        R$pX_D0=pX_D0
        R$pX_D0w=pX_D0w

        return(R)
} ## end of getdist


#' One-Hot Encoding for Categorical Data
#' @description Converts raw categorical string/factor data matrix into numeric one-hot encoded data matrix. Intended to help prepare data to be passed to \code{kbal} argument \code{allx} when categorical data is used. 
#' @param data a dataframe or matrix where columns are string or factor type covariates
#' @return \item{onehot_data}{a matrix of combined sample and population data with rows corresponding to units and columns one-hot encoded categorical covariates}
#' @examples
#' \donttest{
#' #Ex 1. Make up some categorical demographic data
#' dat = data.frame(pid = c(rep("Rep", 20),
#'                          rep("Dem", 20), 
#'                          rep("Ind", 20)), 
#'                  gender = c(rep("female", 35),
#'                             rep("male", 25)))
#' #Convert to one-hot encoded data matrix:
#' onehot_dat = one_hot(dat)
#' }
#' #Ex 2. lalonde data
#' data(lalonde)
#' cat_vars=c("black","hisp","married","nodegr","u74","u75")
#' onehot_lalonde = one_hot(lalonde[, cat_vars])
#' @importFrom stats model.matrix contrasts
#' @export
one_hot <- function(data) {
    onehot_data <- data.frame(lapply(data.frame(data),as.factor))
    onehot_data <- model.matrix(~ ., onehot_data,
                                contrasts.arg = lapply(onehot_data[,,drop = F], 
                                                       contrasts, 
                                                       contrasts=FALSE))
    onehot_data <- onehot_data[, -1]
    
    return(onehot_data)
}


#' Maximum Variance of Gaussian Kernel Matrix
#' @description Searches for the argmax of the variance of the Kernel matrix
#' @param data a matrix of data where rows are all units and columns are covariates. Where all covariates are categorical, this matrix should be one-hot encoded (refer to \code{\link{one_hot}} to produce) with \code{cat_data} argument true.
#' @param useasbases binary vector specifying what observations are to be used in forming bases (columns) of the kernel matrix. Suggested default is: if the number of observations is under 4000, use all observations; when the number of observations is over 4000, use the sampled (control) units only.
#' @param cat_data logical for whether kernel contains only categorical data or not
#' @param maxsearch_b the maximum value of \eqn{b}, the denominator of the Gaussian, searched during maximization.
#' @return \item{b_maxvar}{numeric \eqn{b} value, the denominator of the Gaussian, which produces the maximum variance of \eqn{K} kernel matrix}
#' \item{var_K}{numeric maximum variance of \eqn{K} kernel matrix found with \eqn{b} as \code{b_maxvar}}
#' @examples
#' \donttest{
#' #lalonde with only categorical data
#' data(lalonde)
#' cat_vars <- c("black","hisp","married","nodegr","u74","u75")
#' #Convert to one-hot encoded data matrix:
#' onehot_lalonde = one_hot(lalonde[, cat_vars])
#' colnames(onehot_lalonde)
#' best_b <- b_maxvarK(data = onehot_lalonde, 
#'                     useasbases = 1-lalonde$nsw) 
#'  }
#' @importFrom dplyr '%>%' pull group_by summarise n
#' @importFrom stats optimize na.omit
#' @export
b_maxvarK <- function(data,
                      useasbases,
                      cat_data = TRUE,
                      maxsearch_b = 2000) {

    # Error handling
    if (!is.matrix(data)) stop("`data` must be a matrix.")
    if (!is.numeric(useasbases) || length(useasbases) != nrow(data)) stop("`useasbases` must be a binary vector with the same length as the number of rows in `data`.")
    if (!is.logical(cat_data)) stop("`cat_data` must be a logical value.")
    if (!is.numeric(maxsearch_b) || length(maxsearch_b) != 1) stop("`maxsearch_b` must be a single numeric value.")
    
    #categorical kernel + b range:
    #get raw counts:
    if(cat_data) {
        K <- makeK(data, b=1,
               useasbases = useasbases,
               linkernel = FALSE, scale = FALSE)
        raw_counts <- -log(K)
        n_d <- data.frame(diff = c(raw_counts)) %>% group_by(diff) %>% summarise(n())
        
        #internal function to get the variance of K 
        var_K= function(b, n_d, diag_length){
            d <- n_d[,1] %>% pull()
            n_d <- as.vector(n_d[,2] %>% pull())
            #REMOVING DIAGONAL 0 COUNTS FROM MAXIMIZATION CONSIDERATION
            n_d[1] <- n_d[1] - diag_length
            p_d <- n_d/ sum(n_d) 
            
            mean_k = sum(exp(-1*d/b)*p_d)
            var_k = sum((exp(-1*d/b)-mean_k)^2 * p_d)
            return(var_k)
        }
        
        #does this diag technique work for all shapes of K?
        res = optimize(var_K, n_d, length(diag(K)),
                       interval=c(0, maxsearch_b), maximum=TRUE)
    } else {
        var_K= function(b, data){
            #makeK
            #option 1: do not divide by 2 and scale X cont BEFORE HAND (outside func)
            K <- makeK(data, b=b,
                       useasbases = useasbases,
                       linkernel = FALSE,
                       scale = FALSE)
            #REMOVING DIAGONAL 0 COUNTS FROM MAXIMIZATION CONSIDERATION 
            diag(K) <- NA
            #old: memory issues with certain installations... R is messed up
            #var_k <- var(na.omit(as.vector(K)))
            #some benchmarking showed this seems to be the fastest method
            n = nrow(K)*ncol(K) - sum(useasbases)
            #to match R's var calc 1st denom here needs to be n-1 but to match above var calc
            #using 1/n
            var_k = (1/(n))*(sum(K^2, na.rm = T) - (1/n)*sum(K, na.rm = T)^2)
            return(var_k)
        }
        res = optimize(var_K, data,
                       interval=c(0, maxsearch_b), maximum=TRUE)
    }
    return(list(b_maxvar = res$maximum, 
                var_K = res$objective))
}




#' Drop Multicollinear Columns
#' @description Drops multicollinear columns in order of highest correlation
#' @param allx a matrix of data to check for multicollinearity
#' @param printprogress logical to indicate if progress should be printed out to command line
#' @return \item{allx_noMC}{resulting data matrix of full rank after multicolliear columns have been dropped}
#' \item{dropped_cols}{column names of those dropped}
#' @examples
#' \donttest{
#' # make data with multicollinearity 
#' data <- data.frame(x = rnorm(100),
#'                    y = sample.int(100,100), 
#'                    z = runif(100, 3,6))
#' test = data.frame(mc_1 = data$x,
#'                   mc_2 = data$x*2 + data$y - data$z)
#' dat = cbind(test, data)
#' #run
#' mc_check= drop_multicollin(dat)
#' mc_check$dropped_cols 
#' }
#' @export
drop_multicollin <- function(allx, printprogress = TRUE) {
    
    qr_X = qr(allx)
    multicollin = FALSE
    if(qr_X$rank < ncol(allx)) {
        multicollin = TRUE
    }
    if(multicollin & printprogress) {
        cat("Dropping detected multicollinear columns\n")
    }
    allx_update = allx
    dropped_cols = NULL
    cor = cor(allx_update)
    diag(cor) = 0
    cor[lower.tri(cor)] = 0
    cor = abs(cor)
    all_cor <- sort(c(cor), decreasing = TRUE)
    i = 1
    rank_target = qr(allx)$rank
    while(multicollin == TRUE){
        drop = which(cor == all_cor[i], arr.ind = T)[1,1]
        
        if(qr(allx_update[,-drop])$rank == rank_target) {
            if(!is.null(rownames(which(cor == all_cor[i],
                                       arr.ind  =TRUE)))) {
                dropped_cols = c(dropped_cols, rownames(which(cor == all_cor[i],
                                                              arr.ind  =TRUE))[1])
            } else {
                dropped_cols = c(dropped_cols, paste("column", (which(cor == all_cor[i],
                                                              arr.ind  =TRUE))[1]))
            }
            
            allx_update <- allx_update[,-drop, drop = F]
        }
        #cat(i, drop, ncol(allx_update),"\n")
        i = i + 1
        if(qr_X$rank == ncol(allx_update)) {multicollin = FALSE}
        
    }
    allx = allx_update
        
    return(list(allx_noMC = allx, 
                dropped_cols = dropped_cols))
}





# The main event: Actual kbal function!
#' Kernel Balancing
#'
#' @description Kernel balancing (\code{kbal}) is non-parametric weighting tool to make two groups have a similar distribution of covariates, not only in terms of means or marginal distributions but also on (i) general smooth functions of the covariates, including on (ii) a smoothing estimator of the joint distribution of the covariates. It was originally designed (Hazlett, 2017) to make control and treated groups look alike, as desired when estimating causal effects under conditional ignorabiity. This package also facilitates use of this approach for more general distribution-alignment tasks, such as making a sampled group have a similar distribution of covariates as a target population, as in survey reweighting. The examples below provide an introduction to both settings.
#' 
#' To proceed in the causal effect setting, kbal assumes that the expectation of the non-treatment potential outcome conditional on the covariates falls in a large, flexible space of functions associated with a kernel. It then constructs linear bases for this function space and achieves approximate balance on these bases. The approximation is one that minimizes the worst-case bias that could persist due to remaining imbalances. 
#' 
#' The \code{kbal} function implements kernel balancing using a gaussian kernel to expand the features of \eqn{X_i} to infinite dimensions.  It finds approximate mean balance for the control or sample group and treated group or target population in this expanded feature space by using the first \code{numdims} dimensions of the singular value decomposition of the gaussian kernel matrix. It employs entropy balancing to find the weights for each unit which produce this approximate balance. When \code{numdims} is not user-specified, it searches through increasing dimensions of the SVD of the kernel matrix to find the number of dimensions which produce weights that minimizes the worst-case bias bound with a given \code{hilbertnorm}. It then returns these optimal weights, along with the minimized bias, the kernel matrix, a record of the number of dimensions used and the corresponding bais, as well as an original bias using naive group size weights for comparison. Note that while kernel balancing goes far beyond simple mean balancing, it may not result in perfect mean balance. Users who wish to require mean balancing can specify \code{meanfirst = T} to require mean balance on as many dimensions of the data as optimally feasible. Alternatively, users can manually specify \code{constraint} to append additional vector constraints to the kernel matrix in the bias bound optimization, requiring mean balance on these columns. Note further that \code{kbal} supports three types of input data: fully categorical, fully continuous, or mixed. When data is only categorical, as is common with demographic variables for survey reweighting, users should use argument \code{cat_data = TRUE} and can input their data as factors, numeric, or characters and \code{kbal} will internally transform the data to a more appropriate one-hot encoding and search for the value of \code{b}, the denominator of the exponent in the Gaussian, which maximizes the variance of the kernel matrix. When data is fully continuous, users should use default settings (\code{cat_data = FALSE} and \code{cont_data = FAlSE}, which will scale all columns and again conduct an internal search for the value of \code{b} which maximizes the variance of \code{K}. Note that with continuous data, this search may take considerably more computational time than the categorical case. When data is a mix of continuous and categorical data, users should use argument \code{mixed_data = TRUE}, specify by name what columns are categorical with \code{cat_columns}, and also set the scaling of the continuous variables with \code{cont_scale}. This will result in a one-hot encoding of categorical columns concatenated with the continous columns scaled in accordance with \code{cont_scale} and again an internal search for the value of \code{b} which maximizes the variance in the kernel matrix. Again note that compared to the categorical case, this search will take more computational time. 
#' 
#'
#' @references Hazlett, C. (2017), "Kernel Balancing: A flexible non-parametric weighting procedure for estimating causal effects." Forthcoming in Statistica Sinica. https://doi.org/10.5705/ss.202017.0555
#' 
#' @param allx a data matrix containing all observations where rows are units and columns are covariates. When using only continuous covariates (\code{cat_data = F} and \code{mixed_data = F}), all columns must be numeric. When using categorical data (either \code{cat_data = T} or \code{mixed_data = T}), categorical columns can be characters or numerics which will be treated as factors. Users should one-hot encoded categorical covariates as this transformation occurs internally.
#' @param useasbases optional binary vector to specify what observations are to be used in forming bases (columns) of the kernel matrix to get balance on.  If the number of observations is under 4000, the default is to use all observations. When the number of observations is over 4000, the default is to use the sampled (control) units only.
#' @param b scaling factor in the calculation of gaussian kernel distance equivalent to the entire denominator \eqn{2\sigma^2} of the exponent. Default is to search for the value which maximizes the variance of the kernel matrix.
#' @param sampled a numeric vector of length equal to the total number of units where sampled units take a value of 1 and population units take a value of 0.
#' @param sampledinpop a logical to be used in combination with input \code{sampled} that, when \code{TRUE}, indicates that sampled units should also be included in the target population when searching for optimal weights.
#' @param treatment an alternative input to \code{sampled} and \code{sampledinpop} that is a numeric vector of length equal to the total number of units. Current version supports the ATT estimand. Accordingly, the treated units are the target population, and the control are equivalent to the sampled. Weights play the role of making the control groups (sampled) look like the target population (treated). When specified, \code{sampledinpop} is forced to be \code{FALSE}.
#' @param population.w optional vector of population weights length equal to the number of population units. Must sum to either 1 or the number of population units.
#' @param K optional matrix input that takes a user-specified kernel matrix and performs SVD on it internally in the search for weights which minimize the bias bound.
#' @param K.svd optional list input that takes a user-specified singular value decomposition of the kernel matrix. This list must include three objects \code{K.svd$u}, a matrix of left-singular vectors, \code{K.svd$v}, a matrix of right-singular vectors, and their corresponding singular values \code{K.svd$d}. 
#' @param cat_data logical argument that when true indicates \code{allx} contains only categorical data. When true, the internal construction of the kernel matrix uses a one-hot encoding of \code{allx} (multiplied by a factor of \eqn{\sqrt{0.5}} to compensate for double counting) and the value of \code{b} which maximizes the variance of this kernel matrix. When true, \code{mixed_data}, \code{scale_data}, \code{linkernel}, and \code{drop_MC} should be \code{FALSE}.
#' @param mixed_data logical argument that when true indicates \code{allx} contains a combination of both continuous and categorical data. When true, the internal construction of the kernel matrix uses a one-hot encoding of the categorical variables in \code{allx} as specified by \code{cat_columns} (multiplied by a factor of \eqn{\sqrt{0.5}} to compensate for double counting) concatenated with the remaining continuous variables scaled to have default standard deviation of 1 or that specified in \code{cont_scale}. When both \code{cat_data} and \code{cat_data} are \code{FALSE}, the kernel matrix assumes all continuous data, does not one-hot encode any part of \code{allx} but still uses the value of \code{b} which produces maximal variance in \code{K}.
#' @param cat_columns optional character argument that must be specified when \code{mixed_data} is \code{TRUE} and that indicates what columns of \code{allx} contain categorical variables. 
#' @param cont_scale optional numeric argument used when \code{mixed_data} is \code{TRUE} which specifies how to scale the standard deviation of continuous variables in \code{allx}. Can be either a a single value or a vector with length equal to the number of continuous variables in \code{allx} (columns not specified in \code{cat_columns}) and ordered accordingly.
#' @param scale_data logical when true scales the columns of \code{allx} (demeans and scales variance to 1) before building the kernel matrix internally. This is appropriate when \code{allx} contains only continuous variables with different scales, but is not recommended when \code{allx} contains any categorical data. Default is \code{TRUE} when both \code{cat_data} and \code{mixed_data} are \code{FALSE} and \code{FALSE} otherwise.
#' @param drop_MC logical for whether or not to drop multicollinear columns in \code{allx} before building \code{K}. When either \code{cat_data} or \code{mixed_data} is \code{TRUE}, forced to be \code{FALSE}. Otherwise, with continuous data only, default is \code{TRUE}.
#' @param linkernel logical if true, uses the linear kernel \eqn{K=XX'} which achieves balance on the first moments of \eqn{X} (mean balance). Note that for computational ease, the code employs \eqn{K=X} and adjusts singular values accordingly.
#' @param meanfirst logical if true, internally searches for the optimal number of dimensions of the svd of \code{allx} to append to \code{K} as additional constraints. This will produce mean balance on as many dimensions of \code{allx} as optimally feasible with specified ebalance convergence and a minimal bias bound on the remaining unbalances columns of the left singular vectors of \code{K}. Note that any scaling specified on \code{allx} will be also be applied in the meanfirst routine.
#' @param mf_columns either character or numeric vector to specify what columns of \code{allx} to perform meanfirst with. If left unspecified, all columns will be used. 
#' @param constraint optional matrix argument of additional constraints which are appended to the front of the left singular vectors of \code{K}. When specified, the code conducts a constrained optimization requiring mean balance on the columns of this matrix throughout the search for the minimum bias bound over the dimensions of the left singular vectors of \code{K}. 
#' @param scale_constraint logical for whether constraints in \code{constraint} should be scaled before they are appended to the svd of \code{K}.
#' @param numdims optional numeric argument specifying the number of dimensions of the left singular vectors of the kernel matrix to find balance bypassing the optimization search for the number of dimensions which minimize the biasbound.
#' @param minnumdims numeric argument to specify the minimum number of the left singular vectors of the kernel matrix to seek balance on in the search for the number of dimensions which minimize the bias. Default minimum is 1.
#' @param maxnumdims numeric argument to specify the maximum number of the left singular vectors of the kernel matrix to seek balance on in the search for the number of dimensions which minimize the bias. For a Gaussian kernel, the default is the minimum between 500 and the number of bases given by \code{useasbases}. With a linear kernel, the default is the minimum between 500 and the number of columns in \code{allx}. 
#' @param fullSVD logical argument for whether the full SVD should be conducted internally. When \code{FALSE}, the code uses truncated svd methods from the \code{Rspectra} package in the interest of improving run time. When \code{FALSE}, the code computes only the SVD up to the either 80 percent of the columns of \code{K} or \code{maxnumdims} singular vectors, whichever is larger. When the number of columns is less thanm 80 percent the  number of rows, defaults to full svd.
#' @param incrementby numeric argument to specify the number of dimensions to increase by from \code{minnumdims} to \code{maxnumdims} in each iteration of the search for the number of dimensions which minimizes the bias. Default is 1.
#' @param ebal.maxit maximum number of iterations used by \code{ebalance_custom()} in optimization in the search for weights \code{w}.
#' @param ebal.tol tolerance level used by \code{ebalance_custom()}. 
#' @param ebal.convergence logical to require ebalance convergence when selecting the optimal \code{numdims} dimensions of \code{K} that minimize the biasbound. When constraints are appended to the left singular vectors of \code{K} via \code{meanfirst=TRUE} or \code{constraints}, forced to be \code{TRUE} and otherwise \code{FALSE}.
#' @param maxsearch_b optional argument to specify the maximum b in search for maximum variance of \code{K} in \code{b_maxvarK()}.
#' @param early.stopping logical argument indicating whether bias balance optimization should stop twenty rounds aftering finding a minimum
#' @param printprogress logical argument to print updates throughout.
#'
#' @return \item{w}{a vector of the weights found using entropy balancing on \code{numdims} dimensions of the SVD of the kernel matrix.}
#' \item{biasbound_opt}{a numeric giving the minimal bias bound found using \code{numdims} as the number of dimesions of the SVD of the kernel matrix. When \code{numdims} is user-specified, the bias bound using this number of dimensions of the kernel matrix.}
#'  \item{biasbound_orig}{a numeric giving the bias bound found when all sampled (control) units have a weight equal to one over the number of sampled (control) units and all target units have a weight equal to one over the number of target units.}
#'  \item{biasbound_ratio}{a numeric giving the ratio of \code{biasbound_orig} to\code{biasbound_opt}. Can be informative when comparing the performance of different \code{b} values.} 
#'  \item{dist_record}{a matrix recording the bias bound corresponding to balance on increasing dimensions of the SVD of the kernel matrix starting from \code{minnumdims} increasing by \code{incrementby} to \code{maxnumdims} or until the bias grows to be 1.25 times the minimal bias found.}
#'  \item{numdims}{a numeric giving the optimal number of dimensions of the SVD of the kernel matrix which minimizes the bias bound.}
#'  \item{L1_orig}{a numeric givingthe L1 distance found when all sampled (control) units have a weight equal to one over the number of sampled (control) units and all target units have a weight equal to one over the number of target units.}
#'  \item{L1_opt}{a numeric giving the L1 distance at the minimum bias bound found using \code{numdims} as the number of dimensions of the SVD of the kernel matrix. When \code{numdims} is user-specified, the L1 distance using this number of dimensions of the kernel matrix.}
#'  \item{K}{the kernel matrix}
#'  \item{onehot_dat}{when categorical data is specified, the resulting one-hot encoded categorical data used in the construction of \code{K}. When mixed data is specified, returns concatenated one-hot encoded categorical data and scaled continuous data used to construct \code{K}.}
#'  \item{linkernel}{logical for whether linear kernel was used}
#'  \item{svdK}{a list giving the SVD of the kernel matrix with left singular vectors \code{svdK$u}, right singular vectors \code{svdK$v}, and singular values \code{svdK$d}}
#'  \item{b}{numeric scaling factor used in the the calculation of gaussian kernel  equivalent to the denominator \eqn{2\sigma^2} of the exponent.}
#'  \item{maxvar_K}{returns the resulting variance of the kernel matrix when the \code{b} determined internally as the argmax of the variance \code{K}}
#'  \item{bases}{numeric vector indicating what bases (rows in \code{allx}) were used to construct kernel matrix (columns of K)}
#'  \item{truncatedSVD.var}{when truncated SVD methods are used on symmetric kernel matrices, a numeric which gives the proportion of the total variance of \code{K} captured by the first \code{maxnumdims} singular values found by the truncated SVD. When the kernel matrix is non-symmetric, this is a worst case approximation of the percent variance explained, assuming the remaining unknown singular values are the same magnitude as the last calculated in the truncated SVD.}
#'  \item{dropped_covariates}{provides a vector of character column names for covariates dropped due to multicollinearity.}
#'  \item{meanfirst_dims}{when \code{meanfirst=TRUE} the optimal number of the singular vectors of \code{allx} selected and appended to the front of the left singular vectors of \code{K}}
#'  \item{meanfirst_cols}{when \code{meanfirst=TRUE} \code{meanfirst_dims} first left singular vectors of \code{allx} selected that are appended to the front of the left singular vectors of \code{K} and balanced on}
#'  \item{ebal_error}{when ebalance is unable to find convergent weights, the associated error message it reports}
#' @examples
#' #----------------------------------------------------------------
#' # Example 1: Reweight a control group to a treated to estimate ATT. 
#' # Benchmark using Lalonde et al.
#' #----------------------------------------------------------------
#' #1. Rerun Lalonde example with settings as in Hazlett, C (2017). Statistica Sinica paper:
#' data(lalonde)
#' xvars=c("age","black","educ","hisp","married","re74","re75","nodegr","u74","u75")
#'  \donttest{
#' 
#' kbalout.full= kbal(allx=lalonde[,xvars],
#'                    b=length(xvars),
#'                    treatment=lalonde$nsw, 
#'                    fullSVD = TRUE)
#' summary(lm(re78~nsw,w=kbalout.full$w, data = lalonde))  
#'  }
#'  
#'  #2. Lalonde with categorical data only: u74, u75, nodegree, race, married
#'  cat_vars=c("race_ethnicity","married","nodegr","u74","u75")
#'  \donttest{
#'  kbalout_cat_only = kbal(allx=lalonde[,cat_vars],
#'                          cat_data = TRUE,
#'                          treatment=lalonde$nsw,
#'                          fullSVD = TRUE)
#'  kbalout_cat_only$b
#'  summary(lm(re78~nsw,w=kbalout_cat_only$w, data = lalonde))
#'  }
#'
#'  #3. Lalonde with mixed categorical and continuous data
#'  cat_vars=c("race_ethnicity", "married")
#'  all_vars= c("age","educ","re74","re75","married", "race_ethnicity")
#'  \donttest{
#'  kbalout_mixed = kbal(allx=lalonde[,all_vars],
#'                       mixed_data = TRUE, 
#'                       cat_columns = cat_vars,
#'                       treatment=lalonde$nsw,
#'                       fullSVD = TRUE)
#'  kbalout_mixed$b
#'  summary(lm(re78~nsw,w=kbalout_mixed$w, data = lalonde))
#'  }
#'  
#' #----------------------------------------------------------------
#' # Example 1B: Reweight a control group to a treated to esimate ATT. 
#' # Benchmark using Lalonde et al. -- but just mean balancing now 
#' # via "linkernel".
#' #----------------------------------------------------------------
#'
#' # Rerun Lalonde example with settings as in Hazlett, C (2017). Statistica paper:
#'kbalout.lin= kbal(allx=lalonde[,xvars],
#'                  b=length(xvars),
#'                  treatment=lalonde$nsw, 
#'                  linkernel=TRUE,
#'                  fullSVD=TRUE)
#' 
#' # Check balance with and without these weights:
#'dimw(X=lalonde[,xvars], w=kbalout.lin$w, target=lalonde$nsw)
#'
#'summary(lm(re78~nsw,w=kbalout.lin$w, data = lalonde))
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
#' # kbal correctly downweights white republicans and non-white non-republicans
#' # and upweights the non-white republicans and white non-republicans
#' unique(round(cbind(samp[,-3], k_bal_weight = kbalout$w[sampled==1]),6))
#' @export
kbal = function(allx, 
                useasbases=NULL,
                b=NULL, 
                sampled=NULL, 
                sampledinpop=NULL,
                treatment=NULL,
                population.w = NULL,
                K=NULL,
                K.svd = NULL,
                cat_data = FALSE,
                mixed_data = FALSE,
                cat_columns = NULL,
                cont_scale = NULL, 
                scale_data = NULL,
                drop_MC = NULL,
                linkernel = FALSE,
                meanfirst = FALSE,
                mf_columns = NULL,
                constraint = NULL,
                scale_constraint = TRUE,
                numdims=NULL,
                minnumdims=NULL, 
                maxnumdims=NULL,
                fullSVD = FALSE,
                incrementby=1,
                ebal.maxit = 500,
                ebal.tol=1e-6,
                ebal.convergence = NULL,
                maxsearch_b = 2000, 
                early.stopping = TRUE,
                printprogress = TRUE) {

    N=nrow(allx)
    if(is.null(N)) {
        stop("Please ensure \"allx\" is a matrix or dataframe. If subsetting to a single column, may need to use argument \"drop = F\" ")
    }
    
    # Set ebal.convergence default according to whether there are constraints or not:
    if(is.null(ebal.convergence)) {
      if(is.null(constraint) & meanfirst == FALSE) {
          ebal.convergence=FALSE
          } else {
              ebal.convergence=TRUE
          }
    }
    
#####start of big error catch series to check if data is passed in correctly and
    #default setting/data set up
    
    if(cat_data & mixed_data) {
        warning("\"cat_data\" and \"mixed_data\" should not both be TRUE. Proceeding with \"mixed_data\" TRUE only.", immediate. = TRUE)
        cat_data = FALSE
    }
    
    ####Setting Defaults: scale_data and drop_MC
    #cont data default = scaled 
    if(is.null(scale_data) & !cat_data & !mixed_data) { 
        scale_data = TRUE
    } else if(is.null(scale_data)) { #if cat data at all, default = don't scale
        scale_data = FALSE
    } 
    if(is.null(drop_MC) & !cat_data & !mixed_data) { #cont data default = drop MC
        drop_MC = TRUE
    } else if(is.null(drop_MC)) { #cat data at all default = dn drop MC
        drop_MC = FALSE
    } else if(drop_MC & (cat_data | mixed_data)) { #if user asks for drop and cat, tell them no
        warning("\"drop_MC\" should be FALSE when using categorical data. Proceeding without dropping multicollinear columns. \n", 
                immediate. = TRUE)
        drop_MC = FALSE
    } 
    if(drop_MC) {
        MC_out <- drop_multicollin(allx, printprogress = printprogress)
        dropped_cols <- MC_out$dropped_cols
        allx <- MC_out$allx_noMC
    } else {
        dropped_cols = NULL
    }
    
    
############### ERROR CHECKS ##################    
    
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
            warning("using default parameter \"sampledinpop\" = FALSE \n", immediate. = TRUE)
            sampledinpop = FALSE
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
        if(!is.null(sampledinpop) && sampledinpop == TRUE) {warning("Targeting ATT, which implies sampledinpop=FALSE.\n", immediate. = TRUE)}
        sampledinpop=FALSE
    }

    #5. checking for covariates with no variance
    if(!cat_data & !mixed_data & is.null(K.svd) & is.null(K)) {
        if(sum(apply(allx, 2, class) != "numeric") != 0) {
            stop("One or more column in \"allx\" is non-numeric while \"cat_data\" and \"mixed_data\" are both FALSE. Expecting continuous numeric data.")
        } else if( 0 %in% apply(allx, 2, sd) ) {
            stop("One or more column in \"allx\" has zero variance")
        }
    }
    #6. error catch for NAs in data
    if(sum(is.na(allx) != 0)) {
        stop("\"allx\" contains missing values")
    }
    
    ##### Setting up data: build observed and target from inputs  #####
    if(!is.null(sampled) & sampledinpop==FALSE) {
        observed = sampled
        target = 1-sampled
    } else if(!is.null(sampled)) { 
        observed = sampled
        target = rep(1,N)
    } else{ #note that we checked above for both sampled null, so this must be
        #treated passed in case
        observed = 1-treatment
        target = treatment
    }
    
    ###### Setting defaults - useasbases #####
    #7. checking useasbases if passed in
    if(!is.null(useasbases)) {
        if( (length(useasbases) != N & !linkernel) )  {
            stop("Dimensions of \"useasbases\" do not match data \"allx\"")
        }  
        if(!(all(useasbases %in% c(0,1)))) {
            stop("\"useasbases\" contains non-binary elements")
        }
        #check do not pass in K
        if(!is.null(K) | !is.null(K.svd)) {
            warning("\"useasbases\" argument only used in the construction of the kernel matrix \"K\" and should not be specified when \"K\" or \"K.svd\" is already user-supplied.\n", immediate. = TRUE)
        }
    }
    
    #Setting defaults - useasbases: If we don't specify which observations to use as bases,
    # use all as default unless K is very large, then use sample set.
    if (is.null(useasbases) & N <= 4000 & !linkernel) {
            useasbases = rep(1,N)
    } else if(is.null(useasbases) & !linkernel) {
        if(is.null(K) & is.null(K.svd) ) {
            warning("Dimensions of K greater than 4000, using sampled as default bases\n",
                    immediate. = TRUE)
        }
        useasbases = as.numeric(observed==1)
    }
    #for a linear kernel, the bases are gonna be defined by the cols of the data
    #we will use all of them 
    if (is.null(useasbases) & linkernel) {
        #this does not get used in makeK just for the rest of the
        useasbases = rep(1,ncol(allx))
    } 
 
    ###### Setting defaults: Population  weights #####
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
            stop("\"population.w\" must be positive")
        }
        #check population weights sum to num of treated/population units
        if(round(sum(population.w)) != sum(target)) {
            #allow user to pass in weights that sum to one and transform them here
            if(round(sum(population.w)) == 1) {
                population.w = population.w/mean(population.w)
            } else { #in this case they don't sum to N_t or 1 so ng
                stop("\"population.w\" must sum to either 1 or the number of treated/population units")
            }
        }
        #adding uniform weights for control/sampled units
        w.pop <- rep(1, N)
        w.pop[target==1] = population.w
        
    }

    ##### Setting Defaults: b (maxvar b) #####
    #adding default b within the kbal function rather than in makeK
    #changing default for not cat data to be 2*ncol to match kbal
    if(!is.null(b) && length(b) != 1) {
        stop("\"b \" must be a scalar.")
    }
    maxvar_K_out = NULL
    onehot = NULL
    #cont data only
    if(!cat_data & !mixed_data) {
        
        if(!is.matrix(allx)) {
            allx = as.matrix(allx)
        }
        #check not cat data: 
        if(((is.null(dim(apply(allx, 2, unique))) && sum(lapply(apply(allx, 2, unique), length) <= 10) != 0) | 
           sum(dim(apply(allx, 2, unique))[1] <= 10) != 0)
           && (is.null(K) && is.null(K.svd))) {
            warning("One or more columns of \"allx\" contain less than 10 unique values, but \"cat_data\" and \"mixed_data\" are set to FALSE. Are you sure \"allx\" contains only continuous data?", immediate. = TRUE)
        }
        #if no K and NOT linkern supplied find bmaxvar b
        if(!linkernel & (is.null(K.svd) & is.null(K) & is.null(b))) {
            if(N > 10000) {
                warning("NB: with continuous data and high dimensional \"K\", internal search for b value which produces maximal variance may be time consuming. Consider atlernative choice of b=2*ncol(allx)", 
                        immediate. = TRUE)
            }
            if(printprogress) {
                cat("Searching for b value which maximizes the variance in K: ")
            }
            if(scale_data) {
                allx = scale(allx)
            }    
            res = b_maxvarK(data = allx, 
                            cat_data = cat_data,
                            useasbases = useasbases, 
                            maxsearch_b = maxsearch_b)
            b = res$b_maxvar
            maxvar_K_out = res$var_K
                
            if(printprogress) {
                cat(round(b, 3) ,"selected \n")
            }
        } else if(!is.null(K.svd) | !is.null(K) & is.null(b)) {
            #for later internal checks, not used ofc bc K or svdK is passed in
            b = 2*ncol(allx)
        }
        
        if(meanfirst == TRUE & !is.null(mf_columns)) {
            #note these will be scaled if allx is also scaled (happens above)
            #colnames conversion for mf_columns
            if((class(mf_columns) == "character" & sum(mf_columns %in% colnames(allx)) != length(mf_columns)) |
               (class(mf_columns) == "numeric" &  sum(mf_columns %in% c(1:ncol(allx))) != length(mf_columns))  ) {
                stop("One or more \"mf_columns\" elements does not match the column names in \"allx\" or exceeds the number of columns in \"allx\" ")
            }
            if(class(mf_columns) == "character") { #switch to numeric for ease if input is colnames
                mf_columns = which(colnames(allx) %in% mf_columns)
            }
            allx_mf = allx[, mf_columns]
            
        } else if(!is.null(mf_columns)) {
            warning(" \"mf_columns\" are specified when \"meanfirst\" is FALSE. ignoring input.", immediate. = T)
        }
        
    } else if(cat_data) { #cat_data = TRUE, onehot encode data and find maxvar b
        if(!is.null(cont_scale)) {
            warning("\"cont_scale\" only used with mixed data. Ignoring.\n",
                    immediate. = TRUE)
        }
        
        #mf cols:
        if(meanfirst == TRUE & !is.null(mf_columns)) {
            
            if((class(mf_columns) == "character" & sum(mf_columns %in% colnames(allx)) != length(mf_columns)) |
               (class(mf_columns) == "numeric" &  sum(mf_columns %in% c(1:ncol(allx))) != length(mf_columns))  ) {
                stop("One or more \"mf_columns\" elements does not match the column names in \"allx\" or exceeds the number of columns in \"allx\" ")
            }
            #colnames conversion for mf_columns
            if(class(mf_columns) == "character") { #switch to numeric for ease if input is colnames
                mf_columns = which(colnames(allx) %in% mf_columns)
            }
            allx_mf = one_hot(allx[, mf_columns])
            
        } else if(!is.null(mf_columns)) {
            warning(" \"mf_columns\" are specified when \"meanfirst\" is FALSE. ignoring input.", immediate. = T)
        }
        
        if(!is.null(K.svd) | !is.null(K)) {
            #if meanfirst is true we'll need this
            onehot = one_hot(allx)
            
            #for later internal checks of specified b + passed in K
            if(is.null(b)){ b = 2*ncol(allx) } 
        } else {
            if((is.null(dim(apply(allx, 2, unique))) && sum(lapply(apply(allx, 2, unique), length) == 1) != 0) | sum(dim(apply(allx, 2, unique))[1] == 1) != 0) {
                stop("One or more column in \"allx\" has zero variance")
            }
            if(!is.null(dim(apply(allx, 2, unique))) && 
               dim(apply(allx, 2, unique))[1] == N) {
                stop("\"cat_data\"=TRUE, but all columns of \"allx\" contain n=nrow(allx) unique values indicating continuous data.")
            } 
            if((is.null(dim(apply(allx, 2, unique))) && sum(lapply(apply(allx, 2, unique), length) == 1) > 10) | sum(dim(apply(allx, 2, unique))[1] > 10) != 0) {
                warning("One or more column in \"allx\" has more than 10 unique values while \"cat_data\"= TRUE. Ensure that all variables are categorical.",
                        immediate. = TRUE)
            }
            #for consistency with mixed data, undoing double counts direcly by * by sqrt(0.5)
            allx = sqrt(0.5)*one_hot(allx)
            #still passout the regular one hot data 
            onehot = allx/sqrt(0.5)
        
            #checks
            if(scale_data) {
                warning("Note that \"scale_data\" should be FALSE when using categorical data. Ignoring. \n", 
                        immediate. = TRUE)
            }
            if(linkernel) {
                warning("\"linkernel\" should be FALSE when \"cat_data\" is TRUE. Proceeding with gaussian kernel.\n", 
                        immediate. = TRUE)
                linkernel = FALSE
                #fix useasbases whose default was messed up by this above:
                if (N <= 4000) {
                    useasbases = rep(1,N)
                } else {
                    warning("Dimensions of K greater than 4000, using sampled as default bases\n",
                            immediate. = TRUE)
                    useasbases = as.numeric(observed==1)
                }
            }
            #go get best b
            if(is.null(b)) {
                if(printprogress) {
                    cat("Searching for b value which maximizes the variance in K: ")
                }
                res = b_maxvarK(data = allx, 
                                cat_data = cat_data,
                                useasbases = useasbases, 
                                maxsearch_b = maxsearch_b)
                b = res$b_maxvar
                maxvar_K_out = res$var_K
                if(printprogress) {
                    cat(round(b, 3) ,"selected \n")
                }
            }
        }
    } else { #mixed data: one hot encoded cat data and adjust for double counting, maxvarK
       
        if((!is.null(K.svd) | !is.null(K)) & !meanfirst) {
            #we only end up here if they pass in a K and dont want mf so we dont need to do all the scaling
            #not relevant now w aMF and specified cols
            # warning("\"mixed_data\" TRUE argument only used in the construction of the kernel matrix \"K\" and is not used when \"K\" or \"K.svd\" is already user-supplied.\n", immediate. = TRUE)
            #don't use this it's internal for a check for later if user passes in a b + k.svd
            if(is.null(b)){ b = 2*ncol(allx) } 
        } else {
            if(is.null(cat_columns)) {
                stop("\"cat_columns\" argument must be specified when \"mixed_data\" is TRUE.")
            } else if(class(cat_columns) == "character" & sum(cat_columns %in% colnames(allx)) != length(cat_columns)) {
                stop("One or more \"cat_columns\" elements does not match the column names in \"allx\".")
            } else if(class(cat_columns) == "character") { #switch to numeric for ease if input is colnames
                cat_columns = which(colnames(allx) %in% cat_columns)
            }
            if((is.null(dim(apply(allx, 2, unique))) && sum(lapply(apply(allx, 2, unique), length) == 1) != 0) | sum(dim(apply(allx, 2, unique))[1] == 1) != 0) {
                stop("One or more column in \"allx\" has zero variance")
            }
            if(!is.null(dim(apply(allx[,cat_columns, drop= F], 2, unique))) && 
               dim(apply(allx[,cat_columns, drop= F], 2, unique))[1] == N) {
                stop("\"mixed_data\"=TRUE, but one or more categorical columns of \"allx\" specified by \"cat_columns\" contain n=nrow(allx) unique values indicating continuous data.")
            } 
            if((is.null(dim(apply(allx[,cat_columns, drop= F], 2, unique))) && sum(lapply(apply(allx[,cat_columns, drop= F], 2, unique), length) >10 ) != 0) ) {
                warning("\"mixed_data\"=TRUE, but one or more column in \"allx\" designated as categorcial by \"cat_columns\" has more than 10 unique values. Ensure that all variables are categorical.", immediate. = TRUE)
            }
            #factor of sqrt(0.5) to adjust for double counting s.t. cont and cat contribute equally in kernel dist
            allx_cat = sqrt(0.5)*one_hot(allx[,cat_columns, drop= F])
            
            if((is.null(dim(apply(allx[, -cat_columns, drop = F], 2, unique))) && sum(lapply(apply(allx[, -cat_columns, drop = F], 2, unique), length) <= 10) != 0) | sum(dim(apply(allx[, -cat_columns, drop = F], 2, unique))[1] <= 10) != 0) {
                warning(" \"mixed_data\"= TRUE, but one or more columns of \"allx\" designated as continuous by omission from \"cat_columns\" contain less than 10 unique values. Are you sure  all categorical variables are specified in \"cat_columns\"?", immediate. = TRUE)
            }
            
            if(is.null(cont_scale)) {
                warning("When combining continuous and categorical data, scaling choices for continuous variables imply different relative variable importance in the kernel distance and affect the performance of kbal. A default scaling of sd=1 for all continuous variables will be used, but users are encouraged to carefully think about the most appropriate scaling.\n", immediate. = T) 
                allx_cont <- scale(allx[, -cat_columns, drop = F])
            } else {
                if(!is.numeric(cont_scale) | !(length(cont_scale) == 1 | length(cont_scale) == ncol(allx[, -cat_columns, drop = F])) ) {
                    stop("\"cont_scale\" must be a numeric vector of either length 1 or length equal to the number of continuous columns in \"allx\".")
                }
                #user specified scaling
                allx_cont <- t(t(allx[, -cat_columns, drop = F])/(apply(allx[, -cat_columns, drop = F], 2, sd)*(1/cont_scale)))
            }
            
            #get mf cols before combining all these columns together
            #mf cols
            if(meanfirst == T & !is.null(mf_columns)) {
                if(class(cat_columns) != class(mf_columns)) {
                    stop("please ensure \"mf_columns\" and \"cat_columns\" are of the same type, either character or numeric")
                }
                #this inherits the scaling decisions about all x
                #colnames conversion for mf_columns
                if((class(mf_columns) == "character" & sum(mf_columns %in% colnames(allx)) != length(mf_columns)) |
                   (class(mf_columns) == "numeric" &  sum(mf_columns %in% c(1:ncol(allx))) != length(mf_columns))  ) {
                    stop("One or more \"mf_columns\" elements does not match the column names in \"allx\" or exceeds the number of columns in \"allx\" ")
                } 
                
                if(sum(mf_columns %in% cat_columns) > 0) {
                    allx_mf_cat = one_hot(allx[, mf_columns[which(mf_columns %in% cat_columns)], drop = F ]) 
                    allx_mf_cont = allx_cont[, mf_columns[-which(mf_columns %in% cat_columns)], drop = F]
                    #scaling? is just inhereited
                    allx_mf = cbind(allx_mf_cont, allx_mf_cat)
                } else {
                    allx_mf = allx_cont[, mf_columns, drop = F]
                }
            } else if(!is.null(mf_columns)) {
                warning(" \"mf_columns\" are specified when \"meanfirst\" is FALSE. ignoring input.", immediate. = T)
            }
            #now we can combine safely
            allx <- cbind(allx_cat, allx_cont)
            onehot = cbind(allx_cat/sqrt(0.5), allx_cont)
            
            #checks
            if(scale_data) {
                warning("Note that when \"mixed_data\" is TRUE, scaling is only performed on the continuous data in accordance with \"cont_scale\" and \"scale_data\"=TRUE is not used.\n", 
                        immediate. = TRUE)
            }
            if(linkernel) {
                warning("\"linkernel\" should be FALSE when \"mixed_data\" is TRUE. Proceeding with gaussian kernel.\n", 
                        immediate. = TRUE)
                #otherwise we get mc issues bc of the onehot encoding
                linkernel = FALSE
                #fix useasbases whose default was messed up by this above:
                if (N <= 4000) {
                    useasbases = rep(1,N)
                } else {
                    warning("Dimensions of K greater than 4000, using sampled as default bases\n",
                            immediate. = TRUE)
                    useasbases = as.numeric(observed==1)
                }
            }
            
            #go get best b
            if(is.null(b)) {
                if(N > 10000) {
                    warning("NB: with continuous data and high dimensional \"K\", internal search for b value which produces maximal variance may be time consuming. Consider atlernative choice of b=2*ncol(allx)", 
                            immediate. = TRUE)
                }
                if(printprogress) {
                    cat("Searching for b value which maximizes the variance in K: ")
                }
                res = b_maxvarK(data = allx, 
                                cat_data = cat_data,
                                useasbases = useasbases)
                b = res$b_maxvar
                maxvar_K_out = res$var_K
                if(printprogress) {
                    cat(round(b, 3) ,"selected \n")
                }
            }
        }
    }
    
    
    ###### Setting defaults: minnumdims, maxnumdims #####
    #8. now checking maxnumdims: the most dims you can use is the number of bases
    if (!is.null(maxnumdims) && maxnumdims>sum(useasbases)) {
        warning("Cannot allow dimensions of K to be greater than the number of bases. Reducing \"maxnumdims\". \n", immediate. = TRUE)
        maxnumdims=sum(useasbases)
    }#make sure don't send in a neg
    if (!is.null(maxnumdims) && maxnumdims<=0) {
        stop("\"maxnumdims\" must be greater than zero")
    } #when linkernel ==TRUE ensure maxnumdims not greater than cols of X
    if(linkernel == TRUE && cat_data == FALSE && !is.null(maxnumdims) && maxnumdims > ncol(allx)) {
        warning("When using a linear kernel, cannot allow dimensions of K to be greater than the number of columns in \"allx\". Reducing \"maxnumdims\" to the number of coumns in \"allx\ \n", 
                immediate. = TRUE )
        maxnumdims = ncol(allx)
    }
    if(linkernel == TRUE && !is.null(minnumdims) && minnumdims > ncol(allx)) {
        warning("When using a linear kernel, cannot allow dimensions of K to be greater than the number of columns in \"allx\". Reducing \"minnumdims\" to 1.\n", immediate. = TRUE )
        minnumdims = 1
    }
    
    #Setting defaults: minnumdims, maxnumdims
    if (is.null(minnumdims)){minnumdims=1}
    
    if (is.null(maxnumdims)){
        if(linkernel == FALSE) {
            maxnumdims= min(500, sum(useasbases))  
            if(!is.null(K)) {maxnumdims = ncol(K)}
            if(!is.null(K.svd)) {maxnumdims = ncol(K.svd$u)}
        } else { maxnumdims = min(ncol(allx), 500) } 
        trunc_svd_dims = round(.8*sum(useasbases))
    } else{#2021: if pass max, want to still get svd out more dims so that bb is correct
        trunc_svd_dims = round(max(.8*sum(useasbases), maxnumdims))
        #in case that's bigger than the columns we have is checked below afer we have have K
    }
    #catch for if user passes in K with more rows than cols
    if(!is.null(K) && trunc_svd_dims >= ncol(K)) {
        warning("The number of columns of user-supplied \"K\" is less than 80% of the number of rows, indicating a full svd would not be overly time consuming. Conducting full svd.\n", 
                immediate. = TRUE)
        fullSVD = TRUE
    }

    #9. now checking numdims if passed in
    if(!is.null(numdims) && numdims>maxnumdims) { #check not over max (mainly for truc SVD)
        warning("\"numdims\" cannot exceed \"maxnumdims\". Reducing to maximum allowed.\n",
                immediate. = TRUE)
        numdims= maxnumdims
    } else if(!is.null(numdims) && numdims <= 0) { #check not leq zero
        stop("Specified \"numdims\" must be greater than zero")
    }
    
    #10. incrementby not null and geq 1
    if(is.null(incrementby) || incrementby < 1) {
        warning(" \"incrementby\" must be greater than or equal to 1. Setting \"incrementby\" to be 1.\n", immediate. = TRUE)
        incrementby = 1
    }
   
#####end of big error catch series and data setup

#if pass maxnumdims = N then they want the full svd, so change this for them to avoid all those if statement checks below just to do the same thing
    if(!linkernel && maxnumdims == nrow(allx)) {
        fullSVD = TRUE
    } else if(linkernel & maxnumdims == ncol(allx)) { #for linear kernel this is maxnumdims = ncol
        fullSVD = TRUE
    } 

    ############ Direct CONSTRAINT #############
   # if user passes in constraint to append, ensure it's scaled and dn have mc issues
    if(!is.null(constraint)) {
        if(scale_constraint) {
            constraint <- scale(constraint)
        }
        qr_constr = qr(constraint)
        multicollin_constr = FALSE
        if(qr_constr$rank < ncol(constraint)) {
            MC_out <- drop_multicollin(constraint)
            warning("\"constraint\" contains collinear columns: ", MC_out$dropped_cols,
                    " which will be dropped\n",
                    immediate. = TRUE)
            constraint <- MC_out$allx_noMC
        }
    }
    
    ################## MEAN FIRST #################
    meanfirst_dims = NULL
    if(meanfirst == TRUE) {
        if(!is.null(constraint)) {
            warning("\"constraint\" argument is not used when \"meanfirst\" is TRUE.\n", immediate. = TRUE)
        }
        #note that b and useasbases are irrelevant here since we're using a linear kernel
        #if the user did not specify the mf cols pass all along
        if(is.null(mf_columns)) {
            #for cat data, we just correct the 0.5 factor implicit in allx, which is not in onehot
            if(cat_data) {
                allx_mf = onehot
            #mixed data 
            } else if(mixed_data) {
                #we want to pass through the scaling decisions made for the cont data and also undo the sqrt(.5) factor on onehot
                allx_mf = cbind(allx_cont, onehot)
            #for continuous data we want to pass through the same scaling choices which have arleady been performed above
            } else {
                allx_mf = allx
            }
        }
        kbalout.mean = suppressWarnings(kbal(allx=allx_mf, 
                           treatment=treatment,
                           sampled = sampled,
                           scale_data = TRUE, 
                           #we don't want to do drop mc cols with cat data
                           #with mf routine we take the svd anyway so even w other data
                           #this issue will resolve itself
                           drop_MC = FALSE, 
                           sampledinpop = sampledinpop,
                           useasbases = useasbases,
                           meanfirst = FALSE,
                           ebal.tol = ebal.tol,
                           ebal.maxit = ebal.maxit,
                           ebal.convergence = TRUE, 
                           linkernel = TRUE,
                           printprogress = FALSE))
        
        constraint_svd_keep = kbalout.mean$svdK$u[, 1:kbalout.mean$numdims,drop = F]
        if(printprogress) {
            cat("Selected", kbalout.mean$numdims,
                "dimensions of \"allx\" to use as mean balance constraints. \n")
        }
        meanfirst_dims = kbalout.mean$numdims
        #for now this allows the manual constraint method and this svd way
        constraint = constraint_svd_keep
    }
    
    
    ########### BUILDING K ################
    #pass out b and bases used if kernel built internally and gaussian  
    b_out = ifelse((linkernel | !is.null(K.svd) | !is.null(K)), NA, b)
    if(is.na(b_out)){ b_out = NULL}
    if(linkernel | !is.null(K.svd) | !is.null(K)) {
        useasbases_out = NULL
    } else{ useasbases_out = useasbases}
    
    
    #Setting Up K: Make the kernel or take in user K or user K.svd and check those work with dims
   #first check didn't pass both
    #CASE 1: if user pases K without svd, always conduct SVD upto maxnumdims or numdims
    if(!is.null(K) & is.null(K.svd)) { 
        #error checks:
        #check maxnumdims
        if(maxnumdims > ncol(K) ) {
          warning("\"maxnumdims\" cannot exceed number of columns of \"K\". Reducing to maximum allowed.\n", immediate. = TRUE)
          maxnumdims = ncol(K)
        }
        #check numdims
        if(!is.null(numdims) && numdims > ncol(K) ) {
          warning("\"numdims\" cannot exceed number of columns of \"K\". Reducing to maximum allowed.\n", immediate. = TRUE)
          numdims = ncol(K)
        }
        #check minnumdims 
        if(minnumdims > maxnumdims) {
          warning("\"minnumdims\" cannot exceed \"maxnumdims\". Reducing \"minnumdims\" to 1.\n", immediate. = TRUE)
          minnumdims = 1
         }
        #warning that linkernel meaningless
        if(!is.null(linkernel) && linkernel) {
          warning("\"linkernel\" argument only used in the construction of the kernel matrix \"K\" and is not used when \"K\" or \"K.svd\" is already user-supplied.\n", immediate. = TRUE)
        }
        if(!is.null(b) && b != 2*ncol(allx) ) {
          warning("\"b\" argument only used in the construction of the kernel matrix \"K\" and is not used when \"K\" or \"K.svd\" is already user-supplied.\n", 
                  immediate. = TRUE)
        }
        #provided we pass all those checks get svd with RSpectra
        if(printprogress == TRUE) {cat("Using user-supplied K \n")}
        #if user does not ask for fullsvd, and does not give numdims, get svd upto maxnumdims
        if(!fullSVD) { 
            if(printprogress) {
              cat("Running trucated SVD on kernel matrix up to",trunc_svd_dims, "dimensions \n")
            }
            #for a symmetric K just do eigs_sym as is:
            if(nrow(K) == ncol(K)) {
              rspec.out= suppressWarnings(RSpectra::eigs_sym(K, trunc_svd_dims))
              #add negative evals catch via bound at e-12 = 0
              #for now just set values to be zero
              rspec.out$values[abs(rspec.out$values) <= 1e-12 ] = 0
              if(sum(sign(rspec.out$values) == -1) != 0) {
                  stop("Trucated SVD produced negative eigenvalues, please rerun using \"fullSVD=TRUE\" ")
              }
              svd.out = list(u = rspec.out$vectors, d = rspec.out$values,
                             v= rspec.out$vectors)
              #tr(K) = sum eigenvalues
              var_explained <- round(sum(svd.out$d)/nrow(K),6)
              if(printprogress) {
                  cat("Truncated SVD with", 
                      trunc_svd_dims,
                      "first singular values accounts for", 
                      round(100*var_explained, 2),
                      "% of the variance of \"K\" \n")
              }
              if(var_explained < .999) {
                  warning("Truncated SVD with only ", trunc_svd_dims,
                          " first singular values only accounts for ", round(100*var_explained, 2) ,
                          " of the variance of \"K\". The biasbound optimization may not perform as expected. You many want to increase \"maxnumdims\" to capture more of the total varince of \"K\".\n", immediate. = TRUE)
              }
          } else { #use svds, suppressing warnings that it prints if uses full size svd
             
              svd.out= RSpectra::svds(K, trunc_svd_dims)
              if(sum(sign(svd.out$d) == -1) != 0) {
                  stop("Trucated SVD produced negative eigenvalues, please rerun using \"fullSVD=TRUE\" ")
              }
              #don't know the total sum of eigenvalues, but sum up to calculated and assume all remaining are = to last
              max_rank = min(nrow(K), ncol(K))
              worst_remaining_var <- sum(svd.out$d) + (max_rank-trunc_svd_dims)*svd.out$d[length(svd.out$d)]
              var_explained <- round(sum(svd.out$d)/worst_remaining_var,6)
              
              if(printprogress) {
                  cat("When bases are chosen such that \"K\" is nonsymmetric, the proportion of total variance in \"K\" accounted for by the truncated SVD with only", 
                      trunc_svd_dims,
                      "first singular values is lower bounded (worst-case) to explain",
                      round(100*var_explained, 2), "% of the variance of \"K\" \n")
              }
              if(var_explained < .999) {
                  warning("Truncated SVD with only ", trunc_svd_dims,
                          " first singular values accounts for, worst-case, approximately only ",
                          round(100*var_explained, 2),
                          "% of the variance of \"K\". The biasbound optimization may not perform as expected. You many want to increase \"maxnumdims\" to capture more of the total varince of \"K\".\n", 
                          immediate. = TRUE)
              }
    
          }
            
            U=svd.out$u
            #if user asked for full SVD, get it
        } else { 
          if(printprogress) {cat("Running full SVD on kernel matrix \n")}
          svd.out = svd(K)
          U = svd.out$u
          var_explained = NULL
        }
    #CASE 2: if user pases in K.svd (with or without K) dn conduct SVD   
    } else if(!is.null(K.svd)) {
        #NB: we require K.svd to have $u and $d $v just as a real svd would
        #error catches
        #check maxnumdims
        if(maxnumdims > ncol(K.svd$u) ) {
            warning("\"maxnumdims\" cannot exceed number of columns of \"K.svd\". Reducing to maximum allowed.\n", immediate. = TRUE)
            maxnumdims = ncol(K.svd$u)
        }
        #check numdims
        if(!is.null(numdims) && numdims > ncol(K.svd$u) ) {
            warning("\"numdims\" cannot exceed number of columns of \"K.svd\". Reducing to maximum allowed.\n", immediate. = TRUE)
            numdims = ncol(K)
         }
        #check minnumdims 
        if(minnumdims > maxnumdims) {
            warning("\"minnumdims\" cannot exceed \"maxnumdims\". Reducing \"minnumdims\" to 1.\n", immediate. = TRUE)
            minnumdims = 1
        }
        if(!is.null(linkernel) && linkernel) { #only if linkernel = TRUE
            warning("\"linkernel\" argument only used in the construction of the kernel matrix \"K\" and is not used when \"K\" or \"K.svd\" is already user-supplied.\n", immediate. = TRUE)
        }
        if(b != 2*ncol(allx)) {
            warning("\"b\" argument only used in the construction of the kernel matrix \"K\" and is not used when \"K\" or \"K.svd\" is already user-supplied.\n", immediate. = TRUE)
        }
        if(sum(c(length(ls(K.svd)) >= 2, c("d", "u") %in% ls(K.svd))) != 3 ) {
            stop("\"K.svd\" must be a list object containing \"u\" the left singular vectors and \"d\" the singular values.")
        } else if(ncol(K.svd$u) != length(K.svd$d)) {
            stop("\"K.svd\" must be a list object containing \"u\" the left singular vectors and \"d\" the singular values. Dimensions of \"u\" do not match dimensions of \"d\".")
        }
        svd.out = K.svd
        U = K.svd$u
        var_explained = NULL
        if(!is.null(K)) {
            warning("\"K\" only used for calculating L1 distance. All balancing and weight construction only relies on \"K.svd\" input\n")
        } else {
            #reconstruct K for the L1 distance
            d_diag <- matrix(0, nrow = ncol(K.svd$u), ncol = ncol(K.svd$v))
            diag(d_diag) <- K.svd$d
            K = K.svd$u %*% d_diag  %*% t(K.svd$v)
        }
    #CASE 3: if user does not specify either K or K.svd, build K and get svd of K
    } else { 
        if(printprogress == TRUE) {cat("Building kernel matrix\n")}
        if(linkernel == FALSE) {
            K = makeK(allx = allx, useasbases = useasbases, b=b, 
                      scale = scale_data)
        } else {
            K = makeK(allx = allx,
                      scale = scale_data,
                      useasbases = useasbases,  #unnecc/bases not used for lin kernel
                      linkernel = TRUE)
        }
        #if user does not ask for full svd, and does not pass in numdims, get svd upto maxnumdim
        if(!fullSVD) {
            if(printprogress) {cat("Running truncated SVD on kernel matrix up to",
                                 trunc_svd_dims, "dimensions \n")}
            #for a symmetric K just do eigs_sym as is:
            if(nrow(K) == ncol(K)) {
                rspec.out= suppressWarnings(RSpectra::eigs_sym(K, trunc_svd_dims))
                rspec.out$values[ abs(rspec.out$values) <= 1e-12 ] = 0
                if(sum(sign(rspec.out$values) == -1) != 0) {
                    stop("Trucated SVD produced negative eigenvalues, please rerun using \"fullSVD=TRUE\" ")
                }
                svd.out = list(u = rspec.out$vectors, d = rspec.out$values,
                               v= rspec.out$vectors)
                var_explained = round(sum(svd.out$d)/nrow(K),6)
                if(printprogress) {cat("Truncated SVD with", trunc_svd_dims,
                                     "first singular values accounts for", var_explained,
                                     "of the variance of \"K\" \n")}
              
                if(var_explained < .999) {
                    warning("Truncated SVD with ", trunc_svd_dims,
                          " first singular values only accounts for ", var_explained,
                          " of the variance of \"K\". The biasbound optimization may not perform as expected. You many want to increase \"maxnumdims\" to capture more of the variance of \"K\" \n", immediate. = TRUE)
                }
            } else { #use truncated svd
                svd.out= RSpectra::svds(K, round(trunc_svd_dims))
                if(sum(sign(svd.out$d) == -1) != 0) {
                    stop("Trucated SVD produced negative eigenvalues, please rerun using \"fullSVD=TRUE\" ")
                }
                #don't know the total sum of eigenvalues, but sum up to calculated and assume all remaining are = to last
                max_rank = min(nrow(K), ncol(K))
                worst_remaining_var <- sum(svd.out$d) + (max_rank-trunc_svd_dims)*svd.out$d[length(svd.out$d)]
                var_explained <- round(sum(svd.out$d)/worst_remaining_var,6)

                if(printprogress) {
                    cat("When bases are chosen such that \"K\" is nonsymmetric, the proportion of total variance in \"K\" accounted for by the truncated SVD with only",
                        trunc_svd_dims,
                        "first singular values is lower bounded (worst-case) to explain",
                        round(100*var_explained, 2), "% of the variance of \"K\" \n")
                }
            }
            U=svd.out$u
        #if user askes for full svd, go get it
        } else {
            if(printprogress) {cat("Running full SVD on kernel matrix \n")}
            svd.out = svd(K)
            U = svd.out$u
            var_explained = NULL
        }
    }
    #adjust singular values for linear kernel since we just do svd(X) instead of svd(XX')
    if(linkernel == TRUE) {
        svd.out$d = svd.out$d^2
    }
    
    ####### Adding Constraint to minimization: paste constraint vector to front of U
    if(!is.null(constraint)) {
        #check dims of constraint
        #this will cause problems if pass in one constraint as a vector, it needs to be a 1 column matrix
        if(!class(constraint)[1] %in% c("matrix", "data.frame")) {
            stop("\"constraint\" must be a matrix")
        }
        if(nrow(constraint) != N) { 
            stop("\"constraint\" must have the same number of rows as \"allx\"")
        } 
        minnumdims <- ncol(constraint) + minnumdims
        maxnumdims <- ncol(constraint) + maxnumdims
        #binds constraint to front of separate object U, leaving original U in svd.out
        #this way we run biasbound() on svd.out to get biasbound on svd(K) only
        #but use U for getw to get weights that balance on constraints as well
        #note also that L1 distance uses K and that remains unchanged so it is ok
        #while we have a K matrix (if pass in svd not as straight forward)
        U = cbind(constraint, U) 
        #if numdims given move it up to accomodate the constraint
        if(!is.null(numdims)) {
            numdims = numdims + ncol(constraint)
        }
    }
    
    ################### BASELINE ##############################
    # Get biasbound with no improvement in balance:
    #recall: w.pop is flat weights for sampled, user specified weights for population
    biasbound_orig=biasbound(w = rep(1,N),
                             w.pop = w.pop, 
                             observed=observed, 
                             target = target,
                             svd.out = svd.out, hilbertnorm = 1)
    #NB: not well defined when pass in K.svd not K (added warning above)
    getdist.orig = getdist(target=target, observed = observed,
                           w = rep(1,N), w.pop = w.pop, K=K)
    L1_orig = getdist.orig$L1
    if(printprogress==TRUE) {
        cat("Without balancing, biasbound (norm=1) is",
            round(biasbound_orig,5), "and the L1 discrepancy is", round(L1_orig,3), "\n")
        }

    ############# NUMDIMS GIVEN  ###################
    # If numdims given, just get the weights in one shot:
    if(!is.null(numdims)) {
        U2=U[,1:numdims, drop=FALSE]
        U2.w.pop <- w.pop*U2
        getw.out= getw(target=target, observed=observed, svd.U=U2.w.pop, 
                       ebal.tol = ebal.tol, ebal.maxit = ebal.maxit)
        converged = getw.out$converged
        biasboundnow=biasbound( w = getw.out$w,
                                observed=observed,  target = target,
                                svd.out = svd.out, 
                                w.pop = w.pop, 
                                hilbertnorm = 1)
        if(!converged) {
            numpass = numdims
            if(!is.null(constraint)) {numpass = numdims - ncol(constraint)}
            warning("With user-specified ", 
                    numpass,
                    " dimensions ebalance did not converge within tolerance. Disregarding ebalance convergence and returining weights, biasbound, and L1 distance for requested dimensions.\n")
            }
        if(printprogress == TRUE & is.null(constraint)) {
            cat("With user-specified", numdims,"dimensions, biasbound (norm=1) of ",
                round(biasboundnow,5), " \n")
        } else if(printprogress) {
            numdims = numdims - ncol(constraint)
            cat("With user-specified",numdims,"dimensions of K, biasbound (norm=1) of ",
                round(biasboundnow,5), " \n")
        }
        
        #stuff to set so we can skip the entire if statement below and just printout
        dist.record = biasboundnow
        dist_pass = rbind(numdims, dist.record, converged)
        biasbound_opt = biasboundnow
        dist.orig= biasbound_orig
        L1_optim = getdist(target=target, observed = observed,
                           w = getw.out$w, w.pop = w.pop, K=K)$L1
    }
 
  
    ############# BIASBOUND OPTIMIZATION ###################
    # If numdims not given, we search to minimize biasbound:
    if(is.null(numdims)) {
        thisnumdims=minnumdims
        dist.record=NULL
        convergence.record = NULL
        keepgoing=TRUE
        wayover=FALSE
        mindistsofar=998
        # keep track of how many iterations have increased distance or diverged as stopping criteria
        dist.incr.rounds = 0
        
        while(keepgoing==TRUE) {
            U_try=U[,1:thisnumdims, drop=FALSE]
            U_try.w.pop <- w.pop*U_try
            getw.out= suppressWarnings(getw(target = target,
                                            observed=observed, svd.U = U_try.w.pop, 
                                            ebal.tol = ebal.tol, ebal.maxit = ebal.maxit))
            convergence.record = c(convergence.record, getw.out$converged)
            
            biasboundnow=biasbound(w = getw.out$w,
                                   observed=observed,  target = target,
                                   svd.out = svd.out, 
                                   w.pop = w.pop, 
                                   hilbertnorm = 1)
            if(printprogress == TRUE & is.null(constraint)) {
                cat("With",thisnumdims,"dimensions of K, ebalance convergence is", 
                    getw.out$converged ,"yielding biasbound (norm=1) of",
                    round(biasboundnow,5), " \n")
            } else if(printprogress == TRUE) {
                cat("With",thisnumdims - ncol(constraint),"dimensions of K, ebalance convergence is",
                    getw.out$converged, "yielding biasbound (norm=1) of",
                    round(biasboundnow,5), " \n")
            }
            
            dist.record=c(dist.record,biasboundnow)
            dist.now=biasboundnow # To make more generic, distance could be any measure.
            dist.orig=biasbound_orig
            
            thisnumdims=thisnumdims+incrementby
            
            if(!is.na(dist.now)) {
                if(dist.now<mindistsofar){mindistsofar=dist.now} 

                # check if optimization rounds are increasing the bias bound
                if (dist.now > mindistsofar) {
                  dist.incr.rounds = dist.incr.rounds  + 1
                } else {
                  dist.incr.rounds = 0
                }

                wayover=(dist.now/mindistsofar)>1.25
                # if early stopping is enabled, we break the loop if the distance has been increasing for 20 rounds
                keepgoing=(thisnumdims<=maxnumdims) & (wayover==FALSE) & (early.stopping == FALSE | dist.incr.rounds <= 20) #& (dist.now<dist.orig)
                #& dist.now!="error"
                # (dist.now>mindistsofar)  # XXX this was in there, but needed?
                # Need to work on "keepgoing" for case where ebal fails.
            }
        } # End of while loop for "keepgoing"
    
        if(is.null(constraint)) {
            dimseq=seq(minnumdims,maxnumdims,incrementby)
        } else {
            dimseq=seq(minnumdims-ncol(constraint),maxnumdims,incrementby)
        }
        #building record of search
        dist_pass = rbind(dimseq[1:length(dist.record)], dist.record, convergence.record)
        rownames(dist_pass) <- c("Dims", "BiasBound", "Ebal Convergence")
        

        ########### CHECKING CONVERGENCE #############
        #1) user requires ebal convergence and we did have some convergence
        if(ebal.convergence & sum(convergence.record) != 0) {
            min_converged = min(dist.record[convergence.record], na.rm=TRUE)
            numdims=dimseq[which(dist.record==min_converged)]
            
            # If nothing improved balance, there will be multiple minima.
            # Throw warning, and choose the fewest numdims.
            if(length(numdims)>1) {
                warning("Lack of improvement in balance; choosing fewest dimensions to balance on among those with the same (lack of) improvement. But beware that balance is likely poor.\n",
                            immediate. = TRUE)
                numdims=min(numdims)
            }
            # Finally, we didn't save weights each time, so go back and re-run
            # at optimal  number of dimensions
            if(printprogress == TRUE) {
                cat("Re-running at optimal choice of numdims,", numdims, "\n")
            }
            if(is.null(constraint)) {
                U_final.w.pop <- w.pop*U[,1:numdims, drop = FALSE]
            } else {
                U_final.w.pop <- w.pop*U[,1:(numdims+ncol(constraint)), drop = FALSE]
            }
            getw.out = getw(target= target, observed=observed, 
                            svd.U=U_final.w.pop, 
                            ebal.tol = ebal.tol, ebal.maxit = ebal.maxit)
            biasbound_opt= biasbound(w = getw.out$w, observed=observed, target = target, 
                                     svd.out = svd.out, 
                                     w.pop = w.pop,
                                     hilbertnorm = 1)
            #NB: not well defined iF K.svd passed in
            L1_optim = getdist(target=target, observed = observed,
                               w = getw.out$w, w.pop = w.pop, K=K)$L1
        #2) user asked for convergence, but NO CONVERGED DIMS  
        } else if(ebal.convergence) {
            #if we sent in constraints let's return weights just on these 
            if(!is.null(constraint)) {
                U_constraint=U[,1:(minnumdims-1), drop=FALSE]
                U_c.w.pop <- w.pop*U_constraint
                getw.out= getw(target = target, observed=observed, svd.U = U_c.w.pop, 
                               ebal.tol = ebal.tol, ebal.maxit = ebal.maxit)
                convergence.record = getw.out$converged
                warning("Ebalance did not converge within tolerance for any ",
                        dist_pass[1,1],"-",
                        dist_pass[1,ncol(dist_pass)],
                        " searched dimensions of K.\nNo optimal numdims to return. Returning biasbound, L1 distance, and weights from balance on constraints only with ebalance convergence,",
                        convergence.record)
                
                biasbound_opt=biasbound(w = getw.out$w,
                                        observed=observed, 
                                        target = target,
                                        svd.out = svd.out, 
                                        w.pop = w.pop, 
                                        hilbertnorm = 1)
                L1_optim = getdist(target=target, observed = observed,
                                   w = getw.out$w, w.pop = w.pop, K=K)$L1
                numdims = NULL
            } else { #no constraint and no convergence
                warning("Ebalance did not converge within tolerance for any ",dist_pass[1,1],"-",
                        dist_pass[1,ncol(dist_pass)],
                        " searched dimensions of K.\nDisregarding ebalance convergence and returning biasbound, L1 distance, and weights that yeild the minimum biasbound.")
                #disregard convergence and pick minnumdims
                numdims=dimseq[which(dist.record==min(dist.record,na.rm=TRUE))]
                U_final.w.pop <- w.pop*U[,1:numdims, drop = FALSE]
                getw.out = getw(target= target, observed=observed,
                                svd.U=U_final.w.pop, 
                                ebal.tol = ebal.tol, ebal.maxit = ebal.maxit)
                biasbound_opt= biasbound(w = getw.out$w, observed=observed, target = target, 
                                         svd.out = svd.out, 
                                         w.pop = w.pop,
                                         hilbertnorm = 1)
                #NB: not well defined iF K.svd passed in
                L1_optim = getdist(target=target, observed = observed,
                                   w = getw.out$w, w.pop = w.pop, K=K)$L1
            }   
        #3) User did not ask for convergence
        } else { #we don't for convergence at all, we just find min biasbound
            if(is.null(constraint)) {
                dimseq=seq(minnumdims,maxnumdims,incrementby)
                numdims=dimseq[which(dist.record==min(dist.record,na.rm=TRUE))]
                if (length(numdims)>1) {
                    warning("Lack of improvement in balance; choosing fewest dimensions to balance on among those with the same (lack of) improvement. But beware that balance is likely poor.\n",
                            immediate. = TRUE)
                    numdims=min(numdims)
                }
                U_final.w.pop <- w.pop*U[,1:numdims, drop = FALSE]
            } else {
                dimseq=seq(minnumdims-ncol(constraint),maxnumdims,incrementby)
                numdims=dimseq[which(dist.record==min(dist.record,na.rm=TRUE))]
                if (length(numdims)>1) {
                    warning("Lack of improvement in balance; choosing fewest dimensions to balance on among those with the same (lack of) improvement. But beware that balance is likely poor.\n",
                            immediate. = TRUE)
                    numdims=min(numdims)
                }
                U_final.w.pop <- w.pop*U[,1:(numdims + ncol(constraint)), drop = FALSE]
            }
            
            if(printprogress == TRUE) {
                cat("Disregarding ebalance convergence and re-running at optimal choice of numdims,",
                    numdims, "\n")
            }
            
            getw.out = getw(target= target, observed=observed, 
                            svd.U=U_final.w.pop, 
                            ebal.tol = ebal.tol, ebal.maxit = ebal.maxit)
            biasbound_opt= biasbound(w = getw.out$w, observed=observed, target = target, 
                                     svd.out = svd.out, 
                                     w.pop = w.pop,
                                     hilbertnorm = 1)
            #NB: not well defined iF K.svd passed in
            L1_optim = getdist(target=target, observed = observed,
                               w = getw.out$w, w.pop = w.pop, K=K)$L1
    
        }
    }
    if(meanfirst) {
        cat("Used", meanfirst_dims, "dimensions of \"allx\" for mean balancing, and an additional",
            numdims, "dimensions of \"K\" from kernel balancing.\n") 
    }
    
    ebal_error = getw.out$ebal_error
        
    R=list()
    R$w= getw.out$w
    R$biasbound_opt=biasbound_opt
    R$biasbound_orig=dist.orig
    R$biasbound_ratio= dist.orig/biasbound_opt
    R$dist_record= dist_pass
    R$numdims=numdims
    R$L1_orig = L1_orig
    R$L1_opt = L1_optim
    R$K = K
    R$onehot_data = onehot
    R$linkernel = linkernel
    R$svdK = svd.out
    R$b = b_out
    R$maxvar_K = maxvar_K_out
    R$bases = useasbases_out
    R$truncatedSVD.var = var_explained
    R$dropped_covariates = dropped_cols
    R$meanfirst_dims = meanfirst_dims
    R$appended_constraint_cols = constraint
    R$ebal_error = ebal_error
  return(R)
} # end kbal main function





