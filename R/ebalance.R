#' Modified version of ebalance (originally from Jens Hainmueller)
#' @description  This is a custom version of the \code{ebal} (entropy balancing) package by Jens Hainmueller.  Chooses weights on controls to make
#' covariate means equal to those of treated. This version differs from \code{ebal} only in that it handles
#' cases where there is only a single unit, which otherwise causes a problem in the original code.
#' @param Treatment a numeric vector of length equal to the total number of units where treated (population) units take a value of 1 and control (sampled) units take a value of 0.
#' @param X matrix of data where rows are observations and columns are covariates.
#' @param base.weight an optional numeric vector argument of length equal to the total number of control units to specify the base weight of each control unit within entropy balancing. Default is even weights (1) for all controls.
#' @param norm.constant an optional numeric argument; users should leave unspecified in most cases.
#' @param coefs an optional vector argument of length equal to one more than the number of covariates in \code{X}; users should leave unspecified in most cases. 
#' @param max.iterations numeric maximum number of iterations to use when searching for weights
#' @param constraint.tolerance numeric tolerance level.
#' @param print.level a numeric argument to specify the amount of information printed out. 0 is silent, 1 prints convergence status, 2 prints maximum deviance per iteration, 3 prints loss and step length.
#' @importFrom stats var optimize
#' @return \item{target.margins}{Column sums of \code{X} among the treated units.}
#' \item{co.xdata}{Covariate matrix for the controls only built from \code{X} with an additional appended column of ones.}
#' \item{w}{weights found using ebalance. Note that treated units all receive flat weights of 1}
#' \item{maxdiff}{absolute value of the largest component of the gradient in the last iteration.}
#' \item{norm.constant}{norm constant used}
#' \item{constraint.tolerance}{tolerance used to evaluate convergence}
#' \item{max.iterations}{max iterations used}
#' \item{base.weight}{base weights used}
#' \item{print.level}{print level used}
#' \item{converged}{Convergence status. If ebalance failed to find weights within the specified \code{constraint.tolerance} after \code{max.iterations} this is \code{FALSE}. Note that even if ebalance does not converge, the last iteration's weights \code{w} are returned.}
#' @export

ebalance_custom <-
    function(
        Treatment,
        X,
        base.weight = NULL,
        norm.constant  = NULL,
        coefs = NULL ,
        max.iterations = 200,
        constraint.tolerance = 1e-3,
        print.level=0
    ){
        
        # Checks
        if (sum(Treatment  != 1 & Treatment  != 0) > 0) {
            stop("Treatment indicator ('Treatment') must be a logical variable, TRUE (1) or FALSE (0)")
        }
        if (var(Treatment) == 0) {
            stop("Treatment indicator ('Treatment') must contain both treatment and control observations")
        }
        
        Treatment <- as.numeric(Treatment)
        X  <- as.matrix(X)
        
        if (sum(is.na(X))>0){
            stop("X contains missing data")
        }
        
        if (sum(is.na(Treatment))>0){
            stop("Treatment contains missing data")
        }
        
        if (length(Treatment) != nrow(X)) {
            stop("length(Treatment) != nrow(X)")
        }
        
        if (length(max.iterations) != 1 ) {
            stop("length(max.iterations) != 1")
        }
        if (length(constraint.tolerance) != 1 ) {
            stop("length(constraint.tolerance) != 1")
        }
        
        # set up elements
        ntreated  <- sum(Treatment==1)
        ncontrols <- sum(Treatment==0)
        
        if (is.null(base.weight)) {
            base.weight = rep(1, ncontrols)
        }
        if ( length(base.weight) !=  ncontrols) {
            stop("length(base.weight) !=  number of controls  sum(Treatment==0)")
        }
        
        co.x <- X[Treatment==0,]
        co.x <- cbind(rep(1,ncontrols),co.x)
        
        if(qr(co.x)$rank != ncol(co.x)){
            stop("collinearity in covariate matrix for controls (remove collinear covariates)")
        }
        
        
        tr.total <- apply(as.matrix(X[Treatment==1,,drop=FALSE]),2,sum)
        
        if (is.null(norm.constant)) {
            norm.constant <- ntreated
        }
        if (length(norm.constant) != 1) {
            stop("length(norm.constant) != 1")
        }
        
        tr.total <- c(norm.constant,tr.total)
        
        if(is.null(coefs)) {
            coefs = c(log(tr.total[1]/sum(base.weight)),rep(0,(ncol(co.x)-1)))
        }
        
        if(length(coefs) != ncol(co.x)) {
            stop("coefs needs to have same length as number of covariates plus one")
        }
        
        ## run algo
        eb.out <- eb(tr.total=tr.total,
                     co.x=co.x,
                     coefs=coefs,
                     base.weight=base.weight,
                     max.iterations=max.iterations,
                     constraint.tolerance=constraint.tolerance,
                     print.level=print.level
        )
        
        #redundant because eb() run above will print this out if print.level >=1
        # if(eb.out$converged == TRUE & print.level>=1) {
        #     cat("Converged within tolerance \n")
        # }
        
        z <- list(
            target.margins = tr.total,
            co.xdata = co.x,
            w=eb.out$Weights.ebal,
            coefs=eb.out$coefs,
            maxdiff=eb.out$maxdiff,
            norm.constant = norm.constant,
            constraint.tolerance=constraint.tolerance,
            max.iterations=max.iterations,
            base.weight=base.weight,
            print.level=print.level,
            converged=eb.out$converged
        )
        
        class(z) <- "ebalance"
        return(z)
        
    }



eb <- function(
    tr.total=tr.total,
    co.x=co.x,
    coefs=coefs,
    base.weight=base.weight,
    max.iterations=max.iterations,
    constraint.tolerance=constraint.tolerance,
    print.level=print.level
) {
    
    converged <- FALSE
    for(iter in 1:max.iterations) {
        weights.temp <-  c(exp(co.x %*% coefs))
        weights.ebal <- weights.temp *  base.weight
        co.x.agg   <- c(weights.ebal %*% co.x)
        gradient   <- co.x.agg - tr.total
        if(max(abs(gradient))<constraint.tolerance){
            converged <- TRUE
            break
        }
        if(print.level>=2){ cat("Iteration",iter,"maximum deviation is =",format(max(abs(gradient)),digits=4),"\n") }
        hessian = t(co.x) %*% (weights.ebal * co.x)
        Coefs <- coefs
        newton <- solve(hessian,gradient)
        coefs  <- coefs - newton
        loss.new <- line.searcher(Base.weight=base.weight,Co.x=co.x,Tr.total=tr.total,coefs=coefs,Newton=newton,ss=1)
        loss.old <- line.searcher(Base.weight=base.weight,Co.x=co.x,Tr.total=tr.total,coefs=Coefs,Newton=newton,ss=0)
        if(print.level>=3){cat("new loss",loss.new,"old loss=",loss.old,"\n")}
        
        if (is.na(loss.new)== FALSE && is.na(loss.old)==FALSE) {
            if(loss.old <= loss.new){
                ss.out <- suppressWarnings(optimize(line.searcher,
                                   lower=.00001,upper=1,maximum=FALSE,
                                   Base.weight=base.weight,Co.x=co.x,
                                   Tr.total=tr.total,coefs=Coefs,Newton=newton))
                
                if(print.level>=3){cat("LS Step Length is ",ss.out$minimum,"\n")}
                if(print.level>=3){cat("Loss is",ss.out$objective,"\n")}
                coefs = Coefs - ss.out$minimum*solve(hessian,gradient)
            }
            
        }
        
        
    }
    if(print.level>=1 && converged){cat("Converged within tolerance \n")}
    return(
        list(
            maxdiff=max(abs(gradient)),
            coefs=coefs,
            Weights.ebal=weights.ebal,
            converged=converged
        )
    )
}




# function to conduct line search for optimal step length
line.searcher <- function(
    Base.weight,
    Co.x,
    Tr.total,
    coefs,
    Newton,
    ss)
{
    weights.temp <- c(exp(Co.x %*% (coefs - (ss * Newton) )))
    #weights.temp[is.infinite(weights.temp)] <- 100
    weights.temp <- weights.temp * Base.weight
    Co.x.agg     <- c(weights.temp %*% Co.x)
    maxdiff      <- max(abs(Co.x.agg-Tr.total))
    return(maxdiff)
}