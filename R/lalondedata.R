#' Data from National Supported Work program and Panel Study in Income Dynamics
#' @description
#' Dehejia and Wahba (1999) sample of ata from Lalonde (1986). 
#' This data set includes 185 treated units from the National 
#' Supported Work (NSW) program, paired with 2490 control units
#' drawn from the Panel Study of Income Dynamics (PSID-1).
#'
#' The treatment variable of interestis \code{nsw}, which indicates that an individual 
#' was in the job training program. The main outcome of interest is
#' real earnings in 1978 (\code{re78}).  The remaining variables are characteristics
#' of the individuals, to be used as controls.
#'
#' @format A data frame with 2675 rows and 12 columns.
#' \describe{
#'   \item{nsw}{treatment indicator: participation in the National Supported Work program.}
#'   \item{re78}{real earnings in 1978 (outcome)}
#'   \item{u78}{unemployed in 1978; actually an indicator for zero income in 1978.}
#'   \item{age}{age in years}
#'   \item{black}{indicator for identifying as black.}
#'   \item{hisp}{indicator for identifying as hispanic}
#'   \item{married}{indicator for being married}
#'   \item{re74}{real income in 1974}
#'   \item{re75}{real income in 1975}
#'   \item{u74}{unemployment in 1974; actually an indicator for zero income in 1974}
#'   \item{u75}{unemployment in 1975; actually an indicator for zero income in 1975}
#' }
#' @references Dehejia, Rajeev H., and Sadek Wahba. "Causal effects in nonexperimental studies: Reevaluating the evaluation of training programs." Journal of the American statistical Association 94.448 (1999): 1053-1062.
#' 
#'LaLonde, Robert J. "Evaluating the econometric evaluations of training programs with experimental data." The American economic review (1986): 604-620.
"lalonde"
