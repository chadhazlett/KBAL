##############################
## SE functions for kpop output
##
## calc_fixedweight_var  -- internal helper
## calc_linearized_var   -- exported
## calc_ses_cat_var      -- exported
## kpop_summary          -- exported
## print.kpop_summary    -- exported
##############################


## Internal: fixed-weight variance for a weighted mean.
##
## Treats weights as fixed/known constants. Understates uncertainty when
## weights were estimated (e.g. via kpop), but fast and simple. Used as a
## comparison column in calc_ses_cat_var and kpop_summary.
##
## @param outcome Numeric vector of outcome values.
## @param weights Numeric vector of weights (same length as outcome).
## @return Scalar variance of the weighted mean.
## @keywords internal
calc_fixedweight_var <- function(outcome, weights) {
  y     <- outcome
  w     <- weights / sum(weights)
  n_eff <- 1 / sum(w^2)
  mu_w  <- sum(w * y)
  (n_eff / (n_eff - 1)) * sum(w^2 * (y - mu_w)^2)
}


#' Linearized variance for a kpop-weighted mean
#'
#' @description Uses the kernel dimensions stored in the \code{kpop()} output
#'   to fit a weighted ridge regression of the outcome, then applies a
#'   sandwich-type linearization that accounts for the fact that the weights
#'   were estimated rather than fixed. This is the preferred variance estimator
#'   for kpop-weighted analyses.
#'
#'   Internally calls \code{glmnet::cv.glmnet} with \code{alpha = 0} (ridge).
#'   Set a random seed before calling for reproducibility.
#'
#' @param kpop_object Object returned by \code{\link{kpop}} or
#'   \code{\link{kbal}}.
#' @param outcome Numeric vector of outcomes for the units of interest.
#'   Must already be filtered to non-NA rows; length must equal
#'   \code{length(row_idx)}.
#' @param row_idx Integer vector of row positions within the original sample
#'   corresponding to \code{outcome}. Used to select the matching rows of the
#'   kpop object (weight vector, SVD matrix, constraint columns). Defaults to
#'   \code{1:length(outcome)} (all sample units). Does \emph{not} re-index
#'   \code{outcome}.
#' @param sample_size Retained for API compatibility; not used in the current
#'   formula (effective \emph{n} is derived from the normalized weights).
#' @param use_1se_lambda Logical. If \code{TRUE} (default), uses the
#'   \code{lambda.1se} rule from cross-validation; if \code{FALSE}, uses
#'   \code{lambda.min}.
#'
#' @return Scalar variance of the kpop-weighted mean, or \code{NA_real_} if
#'   the ridge regression cannot be fit (e.g. near-constant outcome after
#'   glmnet's internal standardization).
#'
#' @seealso \code{\link{kpop_summary}}, \code{\link{calc_ses_cat_var}}
#'
#' @importFrom glmnet cv.glmnet
#' @importFrom stats predict
#' @export
calc_linearized_var <- function(kpop_object,
                                outcome,
                                row_idx        = NULL,
                                sample_size    = NULL,
                                use_1se_lambda = TRUE) {

  ## outcome must already be the filtered subset corresponding to row_idx.
  ## row_idx is used only to select matching rows of the kpop object
  ## (weights, SVD matrix, constraint columns); it does NOT re-index outcome.
  n <- length(outcome)
  if (is.null(row_idx)) row_idx <- seq_len(n)

  ## Trivial constant outcome: SE is 0
  if (length(unique(outcome)) < 2) return(0)

  ## Weights for the subset; normalize to sum to 1
  w_raw <- as.numeric(kpop_object$w[row_idx])
  if (any(!is.finite(w_raw)) || any(w_raw < 0)) return(NA_real_)
  if (sum(w_raw) <= 0) return(NA_real_)
  w <- w_raw / sum(w_raw)

  ## Kernel basis (right singular vectors) for selected rows
  Kdims <- kpop_object$numdims
  V     <- as.matrix(kpop_object$svdK$v[row_idx, seq_len(Kdims), drop = FALSE])

  ## Design matrix: unpenalized appended constraints (if any) + penalized kernel dims
  if (!is.null(kpop_object$appended_constraint_cols)) {
    A       <- as.matrix(kpop_object$appended_constraint_cols[row_idx, , drop = FALSE])
    X       <- cbind(A, V)
    penalty <- c(rep(0, ncol(A)), rep(1, Kdims))
  } else {
    X       <- V
    penalty <- rep(1, Kdims)
  }

  ## Weighted ridge regression with CV-chosen lambda
  cv_fit <- tryCatch(
    glmnet::cv.glmnet(X, outcome, alpha = 0, penalty.factor = penalty, weights = w),
    error = function(e) NULL
  )
  if (is.null(cv_fit)) return(NA_real_)

  lambda_chosen <- if (use_1se_lambda) cv_fit$lambda.1se else cv_fit$lambda.min
  yhat <- as.numeric(stats::predict(cv_fit$glmnet.fit, s = lambda_chosen, newx = X))
  e    <- outcome - yhat

  ## Influence contributions: z_i = w_i * e_i
  z      <- w * e
  zbar   <- mean(z)
  var_mu <- sum((z - zbar)^2)

  ## Small-sample correction via Kish effective sample size
  n_eff <- 1 / sum(w^2)
  if (is.finite(n_eff) && n_eff > 1) var_mu <- var_mu * n_eff / (n_eff - 1)

  var_mu
}


#' Weighted proportions and SEs for a categorical outcome
#'
#' @description One-hot encodes a categorical outcome and applies both the
#'   linearized and fixed-weight SE estimators to each level. Rows with
#'   \code{NA} in the outcome are dropped automatically; the correct rows of
#'   the kernel SVD and constraint matrices are selected via \code{row_idx}.
#'
#' @param kpop_object Object returned by \code{\link{kpop}} or
#'   \code{\link{kbal}}.
#' @param poll Data frame of respondent data containing \code{outcome_var}.
#'   \code{nrow(poll)} must equal the number of sample units (the full sample
#'   data frame, not pre-filtered -- NA filtering is handled internally).
#' @param outcome_var Quoted column name of the categorical outcome variable.
#' @param sample_size Optional override passed to
#'   \code{\link{calc_linearized_var}}. Defaults to the number of non-NA
#'   respondents for this variable.
#'
#' @return Data frame with columns: \code{level}, \code{weighted_mean},
#'   \code{linearized_SE}, \code{fixed_SE}.
#'
#' @seealso \code{\link{calc_linearized_var}}, \code{\link{kpop_summary}}
#'
#' @importFrom stats as.formula model.matrix weighted.mean
#' @importFrom dplyr bind_rows bind_cols
#' @export
calc_ses_cat_var <- function(kpop_object,
                             poll,
                             outcome_var,
                             sample_size = NULL) {

  keep     <- !is.na(poll[[outcome_var]])
  row_idx  <- which(keep)
  poll_sub <- poll[keep, , drop = FALSE]
  weights  <- kpop_object$w[row_idx]

  if (is.null(sample_size)) sample_size <- length(row_idx)

  one_hot <- stats::model.matrix(
    stats::as.formula(paste0("~ as.factor(", outcome_var, ") - 1")),
    data = poll_sub
  )

  res <- lapply(seq_len(ncol(one_hot)), function(j) {
    outcome       <- one_hot[, j]
    weighted_mean <- stats::weighted.mean(outcome, weights)
    linearized_SE <- sqrt(calc_linearized_var(kpop_object,
                                              outcome     = outcome,
                                              row_idx     = row_idx,
                                              sample_size = sample_size))
    fixed_SE      <- sqrt(calc_fixedweight_var(outcome, weights))
    data.frame(weighted_mean = weighted_mean,
               linearized_SE = linearized_SE,
               fixed_SE      = fixed_SE)
  })
  res <- dplyr::bind_rows(res)

  level_names <- gsub(paste0("as.factor[(]", outcome_var, "[)]"), "",
                      colnames(one_hot))
  dplyr::bind_cols(level = level_names, res)
}


#' Summarise outcomes from a kpop-weighted analysis
#'
#' @description For each outcome variable supplied, computes the unweighted
#'   mean (or proportion), kpop-weighted mean (or proportion), linearized SE,
#'   and fixed-weight SE. Continuous variables produce one row; categorical
#'   variables produce one row per level. Weight-quality diagnostics (ESS,
#'   bias bound) from the kpop object are also returned.
#'
#'   A variable is treated as categorical if it is a factor or character, or
#'   if it is numeric with 10 or fewer unique non-missing values, or if its
#'   name is listed in \code{cat_outcomes}.
#'
#'   Internally calls \code{glmnet::cv.glmnet} for the linearized SE.
#'   Pass \code{seed} for reproducible cross-validation.
#'
#' @param kpop_object Object returned by \code{\link{kpop}}.
#' @param data Sample-only data frame passed to \code{kpop()} as
#'   \code{sample_data}. \code{nrow(data)} must equal the number of sample
#'   units (\code{kpop_object$w[1:nrow(data)]} are the sample weights).
#' @param outcomes Character vector of column names in \code{data} to
#'   summarise.
#' @param cat_outcomes Optional character vector of column names to treat as
#'   categorical regardless of the auto-detection rule.
#' @param sample_size Optional override passed to
#'   \code{\link{calc_linearized_var}}.
#' @param seed Integer seed passed to \code{set.seed} before glmnet CV calls.
#'   \code{NULL} (default) leaves the RNG state unchanged.
#' @param use_1se_lambda Logical passed to \code{\link{calc_linearized_var}}.
#'   Default \code{TRUE}.
#'
#' @return A list of class \code{"kpop_summary"} with components:
#'   \item{estimates}{Data frame with columns \code{outcome}, \code{level},
#'     \code{n}, \code{unweighted}, \code{kpop_weighted}, \code{linearized_SE},
#'     \code{fixed_SE}. Continuous outcomes have \code{level = NA}.}
#'   \item{diagnostics}{Single-row data frame with \code{n_sample}, \code{ESS},
#'     \code{ESS_pct}, \code{biasbound_orig}, \code{biasbound_opt},
#'     \code{biasbound_ratio}, \code{numdims}.}
#'
#' @examples
#' \donttest{
#' data("lalonde")
#' set.seed(1)
#' sample_rows <- lalonde$nsw == 0
#' pop_rows    <- lalonde$nsw == 1
#' xvars <- c("age", "black", "educ", "hisp", "married", "re74", "re75")
#' out <- kpop(sample_data     = lalonde[sample_rows, xvars],
#'             population_data = lalonde[pop_rows,    xvars],
#'             mixed_data      = FALSE,
#'             minnumdims      = 1,
#'             maxnumdims      = 10)
#' kpop_summary(out,
#'              data     = lalonde[sample_rows, ],
#'              outcomes = c("re78", "black"),
#'              seed     = 2)
#' }
#'
#' @seealso \code{\link{calc_linearized_var}}, \code{\link{calc_ses_cat_var}}
#' @importFrom stats weighted.mean
#' @export
kpop_summary <- function(kpop_object,
                         data,
                         outcomes,
                         cat_outcomes   = NULL,
                         sample_size    = NULL,
                         seed           = NULL,
                         use_1se_lambda = TRUE) {

  if (!is.data.frame(data))
    stop("`data` must be a data frame.")
  if (!is.character(outcomes) || length(outcomes) == 0)
    stop("`outcomes` must be a non-empty character vector of column names.")
  missing_cols <- setdiff(outcomes, names(data))
  if (length(missing_cols) > 0)
    stop("Columns not found in `data`: ", paste(missing_cols, collapse = ", "))

  if (!is.null(seed)) set.seed(seed)

  n_sample <- nrow(data)

  is_categorical <- function(var) {
    if (!is.null(cat_outcomes) && var %in% cat_outcomes) return(TRUE)
    x <- data[[var]]
    if (is.factor(x) || is.character(x)) return(TRUE)
    if (is.numeric(x) && length(unique(x[!is.na(x)])) <= 10) return(TRUE)
    FALSE
  }

  est_list <- lapply(outcomes, function(var) {

    if (is_categorical(var)) {
      ## --- Categorical: one row per level ---
      keep  <- !is.na(data[[var]])
      n_obs <- sum(keep)

      ## Unweighted proportions
      tbl <- prop.table(table(data[[var]][keep]))
      unweighted_df <- data.frame(
        level      = as.character(names(tbl)),
        unweighted = as.numeric(tbl),
        stringsAsFactors = FALSE
      )

      se_res       <- calc_ses_cat_var(kpop_object, data, var,
                                       sample_size = sample_size)
      se_res$level <- as.character(se_res$level)

      res <- merge(unweighted_df, se_res, by = "level", sort = FALSE)
      res <- res[order(match(res$level, unweighted_df$level)), ]
      res$n <- n_obs
      names(res)[names(res) == "weighted_mean"] <- "kpop_weighted"
      res <- res[, c("level", "n", "unweighted", "kpop_weighted",
                     "linearized_SE", "fixed_SE")]
      data.frame(outcome = var, res, stringsAsFactors = FALSE,
                 row.names = NULL)

    } else {
      ## --- Continuous: one row ---
      keep    <- !is.na(data[[var]])
      row_idx <- which(keep)
      y       <- data[[var]][keep]
      w_raw   <- kpop_object$w[row_idx]

      unwtd    <- mean(y)
      kpop_wtd <- stats::weighted.mean(y, w_raw)

      lin_var <- calc_linearized_var(kpop_object, y,
                                     row_idx        = row_idx,
                                     sample_size    = sample_size,
                                     use_1se_lambda = use_1se_lambda)
      fix_var <- calc_fixedweight_var(y, w_raw)

      lin_se <- if (!is.na(lin_var) && lin_var >= 0) sqrt(lin_var) else NA_real_
      fix_se <- if (!is.na(fix_var) && fix_var >= 0) sqrt(fix_var) else NA_real_

      data.frame(
        outcome       = var,
        level         = "(cont.)",
        n             = sum(keep),
        unweighted    = unwtd,
        kpop_weighted = kpop_wtd,
        linearized_SE = lin_se,
        fixed_SE      = fix_se,
        stringsAsFactors = FALSE
      )
    }
  })

  estimates        <- do.call(rbind, est_list)
  rownames(estimates) <- NULL

  ## Diagnostics from kpop object
  w_sample <- kpop_object$w[seq_len(n_sample)]
  ess_val  <- sum(w_sample)^2 / sum(w_sample^2)

  diagnostics <- data.frame(
    n_sample        = n_sample,
    ESS             = ess_val,
    ESS_pct         = ess_val / n_sample,
    biasbound_orig  = kpop_object$biasbound_orig,
    biasbound_opt   = kpop_object$biasbound_opt,
    biasbound_ratio = kpop_object$biasbound_ratio,
    numdims         = kpop_object$numdims,
    stringsAsFactors = FALSE
  )

  structure(
    list(estimates = estimates, diagnostics = diagnostics),
    class = "kpop_summary"
  )
}


#' Print method for kpop_summary objects
#'
#' @param x Object of class \code{"kpop_summary"}.
#' @param digits Number of significant digits for the estimates table.
#'   Default \code{3}.
#' @param ... Further arguments (currently unused).
#' @return Invisibly returns \code{x}.
#' @export
print.kpop_summary <- function(x, digits = 3, ...) {
  d <- x$diagnostics
  cat(sprintf(
    "-- kpop diagnostics ------------------------------------------\n  n = %d  |  ESS = %.1f (%.1f%%)  |  bias-bound ratio = %.4f  |  dims = %d\n",
    d$n_sample, d$ESS, d$ESS_pct * 100, d$biasbound_ratio, d$numdims
  ))
  cat("-- estimates -------------------------------------------------\n")
  est <- x$estimates
  num_cols <- c("unweighted", "kpop_weighted", "linearized_SE", "fixed_SE")
  for (col in intersect(num_cols, names(est))) {
    est[[col]] <- signif(est[[col]], digits)
  }
  print(est, row.names = FALSE)
  invisible(x)
}
