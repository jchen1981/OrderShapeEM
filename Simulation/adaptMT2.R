

#===============================================================
# adapt_model class
#===============================================================

#' adapt_model Objects for M-steps
#'
#'
#' \code{adapt_model} objects provide the functions and their arguments in computing the M-steps.
#' Each object can be passed to \code{\link{adapt}} as a candidate model.
#'
#' \code{pifun} should be in the form of \code{pifun(formula, data, family, weights, ...)} or \code{pifun(x, y, family, ...)}.
#' The former includes \code{\link[stats]{glm}} and \code{\link[mgcv]{gam}} and the latter includes \code{\link[glmnet]{glmnet}}.
#' The outputs should be in the form of \code{list(fitv = , info = , ...)} where \code{fitv} gives the estimate of pi(x),
#' as a vector with the same order of \code{x}, and \code{info} should at least contain a key \code{df} if model selection is used, i.e. \code{info = list(df = , ...)}
#'
#' \code{mufun} should be in the form of \code{pifun(formula, data, family, weights, ...)} or \code{pifun(x, y, family, weights, ...)}.
#' Note that \code{mufun} must take \code{weights} as an input. The outputs should be in the same form as \code{pifun} except that \code{fitv} should give the estimate of mu(x).
#'
#' When \code{pifun} / \code{mufun} takes the form of \code{(formula, family, ...)}, \code{piargs} / \code{muargs} should at least contain a key \code{formula}; when \code{pifun} / \code{mufun} takes the form of \code{(x, y, family, ...)}, \code{piargs} / \code{muargs} can be empty.
#'
#' For glm/gam/glmnet, one can use the shortcut by running \code{\link{gen_adapt_model}} with name = "glm" or "gam" or "glmnet" but without specifying \code{pifun}, \code{mufun}, \code{pifun_init} and \code{mufun_init}. See examples below.
#'
#' @param pifun a function to fit pi(x). See Details
#' @param mufun a function to fit mu(x). See Details
#' @param pifun_init a function to fit pi(x) at the initial step
#' @param mufun_init a function to fit mu(x) at the initial step
#' @param piargs a list. Arguments for "pifun". An empty list as default
#' @param muargs a list. Arguments for "mufun". An empty list as default
#' @param piargs_init a list. Arguments for piargs_init. An empty list as default
#' @param muargs_init a list. Arguments for muargs_init. An empty list as default
#' @param name a string. An optional argument for the user-specified name of the model. An empty string as default.
#'
#' @return
#' \item{name}{same as the input \code{name}}
#' \item{algo}{a list recording \code{pifun}, \code{mufun}, \code{pifun_init} and \code{mufun_init}}
#' \item{args}{a list recording \code{piargs}, \code{muargs}, \code{piargs_init} and \code{muargs_init}}
#'
#' @examples
#' \donttest{
#' # Exemplary code to generate 'adapt_model' for logistic-Gamma glm  with naive initialization.
#' # The real implementation in the package is much more complicated.
#'
#' # pifun as a logistic regression
#' pifun <- function(formula, data, weights, ...){
#'   glm(formula, data, weights = weights, family = binomial(),  ...)
#' }
#' # pifun_init as a constant
#' pifun_init <- function(x, pvals, s, ...){
#'   rep(0.1, length(pvals))
#' }
#' # mufun as a Gamma GLM
#' mufun <- function(formula, data, weights, ...){
#'   glm(formula, data, weights = weights, family = Gamma(), ...)
#' }
#' # mufun_init as a constant
#' mufun_init <- function(x, pvals, s, ...){
#'   rep(1.5, length(pvals))
#' }
#'
#' library("splines") # for using ns() in the formula
#' piargs <- list(formula = "ns(x, df = 8)")
#' muargs <- list(formula = "ns(x, df = 8)")
#' name <- "glm"
#'
#' mod <- gen_adapt_model(pifun, mufun, pifun_init, mufun_init,
#'                        piargs, muargs, name = name)
#' mod
#'
#' # Using shortcut for GLM. See the last paragraph of Details.
#' mod2 <- gen_adapt_model(name = "glm", piargs = piargs, muargs = muargs)
#' mod2
#' }
#'
#' @export
gen_adapt_model <- function(pifun = NULL,
		mufun = NULL,
		pifun_init = NULL,
		mufun_init = NULL,
		piargs = list(),
		muargs = list(),
		piargs_init = list(),
		muargs_init = list(),
		name = ""){
	args <- list(piargs = piargs, muargs = muargs,
			piargs_init = piargs_init, muargs_init = muargs_init)
	
	if (is.null(pifun) && is.null(mufun) &&
			is.null(pifun_init) && is.null(mufun_init)){
		model <- structure(
				list(name = name, args = args),
				class = "adapt_model"
		)
		return(model)
	}
	
	if (!is.function(pifun)){
		stop("\"pifun\" must be a function.")
	}
	if (!is.function(mufun)){
		stop("\"mufun\" must be a function.")
	}
	if (!is.function(pifun_init)){
		stop("\"pifun_init\" must be a function.")
	}
	if (!is.function(mufun_init)){
		stop("\"mufun_init\" must be a function.")
	}
	
	algo <- list(pifun = pifun, mufun = mufun,
			pifun_init = pifun_init, mufun_init = mufun_init)
	
	model <- structure(
			list(name = name, algo = algo, args = args),
			class = "adapt_model"
	)
	return(model)
}

gen_adapt_model_glm <- function(dist,
		piargs = list(),
		muargs = list()){
	pifun <- function(formula, data, weights, ...){
		safe_glm(formula, data, weights = weights,
				family = quasibinomial(), ...)
	}
	
	mufun <- function(formula, data, weights, ...){
		safe_glm(formula, data, weights = weights,
				family = dist$family, ...)
	}
	
	pifun_init <- function(x, pvals, s, ...){
		J <- ifelse(
				pvals < s | pvals > 1 - s, 1,
				2 * s / (2 * s - 1)
		)
		if (any(s >= 0.49)){
			J[s >= 0.49] <- 0
		}
		
		fun <- function(formula, data, ...){
			safe_glm(formula, data, family = gaussian(), ...)
		}
		
		args <- list(...)
		args <- complete_args(x, J, fun, args)
		
		fit_pi(fun, args, type = "init")
	}
	
	mufun_init <- function(x, pvals, s, ...){
		phat <- ifelse(
				pvals < s | pvals > 1 - s,
				pmin(pvals, 1 - pvals),
				pvals
		)
		phat <- pminmax(phat, 1e-15, 1-1e-15)
		yhat <- dist$g(phat)
		
		fun <- function(formula, data, ...){
			safe_glm(formula, data, family = dist$family, ...)
		}
		
		args <- list(...)
		args <- complete_args(x, yhat, fun, args)
		
		fit_mu(fun, args, dist, type = "init")
	}
	
	if (is.null(piargs$formula) || is.null(muargs$formula)){
		stop("Argument \"formula\" is missing from \"piargs\" or \"muargs\".")
	}
	
	piargs_init <- piargs
	muargs_init <- muargs
	
	gen_adapt_model(pifun, mufun, pifun_init, mufun_init,
			piargs, muargs, piargs_init, muargs_init,
			name = "glm")
}

gen_adapt_model_gam <- function(dist,
		piargs = list(),
		muargs = list()){
	pifun <- function(formula, data, weights, ...){
		safe_gam(formula, data, weights = weights,
				family = quasibinomial(), ...)
	}
	
	mufun <- function(formula, data, weights, ...){
		safe_gam(formula, data, weights = weights,
				family = dist$family, ...)
	}
	
	pifun_init <- function(x, pvals, s, ...){
		J <- ifelse(
				pvals < s | pvals > 1 - s, 1,
				2 * s / (2 * s - 1)
		)
		if (any(s >= 0.49)){
			J[s >= 0.49] <- 0
		}
		
		fun <- function(formula, data, ...){
			safe_gam(formula, data, family = gaussian(), ...)
		}
		
		args <- list(...)
		args <- complete_args(x, J, fun, args)
		
		fit_pi(fun, args, type = "init")
	}
	
	mufun_init <- function(x, pvals, s, ...){
		phat <- ifelse(
				pvals < s | pvals > 1 - s,
				pmin(pvals, 1 - pvals),
				pvals
		)
		phat <- pminmax(phat, 1e-15, 1-1e-15)
		yhat <- dist$g(phat)
		
		fun <- function(formula, data, ...){
			safe_gam(formula, data, family = dist$family, ...)
		}
		
		args <- list(...)
		args <- complete_args(x, yhat, fun, args)
		
		fit_mu(fun, args, dist, type = "init")
	}
	
	if (is.null(piargs$formula) || is.null(muargs$formula)){
		stop("Argument \"formula\" is missing from \"piargs\" or \"muargs\".")
	}
	
	piargs_init <- piargs
	muargs_init <- muargs
	
	gen_adapt_model(pifun, mufun, pifun_init, mufun_init,
			piargs, muargs, piargs_init, muargs_init,
			name = "gam")
}

gen_adapt_model_glmnet <- function(dist,
		piargs = list(),
		muargs = list()){
	pifun <- function(x, y, weights, ...){
		safe_glmnet(x, y, weights = weights,
				family = "binomial", ...)
	}
	
	
	mufun <- function(x, y, weights, ...){
		safe_glmnet(x, y, weights = weights,
				family = dist$family, ...)
	}
	
	pifun_init <- function(x, pvals, s, ...){
		J <- ifelse(
				pvals < s | pvals > 1 - s, 1,
				2 * s / (2 * s - 1)
		)
		if (any(s >= 0.49)){
			J[s >= 0.49] <- 0
		}
		
		fun <- function(x, y, ...){
			safe_glmnet(x, y, family = "gaussian", ...)
		}
		
		args <- list(...)
		args <- complete_args(x, J, fun, args)
		
		fit_pi(fun, args, type = "init")
	}
	
	mufun_init <- function(x, pvals, s, ...){
		phat <- ifelse(
				pvals < s | pvals > 1 - s,
				pmin(pvals, 1 - pvals),
				pvals
		)
		phat <- pminmax(phat, 1e-15, 1-1e-15)
		yhat <- dist$g(phat)
		
		fun <- function(x, y, ...){
			safe_glmnet(x, y, family = dist$family$family, ...)
		}
		
		args <- list(...)
		args <- complete_args(x, yhat, fun, args)
		
		fit_mu(fun, args, dist, type = "init")
	}
	
	piargs_init <- piargs
	muargs_init <- muargs
	
	gen_adapt_model(pifun, mufun, pifun_init, mufun_init,
			piargs, muargs, piargs_init, muargs_init,
			name = "glmnet")
}

#---------------------------------------------------------------
# Wrappers of AdaPT
#---------------------------------------------------------------

check_formulas <- function(formulas){
	err <- FALSE
	if (is.character(formulas)){
		formulas <- as.list(formulas)
	} else if (is.list(formulas)){
		err <- any(!sapply(formulas, class) %in% c("character", "formula"))
	} else {
		err <- TRUE
	}
	
	if (err){
		stop("Invalid formulas")
	}
	
	return(formulas)
}

#' Adaptive P-value Thresholding with Generalized Linear Models
#'
#' \code{adapt_glm} is a wrapper of \code{\link{adapt}} that fits pi(x) and mu(x) by \code{\link[stats]{glm}}.
#'
#' \code{pi_formulas} and \code{mu_formulas} can either be a list or a vector with each element being a string or a formula. For instance, suppose \code{x} has a single column with name \code{x1}, the following five options are valid for the same inputs (\code{\link[splines]{ns}} forms a spline basis with \code{df} knots):
#' \enumerate{
#' \item{c("x1", "ns(x1, df = 8)");}
#' \item{c("~ x1", "~ ns(x1, df = 8)");}
#' \item{list("x1", "ns(x1, df = 8)");}
#' \item{list("~ x1", "~ ns(x1, df = 8)");}
#' \item{list(~ x1, ~ ns(x1, df = 8))}
#' }
#' There is no need to specify the name of the response variable, as this is handled in the function.
#'
#' When \code{x} has a few variables, it is common to use non-parametric GLM by replacing \code{x} by a spline basis of \code{x}. In this case, \code{\link[splines]{ns}} from \code{library(splines)} package is suggested.
#' 
#' @param pi_formulas a vector/list of strings/formulas. Formulas for fitting pi(x) by glm. See Details
#' @param mu_formulas a vector/list of strings/formulas. Formulas for fitting mu(x) by glm. See Details
#' @param piargs a list. Other arguments passed to glm for fitting pi(x)
#' @param muargs a list. Other arguments passed to glm for fitting mu(x)
#' @param ... other arguments passed to \code{\link{adapt}} (except \code{models})
#' @inheritParams adapt
#'
#' @examples
#' \donttest{
#' # Load estrogen data
#' data(estrogen)
#' pvals <- as.numeric(estrogen$pvals)
#' x <- data.frame(x = as.numeric(estrogen$ord_high))
#' dist <- beta_family()
#'
#' # Subsample the data for convenience
#' inds <- (x$x <= 5000)
#' pvals <- pvals[inds]
#' x <- x[inds,,drop = FALSE]
#' 
#' # Run adapt_glm
#' library("splines")
#' formulas <- paste0("ns(x, df = ", 6:10, ")")
#' res <- adapt_glm(x = x, pvals = pvals, pi_formulas = formulas,
#'                  mu_formulas = formulas, dist = dist, nfits = 10)
#'
#' # Run adapt by manually setting models for glm
#' models <- lapply(formulas, function(formula){
#'     piargs <- muargs <- list(formula = formula)
#'     gen_adapt_model(name = "glm", piargs = piargs, muargs = muargs)
#' })
#' res2 <- adapt(x = x, pvals = pvals, models = models,
#'               dist = dist, nfits = 10)
#' 
#' # Check equivalence
#' identical(res, res2)
#' }
#'
#' @seealso
#' \code{\link{adapt}}, \code{\link{adapt_gam}}, \code{\link{adapt_glmnet}}, \code{\link[stats]{glm}}, \code{\link[splines]{ns}}
#' 
#' @export
adapt_glm2 <- function(x, pvals, pi_formulas, mu_formulas,
		dist = beta_family(),
		s0 = rep(0.45, length(pvals)),
		alphas = seq(0.01, 1, 0.01),
		piargs = list(), muargs = list(),
		...){
	if (!is.data.frame(x)){
		stop("\'x\' must be a data.frame")
	}
	
	pi_formulas <- check_formulas(pi_formulas)
	mu_formulas <- check_formulas(mu_formulas)
	stopifnot(length(pi_formulas) == length(mu_formulas))
	
	models <- lapply(1:length(pi_formulas), function(i){
				piargs <- c(list(formula = pi_formulas[[i]]), piargs)
				muargs <- c(list(formula = mu_formulas[[i]]), muargs)
				gen_adapt_model(name = "glm", piargs = piargs, muargs = muargs)
			})
	
	adapt2(x, pvals, models, dist, s0, alphas, ...)
}

#' Adaptive P-value Thresholding with Generalized Additive Models
#'
#' \code{adapt_gam} is a wrapper of \code{\link{adapt}} that fits pi(x) and mu(x) by \code{\link[mgcv]{gam}} from \code{mgcv} package.
#'
#' \code{pi_formulas} and \code{mu_formulas} can either be a list or a vector with each element being a string or a formula. For instance, suppose \code{x} has a single column with name \code{x1}, the following five options are valid for the same inputs (\code{\link[splines]{ns}} forms a spline basis with \code{df} knots and \code{\link[mgcv]{s}} forms a spline basis with knots automatically selected by generalized cross-validation):
#' \enumerate{
#' \item{c("x1", "ns(x1, df = 8)", "s(x1)");}
#' \item{c("~ x1", "~ ns(x1, df = 8)", "s(x1)");}
#' \item{list("x1", "ns(x1, df = 8)", "s(x1)");}
#' \item{list("~ x1", "~ ns(x1, df = 8)", "s(x1)");}
#' \item{list(~ x1, ~ ns(x1, df = 8), s(x1))}
#' }
#' There is no need to specify the name of the response variable, as this is handled in the function.
#'
#' When \code{x} has a few variables, it is common to use non-parametric GLM by replacing \code{x} by a spline basis of \code{x}. In this case, \code{\link[splines]{ns}} from \code{library(splines)} package or \code{\link[mgcv]{s}} from \code{mgcv} package are suggested. When \code{\link[mgcv]{s}} (from \code{mgcv} package) is used, it is treated as a single model because the knots will be selected automatically.
#' 
#' @param pi_formulas a vector/list of strings/formulas. Formulas for fitting pi(x) by gam. See Details
#' @param mu_formulas a vector/list of strings/formulas. Formulas for fitting mu(x) by gam. See Details
#' @param piargs a list. Other arguments passed to gam for fitting pi(x)
#' @param muargs a list. Other arguments passed to gam for fitting mu(x)
#' @param ... other arguments passed to \code{\link{adapt}} (except \code{models})
#' @inheritParams adapt
#'
#' @seealso
#' \code{\link{adapt}}, \code{\link{adapt_glm}}, \code{\link{adapt_glmnet}}, \code{\link[mgcv]{gam}}, \code{\link[splines]{ns}}, \code{\link[mgcv]{s}}
#'
#' @examples
#' \donttest{
#' # Generate a 2-dim x
#' n <- 400
#' x1 <- x2 <- seq(-100, 100, length.out = 20)
#' x <- expand.grid(x1, x2)
#' colnames(x) <- c("x1", "x2")
#'
#' # Generate p-values (one-sided z test)
#' # Set all hypotheses in the central circle with radius 30 to be
#' # non-nulls. For non-nulls, z~N(2,1) and for nulls, z~N(0,1).
#' H0 <- apply(x, 1, function(coord){sum(coord^2) < 900})
#' mu <- ifelse(H0, 2, 0)
#' set.seed(0)
#' zvals <- rnorm(n) + mu
#' pvals <- 1 - pnorm(zvals)
#'
#' # Run adapt_gam with a 2d spline basis
#' library("mgcv")
#' formula <- "s(x1, x2)"
#' dist <- beta_family()
#' res <- adapt_gam(x = x, pvals = pvals, pi_formulas = formula,
#'                  mu_formulas = formula, dist = dist, nfits = 5)
#' }
#' 
#' 
#' @export
adapt_gam <- function(x, pvals, pi_formulas, mu_formulas,
		piargs = list(), muargs = list(),
		dist = beta_family(),
		s0 = rep(0.45, length(pvals)),
		alphas = seq(0.01, 1, 0.01),
		...){
	if (!is.data.frame(x)){
		stop("\'x\' must be a data.frame")
	}
	
	if (!requireNamespace("mgcv", quietly = TRUE)){
		stop("'mgcv' package is required for 'adapt_gam'. Please intall.")
	}
	
	pi_formulas <- check_formulas(pi_formulas)
	mu_formulas <- check_formulas(mu_formulas)
	stopifnot(length(pi_formulas) == length(mu_formulas))
	
	models <- lapply(1:length(pi_formulas), function(i){
				piargs <- c(list(formula = pi_formulas[[i]]), piargs)
				muargs <- c(list(formula = mu_formulas[[i]]), muargs)
				gen_adapt_model(name = "gam", piargs = piargs, muargs = muargs)
			})
	
	adapt(x, pvals, models, dist, s0, alphas, ...)
}

#' Adaptive P-value Thresholding with L1/L2 Penalized Generalized Linear Models
#'
#' \code{adapt_glmnet} is a wrapper of \code{\link{adapt}} that fits pi(x) and mu(x) by \code{\link[glmnet]{glmnet}} from \code{glmnet} package.
#'
#' \code{adapt_glmnet} by default implements LASSO on \code{x} with lambda selected by cross-validation. Specify in \code{piargs} and \code{muargs} if ridge or elastic-net penalty is needed.
#' 
#' @param piargs a list. Other arguments passed to glmnet for fitting pi(x)
#' @param muargs a list. Other arguments passed to glmnet for fitting mu(x)
#' @param ... other arguments passed to \code{\link{adapt}} (except \code{models})
#' @inheritParams adapt
#'
#' 
#' @seealso
#' \code{\link{adapt}}, \code{\link{adapt_glm}}, \code{\link{adapt_gam}}, \code{\link[glmnet]{glmnet}}
#'
#' @examples
#' \donttest{
#' # Generate a 100-dim covariate x
#' set.seed(0)
#' m <- 100
#' n <- 1000
#' x <- matrix(runif(n * m), n, m)
#'
#' # Generate the parameters from a conditional two-group
#' # logistic-Gamma GLM  where pi(x) and mu(x) are both
#' # linear in x. pi(x) has an intercept so that the average
#' # of pi(x) is 0.3
#' inv_logit <- function(x) {exp(x) / (1 + exp(x))}
#' pi1 <- 0.3
#' beta.pi <- c(3, 3, rep(0, m-2))
#' beta0.pi <- uniroot(function(b){
#'     mean(inv_logit(x %*% beta.pi + b)) - pi1
#' }, c(-100, 100))$root
#' pi <- inv_logit(x %*% beta.pi + beta0.pi)
#' beta.mu <- c(2, 2, rep(0, m-2))
#' beta0.mu <- 0
#' mu <- pmax(1, x %*% beta.mu + beta0.mu)
#'
#' # Generate p-values
#' H0 <- as.logical(ifelse(runif(n) < pi, 1, 0))
#' y <- ifelse(H0, rexp(n, 1/mu), rexp(n, 1))
#' pvals <- exp(-y)
#'
#' # Run adapt_glmnet
#' res <- adapt_glmnet(x, pvals, s0 = rep(0.15, n), nfits = 5)
#' }
#' @export
adapt_glmnet <- function(x, pvals,
		piargs = list(), muargs = list(),
		dist = beta_family(),
		s0 = rep(0.45, length(pvals)),
		alphas = seq(0.01, 1, 0.01),
		...){
	if (!is.matrix(x) && !inherits(x, "sparseMatrix")){
		stop("Invalid \'x\'. See \'?glmnet\' for details.")
	}
	
	if (!requireNamespace("glmnet", quietly = TRUE)){
		stop("'glmnet' package is required for 'adapt_glmnet'. Please intall.")
	}
	
	models <- gen_adapt_model(name = "glmnet", piargs = piargs, muargs = muargs)
	
	adapt(x, pvals, models, dist, s0, alphas, ...)
}

#---------------------------------------------------------------
# Compute oracle local FDR estimate by revealing all
# p-values and compute the correlation between the oracle
# estimate and estimate within each step
#---------------------------------------------------------------

#' Quantifying Information Loss of Adaptive P-Value Thresholding
#'
#' \code{corr_lfdr} computes the oracle local FDR estimate, by using revealing all p-values, and computes the Pearson correlation between it and the estimate within each step of \code{adapt}.
#'
#' @param obj an 'adapt' object. Output of \code{\link{adapt}} function
#' @param x covariates (i.e. side-information). Should be compatible to \code{models}.
#' @param pvals a vector of values in [0, 1]. P-values
#' @param model an optional argument. If \code{model = NULL} then the last model in \code{obj$models} is used for fitting the oracle model (i.e. with all p-values revealed). Otherwise it should be an 'adapt_model' object
#' @param niter_oracle an positive integer. Number of iterations in EM algorithm
#'
#' @return
#' \itemize{
#' \item{corr}{a vector of values in [0, 1]. Pearson correlation of oracle local FDR estimate and the estimates within each step. Each value corresponds to an entry of \code{obj$params}}
#' \item{oracle_lfdr}{a vector of values in [0, 1]. Oracle local FDR estimate}
#' \item{lfdr}{a matrix of values in [0, 1]. Local FDR estimates within each step.}
#' \item{alphas}{a vector of values in [0, 1]. The target FDR levels corresponding to each local FDR estimate}
#' \item{nmasks}{a vector of integers. The number of masked p-values corresponding to each local FDR estimate}
#' }
#'
#' @examples
#' \donttest{
#' # Load estrogen data
#' data(estrogen)
#' pvals <- as.numeric(estrogen$pvals)
#' x <- data.frame(x = as.numeric(estrogen$ord_high))
#' dist <- beta_family()
#'
#' # Subsample the data for convenience
#' inds <- (x$x <= 5000)
#' pvals <- pvals[inds]
#' x <- x[inds,,drop = FALSE]
#'
#' # Run adapt_glm
#' library("splines")
#' formulas <- paste0("ns(x, df = ", 6:10, ")")
#' res <- adapt_glm(x = x, pvals = pvals, pi_formulas = formulas,
#'                  mu_formulas = formulas, dist = dist, nfits = 10)
#'
#' # Run corr_lfdr
#' obj <- corr_lfdr(res, x, pvals)
#' obj$corr
#' }
#'
#' @export
corr_lfdr <- function(obj, x, pvals, model = NULL,
		niter_oracle = 100){
	if (class(obj)[1] != "adapt"){
		stop("\'obj\' is not of class \'adapt\"")
	}
	n <- length(pvals)
	
	params_list <- obj$params
	dist <- obj$dist
	m <- length(params_list)
	alphas <- rep(0, m)
	nmasks <- rep(0, m)
	lfdr <- matrix(0, n, m)
	for (i in 1:m){
		params <- list(pix = params_list[[i]]$pix,
				mux = params_list[[i]]$mux)
		alphas[i] <- params_list[[i]]$alpha
		nmasks[i] <- params_list[[i]]$nmasks
		lfdr[, i] <- compute_lfdr_mix(pvals, dist, params)
	}
	
	## Oracle lfdr
	if (is.null(model)){
		model <- obj$models[[length(obj$models)]]
	}
	oracle_params <- EM_mix(x, pvals, rep(0, n), dist, model,
			params0 = params,
			niter = niter_oracle)$params
	oracle_lfdr <- compute_lfdr_mix(pvals, dist, oracle_params)
	
	## Correlation
	corr <- apply(lfdr, 2, function(x){cor(x, oracle_lfdr)})
	
	return(list(corr = corr,
					oracle_lfdr = oracle_lfdr,
					lfdr = lfdr,
					alphas = alphas, nmasks = nmasks))
}

#---------------------------------------------------------------
# Compute oracle local FDR estimate by revealing all p-values
#---------------------------------------------------------------

#' Fitting Conditional Two-Groups Models on Unmasked P-Values
#'
#' \code{ctgm_lfdr} computes the oracle local FDR estimate, by using all p-values without masking.
#'
#' \code{ctgm_lfdr} implements the EM algorithm to fit pi(x) and mu(x) on unmasked p-values. Although it is not related to FDR control of AdaPT, it provides useful measures for post-hoc justification and other purposes.
#' For instance, one can use these local FDR estimates for prioritizing the hypotheses if strict FDR control is not required.
#'
#' In contrast to \code{adapt}, \code{cytm_lfdr} does not guarantee FDR control unless the model is correctly specified. It is recommended to use \code{ctgm_lfdr} only when FDR control is not required.
#'
#' \code{x} should have a type compatible to the fitting functions in \code{models}. For GLM and GAM, \code{x} should be a data.frame. For glmnet, \code{x} should be a matrix.
#'
#' \code{models} could either be an \code{adapt_model} object, if a single model is used, or a list of \code{adapt_model} objects, each of which corresponding to a model. Each element should be generated by \code{\link{gen_adapt_model}}. For glm/gam/glmnet, one can use the shortcut by running \code{\link{gen_adapt_model}} with name = "glm" or "gam" or "glmnet" but without specifying \code{pifun}, \code{mufun}, \code{pifun_init} and \code{mufun_init}. See examples below.
#'
#' When \code{type = "over"}, it yields a conservative estimate of local FDR
#' \deqn{lfdr(p) = (1 - \pi_{1} + \pi_{1}f_{1}(1)) / (1 - \pi_{1} + \pi_{1}f_{1}(p)).}
#' When \code{type = "raw"}, it yields the original local FDR.
#' \deqn{lfdr(p) = (1 - \pi_{1}) / (1 - \pi_{1} + \pi_{1}f_{1}(p)).}
#' The former is shown to be more stable and reliable because the weak identifiability in conditional mixture models.
#'
#' @param x covariates (i.e. side-information). Should be compatible to \code{models}. See Details
#' @param pvals a vector of values in [0, 1]. P-values
#' @param models an object of class "\code{adapt_model}" or a list of objects of class "adapt_model". See Details
#' @param dist an object of class "\code{\link{gen_exp_family}}". \code{\link{beta_family}()} as default
#' @param type a character. Either "over" or "raw" indicating the type of local FDR estimates. See Details
#' @param params0 a list in the form of list(pix = , mux = ). Initial values of pi(x) and mu(x). Both can be set as NULL
#' @param niter a positive integer. Number of EM iterations.
#' @param cr a string. The criterion for model selection with BIC as default. Also support AIC, AICC and HIC
#' @param verbose a logical values in the form of list(fit = , ms = ). Indicate whether the progress of model fitting and model selection is displayed
#'
#' @return
#' \itemize{
#' \item{lfdr}{a vector of values in [0, 1]. Local FDR estimates of each hypothesis.}
#' \item{model}{an \code{adapt_model} object. The selected model if multiple models are provided.}
#' }
#'
#' @examples
#' \donttest{
#' # Load estrogen data
#' data(estrogen)
#' pvals <- as.numeric(estrogen$pvals)
#' x <- data.frame(x = as.numeric(estrogen$ord_high))
#' dist <- beta_family()
#'
#' # Subsample the data for convenience
#' inds <- (x$x <= 5000)
#' pvals <- pvals[inds]
#' x <- x[inds,,drop = FALSE]
#'
#' # Generate models for function adapt
#' library("splines")
#' formulas <- paste0("ns(x, df = ", 6:10, ")")
#' models <- lapply(formulas, function(formula){
#'     piargs <- muargs <- list(formula = formula)
#'     gen_adapt_model(name = "glm", piargs = piargs, muargs = muargs)
#' })
#'
#' # Run ctgm_lfdr with two types of lfdr estimates
#' res_over <- ctgm_lfdr(x, pvals, models, type = "over")
#' res_raw <- ctgm_lfdr(x, pvals, models, type = "raw")
#'
#' # Compare two estimates
#' par(mfrow = c(2, 1))
#' hist(res_over$lfdr)
#' hist(res_raw$lfdr)
#' }
#'
#' @export
ctgm_lfdr <- function(x, pvals, models,
		dist = beta_family(),
		type = c("over", "raw"),
		params0 = list(pix = NULL, mux = NULL),
		niter = 50,
		cr = "BIC",
		verbose = TRUE
){
	Mstep_type <- "weighted"
	lfdr_type <- type[1]
	tol <- 0
	
	res <- EM_mix_ms(x, pvals, rep(0, length(pvals)), dist, models, cr,
			params0, niter, tol, verbose, Mstep_type)
	lfdr <- compute_lfdr_mix(pvals, dist, res$params, lfdr_type)
	return(list(lfdr = lfdr, model = res$model))
}

#---------------------------------------------------------------
# Model selection based on partially masked data
#---------------------------------------------------------------

info_cr <- function(loglik, cr, df, n){
	switch(cr,
			"AIC" = 2 * df - 2 * loglik,
			"AICC" = 2 * df * n / (n - df - 1),
			"BIC" = log(n) * df - 2 * loglik,
			"HIC" = 2 * log(log(n)) * df - 2 * loglik)
}

EM_mix_ms <- function(x, pvals, s, dist, models,
		cr = "BIC",
		params0 = list(pix = NULL, mux = NULL),
		niter = 20, tol = 1e-4,
		verbose = TRUE,
		type = "unweighted"){
	n <- length(pvals)
	info_cr_val <- -Inf
	
	m <- length(models)
	if (verbose){
		cat("Model selection starts!\n")
		cat("Shrink the set of candidate models if it is too time-consuming.")
		cat("\n")
		pb <- txtProgressBar(min = 0, max = m, style = 3, width = 50)
	}
	for (i in 1:m){
		model <- complete_model(models[[i]], dist)
		fit <- try(
				EM_mix(x, pvals, s, dist, model, params0, niter, tol,
						type = type),
				silent = TRUE
		)
		if (class(fit)[1] == "try-error"){
			warning(paste0("Model ", i, " fails."))
			next
		}
		
		loglik <- fit$loglik
		df <- fit$info$pi$df + fit$info$mu$df
		val <- info_cr(loglik, cr, df, n)
		if (val > info_cr_val){
			params <- fit$params
			info_cr_val <- val
			best_model <- models[[i]]
			best_model_info <- fit$info
		}
		
		if (verbose){
			setTxtProgressBar(pb, i)
		}
	}
	if (verbose){    
		cat("\n")
	}
	if (info_cr_val == -Inf){
		stop("All models fail.")
	}
	return(list(model = best_model,
					params = params,
					info = best_model_info))
}
#---------------------------------------------------------------
# EM algorithm to fit a mixture model.
#---------------------------------------------------------------

EM_loglik <- function(pvals, dist, pix, mux, Hhat, bhat){
	loglik1 <- sum(Hhat * log(pix) + (1 - Hhat) * log(1 - pix))
	loglik2 <- sum(Hhat * bhat * log(dist$h(pvals, mux)) +
					Hhat * (1 - bhat) * log(dist$h(1 - pvals, mux)))
	return(loglik1 + loglik2)
}

EM_mix <- function(x, pvals, s, dist, model,
		params0 = list(pix = NULL, mux = NULL),
		niter = 10, tol = 1e-4,
		verbose = FALSE,
		type = "unweighted"){
	model <- complete_model(model, dist)
	if (verbose){
		cat("Model fitting starts!\n")
		cat("Reduce niter if it is too time-consuming.\n")
		pb <- utils::txtProgressBar(min = 0, max = niter, style = 3, width = 50)
	}
	
	if (is.null(params0$pix) || is.null(params0$mux)){
		piargs_init <- c(list(x = x, pvals = pvals, s = s),
				model$args$piargs_init)
		pix <- do.call(model$algo$pifun_init, piargs_init)$fitv
		muargs_init <- c(list(x = x, pvals = pvals, s = s),
				model$args$muargs_init)
		mux <- do.call(model$algo$mufun_init, muargs_init)$fitv
		
		## init_res <- init_mix(
		## x, pvals, s, dist,
		## model$algo$pifun_init, model$algo$mufun_init,
		## model$args$piargs_init, model$args$muargs_init)
		
		old_pix <- pix ## <- init_res$pix
		old_mux <- mux ## <- init_res$mux
	} else {
		old_pix <- pix <- params0$pix
		old_mux <- mux <- params0$mux
	}
	
	for (step in 1:niter){
		Estep_res <-
				Estep_mix(pvals, s, dist, pix, mux)
		Mstep_res <-
				Mstep_mix(x, pvals, dist,
						Estep_res$Hhat, Estep_res$bhat, 
						model$algo$pifun, model$algo$mufun,
						model$args$piargs, model$args$muargs,
						type = type[1])
		pix <- Mstep_res$pix
		mux <- Mstep_res$mux
		if (max(abs(mux - old_mux)) < tol &&
				max(abs(pix - old_pix)) < tol){
			break
		}
		old_pix <- pix
		old_mux <- mux
		if (verbose){
			utils::setTxtProgressBar(pb, step)
		}        
	}
	if (verbose){    
		cat("\n")
	}
	params <- list(pix = pix, mux = mux)
	loglik <- EM_loglik(pvals, dist, params$pix, params$mux,
			Estep_res$Hhat, Estep_res$bhat)
	info <- list(pi = Mstep_res$pi_info, mu = Mstep_res$mu_info)
	
	return(list(params = params, loglik = loglik, info = info))
}
#===============================================================
# Compute E-step for the mixture model.
#===============================================================

## ' Computing E-step for Mixture Models
## '
## ' \code{Estep_mix} computes the E-step, namely imputation of missing values, including the p-values and indicators of null/non-null hypotheses.
## '
## ' The p-values are assumed to be generated from an covariate-varying mixture model
## ' \deqn{H_{i}\sim Ber(\pi(x_{i})), p_{i} \mid (H_{i} = 0) \sim U([0, 1]), p_{i} \mid (H_{i} = 1) \sim h(\cdot; \mu(x_{i}))}{Hi ~ Ber(\pi(xi)), pi | (Hi = 0) ~ U([0, 1]), pi | (Hi = 1) ~ h(p; \mu(xi))}
## ' where \eqn{h(p; \mu)} is the density of an exponential family (see \code{\link{exp_family}}). Given a threshold curve s(x), the partially i-th masked p-value is defined as \eqn{p_{i}}{pi} if \eqn{p_{i}\in [s(x_{i}), 1-s(x_{i})]}{pi in [s(xi), 1-s(xi)]}. The E-step computes the expectation of \eqn{H_{i}}{Hi} and \eqn{g(p_{i})H_{i}}{g(pi)Hi} given the parameters \eqn{\pi(x_{i})}{\pi(xi)} and \eqn{\mu(x_{i})}{\mu(xi)}, based on which computes the imputed values for \eqn{H_{i}}{Hi} and \eqn{p_{i}}{pi}. See ...
## '
## ' @param pvals a vector of values in [0, 1]. P-values
## ' @param s a vector of values in [0, 1]. Threshold curve
## ' @param dist an object of class "\code{\link{exp_family}}"
## ' @param pix a vector of values in [0, 1]. \eqn{\pi(x_{i})}{\pi(xi)}
## ' @param mux a vector of values. \eqn{\mu(x_{i})}{\mu(xi)}
## ' @return Imputed values, including
## ' \item{Hhat}{a vector of values in [0, 1]. Imputed values for \eqn{H_{i}}{Hi}'s}
## ' \item{phat}{a vector of values in [0, 1]. Imputed values for \eqn{p_{i}}{pi}'s}
Estep_mix <- function(pvals, s, dist, pix, mux){
	hp <- dist$h(pvals, mux)
	hp_mir <- dist$h(1 - pvals, mux)
	Hhat <- ifelse(
			pvals < s | pvals > 1 - s,
			1 / (1 + 2 * (1 - pix) / pix / (hp + hp_mir)),
			1 / (1 + (1 - pix) / pix / hp)
	)
	Hhat <- pminmax(Hhat, 1e-5, 1-1e-5)
	bhat <- ifelse(
			pvals < s | pvals > 1 - s,
			hp / (hp + hp_mir),
			1
	)
	
	if (any(is.na(Hhat))){
		stop("Hhat in the E-step has NAs.")
	}
	if (any(is.na(bhat))){
		stop("bhat in the E-step has NAs.")
	}
	
	return(list(Hhat = Hhat, bhat = bhat))
}
#===============================================================
# exp_family class
#===============================================================

#' Generate exp_family Objects for Exponential Families
#'
#' \code{exp_family} objects contain all required information in an exponential family to perform the E-step. The exponential function is encoded by
#' \deqn{h(p; \mu) = \exp\{(\eta(\mu) - \eta(\mu^{*})) g(p) - (A(\mu) - A(\mu^{*}))\}}{h(p; \eta) = exp{(\eta(\mu) - \eta(\mu*)) g(p) - (A(\mu) - A(\mu*))}}
#' where \eqn{g(p)} is an arbitrary transformation, \eqn{\mu} is the
#' \emph{mean parameter}, \eqn{\eta} is the natural parameter,
#' and \eqn{A(\mu)} is the partition function. The extra redundant
#' parameter \eqn{\mu^{*}}{\mu*} is to guarantee that \eqn{U([0, 1])}
#' belongs to the class.
#'
#' Beta family (\code{beta_family()}): modeling p-values as Beta-distributed random variables, i.e. \eqn{g(p) = -log(p)}, \eqn{\eta(\mu) = -1 / \mu}, \eqn{\mu* = 1}, \eqn{A(\mu) = log(\mu)}, name = "beta" and family = Gamma(). Beta-family is highly recommended for general problems and used as default.
#'
#' Inverse-gaussian family (\code{inv_gaussian_family()}): modeling p-values as transformed z-scores, i.e. \eqn{g(p) = \Phi^{-1}(p) (\Phi is the c.d.f. of a standard normal random variable)}, \eqn{\eta(\mu) = \mu}, \eqn{\mu* = 0}, \eqn{A(\mu) = \mu^2 / 2}, name = "inv_gaussian" and family = gaussian().
#'
#' @param g a function. An transformation of p-values
#' @param ginv a function. The inverse function of \code{g}
#' @param eta a function. The natural parameter as a function of the mean parameter \code{mu}
#' @param mustar a scalar. The mean parameter that gives \eqn{U([0, 1])}
#' @param A a function. The partition function
#' @param name a string. A name for the family. NULL by default
#' @param family an object of class "\code{\link[stats]{family}}" from \code{stats} package. The family used for model fitting in \code{\link[stats]{glm}}, \code{\link[mgcv]{gam}}, \code{\link[glmnet]{glmnet}}, etc..
#' @return an object of class "exp_family". This includes all inputs and  \code{h}, the density function.
#'
#' @export
#'
gen_exp_family <- function(g, ginv, eta, mustar, A,
		name = NULL, family = NULL){
	h <- function(p, mu){
		ifelse(mu == mustar,
				rep(1, length(p)),
				exp(
						g(p) * (eta(mu) - eta(mustar)) -
								(A(mu) - A(mustar))
				)
		)
	}
	result <- structure(
			list(h = h,
					g = g,
					ginv = ginv,
					eta = eta,
					mustar = mustar,
					A = A,
					name = name,
					family = family),
			class = "exp_family"
	)
	return(result)
}

#' @rdname gen_exp_family
#'
#' @export
#'
beta_family <- function(){
	g <- function(x){
		tmp <- -log(x)
		pmax(pmin(tmp, -log(10^-15)), -log(1-10^-15))
	}
	
	ginv <- function(x){
		tmp <- exp(-x)
		pmax(pmin(tmp, 1 - 10^-15), 10^-15)
	}
	eta <- function(mu){-1/mu}
	mustar <- 1
	A <- log
	name <- "beta"
	family <- Gamma()
	
	gen_exp_family(g, ginv, eta, mustar, A, name, family)
}

#' @rdname gen_exp_family
#'
#' @export
#'
inv_gaussian_family <- function(){
	g <- function(x){
		tmp <- qnorm(1 - x)
		pmax(pmin(tmp, qnorm(1 - 10^-15)), qnorm(10^-15))
	}
	ginv <- function(x){
		tmp <- 1 - pnorm(x)
		pmax(pmin(tmp, 1 - 10^-15), 10^-15)
	}
	eta <- function(mu){mu}
	mustar <- 0
	A <- function(mu){mu^2/2}
	name <- "inv_gaussian"
	family <- gaussian()
	
	gen_exp_family(g, ginv, eta, mustar, A, name, family)
}
#---------------------------------------------------------------
# Helpers for fitting, including M-steps and initialization
#---------------------------------------------------------------
fit_pi <- function(fun, args, type = c("Mstep", "init")){
	type <- type[1]
	if (type == "Mstep"){
		fun_name <- "\'pifun\'"
	} else if (type == "init"){
		fun_name <- "\'pifun_init\'"
	}
	
	if (is.null(args)) {
		stop(paste0(fun_name, " has irregular input types. Replace it with another function or writing a wrapper of ", fun_name, " with regular types of input (formula = , data = ) or (x = x, y = y) or (X = x, y = y)"))
	}
	
	fit <- try(do.call(fun, args))
	if (class(fit) == "try-error"){
		warning("Initialization of pi(x) by \'pifun_init\' fails. Initialize pi(x) as a constant 0.1 by default. Please check if too few p-values lying in [s0, 1-s0] and decrease s0 or change \'pifun_init\' if so. ")
	}
	if (is.null(fit$fitv)){
		stop(paste0(fun_name, " does not output \'fitv\'. Replace it with another function or change the name for fitted value to \'fitv\'"))
	}
	
	pix <- as.numeric(fit$fitv)
	if (any(is.na(pix))){
		stop("Initialization of \'pix\' has NAs")
	}
	pix <- pminmax(pix, 0, 1)
	
	return(list(fitv = pix, info = fit$info))
}

fit_mu <- function(fun, args, dist, type = c("Mstep", "init")){
	type <- type[1]
	if (type == "Mstep"){
		fun_name <- "\'pifun\'"
	} else if (type == "init"){
		fun_name <- "\'pifun_init\'"
	}
	
	if (is.null(args)) {
		stop(paste0(fun_name, " has irregular input types. Replace it with another function or writing a wrapper of ", fun_name, " with regular types of input (formula = , data = ) or (x = x, y = y) or (X = x, y = y)"))
	}
	
	fit <- do.call(fun, args)
	if (is.null(fit$fitv)){
		stop(paste0(fun_name, " does not output \'fitv\'. Replace it with another function or change the name for fitted value to \'fitv\'"))
	}
	
	mux <- as.numeric(fit$fitv)
	if (any(is.na(mux))){
		stop("Initialization of \'mux\' has NAs")
	}
	if (dist$family$family == "Gamma"){
		mux <- pmax(mux, 1)
	} else if (dist$family$family == "gaussian"){
		mux <- pmax(mux, 0)
	}
	
	return(list(fitv = mux, info = fit$info))
}
#---------------------------------------------------------------
# Helpers
#---------------------------------------------------------------

pminmax <- function(x, low, up){
	pmin(pmax(x, low), up)
}

logit <- function(x){
	log(x / (1 - x))
}

inv_logit <- function(x){
	exp(x) / (1 + exp(x))
}

func_input_type <- function(fun){
	argnames <- formalArgs(fun)
	if ("formula" %in% argnames){
		return("formula")
	} else if ("x" %in% argnames){
		return("xy")
	} else if ("X" %in% argnames){
		return("Xy")
	}
}

find_newname <- function(names_vec){
	name <- "aaa"
	while (name %in% names_vec){
		name <- paste0(name, "a")
	}
	return(name)
}

complete_pkg <- function(formula){
	formula <- as.character(formula)
	formula <- tail(formula, 1)
	formula <- tail(strsplit(formula, "~")[[1]], 1)    
	formula <- paste0(" ", formula)    
	if (grepl("ns\\(", formula)){
		if (!requireNamespace("splines", quietly = TRUE)){
			stop("package \'splines\' not found. Please install.")
		}
		formula <- gsub("([^:])ns\\(", "\\1splines::ns\\(", formula)
	}
	if (grepl("[^a-z]s\\(", formula)){
		if (!requireNamespace("mgcv", quietly = TRUE)){
			stop("package \'mgcv\' not found. Please install.")
		}
		formula <- gsub("([^:a-z])s\\(", "\\1mgcv::s\\(", formula)
	}
	return(formula)
}


complete_formula <- function(formula, response_name){
	if (is.null(formula)){
		stop("No formula is found. Please specify a formula ")
	}
	formula <- as.character(formula)
	formula <- tail(formula, 1)
	formula <- tail(strsplit(formula, "~")[[1]], 1)
	formula <- paste0(" ", formula)
	## completed_formula <- as.formula(
	##     paste(response_name, "~", formula),
	##     env = environment(args$formula))
	completed_formula <- paste0(response_name, " ~", formula)
	
	return(completed_formula)
}

complete_args <- function(x, response, fun,
		args = NULL,
		weights = NULL,
		force_integer = FALSE){
	input_type <- func_input_type(fun)
	if (!input_type %in% c("formula", "xy", "Xy")){
		stop("Wrong input type.")
	}
	
	response_name <- find_newname(colnames(x))
	
	if (input_type == "formula"){
		if (is.null(args) || !"formula" %in% names(args)){
			stop("Formula is not found. Please specify a formula for the fitting function.")
		}
		data <- cbind(data.frame(response), x)        
		colnames(data)[1] <- response_name
		args$formula <-  complete_formula(args$formula, response_name)
		data_args <- c(list(data = data), args)
	} else if (input_type == "xy"){
		data_args <- c(
				list(x = x, y = response),
				args)
	} else if (input_type == "Xy"){
		data_args <- c(
				list(X = x, y = response),
				args)
	} 
	
	data_args <- c(data_args, list(weights = weights))
	
	return(data_args)
}

complete_model <- function(model, dist){
	if (is.null(model$algo)){
		switch(model$name,
				"glm" = gen_adapt_model_glm(
						dist, model$args$piargs, model$args$muargs
				),
				"gam" = gen_adapt_model_gam(
						dist, model$args$piargs, model$args$muargs
				),
				"glmnet" = gen_adapt_model_glmnet(
						dist, model$args$piargs, model$args$muargs
				),
				stop("\'model$name\' not found in the library")
		)
	} else {
		model
	}
}
#===============================================================
# Compute M-step for the mixture model.
#===============================================================

Mstep_mix <- function(x, pvals, dist,
		Hhat, bhat,
		pifun, mufun,
		piargs = NULL, muargs = NULL,
		type = "unweighted"){
	if (!"weights" %in% formalArgs(mufun)){
		stop("'mufun' does not have input 'weights'")
	}
	
	n <- length(Hhat)
	x_aug <- rbind(x, x)
	H_aug <- c(rep(1, n), rep(0, n))
	weights <- c(Hhat, 1 - Hhat)
	piargs <- complete_args(x_aug, H_aug, pifun, piargs, weights)
	pi_res <- fit_pi(pifun, piargs, type = "Mstep")
	pi_res$fitv <- pi_res$fitv[1:n]
	
	y_aug <- c(dist$g(pvals), dist$g(1 - pvals))
	
	if (type == "weighted"){
		weights <- c(Hhat * bhat, Hhat * (1 - bhat))
	} else if (type == "unweighted"){
		weights <- c(bhat, 1 - bhat)
	}
	muargs <- complete_args(x_aug, y_aug, mufun, muargs, weights)
	mu_res <- fit_mu(mufun, muargs, dist, type = "Mstep")
	mu_res$fitv <- mu_res$fitv[1:n]
	
	res <- list(pix = pi_res$fitv,
			mux = mu_res$fitv,
			pi_info = pi_res$info,
			mu_info = mu_res$info)
	
	return(res)
}
#---------------------------------------------------------------
# Functions to fit models in adapt
#---------------------------------------------------------------

safe_glm <- function(formula, family, data, weights = NULL,
		...){
	options(warn = -1)
	
	formula <- as.formula(formula)
	if (family$link %in% c("inverse", "log")){
		fit <- try(glm(formula, family, data, weights, ...),
				silent = TRUE)
		if (class(fit)[1] == "try-error"){
			mod_mat <- model.matrix(formula, data = data)
			p <- ncol(mod_mat) - 1
			start <- c(1, rep(0, p))
			fit <- glm(formula, family, data, weights,
					start = start, ...)
		}
	} else {
		fit <- glm(formula, family, data, weights, ...)
	}
	
	fitv <- as.numeric(
			predict(fit, type = "response")
	)
	
	df <- fit$rank
	info <- list(df = df)
	
	options(warn = 0)
	
	return(list(fitv = fitv, info = info))
}

safe_gam <- function(formula, family, data, weights = NULL,
		...){
	options(warn = -1)
	
	formula <- as.formula(formula)
	if (family$link %in% c("inverse", "log")){
		fit <- try(mgcv::gam(formula, family, data, weights, ...),
				silent = TRUE)
		if (class(fit)[1] == "try-error"){
			mod_mat <- model.matrix(formula, data = data)
			p <- ncol(mod_mat) - 1
			start <- c(1, rep(0, p))
			fit <- mgcv::gam(formula, family, data, weights,
					start = start, ...)
		}
	} else {
		fit <- mgcv::gam(formula, family, data, weights, ...)
	}
	
	fitv <- as.numeric(
			predict(fit, type = "response")
	)
	
	df <- fit$rank
	info <- list(df = df)
	
	options(warn = 0)
	
	return(list(fitv = fitv, info = info))
}

safe_glmnet <- function(x, y, family, weights = NULL,
		...){
	options(warn = -1)
	
	if (class(family)[1] == "family"){
		family <- family$family
	}
	
	if (family %in% c("gaussian", "binomial", "poisson", "multinomial", "cox", "mgaussian")){
		if (is.null(weights)){
			fit <- glmnet::cv.glmnet(x, y, 
					family = family, ...)
		} else {
			weights <- pminmax(weights, 1e-5, 1-1e-5)
			fit <- glmnet::cv.glmnet(x, y, weights,
					family = family, ...)
		}
	} else if (family == "Gamma"){
		if (is.null(weights)){
			fit <- HDtweedie::cv.HDtweedie(x, y, p = 2,
					standardize = TRUE,
					...)
		} else {
			weights <- pminmax(weights, 1e-5, 1-1e-5)
			fit <- HDtweedie::cv.HDtweedie(x, y, p = 2,
					weights = weights,
					standardize = TRUE,
					...)
		}
	}
	
	fitv <- as.numeric(
			predict(fit, newx = x, s = "lambda.min",
					type = "response")
	)
	
	beta <- coef(fit, s = "lambda.min")
	vi <- as.numeric(beta != 0)[-1]    
	df <- sum(vi) + 1
	info <- list(df = df, vi = vi)
	
	options(warn = 0)
	
	return(list(fitv = fitv, info = info))
}

#---------------------------------------------------------------
# Adaptive P-Value Thresholding For Multiple Testing With Side Information
# Lihua Lei & William Fithian,
#   "AdaPT: An interactive procedure for multiple testing with side information"
# Available from http://arxiv.org/abs/1609.06035
#---------------------------------------------------------------

fdp_hat2 <- function(A, R){
	(0 + A) / pmax(1, R)
}


compute_lfdr_mix <- function(pvals, dist, params,
                             type = "over"){
    pix <- pmax(params$pix, 1e-5)
    mux <- params$mux
    if (dist$family$family == "Gamma"){
        mux <- pmax(mux, 1 + 1e-5)
    }
    type <- type[1]
    if (type == "over"){
        lfdr <- (pix * dist$h(1, mux) + 1 - pix) /
            (pix * dist$h(pvals, mux) + 1 - pix)
    } else if (type == "raw"){
        lfdr <- (1 - pix) / (pix * dist$h(pvals, mux) + 1 - pix)
    }
    return(lfdr)
}

compute_threshold_mix <- function(dist, params, lfdr_lev,
                                  type = "over"){
    pix <- pmax(params$pix, 1e-5)
    mux <- params$mux
    if (dist$family$family == "Gamma"){
        mux <- pmax(mux, 1 + 1e-5)
    }
    if (lfdr_lev == 0 || lfdr_lev == 1){
        return(rep(lfdr_lev, length(pix)))
    }

    type <- type[1]
    if (type == "over"){
        val1 <- dist$h(1, mux) / lfdr_lev +
            (1 - pix) / pix * (1 - lfdr_lev) / lfdr_lev
    } else if (type == "raw"){
        val1 <- (1 - pix) / pix * (1 - lfdr_lev) / lfdr_lev
    }

    val2 <- (log(val1) + dist$A(mux) - dist$A(dist$mustar)) /
        (dist$eta(mux) - dist$eta(dist$mustar))
    dist$ginv(val2)
}

create_stamps <- function(nmasks, nfits, nms){
    fit_stamps <- c(seq(0, nmasks, floor(nmasks / nfits)))[1:nfits]
    stamps <- data.frame(stamp = fit_stamps, type = "fit",
                         stringsAsFactors = FALSE)
    if (!is.null(nms)){
        ms_inds <- seq(1, nfits, floor(nfits / nms))[1:nms]
        stamps[ms_inds, 2] <- "ms"
    }
    stamps <- rbind(stamps, data.frame(stamp = nmasks, type = "end"))
    return(stamps)
}

check_pkgs <- function(models){
    if (class(models) == "adapt_model"){
        models <- list(models)
    }

    models_names <- sapply(models, function(model){model$name})

    if ("gam" %in% models_names){
        ind <- which(models_names == "gam")
        tmp <- sapply(models, function(model){is.null(model$algo)})
        if (any(tmp) && !requireNamespace("mgcv", quietly = TRUE)){
            stop("'mgcv' package not found. Please install.")
        }
    }
    if ("glmnet" %in% models_names){
        ind <- which(models_names == "glmnet")
        tmp <- sapply(models, function(model){is.null(model$algo)})
        if (any(tmp) && !requireNamespace("glmnet", quietly = TRUE)){
            stop("'glmnet' package not found. Please install.")
        }
    }
    return(invisible())
}

#' Adaptive P-value Thresholding
#'
#' \code{adapt} is a framework allowing for arbitrary exponential families for computing E-steps and arbitrary algorithms for fitting M-steps.
#'
#' \code{x} should have a type compatible to the fitting functions in \code{models}. For GLM and GAM, \code{x} should be a data.frame. For glmnet, \code{x} should be a matrix.
#'
#' \code{models} could either be an \code{adapt_model} object, if a single model is used, or a list of \code{adapt_model} objects, each of which corresponding to a model. Each element should be generated by \code{\link{gen_adapt_model}}. For glm/gam/glmnet, one can use the shortcut by running \code{\link{gen_adapt_model}} with name = "glm" or "gam" or "glmnet" but without specifying \code{pifun}, \code{mufun}, \code{pifun_init} and \code{mufun_init}. See examples below.
#'
#' \code{nfits} is the number of model fitting steps plus \code{nms}, the model selection steps, if \code{models} contains multiple \code{adapt_model} objects. Suppose M is the number of masked p-values at the initial step, then the model is updated at the initial step and at every time when [M/\code{nfits}] more p-values are revealed. If \code{nms > 0}, model selection is performed at the initial step an at every time when [M/\code{nms}] more p-values are revealed. Between two consecutive model selection steps, the model selected from the last step is used for model fitting. For example, when M = 10000, nfits = 10 and nms = 2, model selection will be performed at the initial step and when 5000 p-values are revealed, while the model fitting will be performed when 1000, 2000, 3000, 4000, 6000, 7000, 8000, 9000 p-values are revealed.
#'
#' \code{verbose} has three elements: \code{print}, \code{fit} and \code{ms}. If \code{print = TRUE}, the progress of the main procedure is outputted to the console, in the form of "alpha = 0.05: FDPhat 0.0333, Number of Rej. 30" (where the numbers are made up for illustration). If \code{fit = TRUE}, a progress bar for the model fitting is outputted to the console. Similarly, if \code{ms = TRUE}, a progress bar for the model selection is outputted to the console.
#'
#' For ultra-large scale problems (n > 10^5), it is recommended to keep \code{alphas} short because the output \code{s} is of size n x \code{length(alphas)}.
#' is \code{length(alphas)}.
#'
#' The output \code{qvals} gives the q-values of each hypothesis. \code{qvals[i]} is defined as the minimum target FDR level such that \code{pvals[i]} is rejected. For hypotheses with p-values above s0, the q-values are set to be Inf because they are never rejected by AdaPT for whatever alpha.
#'
#' The output \code{order} gives the order of (the indices of) p-values being revealed, i.e. being in the region (s, 1-s). The latter hypotheses appeared in \code{order} have smaller q-values (i.e. are more likely to be rejected).
#'
#' @param x covariates (i.e. side-information). Should be compatible to \code{models}. See Details
#' @param pvals a vector of values in [0, 1]. P-values
#' @param models an object of class "\code{adapt_model}" or a list of objects of class "adapt_model". See Details
#' @param dist an object of class "\code{\link{gen_exp_family}}". \code{\link{beta_family}()} as default
#' @param s0 a vector of values in [0, 0.5). Initial threshold.
#' @param alphas a vector of values in (0, 1). Target FDR levels.
#' @param params0 a list in the form of list(pix = , mux = ). Initial guess of pi(x) and mu(x). NULL as default
#' @param nfits a positive integer. Number of model-fitting steps. See Details
#' @param nms a non-negative integer. Number of model selection steps. See Details
#' @param niter_fit a positive integer. Number of EM iterations in model fitting
#' @param tol a positive scalar. EM algorithm stops when pi(x) and mu(x) in consecutive steps differ by at most 'tol' for each element
#' @param niter_ms a positive integer. Number of EM iterations in model selection
#' @param cr a string. The criterion for model selection with BIC as default. Also support AIC, AICC and HIC
#' @param verbose a list of logical values in the form of list(print = , fit = , ms = ). Each element indicates whether the relevant information is outputted to the console. See Details
#'
#' @return
#' \item{nrejs}{a vector of integers. Number of rejections for each alpha}
#' \item{rejs}{a list of vector of integers. The set of indices of rejections for each alpha}
#' \item{s}{a matrix of size \code{length(pvals) X length(alphas)}. Threshold curves for each alpha}
#' \item{params}{a list. Each element is a list in the form of \code{list(pix = , mux = , alpha = , nmasks =)}, recording the parameter estimates, the achieved alpha and the number of masked p-values. To avoid massive storage cost, it only contains the information when a new target FDR level is achieved. As a result, it might be  shorter than \code{nfits}.}
#' \item{qvals}{a vector of values in [0, 1]U{Inf}. Q-values. See Details}
#' \item{order}{a permutation of \code{1:length(pvals)}. Indices of hypotheses arranged in the order of reveal. See Details}
#' \item{alphas}{same as the input \code{alphas}}
#' \item{dist}{same as the input \code{dist}}
#' \item{models}{a list of \code{adapt_model} objects of length \code{params}. The model used in each fitting step. As in \code{params}, it only contains the model when a new target FDR level is achieved and each element corresponds to an element of \code{params}.}
#' \item{info}{a list of length \code{nfits}. Each element is a list recording extra information in each fitting step, e.g. degree of freedom (df) and variable importance (vi). As in \code{params}, it only contains the model information when a new target FDR level is achieved and each element corresponds to an element of \code{params}.}
#' \item{args}{a list including the other inputs \code{nfits}, \code{nms}, \code{niter_fit}, \code{niter_ms}, \code{tol}, \code{cr}}.
#'
#' @examples
#' \donttest{
#' # Load estrogen data
#' data(estrogen)
#' pvals <- as.numeric(estrogen$pvals)
#' x <- data.frame(x = as.numeric(estrogen$ord_high))
#' dist <- beta_family()
#'
#' # Subsample the data for convenience
#' inds <- (x$x <= 5000)
#' pvals <- pvals[inds]
#' x <- x[inds,,drop = FALSE]
#'
#' # Generate models for function adapt
#' library("splines")
#' formulas <- paste0("ns(x, df = ", 6:10, ")")
#' models <- lapply(formulas, function(formula){
#'     piargs <- muargs <- list(formula = formula)
#'     gen_adapt_model(name = "glm", piargs = piargs, muargs = muargs)
#' })
#'
#' # Run adapt
#' res <- adapt(x = x, pvals = pvals, models = models,
#'              dist = dist, nfits = 10)
#' }
#'
#' @export
adapt2 <- function(x, pvals, models,
                  dist = beta_family(),
                  s0 = rep(0.45, length(pvals)),
                  alphas = seq(0.01, 1, 0.01),
                  params0 = list(pix = NULL, mux = NULL),
                  nfits = 20, nms = 1,
                  niter_fit = 10, tol = 1e-4,
                  niter_ms = 20, cr = "BIC",
                  verbose = list(print = TRUE, fit = FALSE, ms = TRUE)
                  ){

    Mstep_type <- "unweighted"
    lfdr_type <- "over"

    ## Check if 'pvals' is a vector of values in [0, 1]
    if (!is.numeric(pvals) || min(pvals) < 0 || max(pvals) > 1){
        stop("Invalid p-values")
    }

    ## Check if the size of 'x' matches that of 'pvals'
    if (nrow(x) != length(pvals)){
        stop("'x' must have the same rows as the length of 'pvals'")
    }

    ## Check if 'dist' is of class 'exp_family'
    if (class(dist)[1] != "exp_family"){
        stop("\'dist\' must be of class \'exp_family\'.")
    }

    ## Check if necessary packages are installed.
    check_pkgs(models)

    ## When a single model is provided, set 'nms' to be NULL
    if (class(models) == "adapt_model"){
        model <- models
        nms <- NULL
    } else if (is.list(models)){
        types <- sapply(models, class)
        if (any(types != "adapt_model")){
            stop("All models should be of class \'adapt_model\'.")
        }
        if (!is.integer(nms) || nms <= 0){
            nms <- 1
        } else if (nms > nfits){
            warning("Model selection cannot be more frequent than model fitting. Set \'nms\' to \'nfits\'")
            nms <- nfits
        }
        if (length(models) == 1){
            model <- models[[1]]
            nms <- NULL
        }
    }

    ## Create time stamps when model is fitted or model selection is performed
    nmasks <- sum(pvals <= s0) + sum(pvals >= 1 - s0)
    stamps <- create_stamps(nmasks, nfits, nms)

    ## Create root arguments to simplify fitting and model selection
    fit_args_root <- list(
        x = x, pvals = pvals, dist = dist,
        niter = niter_fit, tol = tol,
        verbose = verbose$fit, type = Mstep_type
        )
    if (any(stamps[, 2] == "ms")){
        ms_args_root <- list(
            x = x, pvals = pvals, dist = dist, models = models,
            cr = cr, niter = niter_ms, tol = tol,
            verbose = verbose$ms, type = Mstep_type
            )
    }

    ## Initialization
    n <- length(pvals)
    params <- params0
    s <- s0
    A <- sum(pvals >= 1 - s)
    R <- sum(pvals <= s)
    minfdp <- fdp_hat2(A, R) # initial FDPhat

    ## Remove the alphas greater than the initial FDPhat, except the smallest one among them
    alphas <- sort(alphas)
    if (min(alphas) >= minfdp){
        warning("Task completed! Initial \'s0\' has guaranteed FDR control for all alpha's in \'alphas\'.")
        alphaind <- 0
    } else if (max(alphas) < minfdp){
        alphaind <- length(alphas)
    } else {
        alphaind <- max(which(alphas <= minfdp))
    }

    m <- length(alphas)
    nrejs_return <- rep(0, m) # number of rejections
    s_return <- matrix(0, n, m) # threshold
    params_return <- list() # parameters (including pix and mux)
    model_list <- list() # all selected models
    info_list <- list() # other information (df, vi, etc.)
    reveal_order <- which((pvals > s) & (pvals < 1 - s)) # the order to be revealed
    if (length(reveal_order) > 0){
        init_pvals <- pvals[reveal_order]
        reveal_order <- reveal_order[order(init_pvals, decreasing = TRUE)]
    }
    fdp_return <- c(rep(Inf, length(reveal_order)), minfdp) # fdphat along the whole path

    if (m > alphaind){
        nrejs_return[(alphaind + 1):m] <- R
        s_return[, (alphaind + 1):m] <- s0
    }
    ## alphas <- alphas[1:m]

    for (i in 1:(nrow(stamps) - 1)){
        if (alphaind == 0){
            if (verbose$print){
                cat("Task completed!")
            }
            break
        }

        alpha <- alphas[alphaind]
        # mask <- (pvals <= s) | (pvals >= 1 - s)
        mask <- rep(TRUE, n)
        mask[reveal_order] <- FALSE
        nmasks <- sum(mask)
        A <- sum(pvals >= 1 - s)
        R <- sum(pvals <= s)
        start <- stamps[i, 1]
        end <- stamps[i + 1, 1]
        nreveals <- min(nmasks, end - start) # number of hypotheses to be revealed
        type <- stamps[i, 2] # type of model update

        ## Model selection or model fitting
        if (type == "ms"){
            ms_args <- c(
                list(s = s, params0 = params),
                ms_args_root
                )
            ## Use "EM_mix_ms" from "EM-mix-ms.R"
            ms_res <- do.call(EM_mix_ms, ms_args)
            params <- ms_res$params
            model <- ms_res$model
            modinfo <- ms_res$info
        } else if (type == "fit"){
            fit_args <- c(
                list(s = s, params0 = params, model = model),
                fit_args_root
                )
            ## Use "EM_mix" from "EM-mix.R"
            fit_res <- do.call(EM_mix, fit_args)
            params <- fit_res$params
            modinfo <- fit_res$info
        }

        if (length(params_return) == 0 ||
            alpha < params_return[[length(params_return)]]$alpha){
            params <- c(params,
                        list(alpha = alpha, nmasks = nmasks))
            params_return <- append(params_return, list(params))
            model_list <- append(model_list, list(model))
            info_list <- append(info_list, modinfo)
        }

        ## Estimate local FDR
        lfdr <- compute_lfdr_mix(
            pmin(pvals, 1 - pvals),
            dist, params, lfdr_type)
        ## Find the top "nreveals" hypotheses with highest lfdr
        lfdr[!mask] <- -Inf
        inds <- order(lfdr, decreasing = TRUE)[1:nreveals]
        reveal_order <- c(reveal_order, inds)
        ## Shortcut to calculate FDPhat after revealing the hypotheses one by one
        Adecre <- cumsum(pvals[inds] >= 1 - s[inds])
        Rdecre <- cumsum(pvals[inds] <= s[inds])
        fdp <- fdp_hat2(A - Adecre, R - Rdecre)
        fdp_return <- c(fdp_return, fdp)
        fdp <- pmin(fdp, minfdp)
        ## Calculate the current minimum FDPhat
        minfdp <- min(fdp)

        while (alphaind > 0){# check if lower FDR level is achieved
            alpha <- alphas[alphaind]
            if (any(fdp <= alpha)){
                breakpoint <- which(fdp <= alpha)[1]
                lfdr_lev <- lfdr[inds[breakpoint]]
                snew <- compute_threshold_mix(dist, params, lfdr_lev, lfdr_type)
                snew <- pmin(s, snew)

                ## Sanity check to avoid rounding errors
                tmp_pvals <- pvals[inds[1:breakpoint]]
                tmp_pvals <- pmin(tmp_pvals, 1 - tmp_pvals)
                tmp_inds <- which(tmp_pvals <= snew[inds[1:breakpoint]])
                if (length(tmp_inds) > 0){
                    snew[inds[tmp_inds]] <- pmin(snew[inds[tmp_inds]], tmp_pvals[tmp_inds] - 1e-15)
                }

                s_return[, alphaind] <- snew

                fdpnew <- fdp[breakpoint]
                Rnew <- sum(pvals <= snew)
                nrejs_return[alphaind] <- Rnew
                if (verbose$print){
                    cat(paste0(
                        "alpha = ", alpha, ": FDPhat ",
                        round(fdpnew, 4), ", Number of Rej. ",
                        Rnew, "\n"))
                }

                alphaind <- alphaind - 1
            } else {
                break
            }
        }

        if (alphaind == 0){ # check again to save computation
            if (verbose$print){
                cat("Task completed!\n")
            }
            break
        }

        ## Update s(x)
        final_lfdr_lev <- lfdr[tail(inds, 1)]
        snew <- compute_threshold_mix(dist, params, final_lfdr_lev, lfdr_type)
        s <- pmin(s, snew)

        ## Sanity check to avoid rounding errors
        tmp_pvals <- pvals[inds]
        tmp_pvals <- pmin(tmp_pvals, 1 - tmp_pvals)
        tmp_inds <- which(tmp_pvals <= snew[inds])
        if (length(tmp_inds) > 0){
            s[inds[tmp_inds]] <- pmin(s[inds[tmp_inds]], tmp_pvals[tmp_inds] - 1e-15)
        }
    }

    remain_inds <- (1:n)[-reveal_order]
    if (length(remain_inds) > 0){
        tmp_pvals <- pvals[remain_inds]
        tmp_pvals <- pmin(tmp_pvals, 1 - tmp_pvals)
        remain_reveal_order <- remain_inds[order(tmp_pvals, decreasing = TRUE)]
        reveal_order <- c(reveal_order, remain_reveal_order)
        fdp_return <- c(fdp_return, rep(minfdp, length(remain_inds)))
    }

    rejs_return <- apply(s_return, 2, function(s){which(pvals <= s)})

    qvals <- rep(1, n)
    qvals[reveal_order] <- ifelse(pvals[reveal_order] <= s0[reveal_order],
                                  cummin(fdp_return[1:n]), rep(Inf, n))

    args <- list(nfits = nfits, nms = nms,
                 niter_fit = niter_fit, niter_ms = niter_ms,
                 tol = tol, cr = cr)

    res <- structure(
        list(nrejs = nrejs_return,
             rejs = rejs_return,
             s = s_return,
             params = params_return,
             qvals = qvals,
             order = reveal_order,
             alphas = alphas,
             dist = dist,
             models = model_list,
             info = info_list,
             args = args),
        class = "adapt")
    return(res)
}
