#' Set the control parameters for OrderShapeEM
#' 
#' The function provides an interface to set the parameters for OrderShapeEM. 
#' 
#' @param maxIter an integer value indicating the maximum number of iterations. Default is 250.
#' @param tol a numeric value giving the tolerance in the relative change in the log likelihood below which the algorithm is considered to be converged.
#' Default is 1e-3.
#' @param trace a logical value indicating whether to print out the EM process. Default is FALSE.
#' @param  pi0.init a scalar giving the initial estimates of pi0. Default is 0.95.
#' @param  k.init a scalar giving the initial estimates of f1. f1 is approximated by a beta distribution with the shape parameter k. 
#' Default is 0.75.
#' @param pvals.cutoff  a numeric value to replace p-values below that value, which is used to increase the stability of the algorithm. Default is 1e-15.
#' @return A list with the set control parameters.
#' 
#' @rdname OrderShapeEM.control
#' @export

OrderShapeEM.control <- function (maxIter = 250, tol = 1e-3,  trace = FALSE,
		pi0.init = 0.95, k.init = 0.75, pvals.cutoff = 1e-15) {
	# Args:
	#   start: a vector of starting values for pi0 (null proportion), mean, sd under H1
	rval <- list(maxIter = maxIter, trace = trace, tol = tol,  pi0.init = pi0.init, 
			k.init = k.init, pvals.cutoff = pvals.cutoff)
	rval
}


#' False discovery rate control using auxiliary (order) information 
#'
#' The function implements a scalable and tuning parameter-free FDR control procedure for large-scale multiple testing exploiting
#' the auxiliary information that reflects the order of the prior null probability. Both the prior null probabilities and the alternative
#' p-value distribution are estimated using isotonic regression (pool-adjacent-violators algorithm).
#'
#' @param pvals a numeric vector of the p-values.
#' @param order.var a numeric vector of covariate values reflective of the order of the prior null probability of the hypotheses.
#' @param control  a list of control arguments for the EM algorithm
#' \itemize{
#' \item{maxIter}{an integer value indicating the maximum number of iterations.}
#' \item{tol}{a numeric value giving the tolerance in the relative change in the log likelihood below which the algorithm is considered to be converged.}
#' \item{trace}{a logical value indicating whether to print out the EM process.}
#' \item{pi0.init, k.init}{two scalars giving the initial estimates of pi0 and f1. f1 is approximated by a beta distribution with the shape parameter k. 
#' Default are 0.95 and 0.75 for \code{pi0.init}, \code{k.init}, respectively.}
#' \item{pvals.cutoff}{a numeric value to replace p-values below that value, which is used to increase the stability of the algorithm.}
#' }
#'
#' @return A list with the elements
#' \item{call}{the call made.}
#' \item{pi0.unadj}{a numeric vector of the estimated null probabilities before calibration.}
#' \item{pi0}{a numeric vector of the estimated null probabilities after calibration.}
#' \item{delta}{a numeric value of the amount of calibration.}
#' \item{lfdr}{a numeric vector of the local fdrs.}
#' \item{fdr}{a numeric vector of the adjusted p-values.}
#' \item{f1}{a vector of the estimated densities for the p-value under the alternative.}
#' \item{pi0.step}{an object of class \code{stepfun} for the pi0.}
#' \item{f1.step}{an object of class \code{stepfun} for the f1.}
#' \item{convergence}{a list containing the convergence information. \code{code}: 1 - converged, 0 - not converged;
#' \code{iter}: the number of iterations performed.}
#' \item{loglik}{a numeric value for the log likelihood.}
#' 
#' @author Jun Chen, Xianyang Zhang
#' @references Hongyuan Cao, Jun Chen, Xianyang Zhang. Optimal false discovery rate control for large-scale multiple testing with
#'  auxiliary information. Submitted.
#' @keywords FDR
#' @importFrom stats stepfun knots sd
#' @importFrom Iso pava
#' @importFrom qvalue qvalue pi0est
#' @examples
#'
#' set.seed(123)
#' data.obj <- SimulateData(prior.strength = 'Moderate', sig.density = 'Low', sig.strength = 'Weak',
#' feature.no = 5000)
#' orderfdr.obj <- OrderShapeEM(data.obj$pvalue, data.obj$prior, OrderShapeEM.control(trace = TRUE))
#' 
#' # Plot the estimated pi0 and f1
#' par(mfrow = c(1, 2))
#' plot(orderfdr.obj$pi0.step, xlab = 'index', ylab = 'pi0', do.points = FALSE,
#'		xlim = c(1, length(data.obj$pvalue)), main = 'Null probability')
#' plot(orderfdr.obj$f1.step, xlab = 'p value', ylab = 'f1', do.points = FALSE, 
#'		xlim = c(min(data.obj$pvalue), 1), log = 'x', main = 'Alternative distribution')
#' 
#' # Calculate the number of true positives and the false discovery proportion
#' sum(orderfdr.obj$fdr <= 0.05 & data.obj$truth)
#' sum(orderfdr.obj$fdr <= 0.05 & !data.obj$truth) / max(sum(orderfdr.obj$fdr <= 0.05), 1)
#' 
#' # Compare to the BH procedure
#' sum(p.adjust(data.obj$p.value, 'fdr') <= 0.05 & data.obj$truth)
#' sum(p.adjust(data.obj$p.value, 'fdr') <= 0.05 & !data.obj$truth) / 
#' max(sum(p.adjust(data.obj$p.value, 'fdr') <= 0.05), 1)
#' 
#' @rdname OrderShapeEM
#' @export

OrderShapeEM <- function (pvals, order.var,  control = OrderShapeEM.control()) {
	
	maxIter <- control$maxIter
	tol <- control$tol
	trace <- control$trace
	pvals.cutoff <- control$pvals.cutoff
	k.init <- control$k.init
	pi0.init <- control$pi0.init
	
	# Global pi0 estimate
	pi0.est <- max(c(pi0est(pvals, pi0.method = 'bootstrap')$pi0,  
					pi0est(pvals, pi0.method = 'smoother')$pi0))
	
	
	pvals[pvals < pvals.cutoff] <- pvals.cutoff
	out <- sort(pvals, index.return=TRUE)
	pvals <- out$x
	index0 <- out$ix
	
	order.var <- order.var[index0]
	index <- sort(order.var, index.return = TRUE)$ix
	
	pvals.diff <- c(pvals[1], diff(pvals))   
	m <- length(pvals)
	
	# Intialization
	pi0 <- rep(pi0.init, m)
	pi1 <- 1 - pi0	
	f0 <- rep(1, m)
	f1 <- (1 - k.init) * pvals  ^ (-k.init)
	
	iter <- 0
	loglik0 <- 100
	converge <- list()
	
	if (trace) cat('Begin iteration ...\n')
	
	while (TRUE) {
		iter <- iter + 1
		
		# E-step
		f01 <- pi1 * f1 + pi0 * f0
		q0 <- pi0 * f0 / f01
		q1 <- pi1 * f1 / f01
		
		# M-step
		# Update the f1 distribution
		y <- - pvals.diff * sum(q1, na.rm = TRUE) / q1
		
		y[is.na(y)] <- -Inf
		y[y == 0] <- max(y[y != 0])
		y[y == -Inf] <- min(y[y != -Inf])
		
		out <- pava(y, q1, decreasing = TRUE, long.out = FALSE, stepfun = FALSE)
		
		f1 <- - 1 / out
		f1 <-  f1 / sum(f1 * pvals.diff)
		
		# Update the pi0s
		pi0 <- pava(q0[index], decreasing = FALSE, long.out = FALSE, stepfun = FALSE)
		pi0[index] <- pi0
		pi1 <- 1 - pi0
		
		loglik <- sum(log(f01))
		loglik.delta <- abs((loglik - loglik0) / loglik0)
		loglik0 <- loglik
		
		if (trace) cat(loglik, '\n')
		if (loglik.delta <= tol & iter <= maxIter){
			if (trace) cat('Converged!\n')
			converge$code <- 1
			converge$iter <- iter
			break
		} 
		if (iter > maxIter) {
			if (trace) cat('Maximum iteration reached!\n')
			converge$code <- 0
			converge$iter <- iter
			break
		}
	}
	
	stepfun1 <- pava(q0[index], decreasing = FALSE, long.out = FALSE, stepfun = TRUE)
	stepfun2 <- pava(y, q1, decreasing = TRUE, long.out = FALSE, stepfun = TRUE)
	stepfun2 <- stepfun(x = pvals[knots(stepfun2 )],  y = c(f1[knots(stepfun2)[1] - 1], f1[knots(stepfun2)]))
	
	
	# pi0 calibration
	delta <- (pi0.est - mean(pi0)) / mean(1 - pi0)
	delta <- ifelse(delta < 0 | is.nan(delta), 0, delta)
	
	
	pi0.unadj <- pi0
	pi0 <- pi0 + delta * (1 - pi0)
	pi1 <- 1 - pi0
	
	
	stepfun1 <- stepfun(x = knots(stepfun1), y = c(min(pi0), sort(pi0)[knots(stepfun1)]))
	
	f01 <- pi1 * f1 + pi0 * f0
	q0 <- pi0 * f0 / f01
	q1 <- pi1 * f1 / f01
	
	index0 <- order(index0)
	lfdr <- q0[index0]
	pi0 <- pi0[index0] 
	
	
	out <- sort(lfdr, index.return = TRUE)
	fdr <- cumsum(out$x) / (1:m)
	fdr <- fdr[order(out$ix)]
	
	# don't allow negative delta
	f1 <- f1[index0]
	pi0.unadj <- pi0.unadj[index0]
	
	
	
	cat('Finished!\n')
	
	# An empirical rule
	if (sd(pi0.unadj) <= 0.025) {
		warning("The order information seems weak. Our procedure could be slightly less powerful than the Storey's procedure.\n")
	}
	
	return(list(call = match.call(), lfdr = lfdr, fdr = fdr,  pi0.unadj = pi0.unadj, pi0 = pi0, f1 = f1, pi0.global = pi0.est,
					pi0.step = stepfun1, f1.step = stepfun2, loglik = loglik,  delta = delta, converge = converge))
	
}


#' Simulate p-values and the auxiliary covariate under various scenarios.
#'
#' The function simulates p-values and the auxiliary covariate under different signal structures (density and strength) and covariate informativeness. 
#'
#' @param prior.strength a character string from \code{'Weak', 'Moderate', 'Strong'} indicating the covariate informativeness.
#' @param feature.no an integer, the number of features to be simulated.
#' @param sig.dist  a character string from \code{'Normal', 'Gamma'} indicating the distribution of the z-value under the alternative.
#' @param sig.density  a character string from \code{'None', 'Lower', 'Low', 'Medium', 'High'} indicating the level of the signal density.
#' @param sig.strength a character string from \code{'Weak', 'Moderate', 'Strong'} indicating the level of the signal strength.
#' @return A list with the elements
#' \item{pvalue}{a numeric vector of p-values.}
#' \item{prior}{a vector of covariate values reflecting the order of the prior null probabilities.}
#' \item{truth}{a vector of logical values indicating H0 (=0) or H1 (=1).}
#' 
#' @author Jun Chen, Xianyang Zhang
#' @references Hongyuan Cao, Jun Chen, Xianyang Zhang. Optimal false discovery rate control for large-scale multiple testing with
#'  auxiliary information. Submitted.
#' @keywords Simulation
#' @importFrom stats rnorm runif rbeta rgamma pnorm rbinom
#' @importFrom Iso pava
#' @importFrom qvalue qvalue pi0est
#' @rdname SimulateData
#' @export
SimulateData <- function (
		prior.strength = c('Weak', 'Moderate', 'Strong'), 
		feature.no     = 10000,
		sig.dist       = c('Normal', 'Gamma'), 
		sig.density    = c('None', 'Lower', 'Low', 'Medium', 'High'), 
		sig.strength   = c('Weak', 'Moderate', 'Strong')) {
	
	prior.strength <- match.arg(prior.strength)
	sig.dist <- match.arg(sig.dist)
	sig.density <- match.arg(sig.density)
	sig.strength <- match.arg(sig.strength)
	
	if (sig.density == 'None') {
		pi <- rnorm(feature.no)
		pvalue <- runif(feature.no)
		truth <- rep(0, feature.no)
		return(list(pvalue = pvalue, prior = pi, truth = truth))
	}
	
	ks <- c(High=0.8, Medium=0.9, Low=0.95, Lower=0.99)
	fs <- c(High=4, Medium=9, Low=19, Lower=99)
	
	if (prior.strength == 'Weak') {
		pi <- rnorm(feature.no, mean = ks[sig.density], sd = 0.005)
		pi[pi > 1] <- 1
		pi[pi < 0] <- 0
		truth <- rbinom(feature.no, prob = 1 - pi, size = 1)
	}
	
	b <- 0.5
	if (prior.strength == 'Moderate') {
		pi <- rbeta(feature.no, shape1 = b * fs[sig.density], shape2 = b)
		pi[pi > 1] <- 1
		pi[pi < 0] <- 0
		truth <- rbinom(feature.no, prob = 1 - pi, size = 1)
	}
	
	if (prior.strength == 'Strong') {
		pi <- c(rnorm(round(feature.no * (1 - ks[sig.density])), mean = 0.2, sd = 0.05), 
				rnorm(feature.no * ks[sig.density], mean = 0.95, sd = 0.005))
		pi[pi > 1] <- 1
		pi[pi < 0] <- 0
		truth <- rbinom(feature.no, prob = 1 - pi, size = 1)
	}
	
	ns <- sum(truth)
	# Generate the signals
	stats <- rnorm(feature.no)
	if (sig.dist == 'Normal') {
		if (sig.strength == 'Strong') stats[truth == 1] <- rnorm(ns, 3, 1)
		if (sig.strength == 'Moderate') stats[truth == 1] <- rnorm(ns, 2.5, 1)
		if (sig.strength == 'Weak') stats[truth == 1] <- rnorm(ns, 2, 1)
	}
	
	if (sig.dist == 'Gamma') {
		if (sig.strength == 'Strong') stats[truth == 1] <- rgamma(ns, shape = 2, scale = 1 / sqrt(2)) - sqrt(2) + 3
		if (sig.strength == 'Moderate') stats[truth == 1] <- rgamma(ns, shape = 2, scale = 1 / sqrt(2)) - sqrt(2) + 2.5
		if (sig.strength == 'Weak') stats[truth == 1] <- rgamma(ns, shape = 2, scale = 1 / sqrt(2)) - sqrt(2) + 2
	}
	
	# Generate p values - one sided
	pvals <- (1 - pnorm(stats))
	
	return(list(pvalue = pvals, prior = pi, truth = truth))
}



