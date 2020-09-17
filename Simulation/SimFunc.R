################################################################################
#																			   #
#         Simulation functions for all scenarios in the manuscript             #
#                                                                              #
################################################################################

# The following simulation function is for figure 2, 6, 9, 10, 13, 14, 15, 16
SimulateData4 <- function (
		prior.strength = c('Weak', 'Moderate', 'Strong'), 
		feature.no     = 10000,
		sig.dist       = c('Normal', 'Gamma'), 
		sig.density    = c('None', 'Lower', 'Low', 'Medium', 'High'), 
		sig.strength   = c('Weak', 'Moderate', 'Strong')) {
	
	
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
		pi <- c(rnorm(round(feature.no * (1 - ks[sig.density])), mean = 0.2, sd = 0.05), rnorm(feature.no * ks[sig.density], mean = 0.95, sd = 0.005))
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
	pvalue <- (1 - pnorm(stats))
	
	return(list(pvalue = pvalue, prior = pi, truth = truth))
}



# The following simulation function is for figure 7 - noisy order
SimulateData5 <- function (
		prior.strength = c('Weak', 'Moderate', 'Strong'), 
		feature.no     = 10000,
		sig.dist       = c('Normal', 'Gamma'), 
		sig.density    = c('Low', 'Medium', 'High'), 
		sig.strength   = c('Weak', 'Moderate', 'Strong')) {
	
	
	ks <- c(High=0.8, Medium=0.9, Low=0.95)
	fs <- c(High=4, Medium=9, Low=19)
	
	pi <- c(rnorm(round(feature.no * (1 - ks[sig.density])), mean = 0.2, sd = 0.05), rnorm(feature.no * ks[sig.density], mean = 0.95, sd = 0.005))
	pi[pi > 1] <- 1
	pi[pi < 0] <- 0
	truth <- rbinom(feature.no, prob = 1 - pi, size = 1)
	
	if (prior.strength == 'Weak') {
		pi <- sample(pi)
	}
	
	
	if (prior.strength == 'Moderate') {
		ind <- sample(1:feature.no, feature.no / 2)
		pi[ind] <- sample(pi[ind])
	}
	
	if (prior.strength == 'Strong') {
		pi <- pi
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
	pvalue <- (1 - pnorm(stats))
	
	return(list(pvalue = pvalue, prior = pi, truth = truth))
}


# The following simulation function is for figure 3 - correlated hypotheses
SimulateData6 <- function (
		prior.strength = c('Weak', 'Moderate', 'Strong'), 
		feature.no     = 10000,
		sig.dist       = c('PosCor', 'PosNegCor'), 
		sig.density    = c('Low', 'Medium', 'High'), 
		sig.strength   = c('Weak', 'Moderate', 'Strong'),
		nb             = 100,
		rho            = 0.5) {
	
	# Prepare transformation matrix
	bs <- feature.no / nb 
	
	mat <- diag(bs)
	mat[, ] <- rho
	diag(mat) <- 1
	obj <- eigen(mat)
	T1 <- obj$vectors %*% diag(sqrt(obj$values)) %*% t(obj$vectors)
	
	mat <- diag(bs)
	mat[, ] <- -rho
	mat[1:(bs/2), 1:(bs/2)] <- rho
	mat[(bs/2+1):bs, (bs/2+1):bs] <- rho
	diag(mat) <- 1
	obj <- eigen(mat)
	T2 <- obj$vectors %*% diag(sqrt(obj$values)) %*% t(obj$vectors)
	
	
	ks <- c(High=0.8, Medium=0.9, Low=0.95)
	fs <- c(High=4, Medium=9, Low=19)
	
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
		pi <- c(rnorm(round(feature.no * (1 - ks[sig.density])), mean = 0.2, sd = 0.05), rnorm(feature.no * ks[sig.density], mean = 0.95, sd = 0.005))
		pi[pi > 1] <- 1
		pi[pi < 0] <- 0
		truth <- rbinom(feature.no, prob = 1 - pi, size = 1)
	}
	
	ns <- sum(truth)
	
	
	# Generate the signals
	#stats <- rnorm(feature.no)
	if (sig.dist == 'PosCor') {
		stats <- as.vector(T1 %*% matrix(rnorm(feature.no), bs, nb))
		if (sig.strength == 'Strong') stats[truth == 1] <- stats[truth == 1] + 3
		if (sig.strength == 'Moderate') stats[truth == 1] <- stats[truth == 1] + 2.5 
		if (sig.strength == 'Weak') stats[truth == 1] <- stats[truth == 1] + 2
	}
	
	if (sig.dist == 'PosNegCor') {
		stats <- as.vector(T2 %*% matrix(rnorm(feature.no), bs, nb))
		if (sig.strength == 'Strong') stats[truth == 1] <- stats[truth == 1] + 3
		if (sig.strength == 'Moderate') stats[truth == 1] <- stats[truth == 1] + 2.5 
		if (sig.strength == 'Weak') stats[truth == 1] <- stats[truth == 1] + 2
	}
	
	# Generate p values - one sided
	pvalue <- (1 - pnorm(stats))
	
	return(list(pvalue = pvalue, prior = pi, truth = truth))
}

# The following simulation is for figure 11 - varying f1
SimulateData8.4 <- function (
		prior.strength = c('Weak', 'Moderate', 'Strong'), 
		feature.no     = 10000,
		sig.dist       = c('Normal'), 
		sig.density    = c('Lower', 'Low', 'Medium', 'High'), 
		sig.strength   = c('Weak', 'Moderate', 'Strong')) {
	
	
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
		pi <- c(rnorm(round(feature.no * (1 - ks[sig.density])), mean = 0.2, sd = 0.05), rnorm(feature.no * ks[sig.density], mean = 0.95, sd = 0.005))
		pi[pi > 1] <- 1
		pi[pi < 0] <- 0
		truth <- rbinom(feature.no, prob = 1 - pi, size = 1)
	}
	
	ns <- sum(truth)
	
	# Generate the signals
	pi.m <- quantile(pi[truth == 1], 0.2)
	stats <- rnorm(feature.no)
	if (sig.dist == 'Normal') {
		if (sig.strength == 'Strong') stats[truth == 1 & pi > pi.m] <- rnorm(ns, 3, 1)[1:sum(truth == 1 & pi > pi.m)]
		if (sig.strength == 'Moderate') stats[truth == 1 & pi > pi.m] <- rnorm(ns, 2.5, 1)[1:sum(truth == 1 & pi > pi.m)]
		if (sig.strength == 'Weak') stats[truth == 1 & pi > pi.m] <- rnorm(ns, 2, 1)[1:sum(truth == 1 & pi > pi.m)]
		#stats[truth == 1 & pi <= pi.m] <- rnorm(ns, 2, 0.5)[1:sum(truth == 1 & pi <= pi.m)]
		
	}
	
	# Generate p values - one sided
	pvalue <- (1 - pnorm(stats))
	
	# Add a different f1
	pvalue [truth == 1 & pi <= pi.m] <- runif(sum(truth == 1 & pi <= pi.m), 0, 0.02)
	
	return(list(pvalue = pvalue, prior = pi, truth = truth))
}


# The following simulation is for figure 12 - varying f0
SimulateData9.4 <- function (
		prior.strength = c('Weak', 'Moderate', 'Strong'), 
		feature.no     = 10000,
		sig.dist       = c('Normal'), 
		sig.density    = c('Lower', 'Low', 'Medium', 'High'), 
		sig.strength   = c('Weak', 'Moderate', 'Strong')
) {
	
	
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
		pi <- c(rnorm(round(feature.no * (1 - ks[sig.density])), mean = 0.2, sd = 0.05), rnorm(feature.no * ks[sig.density], mean = 0.95, sd = 0.005))
		pi[pi > 1] <- 1
		pi[pi < 0] <- 0
		truth <- rbinom(feature.no, prob = 1 - pi, size = 1)
	}
	
	ns <- sum(truth)
	# Generate the signals
	
	stats <- rnorm(feature.no, 0, 1)
	pi.m <- quantile(pi[truth == 0], 0.8)
#	stats[pi >= pi.m] <- rnorm(sum(pi >= pi.m), -0.5, 1)
	
	if (sig.dist == 'Normal') {
		if (sig.strength == 'Strong') stats[truth == 1] <- rnorm(ns, 3, 1)
		if (sig.strength == 'Moderate') stats[truth == 1] <- rnorm(ns, 2.5, 1)
		if (sig.strength == 'Weak') stats[truth == 1] <- rnorm(ns, 2, 1)
	}
	
	# Generate p values - one sided
	pvalue <- (1 - pnorm(stats))
	
	# Add a different f0 (conservative)
	pvalue [truth == 0 & pi >= pi.m] <- runif(sum(truth == 0 & pi >= pi.m), 0.5, 1)
	
	return(list(pvalue = pvalue, prior = pi, truth = truth))
}


