# For figure 13 - different numbers of features
# Signal density is fixed to be 'Medium'
func <- function(part, paras) {
	require(qvalue)
	require(adaptMT)
	require(splines)
	require(OrderShapeEM)
	
	set.seed(part)
	feature.nos <- paras$feature.nos
	sig.dists <- paras$sig.dists
	prior.strengths <- paras$prior.strengths
	sig.densities <- paras$sig.densities
	sig.strengths <- paras$sig.strengths
	fdr.cutoff <- paras$fdr.cutoff
	resdir <- paras$resdir
	
	resdir <- gsub('/$', '', resdir)
	source(file.path(resdir, "SimFunc.R"))
	source(file.path(resdir, "accumulation_test_functions.R"))
	source(file.path(resdir, 'All_q_est_functions.R'))
	
	prefix <- 'Sim4.6'
	sink(file.path(resdir, paste(prefix, "_",  part, ".log", sep="")))
	cat(date(), '\n')
	
	methods <- c('OrderShapeEM', 'OrderShapeEM+', 'AdaPT', 'AdaPT+', 'BH', 'ST', 'HingeExp', 'SeqStep', 'AdaptiveSeqStep', 'ForwardStop', 'SABHA')
	measures <- c('FDR', 'Power')
	res <- array(NA, c(length(prior.strengths), length(feature.nos), length(sig.dists), length(sig.densities), length(sig.strengths),
					length(methods), length(measures)), 
			dimnames=list(prior.strengths, feature.nos, sig.dists, sig.densities, sig.strengths, methods, measures))	
	
	for (prior.strength in prior.strengths) {
		cat('$')
		for (feature.no in feature.nos) {
			cat('*')
			for (sig.dist in sig.dists) {
				cat('!')
				for (sig.density in sig.densities) {
					cat('%')
					for (sig.strength in sig.strengths) {
						cat('#')
						
						# Different numbers of features 
						m <- switch(sig.density,
								Low = 500, 
								Medium = 1000,
								High = 2000)
						
						data <- SimulateData4(prior.strength, feature.no = m, sig.dist, sig.density = 'Medium', sig.strength)
						prior <- data$prior
						pvalue <- data$pvalue
						truth <- data$truth
						pi0.true <- mean(!truth)
						
						# Call methods
						# OrderShapeEM - the original algorithm
						# OrderShapeEM+ - only use the OrderShapeEM if the order is informative (sd(pi0) > 0.025), else use qvalue
						# AdaPT - the original algorithm
						# AdaPT+ -  without '+1' correction term 
						res.obj <- list()
						try({
									res.obj[['OrderShapeEM']] <- res.obj[['OrderShapeEM+']] <- 
											OrderShapeEM(pvalue, prior)
									if (sd(res.obj[['OrderShapeEM']]$pi0.unadj) < 0.025) res.obj[['OrderShapeEM+']]$fdr <- qvalue(pvalue)$qvalue
								})
						
						try({
									cat('.')
									fdr <- rep(1, feature.no)
									x.df <- data.frame(x = prior)
									formulas <- paste0("ns(x, df=6)")
									adapt1.res <- adapt_glm(x = x.df, pvals = pvalue, pi_formulas = formulas, mu_formulas = formulas, 
											verbose = list(print = FALSE, fit = FALSE, ms = FALSE))
									fdr[pvalue <= adapt1.res$s[, floor(fdr.cutoff * 100)]]  <- fdr.cutoff / 2
									res.obj[['AdaPT']] <- list(fdr=fdr)
									
								})
						
						try({
									source(file.path(resdir, 'adaptMT2.R'))
									cat('.')
									fdr <- rep(1, feature.no)
									x.df <- data.frame(x = prior)
									formulas <- paste0("ns(x, df=6)")
									adapt2.res <- adapt_glm2(x = x.df, pvals = pvalue, pi_formulas = formulas, mu_formulas = formulas, 
											verbose = list(print = FALSE, fit = FALSE, ms = FALSE))
									fdr[pvalue <= adapt2.res$s[, floor(fdr.cutoff * 100)]]  <- fdr.cutoff / 2
									res.obj[['AdaPT+']] <- list(fdr=fdr)
									
								})
						
						res.obj[['BH']] <- list(fdr=p.adjust(pvalue, 'fdr'))
						res.obj[['ST']] <- list(fdr=qvalue(pvalue)$qvalue)
						
						# Li & Barber
						tau = 0.5; eps = 0.1 # parameters for SABHA
						thr = 0.5 # parameter for Storey-BH
						thr1 = 0.1; thr2 = 0.5 # parameters for adaptive SeqStep
						
						
						index0 <- order(prior)
						index1 <- order(index0)
						
						khats <- HingeExp(pvalue[index0] * (1-1/(1+choose(10,5))), alpha = fdr.cutoff, C = 2)
						fdr <- rep(1, feature.no)
						fdr[1:khats] <- fdr.cutoff / 2
						res.obj[['HingeExp']] <- list(fdr=fdr[index1])
						
						khats <- SeqStep(pvalue[index0], alpha = fdr.cutoff, C = 2)
						fdr <- rep(1, feature.no)
						fdr[1:khats] <- fdr.cutoff / 2
						res.obj[['SeqStep']] <- list(fdr=fdr[index1])
						
						
						khats <- Adaptive_SeqStep_method(pvalue[index0], fdr.cutoff, thr1, thr2)
						fdr <- rep(1, feature.no)
						fdr[khats] <- fdr.cutoff / 2
						res.obj[['AdaptiveSeqStep']] <- list(fdr=fdr[index1])
						
						
						khats <- ForwardStop(pvalue[index0] * (1-1/(1+choose(10,5))), alpha = fdr.cutoff)
						fdr <- rep(1, feature.no)
						fdr[1:khats] <- fdr.cutoff / 2
						res.obj[['ForwardStop']] <- list(fdr=fdr[index1])
						
						try({ 
									fdr <- rep(1, feature.no)
									pvals <- pvalue[index0]
									qhat <- Solve_q_step(pvals, 0.5, 0.1)
									fdr[SABHA_method(pvals, qhat, fdr.cutoff, 0.5)]  <- fdr.cutoff / 2
									res.obj[['SABHA']] <- list(fdr=fdr[index1])
								})
						
				
						for (method in methods) {
							
							# Evaluation
							if (is.null(res.obj[[method]])) {
								res[prior.strength, paste(feature.no), sig.dist, sig.density, sig.strength, method, c('FDR', 'Power')] <- NA
								next
							} 
							
							fdr <- res.obj[[method]]$fdr
							
							if (sum(fdr <= fdr.cutoff) == 0) {
								res[prior.strength, paste(feature.no), sig.dist, sig.density, sig.strength, method, 'FDR'] <- 0
								res[prior.strength, paste(feature.no), sig.dist, sig.density, sig.strength, method, 'Power'] <- 0
							} else {
								
								res[prior.strength, paste(feature.no), sig.dist, sig.density, sig.strength, method, 'FDR'] <- mean(truth[fdr <= fdr.cutoff] == 0)
								res[prior.strength, paste(feature.no), sig.dist, sig.density, sig.strength, method, 'Power'] <- sum(truth[fdr <= fdr.cutoff] != 0) / sum(truth != 0)
							}
						}
						
					}
				}
			}
		}
	}
	
	cat('\n', date(), '\n')
	sink()
	save(res, file=file.path(resdir, paste(prefix, "_res",  part, ".Rdata", sep="")))
	return(res)
}

paras <- list()
paras$resdir <- "~/project/orderfdr/mforge"
paras$fdr.cutoff <- 0.05
paras$feature.nos <- 10000
paras$sig.dists <- c('Normal')
paras$prior.strengths <- c('Weak', 'Moderate', 'Strong')
paras$sig.densities <- c('Low', 'Medium', 'High')
paras$sig.strengths <- c('Weak', 'Moderate', 'Strong')

resdir <- paras$resdir
setwd(resdir)
res <- clsapply(1:100, func, paras, queque="1-day",  tempdir="~/project/orderfdr/mforge/Sim4.6")
setwd(resdir)

