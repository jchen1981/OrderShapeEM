# c4d GWAS p-values as auxiliary data

require(qvalue)
require(ggplot2)
require(reshape2)
require(OrderShapeEM)

source('All_q_est_functions.R')
source('accumulation_test_functions.R')

load('exampledata')

p2.2 <- cad$P_c4d
p1.1 <- cad$P_cardiogram

##################################################
# Run OrderShapeEM
sink('Figure5b_time.txt')
cat('Start ...\n')
cat(date())
obj <- OrderShapeEM(p1.1, p2.2)
cat(date())
cat('End ...\n')
sink()
save(obj, file = 'Figure5b_igwas_data.RData')


pdf('Figure5b_pi0_f1_stepfun.pdf',
		height = 5, width = 8)
par(mfrow = c(1, 2))
plot(obj$pi0.step, xlab = 'index', ylab = 'Probability', do.points = FALSE,
		xlim = c(1, length(p1.1)), main = 'pi0')
plot(obj$f1.step, xlab = 'p-value', ylab = 'Density', do.points = FALSE, 
		xlim = c(1e-15, 1), log = 'x', main = 'f1')
dev.off()

##################################################
# Run other methods
alphalist <- 10^seq(-5, -2, len = 50)
tau = 0.5; eps = 0.1 # parameters for SABHA
thr = 0.5 # parameter for Storey-BH
thr1 = 0.1; thr2 = 0.5 # parameters for adaptive SeqStep

num_alpha <- length(alphalist)
methods <- c('SeqStep', 'HingeExp', 'ForwardStop', 'AdaptiveSeqStep', 'SABHA', 'OrderShapeEM', 'BH', 'ST')

pvals <- p1.1[order(p2.2)]
qhat <- Solve_q_step(pvals, tau, eps)
res <- array(NA, c(num_alpha, length(methods)), dimnames = list(paste(alphalist), methods))

load(file = 'Figure5b_igwas_data.RData')
q1 <- obj$fdr
q2 <- p.adjust(p1.1, 'fdr')
q3 <- qvalue(p1.1)$qvalue

for(i in 1:num_alpha){
	cat('.')
	res[i, 'SeqStep'] <-  SeqStep(pvals,alpha=alphalist[i],C=2)
	res[i, 'HingeExp'] <-  HingeExp(pvals*(1-1/(1+choose(10,5))),alpha=alphalist[i],C=2)
	res[i, 'ForwardStop'] <-  ForwardStop(pvals*(1-1/(1+choose(10,5))),alpha=alphalist[i])
	res[i, 'AdaptiveSeqStep'] <-  length(Adaptive_SeqStep_method(pvals, alphalist[i], thr1, thr2))
	res[i, 'SABHA'] <-  length(SABHA_method(pvals, qhat, alphalist[i], tau))
	res[i, 'OrderShapeEM'] <-  sum(q1 <= alphalist[i])
	res[i, 'BH'] <-  sum(q2 <= alphalist[i])
	res[i, 'ST'] <-  sum(q3 <= alphalist[i])
	
}
res <- melt(res)
colnames(res) <- c('Cutoff', 'Method', 'value')
save(res, file = 'Figure5b_plot_data.RData')
##################################################
# Plot
temp <- load(file = 'Figure5b_plot_data.RData')
res <- subset(res, !(Method %in% c('Hybrid', 'SeqStep', 'HingeExp', 'ForwardStop', 'AdaptiveSeqStep', 'BH')))

col <- c(AdaptiveSeqStep="#009E73", SABHA="#56B4E9", OrderShapeEM="#D55E00", ST="#CC79A7")
lty <- c(AdaptiveSeqStep=3, SABHA=3, OrderShapeEM=1,  ST=1)
shape <- c( AdaptiveSeqStep=23, SABHA=22, OrderShapeEM=21,  ST=24)

res <- subset(res, Cutoff <= 0.005)
res2 <- res
res2$Cutoff <- factor(res2$Cutoff)
res2 <- res2[as.numeric(res2$Cutoff) %in% c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45), ]
res2$Cutoff <- as.numeric(as.character(res2$Cutoff))

res2$Method <- factor(res2$Method, levels=c('OrderShapeEM', 'SABHA',   'ST'))
res$Method <- factor(res$Method, levels=c('OrderShapeEM', 'SABHA',   'ST'))


pdf('Figure5b_number_of_hits.pdf', height = 5, width = 6)
ggplot(res, aes(x = Cutoff, y = value, group = Method,  colour = Method, fill = Method, linetype = Method)) +
		geom_line(size=0.25) +
		geom_point(data=res2, aes(x = Cutoff, y = value, shape = Method)) +
		ylab('Number of hits') +
		xlab('Target FDR level') +
		theme_bw() +
		ylim(c(0, 100)) +
		scale_linetype_manual(values = lty) +
		scale_colour_manual(values = col) +
		scale_fill_manual(values = col) +
		scale_shape_manual(values = shape)
dev.off()	

##################################################


