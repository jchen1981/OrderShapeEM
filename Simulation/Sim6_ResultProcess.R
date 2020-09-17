setwd('~/Dropbox/Workspace/MayoClinic/Methodology/2016_11_14_CovariateFDR/Code/mforge/')
scene <- 'Sim6'
load(paste0('res/', scene, '_res.Rdata'))
setwd('pic')

require(reshape)
require(ggplot2)
###################################################################################################
# Combine the results
feature.nos <- '10000'
sig.dists <- c('PosNegCor')
prior.strengths <- c('Weak', 'Moderate', 'Strong')
sig.densities <- c('Low', 'Medium', 'High')
sig.strengths <- c('Weak', 'Moderate', 'Strong')
methods <- c( 'OrderShapeEM', 'OrderShapeEM+', 'AdaPT', 'AdaPT+', 'BH', 'ST', 'HingeExp', 'SeqStep', 'AdaptiveSeqStep', 'ForwardStop', 'SABHA')
measures <- c('FDR', 'Power')
res.a <- array(NA, c(length(prior.strengths), length(feature.nos), length(sig.dists), length(sig.densities), length(sig.strengths),
				length(methods), length(measures), length(res)), 
		dimnames=list(OrderPrior=prior.strengths, FeatureNumber=feature.nos, SignalDist=sig.dists,
				SignalDensity=sig.densities, SignalStrength=sig.strengths, Method=methods, 
				Measure=measures, Iteration=paste(1:length(res))))	

for (i in 1:length(res)) {
	res.a[, , , , , , , i] <- res[[i]]
}

###################################################################################################
# Calculate the standard error
res.df <- melt(res.a)
colnames(res.df)[9] <- 'Value'
# Error counts
error.info <- aggregate(Value ~ OrderPrior + FeatureNumber + SignalDist + SignalDensity + SignalStrength + 
				Method + Measure, res.df, function(x) mean(is.na(x)))
write.csv(error.info, 'NumberOfAlgFailure.csv')

m  <- aggregate(Value ~ OrderPrior + FeatureNumber + SignalDist + SignalDensity + SignalStrength + Method + 
				Measure, res.df, function(x) mean(x[!is.na(x)]))
se <- aggregate(Value ~ OrderPrior + FeatureNumber + SignalDist + SignalDensity + SignalStrength + Method + 
				Measure, res.df, function(x) {
			ind <- !is.na(x)
			sd(x[ind]) / sqrt(length(x[ind]))})
sd <- aggregate(Value ~ OrderPrior + FeatureNumber + SignalDist + SignalDensity + SignalStrength + Method + 
				Measure, res.df, function(x) {
			ind <- !is.na(x)
			sd(x[ind])})

ymin <- m[, ncol(m)] - 1.96 * se[, ncol(m)]
ymin <- ifelse(ymin > 0, ymin, 0)
res.df2 <- cbind(m, SD = sd[, ncol(sd)], ymax=m[, ncol(m)] + 1.96 * se[, ncol(m)], ymin=ymin,  SE=se[, ncol(se)], ErrRate=error.info[, ncol(m)])

res.df2$FeatureNumber <- factor(res.df2$FeatureNumber, levels=feature.nos)
res.df2$SignalDensity <- factor(res.df2$SignalDensity, levels=sig.densities)
res.df2$SignalStrength <- factor(res.df2$SignalStrength, levels=sig.strengths)
res.df2$OrderPrior <- factor(res.df2$OrderPrior, levels=prior.strengths)
res.df2$Method <- factor(res.df2$Method, levels=c( 'OrderShapeEM', 'OrderShapeEM+', 'AdaPT', 'AdaPT+', 'SABHA',  
				'HingeExp', 'SeqStep', 'ForwardStop',  'AdaptiveSeqStep',  'BH', 'ST'))
###################################################################################################
# Plot - Annonation, theme
algs <- c(  'OrderShapeEM', 'OrderShapeEM+', 'AdaPT', 'AdaPT+','SABHA', 
		'AdaptiveSeqStep', 'BH', 'ST', 'HingeExp', 'SeqStep',  'ForwardStop')
algs.ann <-c( 'OrderShapeEM', 'OrderShapeEM', 'AdaPT', 'AdaPT+', 'SABHA', 
		'AdaptiveSeqStep', 'BH', 'ST', 'HingeExp', 'SeqStep',  'ForwardStop')

names(algs.ann) <- algs

nMeth <- length(algs)
cols <- scales::hue_pal()(nMeth)
cols <- c("#D55E00", "#D55E00", "#E69F00", "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#CC79A7",  "#F0E442",   
		"steelblue", "orange", 'darkseagreen', 'firebrick2', 'gold2')[1:nMeth]  # Colorblind friendly
shapes <- c(rep(21:25, 3))[1:nMeth]
ltys <- rep(c(1, 2), ceiling(nMeth/2))[1:nMeth]
names(cols) <- names(shapes) <- names(ltys) <- algs.ann
sig.density.ann <- c('Low'='5% Signal', 'Medium'='10% Signal', 'High'='20% Signal')
sig.strength.ann <- c('Weak'='Weak Signal', 'Moderate'='Moderate Signal', 'Strong'='Strong Signal')

###################################################################################################
# Methods to be plotted and their orders
algs2 <-  c('OrderShapeEM+', 'AdaPT', 'SABHA', 'AdaptiveSeqStep', 'BH', 'ST')
Mea.ann <- c(FDR = 'Empirical False Discovery Rate', Power = 'True Positive Rate')


dodge <- position_dodge(width=0.9)
for (Dist in c('PosNegCor')) {
	cat("*")
	for (Mea in c('FDR', 'Power')) {
		cat(".")
		
		res3 <- subset(res.df2, Measure %in% Mea &  SignalDist %in% Dist & Method %in% algs2, drop=TRUE)	
		
		levels(res3$Method) <- algs.ann[levels(res3$Method)]
		levels(res3$SignalDensity) <- sig.density.ann[levels(res3$SignalDensity)]
		levels(res3$SignalStrength) <- sig.strength.ann[levels(res3$SignalStrength)]
		
		pdf(paste0("Figure2_", Dist, "_", Mea, "_grid_bar+.pdf"), width=8, height=5)
		obj <- ggplot(res3, aes(x=OrderPrior, y=Value,  fill=Method)) +
				geom_bar(stat='identity', position=dodge) +
				geom_errorbar(aes(ymax=ymax, ymin=ymin), position=dodge, width=0.25, size = 0.25)
		
		if (Mea == 'FDR') {
			obj <- obj + geom_hline(aes(yintercept=0.05), linetype=2) + 
					facet_grid(SignalDensity ~ SignalStrength)
		} else {
			obj <- obj + facet_grid(SignalDensity ~ SignalStrength, scales='free')
		}
		
		obj <- obj +
				xlab("Order Informativeness") +
				ylab(Mea.ann[Mea]) +
				scale_fill_manual(values=cols[algs.ann[algs2]]) +
				theme_bw() 
		#			theme(axis.text.x = element_text(angle=90, vjust=1)) 
		print(obj)
		dev.off()
	}
}
