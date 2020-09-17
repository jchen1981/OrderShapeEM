# OrderShapeEM
Optimal false discovery rate control using auxiliary (order) information based on a shape-constrained EM algorithm

## Overview
OrderShapeEM implements the optimal false discovery rate (FDR) control procedure with auxiliary information, particularly for prior ordering information. The framework is based on local FDR with hypothesis-specific null probability. The prior null proabilities are estimated using isotonic regression (PAVA algorithm) with respect to the prior ordering information. The inputs of our OrderShapeEM are simply P-values and their prior ordering. 

## Installation 

```
install.packages(c("Iso"))
if (!requireNamespace("BiocManager", quietly = TRUE))
   install.packages("BiocManager")
BiocManager::install(c("qvalue"))
install.packages("devtools")
devtools::install_github("jchen1981/OrderShapeEM")
```


### An Example
We illustrate the usage of tdfdr package using simulated data.

```
require(OrderShapeEM)
set.seed(123)

data.obj <- SimulateData(prior.strength = 'Moderate', sig.density = 'Low', sig.strength = 'Weak', 
feature.no = 5000)
orderfdr.obj <- OrderShapeEM(data.obj$pvalue, data.obj$prior, OrderShapeEM.control(trace = TRUE))

# Plot the estimated pi0 and f1
par(mfrow = c(1, 2))
plot(orderfdr.obj$pi0.step, xlab = 'index', ylab = 'pi0', do.points = FALSE,
          xlim = c(1, length(data.obj$pvalue)), main = 'Null probability')
plot(orderfdr.obj$f1.step, xlab = 'p value', ylab = 'f1', do.points = FALSE, 
          xlim = c(min(data.obj$pvalue), 1), log = 'x', main = 'Alternative distribution')

# Calculate the number of true positives and the false discovery proportion
sum(orderfdr.obj$fdr <= 0.05 & data.obj$truth)
sum(orderfdr.obj$fdr <= 0.05 & !data.obj$truth) / max(sum(orderfdr.obj$fdr <= 0.05), 1)

# Compare to the BH procedure
sum(p.adjust(data.obj$p.value, 'fdr') <= 0.05 & data.obj$truth)
sum(p.adjust(data.obj$p.value, 'fdr') <= 0.05 & !data.obj$truth) / max(sum(p.adjust(
	data.obj$p.value, 'fdr') <= 0.05), 1)
  
```

### Reproducibility
All the codes to reproduce the results in the manuscript are contained in 'Simulation' and 'RealData' folders. 

