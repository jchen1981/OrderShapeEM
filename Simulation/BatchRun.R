
setwd("~/project/orderfdr/mforge/")
source("Cluster_mayo.R")

scenes <- c('Sim4.2', 'Sim4.3', 'Sim4.4.new',  'Sim4.6', 'Sim5.2', 'Sim6', 'Sim8.4', 'Sim9.4')

for (scene in scenes) {
	source(paste0(scene, '.R'))
}

cat('Submitted!\n')

