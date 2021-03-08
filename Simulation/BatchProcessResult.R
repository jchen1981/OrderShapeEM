scenes <- c('Sim8.4', 'Sim9.4', 'Sim4.1', 'Sim4.2', 'Sim4.3', 'Sim4.4.new', 'Sim4.6', 'Sim5.2', 'Sim6')


for (scene in scenes) {
	cat('.')
	source(paste0('~/Dropbox/Workspace/MayoClinic/Methodology/2016_11_14_CovariateFDR/Code/mforge/code/', scene, '_ResultProcess.R'))
}

cat('Finished!\n')
