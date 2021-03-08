
scenes <- c('Sim8.4', 'Sim9.4', 'Sim4.1', 'Sim4.2', 'Sim4.3', 'Sim4.4.new', 'Sim4.6', 'Sim5.2', 'Sim6')
setwd("~/project/orderfdr/mforge/")

for (scene in scenes) {
	missingfile <- NULL
	res <- NULL
	nJob <- 100
	for (part in 1:nJob) {
		resfile <- paste0('./', scene, '/', part, ".res", sep="")
		if (file.exists(resfile)) {
			load(resfile)
			res[[paste(part)]] <-  res0
			rm(res0)
			cat(".")
		} else {
			# cat("\n", resfile, " not found!\n")
			missingfile <- c(missingfile, resfile)
		}
	}
	
	cat("Missing files:\n")
	if (length(missingfile) == 0) {
		cat("No jobs dropped!\n")
	} else{
		cat(paste(missingfile, collapse="\n"))
	} 
	
	save(res, file=file.path('res', paste(scene, "_res.Rdata", sep="")))
}
