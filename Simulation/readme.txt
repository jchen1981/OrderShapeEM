The simulations are designed to run on the cluster ("mforge" at Mayo-UIUC). Sim4.2, Sim4.3, Sim4.4.new, Sim4.6, Sim5.2, Sim6, Sim8.4 and Sim9.4 represent
different settings described in the text:

#######################
Sim4.2 - figure 2, 6, 16
Sim4.3 - figure 9, 14 
Sim4.4.new - figure 10
Sim4.6 - figure 13
Sim5.2 - figure 7
Sim6 - figure 3
Sim8.4 - figure 11
Sim9.4 - figure 12
#######################

#######################
SimFunc.R - simulation functions for different settings

All_q_est_functions.R - accumulation tests of Ang Li & Rina Foygel Barber
accumulation_test_functions.R - accumulation tests of Ang Li & Rina Foygel Barber
adaptMT2.R - AdaPT without the correction term

Cluster_mayo.R - the clsapply function to run simulations on the Mayo-UIUC mforge cluster

BatchRun.R  - calls "Sim*.R" to submit all the simulation jobs to the cluster
BatchProcess.R - Process the results from the cluster runs
BatchProcessResult.R - calls "Sim*_ResultProcess.R" to generate pictures for the manuscript
#######################

#######################
# Basic process (after creating/setting relevant directories)
BatchRun.R
BatchProcess.R
BatchProcessResult.R
#######################





