
# Extra simulations in the case for d=3 to investigate the performance of SLIDE when signal to noise ratio is varied

# Source the needed files
source("Models_paper.R")
source("AuxillaryFunctions.R")
source("SLIDEfunctions.R")


Nrep=100
library(doParallel)
library(doRNG)
nworkers <- detectCores()
cl <- makePSOCKcluster(nworkers)
registerDoParallel(cl)



#######################################
# Varied signal to noise ratio
#######################################

# specify signal ratios
signal_ratios <- c(0.5, 0.7, 1, 2, 3)

# total number of loops
n_total = Nrep*length(signal_ratios)
signal_all = rep(signal_ratios, Nrep)


# Perform simulations
set.seed(2940765)
results_d3_snr <- foreach(i=1:n_total) %dorng% {
  # Generate data
  data <- generateModel2(n = 100, pvec = c(100,100,100), snr = signal_all[i], orthogonalV = T)

  # Apply SLIDE
  out = compareBCV_vs_Best_structure_refit_center(data$X, pvec = data$pvec, trueU = data$U, trueV = data$V, nl = 50, center = T)

  # Save results
  list(snr = signal_all[i], output = out)
}
save(results_d3_snr, file = "Results_d3_snr.Rda")


