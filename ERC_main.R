# Codes accompanying "Entropy Regularization in Probabilistic Clustering"
# Load relevant libraries, functions and data ----------------------------------
rm(list=ls())
# Set the working directory to the current folder 
# Code to set the working directory to the current folder from RStudio
library(rstudioapi) # version 0.14
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(viridis) # version 0.6.3
library(salso)   # version 0.3.29
library(ggplot2) # version 3.4.2
library(Cairo)   # version 1.6-0

# Load functions
source("ERC_fcts.R")
Save_Plot = F
# Simulation study univariate Gaussian -----------------------------------------
set.seed(0)
# Simulate the data from Gaussian mixture
n           = 1e3             # sample size
Kn_truth    = 3               # true number of components
mu_truth    = c(-4, 0, 4)     # true means
sigma_truth = rep(1,Kn_truth) # true sigma

data        = simdat(n = n, Kn=Kn_truth, mu = mu_truth, sigma=sigma_truth)
y           = as.matrix(data$y) - mean(as.matrix(data$y)) # center the data

# MCMC
# hyperparameters settings
# (alpha (concentration DP), m (mean base), s2 (), sigma2 ())
hyper        = c(1, 0, 1, 1) 
# MCMC quantities
Nsamples     = 2e4

# If you want to run the MCMC otherwise load the results
run_MCMC = F
# MCMC
if(run_MCMC){
  chain = mcmc_DPM_norm_norm(y, 
                             hyper   = hyper, 
                             clus    = FALSE, 
                             totiter = Nsamples)
  # save(chain,"Data-and-Results/chain_sim.Rdata")
  
} else {
  # Load results
  load("Data-and-Results/chain_sim.Rdata")
}
burnin     = 5000
thin       = 1
ps_samples = chain[seq(from=burnin+1, to=Nsamples, by=thin),]
N_ps       = nrow(ps_samples)

# Compute noisy clusters
noisyclus = double(N_ps)
for (iter in 1:N_ps){
  # compute frequencies of clusters in each iteration
  nc              = table(ps_samples[iter,]) 
  # number of cluster with small (<.1) relative freq
  nc_noisy        = nc[nc<(n/10)] 
  # proportion of obs assigned to small clusters in each iteration
  noisyclus[iter] = sum(nc_noisy) / n
}

# Number of iterations 
N_ps
# Number of iterations with more than .1 obs assigned to small clusters
sum(noisyclus>0.1)
# Number of iterations with more than .05 obs assigned to small clusters
sum(noisyclus>0.05)

# Histogram of percentages of obs assigned to small clusters (rel freq <.1)
# Figure 4 (a) of the manuscript
if(Save_Plot){CairoPNG(filename ='Image/hist.png', 
                       width = 800, height = 600)}
par(mar = c(5,5,2,2))
hist(noisyclus, xlab = "% of observations", ylab = "MCMC iterations", 
     col = viridis(3)[3],  main ="", family = "serif", font = 1, 
     font.lab = 2, cex.lab = 2.5, cex.axis = 2.2)
if(Save_Plot){invisible(dev.off())}

# Importance re-sampling for regularized posterior summaries 
# lambda = 10
reg_chain_10 = entropy_regularization(ps_samples, 10)
# lambda = 20
reg_chain_20 = entropy_regularization(ps_samples, 20)

# Compute noisy clusters in regularized clusters with lambda=10
noisyclus_reg_10 = double(N_ps)
for (iter in 1:N_ps){
  # compute frequencies of clusters in each iteration
  nc = table(reg_chain_10[iter,])
  # number of cluster with small (<.1) relative freq
  nc_noisy = nc[nc<(n/10)]
  noisyclus_reg_10[iter] = sum(nc_noisy)/n
}

# Histogram of percentages of obs assigned to small clusters (rel freq <.1)
# after regularization lambda=10
# Figure 4 (b) of the manuscript
if(Save_Plot){CairoPNG(filename ='Image/histreg10.png', 
                       width = 800, height = 600)}
par(mar = c(5,5,2,2))
hist(noisyclus_reg_10, xlab = "% of observations", ylab = "MCMC iterations", 
     col = viridis(3)[3],  main ="", family = "serif", font = 1, font.lab = 2, 
     cex.lab = 2.5, cex.axis = 2.2)
if(Save_Plot){invisible(dev.off())}

# Compute noisy clusters in regularized clusters with lambda=20
noisyclus_reg_20 = double(N_ps)
for (iter in 1:N_ps){
  nc = table(reg_chain_20[iter,])
  nc_noisy = nc[nc<(n/10)]
  noisyclus_reg_20[iter] = sum(nc_noisy)/n
}

# Histogram of percentages of obs assigned to small clusters (rel freq <.1)
# after regularization lambda=20
# Figure 4 (c) of the manuscript
if(Save_Plot){CairoPNG(filename ='Image/histafter.png', 
                       width = 800, height = 600)}
par(mar = c(5,5,2,2))
hist(noisyclus_reg_20, xlab = "% of observations", ylab = "MCMC iterations",
     col = viridis(3)[3],  main ="", family = "serif", font = 1, font.lab = 2, 
     cex.lab = 2.5, cex.axis = 2.2)
if(Save_Plot){invisible(dev.off())}

sum(noisyclus_reg_10>0.1)
sum(noisyclus_reg_10>0.05)

sum(noisyclus_reg_20>0.1)
sum(noisyclus_reg_20>0.05)

# Heatmap of the true clustering (with package "salso")
# Figure 5 (a) of the manuscript
pairwise_truth = psm(data$c_truth$V1, nCores = 1)
if(Save_Plot){CairoPNG(filename ='Image/truth_sim.png', 
                       width = 500, height = 500)}
heatmap(pairwise_truth, 
        Rowv = NA, Colv = NA, scale='none',
        col = viridis(2)[2:1])
if(Save_Plot){invisible(dev.off())}

# Compute the optimal partition under Binder loss (default a=1, lambda=0)
# If you want to run salso otherwise load the results
run_salso = F
if(run_salso){
  Binder_sim = salso(ps_samples, loss=binder(), 
                     maxZealousAttempts=20,  maxNClusters=n, nCores = 4)
} else {
  # Load results
  load("Data-and-Results/Binder_sim_all.Rdata")
  Binder_sim = Binder_sim_all[,"Binder_sim"]
}
# freq Binder point estimate (without regularization)
table(Binder_sim)


# Relative freq of clusters with Binder point estimate
table_Binder_sim = prop.table(table(Binder_sim))[order(table(Binder_sim), 
                                                       decreasing=TRUE)]
# Figure 6 (a) of the manuscript
if(Save_Plot){CairoPNG(filename ='Image/freqsim.png', 
                       width = 800, height = 600)}
par(mar = c(5,5,2,2))
barplot(table_Binder_sim, xlab="clusters", ylab="frequency", col=viridis(4)[3],
        main ="", family = "serif", font = 1, font.lab = 2,
        cex.lab = 2.5, cex.axis = 2.2)
if(Save_Plot){invisible(dev.off())}

# Heatmap of the optimal Binder clustering (with package "salso")
# Compute the adjacency matrix associated with the Binder point estimate
Binder_opt_adj = psm(Binder_sim, nCores = 1)
# Figure 5 (b) of the manuscript
if(Save_Plot){CairoPNG(filename ='Image/bindersim.png', 
                       width = 500, height = 500)}
heatmap(Binder_opt_adj, Rowv=NA, Colv=NA, scale='none', col=viridis(2)[2:1])
if(Save_Plot){invisible(dev.off())}

# Compute the optimal regularized (lambda=10) partition 
# under Binder loss (default a=1) 
# run_salso==T if you want to run salso otherwise load the results
if(run_salso){
  Binder_sim_10 = salso(reg_chain_10, loss=binder(), 
                        maxZealousAttempts=20,  maxNClusters=n, nCores = 4)
} else {
  Binder_sim_10 = Binder_sim_all[,"Binder_sim_10"]
}
# Relative freq of clusters with regularized (lambda=10) Binder point estimate
table_Binder_sim_10 = prop.table(table(Binder_sim_10))[
  order(table(Binder_sim_10), decreasing=TRUE)]
# Figure 6 (b) of the manuscript
if(Save_Plot){CairoPNG(filename ='Image/freq10.png', 
                       width = 800, height = 600)}
par(mar = c(5,5,2,2))
barplot(table_Binder_sim_10, xlab="clusters", ylab="frequency", 
        col=viridis(4)[3], main ="", family = "serif", font = 1, font.lab = 2,
        cex.lab = 2.5, cex.axis = 2.2)
if(Save_Plot){invisible(dev.off())}
# Heatmap of the optimal regularized Binder clustering (lambda=10)
# Compute the adjacency matrix of the reg (lambda=10) Binder point estimate
Binder_opt_adj_10 = psm(Binder_sim_10, nCores = 1)
# Figure 5 (c) of the manuscript
if(Save_Plot){CairoPNG(filename ='Image/binder10.png', 
                       width = 500, height = 500)}
heatmap(Binder_opt_adj_10, Rowv=NA, Colv=NA, scale='none', col=viridis(2)[2:1])
if(Save_Plot){invisible(dev.off())}

# Compute the optimal regularized (lambda=20) partition 
# under Binder loss (default a=1) 
# run_salso==T if you want to run salso otherwise load the results
if(run_salso){
  Binder_sim_20 = salso(reg_chain_20, loss=binder(), 
                        maxZealousAttempts=20,  maxNClusters=n, nCores = 4)
} else {
  Binder_sim_20 = Binder_sim_all[,"Binder_sim_20"]
}

## Compute miss-classification rates
Binder_table = table(data$c_truth$V1,Binder_sim)
sum(Binder_table)-sum(sort(Binder_table,decreasing = T)[1:3])

Binder_10_table = table(data$c_truth$V1,Binder_sim_10)
sum(Binder_10_table)-sum(sort(Binder_10_table,decreasing = T)[1:3])

Binder_20_table = table(data$c_truth$V1,Binder_sim_20)
sum(Binder_20_table)-sum(sort(Binder_20_table,decreasing = T)[1:3])

# Relative freq of clusters with regularized (lambda=20) Binder point estimate
table_Binder_sim_20 = prop.table(table(Binder_sim_20))[
  order(table(Binder_sim_20), decreasing=TRUE)]
# Figure 6 (c) of the manuscript
if(Save_Plot){CairoPNG(filename ='Image/freq20.png', 
                       width = 800, height = 600)}
par(mar = c(5,5,2,2))
barplot(table_Binder_sim_20, xlab="clusters", ylab="frequency", 
        col=viridis(4)[3], main ="", family = "serif", font = 1, font.lab = 2,
        cex.lab = 2.5, cex.axis = 2.2)
if(Save_Plot){invisible(dev.off())}

# Heatmap of the optimal regularized Binder clustering (lambda=20)
# Compute the adjacency matrix of the reg (lambda=20) Binder point estimate
Binder_opt_adj_20 = psm(Binder_sim_20, nCores = 1)
# Figure 5 (d) of the manuscript
if(Save_Plot){CairoPNG(filename ='Image/binder20.png', 
                       width = 500, height = 500)}
heatmap(Binder_opt_adj_20, Rowv=NA, Colv=NA, scale='none', col=viridis(2)[2:1])
if(Save_Plot){invisible(dev.off())}
# Compute the optimal partition under VI loss (default a=1)
# run_salso==T if you want to run salso otherwise load the results
if(run_salso){
  VI_sim     = salso(ps_samples, loss=VI(), 
                     maxZealousAttempts=20,  maxNClusters=n, nCores = 4)
  # save(VI_sim,file="Data-and-Results/VI_sim.Rdata")
} else {
  load("Data-and-Results/VI_sim.Rdata")
}
# freq VI point estimate (without regularization)
table(VI_sim)



# Wine data set ----------------------------------------------------------------
rm(list=ls())
library(viridis) # version 0.6.3
library(salso)   # version 0.3.29
library(ggplot2) # version 3.4.2
library(Cairo)   # version 1.6-0
# Load functions
source("ERC_fcts.R")
Save_Plot = F
# Data available at https://archive.ics.uci.edu/ml/datasets/wine
wine <- read.csv("Data-and-Results/wine.data", header=FALSE)
y = wine[,2:14] # we do multivariate
y = scale(y)
n = nrow(y)

# Heatmap of the ground truth of the clustering in the wine data set 
wine_types = psm(wine[,1], nCores = 1)
table(wine[,1])
# Figure 9 (a) of the manuscript
if(Save_Plot){CairoPNG(filename ='Image/truewine.png', 
                       width = 500, height = 500)}
heatmap(wine_types, Rowv = NA, Colv = NA, scale='none',
        col = viridis(2)[2:1])
if(Save_Plot){invisible(dev.off())}

# MCMC
# hyperparameters settings
# (alpha (concentration DP), m (mean base), s2 (), sigma2 ())
hyper        = c(0.1, 0, 1, 1)
# MCMC quantities
Nsamples     = 2e4

# If you want to run the MCMC otherwise load the results
run_MCMC = F
set.seed(0)
# MCMC
if(run_MCMC){
  chain_wine = mcmc_DPM_norm_norm_multi(y, 
                                        hyper   = hyper, 
                                        clus    = FALSE, 
                                        totiter = Nsamples)
} else {
  # Load results
  load("Data-and-Results/chain_wine.Rdata")
}

burnin     = 5000
thin       = 1
ps_samples = chain_wine[seq(from=burnin+1, to=Nsamples, by=thin),]
N_ps       = nrow(ps_samples)

# Compute the optimal partition under VI loss (default a=1)
VI_wine     = salso(ps_samples, loss=VI(), 
                      maxZealousAttempts=20,  maxNClusters=n, nCores = 4)
table(VI_wine)

# Compute the optimal partition under Binder loss (default a=1)
Binder_wine     = salso(ps_samples, loss=binder(), 
                        maxZealousAttempts=20,  maxNClusters=n, nCores = 4)
table(Binder_wine)

# Binder and VI produce different results
sum(as.integer(VI_wine)!=as.integer(Binder_wine))

# Heatmap of the optimal Binder clustering (with package "salso")
# Compute the adjacency matrix associated with the Binder point estimate
Binder_opt_adj = psm(Binder_wine, nCores = 1)
# Figure 9 (b) of the manuscript
if(Save_Plot){CairoPNG(filename ='Image/binderwine.png', 
                       width = 500, height = 500)}
heatmap(Binder_opt_adj, Rowv=NA, Colv=NA, scale='none', col=viridis(2)[2:1])
if(Save_Plot){invisible(dev.off())}
# Compute the adjacency matrix associated with the VI point estimate
VI_opt_adj = psm(VI_wine, nCores = 1)
# Figure 9 (c) of the manuscript
if(Save_Plot){CairoPNG(filename ='Image/VI_wine.png', 
                       width = 500, height = 500)}
heatmap(VI_opt_adj, Rowv=NA, Colv=NA, scale='none', col=viridis(2)[2:1])
if(Save_Plot){invisible(dev.off())}

# Importance re-sampling for regularized posterior summaries lambda = 50
reg_chain_50 = entropy_regularization(ps_samples, 50)

# Compute the optimal regularized (lambda=50) partition 
# under VI loss (default a=1)
VI_wine_50     = salso(reg_chain_50, loss=VI(), 
                    maxZealousAttempts=20,  maxNClusters=n, nCores = 4)
table(VI_wine_50)

# Compute the optimal regularized (lambda=50) partition 
# under Binder loss (default a=1)
Binder_wine_50     = salso(reg_chain_50, loss=binder(), 
                        maxZealousAttempts=20,  maxNClusters=n, nCores = 4)
table(Binder_wine_50)

# Heatmap of the optimal regularized Binder clustering (lambda=50)
# Compute the adjacency matrix of the reg (lambda=50) Binder point estimate
Binder_wine_opt_adj_50 = psm(Binder_wine_50, nCores = 1)
# Figure 9 (d) of the manuscript
if(Save_Plot){CairoPNG(filename ='Image/correctedbinder.png', 
                       width = 500, height = 500)}
heatmap(Binder_wine_opt_adj_50, Rowv=NA, Colv=NA, 
        scale='none', col=viridis(2)[2:1])
if(Save_Plot){invisible(dev.off())}

# Optimal regularized (lambda=50) VI and Binder partition are equal
sum(as.integer(VI_wine_50)!=as.integer(Binder_wine_50))
# Compute the adjacency matrix of the reg (lambda=50) Binder point estimate

# Figure 7 (a) of the manuscript
if(Save_Plot){CairoPNG(filename ='Image/binderwine_reorder.png', 
                       width = 500, height = 500)}
heatmap(psm(sort(Binder_wine)), Rowv = NA, Colv = NA, scale='none',
        col = viridis(2)[2:1])
if(Save_Plot){invisible(dev.off())}

# Figure 7 (c) of the manuscript
if(Save_Plot){CairoPNG(filename ='Image/correctedbinder_reorder.png', 
                       width = 500, height = 500)}
heatmap(psm(sort(Binder_wine_50)), Rowv = NA, Colv = NA, scale='none',
        col = viridis(2)[2:1])
if(Save_Plot){invisible(dev.off())}

# TBD Figure 7 (b) of the manuscript
if(Save_Plot){CairoPNG(filename ='Image/VI_wine_reorder.png', 
                       width = 500, height = 500)}
heatmap(psm(sort(VI_wine)), Rowv = NA, Colv = NA, scale='none',
        col = viridis(2)[2:1])
if(Save_Plot){invisible(dev.off())}

tab_VI     = table(wine[,1], VI_wine)        # 8 mistakes
tab_Binder = table(wine[,1], Binder_wine)    # 9 mistakes
tab_reg    = table(wine[,1], Binder_wine_50) # 6 mistakes

# Figure 8 (a) of the manuscript
table_Binder = prop.table(table(Binder_wine))[order(table(Binder_wine), 
                                                    decreasing = TRUE)]
names(table_Binder) = sort(names(table_Binder))

if(Save_Plot){CairoPNG(filename ='Image/clustnocorr.png', 
                       width = 800, height = 600)}
par(mar = c(5,5,2,2))
barplot(table_Binder, xlab = "clusters", ylab = "frequency", col =viridis(4)[3],
        main ="",family = "serif", font = 1, font.lab = 2, cex.lab = 2.5, 
        cex.axis = 2.2)
if(Save_Plot){invisible(dev.off())}

# Figure 8 (c) of the manuscript
table_Binder_50 = prop.table(table(Binder_wine_50))[
  order(table(table(Binder_wine_50)), decreasing = TRUE)]
names(table_Binder_50) = sort(names(table_Binder_50))
if(Save_Plot){CairoPNG(filename ='Image/cluster_corr.png', 
                       width = 800, height = 600)}
par(mar = c(5,5,2,2))
barplot(table_Binder_50, xlab = "clusters", ylab = "frequency", 
        col=viridis(4)[3], main ="", family = "serif", font = 1, font.lab = 2, 
        cex.lab = 2.5, cex.axis = 2.2)
if(Save_Plot){invisible(dev.off())}

# Figure 8 (b) of the manuscript
table_VI = prop.table(table(VI_wine))[order(table(VI_wine), 
                                                    decreasing = TRUE)]
names(table_VI) = sort(names(table_VI))

if(Save_Plot){CairoPNG(filename ='Image/VI_clustnocorr.png', 
                       width = 800, height = 600)}
par(mar = c(5,5,2,2))
barplot(table_VI, xlab = "clusters", ylab = "frequency", col =viridis(4)[3],
        main ="",family = "serif", font = 1, font.lab = 2, cex.lab = 2.5, 
        cex.axis = 2.2)
if(Save_Plot){invisible(dev.off())}

# Simulation study multivariate Bernoulli --------------------------------------
rm(list=ls())
library(viridis) # version 0.6.2
library(salso)   # version 0.3.29
library(ggplot2) # version 3.4.2
library(Cairo)   # version 1.6-0
library(gplots)
# Load functions
source("ERC_fcts.R")
Save_Plot = F
set.seed(0)
p = 50
Sigma = matrix(0.5, nrow = p, ncol = p)
diag(Sigma) = 1
normlatent = MASS::mvrnorm(n = 250, mu = rep(0, p), Sigma = Sigma)
pvars = pnorm(normlatent)
binvars = qbinom(pvars, 1, .3)
true_cor = as.vector(cor(binvars))
y = binvars

# MCMC
# MCMC quantities
Nsamples  = 2e4
# If you want to run the MCMC otherwise load the results
run_MCMC  = FALSE
Save_Plot = TRUE
# MCMC #
if(run_MCMC){
  set.seed(0)
  chain_sim_bern = mcmc_DPM_bern_beta_edited(y,  
                                      hyper=c(1, 0.2, 0.2), hprior = c(1,1),
                                      totiter = Nsamples, verbose = 10)
  # save(chain_sim_bern, file="./Data-and-Results/chain_sim_bern.RData")
} else {
  # Load results
  load("Data-and-Results/chain_sim_bern.Rdata")
}
burnin     = 5000
thin       = 1
ps_samples = chain_sim_bern[seq(from=burnin+1, to=Nsamples, by=thin),]
N_ps       = nrow(ps_samples)

# Compute noisy clusters
noisyclus = double(N_ps)
for (iter in 1:N_ps){
  # compute frequencies of clusters in each iteration
  nc              = table(ps_samples[iter,]) 
  # number of cluster with small (<.1) relative freq
  nc_noisy        = nc[nc<(n/10)] 
  # proportion of obs assigned to small clusters in each iteration
  noisyclus[iter] = sum(nc_noisy) / n
}

# Number of iterations 
N_ps
# Number of iterations with more than .1 obs assigned to small clusters
sum(noisyclus>0.1)
# Number of iterations with more than .05 obs assigned to small clusters
sum(noisyclus>0.05)

# Histogram of percentages of obs assigned to small clusters (rel freq <.1)
# Figure 7 (a) of the manuscript
if(Save_Plot){CairoPNG(filename ='Image/hist_bern.png', 
                       width = 800, height = 600)}
par(mar = c(5,5,2,2))
hist(noisyclus, xlab = "% of observations", ylab = "MCMC iterations", 
     col = viridis(3)[3],  main ="", family = "serif", font = 1, 
     font.lab = 2, cex.lab = 2.5, cex.axis = 2.2)
if(Save_Plot){invisible(dev.off())}

# Importance re-sampling for regularized posterior summaries 
# lambda = 10
reg_chain_10 = entropy_regularization(ps_samples, 10)
# lambda = 20
reg_chain_20 = entropy_regularization(ps_samples, 20)

# Compute noisy clusters in regularized clusters with lambda=10
noisyclus_reg_10 = double(N_ps)
for (iter in 1:N_ps){
  # compute frequencies of clusters in each iteration
  nc = table(reg_chain_10[iter,])
  # number of cluster with small (<.1) relative freq
  nc_noisy = nc[nc<(n/10)]
  noisyclus_reg_10[iter] = sum(nc_noisy)/n
}

# Number of iterations with more than .1 obs assigned to small clusters
sum(noisyclus_reg_10>0.1)
# Number of iterations with more than .05 obs assigned to small clusters
sum(noisyclus_reg_10>0.05)


# Histogram of percentages of obs assigned to small clusters (rel freq <.1)
# after regularization lambda=10
# Figure 7 (b) of the manuscript
if(Save_Plot){CairoPNG(filename ='Image/histreg10_bern.png', 
                       width = 800, height = 600)}
par(mar = c(5,5,2,2))
hist(noisyclus_reg_10, xlab = "% of observations", ylab = "MCMC iterations", 
     col = viridis(3)[3],  main ="", family = "serif", font = 1, font.lab = 2, 
     cex.lab = 2.5, cex.axis = 2.2)
if(Save_Plot){invisible(dev.off())}

# Compute noisy clusters in regularized clusters with lambda=100
noisyclus_reg_20 = double(N_ps)
for (iter in 1:N_ps){
  nc = table(reg_chain_20[iter,])
  nc_noisy = nc[nc<(n/10)]
  noisyclus_reg_20[iter] = sum(nc_noisy)/n
}

# Number of iterations with more than .1 obs assigned to small clusters
sum(noisyclus_reg_20>0.1)
# Number of iterations with more than .05 obs assigned to small clusters
sum(noisyclus_reg_20>0.05)

# Histogram of percentages of obs assigned to small clusters (rel freq <.1)
# after regularization lambda=20
# Figure 7 (c) of the manuscript
if(Save_Plot){CairoPNG(filename ='Image/histafter_bern.png', 
                       width = 800, height = 600)}
par(mar = c(5,5,2,2))
hist(noisyclus_reg_20, xlab = "% of observations", ylab = "MCMC iterations",
     col = viridis(3)[3],  main ="", family = "serif", font = 1, font.lab = 2, 
     cex.lab = 2.5, cex.axis = 2.2)
if(Save_Plot){invisible(dev.off())}

# Compute the optimal partition under Binder loss (default a=1, lambda=0)
# If you want to run salso otherwise load the results
run_salso = F
if(run_salso){
  Binder_sim = salso(ps_samples, loss=binder(), nCores = 4)
} else {
  # Load results
  load("Data-and-Results/Binder_bern_all.Rdata")
  Binder_sim = Binder_bern_all[,"Binder_sim"]
}
# freq Binder point estimate (without regularization)
table(Binder_sim)

# Relative freq of clusters with Binder point estimate
table_Binder_sim = prop.table(table(Binder_sim))[order(table(Binder_sim), 
                                                       decreasing=TRUE)]
# Figure 6 (a) of the manuscript
if(Save_Plot){CairoPNG(filename ='Image/freqsim_bern.png', 
                       width = 800, height = 600)}
par(mar = c(5,5,2,2))
barplot(table_Binder_sim, xlab="clusters", ylab="frequency", col=viridis(4)[3],
        main ="", family = "serif", font = 1, font.lab = 2,
        cex.lab = 2.5, cex.axis = 2.2)
if(Save_Plot){invisible(dev.off())}

# # Heatmap of the optimal Binder clustering (with package "salso")
# # Compute the adjacency matrix associated with the Binder point estimate
# Binder_opt_adj = psm(Binder_sim, nCores = 1)
# # Figure 5 (b) of the manuscript
# if(Save_Plot){CairoPNG(filename ='Image/binder_bern.png', 
#                        width = 500, height = 500)}
# heatmap.2(Binder_opt_adj,dendrogram='none', Rowv=TRUE, Colv=TRUE,
#           trace='none',col=viridis(2)[2:1], key=FALSE)
# if(Save_Plot){invisible(dev.off())}

# Compute the optimal regularized (lambda=10) partition 
# under Binder loss (default a=1) 
# run_salso==T if you want to run salso otherwise load the results
if(run_salso){
  Binder_sim_10 = salso(reg_chain_10, loss=binder(), nCores = 4)
} else {
  Binder_sim_10 = Binder_bern_all[,"Binder_sim_10"]
}
# Relative freq of clusters with regularized (lambda=10) Binder point estimate
table_Binder_sim_10 = prop.table(table(Binder_sim_10))[
  order(table(Binder_sim_10), decreasing=TRUE)]
# Figure 6 (b) of the manuscript
if(Save_Plot){CairoPNG(filename ='Image/freq10_bern.png', 
                       width = 800, height = 600)}
par(mar = c(5,5,2,2))
barplot(table_Binder_sim_10, xlab="clusters", ylab="frequency", 
        col=viridis(4)[3], main ="", family = "serif", font = 1, font.lab = 2,
        cex.lab = 2.5, cex.axis = 2.2)
if(Save_Plot){invisible(dev.off())}
# Heatmap of the optimal regularized Binder clustering (lambda=10)
# Compute the adjacency matrix of the reg (lambda=10) Binder point estimate
# Binder_opt_adj_10 = psm(Binder_sim_10, nCores = 1)
# # Figure 5 (c) of the manuscript
# if(Save_Plot){CairoPNG(filename ='Image/binder10_bern.png', 
#                        width = 500, height = 500)}
# heatmap.2(Binder_opt_adj_10,dendrogram='none', Rowv=TRUE, Colv=TRUE,
#           trace='none',col=viridis(2)[2:1], key=FALSE)
# if(Save_Plot){invisible(dev.off())}

# Compute the optimal regularized (lambda=20) partition 
# under Binder loss (default a=1) 
# run_salso==T if you want to run salso otherwise load the results
if(run_salso){
  Binder_sim_20 = salso(reg_chain_20, loss=binder(), 
                        maxZealousAttempts=20,  maxNClusters=n, nCores = 4)
} else {
  Binder_sim_20 = Binder_bern_all[,"Binder_sim_20"]
}

# Relative freq of clusters with regularized (lambda=20) Binder point estimate
table_Binder_sim_20 = prop.table(table(Binder_sim_20))[
  order(table(Binder_sim_20), decreasing=TRUE)]
# Figure 6 (c) of the manuscript
if(Save_Plot){CairoPNG(filename ='Image/freq20_bern.png', 
                       width = 800, height = 600)}
par(mar = c(5,5,2,2))
barplot(table_Binder_sim_20, xlab="clusters", ylab="frequency", 
        col=viridis(4)[3], main ="", family = "serif", font = 1, font.lab = 2,
        cex.lab = 2.5, cex.axis = 2.2)
if(Save_Plot){invisible(dev.off())}

# Heatmap of the optimal regularized Binder clustering (lambda=20)
# Compute the adjacency matrix of the reg (lambda=20) Binder point estimate
# Binder_opt_adj_20 = psm(Binder_sim_20, nCores = 1)
# # Figure 5 (d) of the manuscript
# if(Save_Plot){CairoPNG(filename ='Image/binder20_bern.png', 
#                        width = 500, height = 500)}
# heatmap.2(Binder_opt_adj_20,dendrogram='none', Rowv=TRUE, Colv=TRUE,
#           trace='none',col=viridis(2)[2:1], key=FALSE)
# if(Save_Plot){invisible(dev.off())}

Binder_bern_all        = cbind(Binder_sim, Binder_sim_10, Binder_sim_20)
names(Binder_bern_all) = c("Binder_sim", "Binder_sim_10", "Binder_sim_20")
# save(Binder_bern_all, file="./Data-and-Results/Binder_bern_all.RData")

# Compute the optimal partition under VI loss (default a=1)
# run_salso==T if you want to run salso otherwise load the results
if(run_salso){
  VI_sim     = salso(ps_samples, loss=VI(), nCores = 4)
} else {
  load("Data-and-Results/VI_bern.Rdata")
  VI_sim = VI_bern_all[,"VI_sim"]
}

# freq VI point estimate (without regularization)
table(VI_sim)

# Relative freq of clusters with Binder point estimate
table_VI_sim = prop.table(table(VI_sim))[order(table(VI_sim), 
                                                       decreasing=TRUE)]
# Figure 6 (a) of the manuscript
if(Save_Plot){CairoPNG(filename ='Image/freqsim_bern_VI.png', 
                       width = 800, height = 600)}
par(mar = c(5,5,2,2))
barplot(table_VI_sim, xlab="clusters", ylab="frequency", col=viridis(4)[3],
        main ="", family = "serif", font = 1, font.lab = 2,
        cex.lab = 2.5, cex.axis = 2.2)
if(Save_Plot){invisible(dev.off())}

# Heatmap of the optimal Binder clustering (with package "salso")
# Compute the adjacency matrix associated with the Binder point estimate
# VI_opt_adj = psm(VI_sim, nCores = 1)
# # Figure 5 (b) of the manuscript
# if(Save_Plot){CairoPNG(filename ='Image/VI_bern.png', 
#                        width = 500, height = 500)}
# heatmap.2(VI_opt_adj,dendrogram='none', Rowv=TRUE, Colv=TRUE,
#           trace='none',col=viridis(2)[2:1], key=FALSE)
# if(Save_Plot){invisible(dev.off())}

# Compute the optimal regularized (lambda=10) partition 
# under Binder loss (default a=1) 
# run_salso==T if you want to run salso otherwise load the results
if(run_salso){
  VI_sim_10 = salso(reg_chain_10, loss=VI(), nCores = 4)
} else {
  VI_sim_10 = VI_bern_all[,"VI_sim_10"]
}
# Relative freq of clusters with regularized (lambda=10) Binder point estimate
table_VI_sim_10 = prop.table(table(VI_sim_10))[
  order(table(VI_sim_10), decreasing=TRUE)]
# Figure 6 (b) of the manuscript
if(Save_Plot){CairoPNG(filename ='Image/freq10_bern_VI.png', 
                       width = 800, height = 600)}
par(mar = c(5,5,2,2))
barplot(table_VI_sim_10, xlab="clusters", ylab="frequency", 
        col=viridis(4)[3], main ="", family = "serif", font = 1, font.lab = 2,
        cex.lab = 2.5, cex.axis = 2.2)
if(Save_Plot){invisible(dev.off())}
# Heatmap of the optimal regularized Binder clustering (lambda=10)
# Compute the adjacency matrix of the reg (lambda=10) Binder point estimate
# VI_opt_adj_10 = psm(VI_sim_10, nCores = 1)
# # Figure 5 (c) of the manuscript
# if(Save_Plot){CairoPNG(filename ='Image/binder10_bern_VI.png', 
#                        width = 500, height = 500)}
# heatmap.2(VI_opt_adj_10,dendrogram='none', Rowv=TRUE, Colv=TRUE,
#           trace='none',col=viridis(2)[2:1], key=FALSE)
# if(Save_Plot){invisible(dev.off())}

# Compute the optimal regularized (lambda=20) partition 
# under Binder loss (default a=1) 
# run_salso==T if you want to run salso otherwise load the results
if(run_salso){
  VI_sim_20 = salso(reg_chain_20, loss=VI(), nCores = 4)
} else {
  VI_sim_20 = VI_bern_all[,"VI_sim_20"]
}

# Relative freq of clusters with regularized (lambda=20) Binder point estimate
table_VI_sim_20 = prop.table(table(VI_sim_20))[
  order(table(VI_sim_20), decreasing=TRUE)]
# Figure 6 (c) of the manuscript
if(Save_Plot){CairoPNG(filename ='Image/freq20_bern_VI.png', 
                       width = 800, height = 600)}
par(mar = c(5,5,2,2))
barplot(table_VI_sim_20, xlab="clusters", ylab="frequency", 
        col=viridis(4)[3], main ="", family = "serif", font = 1, font.lab = 2,
        cex.lab = 2.5, cex.axis = 2.2)
if(Save_Plot){invisible(dev.off())}

# Heatmap of the optimal regularized Binder clustering (lambda=20)
# Compute the adjacency matrix of the reg (lambda=20) Binder point estimate
# VI_opt_adj_20 = psm(VI_sim_20, nCores = 1)
# # Figure 5 (d) of the manuscript
# if(Save_Plot){CairoPNG(filename ='Image/binder20_bern_VI.png', 
#                        width = 500, height = 500)}
# heatmap.2(VI_opt_adj_20,dendrogram='none', Rowv=TRUE, Colv=TRUE,
#           trace='none',col=viridis(2)[2:1], key=FALSE)
# if(Save_Plot){invisible(dev.off())}

VI_bern_all        = cbind( VI_sim, VI_sim_10, VI_sim_20)
names(VI_bern_all) = c( "VI_sim", "VI_sim_10", "VI_sim_20")
# save(VI_bern_all, file="./Data-and-Results/VI_bern_all.RData")

