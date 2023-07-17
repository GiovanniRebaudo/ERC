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
Save_Plot = FALSE
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
# (alpha initialization (concentration DP), m (mean base), s2 (), sigma2 ())
hyper        = c(1, 0, 1, 1) 
# alpha (concentration DP) is gamma distributed with parameters 
hprior = c(1,1)
# MCMC quantities
Nsamples     = 2e4

# If you want to run the MCMC otherwise load the results
run_MCMC = T
# MCMC
if(run_MCMC){
  set.seed(0)
  chain = mcmc_DPM_norm_norm(y, 
                             hyper   = hyper, 
                             clus    = FALSE,
                             hprior  = hprior,
                             totiter = Nsamples)
  # save(chain, file="Data-and-Results/chain_sim_app.Rdata")
} else {
  # Load results
  load("Data-and-Results/chain_sim_app.Rdata")
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
if(Save_Plot){CairoPNG(filename ='Image/hist_app.png', 
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
if(Save_Plot){CairoPNG(filename ='Image/histreg10_app.png', 
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
if(Save_Plot){CairoPNG(filename ='Image/histafter_app.png', 
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
pairwise_truth = psm(data$c_truth$V1, nCores = 1)
if(Save_Plot){CairoPNG(filename ='Image/truth_sim_app.png', 
                       width = 500, height = 500)}
heatmap(pairwise_truth, 
        Rowv = NA, Colv = NA, scale='none',
        col = viridis(2)[2:1])
if(Save_Plot){invisible(dev.off())}

# Compute the optimal partition under Binder loss (default a=1, lambda=0)
# If you want to run salso otherwise load the results
run_salso = T
if(run_salso){
  Binder_sim = salso(ps_samples, loss=binder(), nCores = 4)
} else {
  # Load results
  load("Data-and-Results/Binder_sim_all_app.Rdata")
  Binder_sim = Binder_sim_all_app[,"Binder_sim"]
}
# freq Binder point estimate (without regularization)
table(Binder_sim)

# Relative freq of clusters with Binder point estimate
table_Binder_sim = prop.table(table(Binder_sim))[order(table(Binder_sim), 
                                                       decreasing=TRUE)]
if(Save_Plot){CairoPNG(filename ='Image/freqsim_app.png', 
                       width = 800, height = 600)}
par(mar = c(5,5,2,2))
barplot(table_Binder_sim, xlab="clusters", ylab="frequency", col=viridis(4)[3],
        main ="", family = "serif", font = 1, font.lab = 2,
        cex.lab = 2.5, cex.axis = 2.2)
if(Save_Plot){invisible(dev.off())}

# Heatmap of the optimal Binder clustering (with package "salso")
# Compute the adjacency matrix associated with the Binder point estimate
Binder_opt_adj = psm(Binder_sim, nCores = 1)
if(Save_Plot){CairoPNG(filename ='Image/bindersim_app.png', 
                       width = 500, height = 500)}
heatmap(Binder_opt_adj, Rowv=NA, Colv=NA, scale='none', col=viridis(2)[2:1])
if(Save_Plot){invisible(dev.off())}

# Compute the optimal regularized (lambda=10) partition 
# under Binder loss (default a=1) 
# run_salso==T if you want to run salso otherwise load the results
if(run_salso){
  Binder_sim_10 = salso(reg_chain_10, loss=binder(), nCores = 4)
} else {
  Binder_sim_10 = Binder_sim_all_app[,"Binder_sim_10"]
}
# Relative freq of clusters with regularized (lambda=10) Binder point estimate
table_Binder_sim_10 = prop.table(table(Binder_sim_10))[
  order(table(Binder_sim_10), decreasing=TRUE)]
if(Save_Plot){CairoPNG(filename ='Image/freq10_app.png', 
                       width = 800, height = 600)}
par(mar = c(5,5,2,2))
barplot(table_Binder_sim_10, xlab="clusters", ylab="frequency", 
        col=viridis(4)[3], main ="", family = "serif", font = 1, font.lab = 2,
        cex.lab = 2.5, cex.axis = 2.2)
if(Save_Plot){invisible(dev.off())}
# Heatmap of the optimal regularized Binder clustering (lambda=10)
# Compute the adjacency matrix of the reg (lambda=10) Binder point estimate
Binder_opt_adj_10 = psm(Binder_sim_10, nCores = 1)
if(Save_Plot){CairoPNG(filename ='Image/binder10_app.png', 
                       width = 500, height = 500)}
heatmap(Binder_opt_adj_10, Rowv=NA, Colv=NA, scale='none', col=viridis(2)[2:1])
if(Save_Plot){invisible(dev.off())}

# Compute the optimal regularized (lambda=20) partition 
# under Binder loss (default a=1) 
# run_salso==T if you want to run salso otherwise load the results
if(run_salso){
  Binder_sim_20 = salso(reg_chain_20, loss=binder(), nCores = 4)
} else {
  Binder_sim_20 = Binder_sim_all_app[,"Binder_sim_20"]
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
if(Save_Plot){CairoPNG(filename ='Image/freq20_app.png', 
                       width = 800, height = 600)}
par(mar = c(5,5,2,2))
barplot(table_Binder_sim_20, xlab="clusters", ylab="frequency", 
        col=viridis(4)[3], main ="", family = "serif", font = 1, font.lab = 2,
        cex.lab = 2.5, cex.axis = 2.2)
if(Save_Plot){invisible(dev.off())}

# Heatmap of the optimal regularized Binder clustering (lambda=20)
# Compute the adjacency matrix of the reg (lambda=20) Binder point estimate
Binder_opt_adj_20 = psm(Binder_sim_20, nCores = 1)
if(Save_Plot){CairoPNG(filename ='Image/binder20_app.png', 
                       width = 500, height = 500)}
heatmap(Binder_opt_adj_20, Rowv=NA, Colv=NA, scale='none', col=viridis(2)[2:1])
if(Save_Plot){invisible(dev.off())}

Binder_sim_all_app        = cbind(Binder_sim, Binder_sim_10, Binder_sim_20)
names(Binder_sim_all_app) = c("Binder_sim", "Binder_sim_10", "Binder_sim_20")
# save(Binder_sim_all_app, file="./Data-and-Results/Binder_sim_all_app.RData")

# Compute the optimal partition under VI loss (default a=1)
# run_salso==T if you want to run salso otherwise load the results
if(run_salso){
  set.seed(0)
  VI_sim     = salso(ps_samples,   loss=VI(), nCores = 4)
  VI_sim_10  = salso(reg_chain_10, loss=VI(), nCores = 4)
  VI_sim_20  = salso(reg_chain_20, loss=VI(), nCores = 4)
} else {
  load("Data-and-Results/VI_sim_all_app.Rdata")
  VI_sim     = VI_sim_all_app[,"VI_sim"]
  VI_sim_10  = VI_sim_all_app[,"VI_sim_10"]
  VI_sim_20  = VI_sim_all_app[,"VI_sim_20"]
}
VI_sim_all_app        = cbind(VI_sim, VI_sim_10, VI_sim_20)
names(VI_sim_all_app) = c("VI_sim", "VI_sim_10", "VI_sim_20")
# save(VI_sim_all_app, file="./Data-and-Results/VI_sim_all_app.RData")

# freq VI point estimate (without regularization)
table(VI_sim)


# Wine data set ----------------------------------------------------------------
rm(list=ls())
# Load functions
source("ERC_fcts.R")
Save_Plot = FALSE
# Data available at https://archive.ics.uci.edu/ml/datasets/wine
wine <- read.csv("Data-and-Results/wine.data", header=FALSE)
y = wine[,2:14] # we do multivariate
y = scale(y)
n = nrow(y)

# Heatmap of the ground truth of the clustering in the wine data set 
wine_types = psm(wine[,1], nCores = 1)
table(wine[,1])
if(Save_Plot){CairoPNG(filename ='Image/truewine_app.png', 
                       width = 500, height = 500)}
heatmap(wine_types, Rowv = NA, Colv = NA, scale='none',
        col = viridis(2)[2:1])
if(Save_Plot){invisible(dev.off())}

# MCMC
# hyperparameters settings
# (alpha (concentration DP) random, m (mean base), s2 (), sigma2 ())
hyper        = c(0.1, 0, 1, 1)
hprior = c(1,1)
# MCMC quantities
Nsamples     = 2e4

# If you want to run the MCMC otherwise load the results
run_MCMC = T
set.seed(0)
# MCMC
if(run_MCMC){
  chain_wine = mcmc_DPM_norm_norm_multi(y, 
                                        hyper   = hyper, 
                                        clus    = FALSE,
                                        hprior  = hprior,
                                        totiter = Nsamples)
  save(chain_wine, file="Data-and-Results/chain_wine_app.Rdata")
} else {
  # Load results
  load("Data-and-Results/chain_wine_app.Rdata")
}
burnin     = 5000
thin       = 1
ps_samples = chain_wine[seq(from=burnin+1, to=Nsamples, by=thin),]
N_ps       = nrow(ps_samples)

# Compute the optimal partition under VI loss (default a=1)
VI_wine     = salso(ps_samples, loss=VI(), nCores = 4)
table(VI_wine)

# Compute the optimal partition under Binder loss (default a=1)
Binder_wine     = salso(ps_samples, loss=binder(), nCores = 4)
table(Binder_wine)

# Binder and VI produce different results
sum(as.integer(VI_wine)!=as.integer(Binder_wine))

# Heatmap of the optimal Binder clustering (with package "salso")
# Compute the adjacency matrix associated with the Binder point estimate
Binder_opt_adj = psm(Binder_wine, nCores = 1)
if(Save_Plot){CairoPNG(filename ='Image/binderwine_app.png', 
                       width = 500, height = 500)}
heatmap(Binder_opt_adj, Rowv=NA, Colv=NA, scale='none', col=viridis(2)[2:1])
if(Save_Plot){invisible(dev.off())}
# Compute the adjacency matrix associated with the VI point estimate
VI_opt_adj = psm(VI_wine, nCores = 1)
if(Save_Plot){CairoPNG(filename ='Image/VI_wine_app.png', 
                       width = 500, height = 500)}
heatmap(VI_opt_adj, Rowv=NA, Colv=NA, scale='none', col=viridis(2)[2:1])
if(Save_Plot){invisible(dev.off())}

# Importance re-sampling for regularized posterior summaries lambda = 50
reg_chain_50 = entropy_regularization(ps_samples, 50)

# Compute the optimal regularized (lambda=50) partition 
# under VI loss (default a=1)
VI_wine_50     = salso(reg_chain_50, loss=VI(), nCores = 4)
table(VI_wine_50)

# Compute the optimal regularized (lambda=50) partition 
# under Binder loss (default a=1)
Binder_wine_50     = salso(reg_chain_50, loss=binder(), nCores = 4)
table(Binder_wine_50)

# Heatmap of the optimal regularized Binder clustering (lambda=50)
# Compute the adjacency matrix of the reg (lambda=50) Binder point estimate
Binder_wine_opt_adj_50 = psm(Binder_wine_50, nCores = 1)
if(Save_Plot){CairoPNG(filename ='Image/correctedbinder_app.png', 
                       width = 500, height = 500)}
heatmap(Binder_wine_opt_adj_50, Rowv=NA, Colv=NA, 
        scale='none', col=viridis(2)[2:1])
if(Save_Plot){invisible(dev.off())}

# Optimal regularized (lambda=50) VI and Binder partition are different
sum(as.integer(VI_wine_50)!=as.integer(Binder_wine_50))
# Compute the adjacency matrix of the reg (lambda=50) Binder point estimate

if(Save_Plot){CairoPNG(filename ='Image/binderwine_reorder_app.png', 
                       width = 500, height = 500)}
heatmap(psm(sort(Binder_wine)), Rowv = NA, Colv = NA, scale='none',
        col = viridis(2)[2:1])
if(Save_Plot){invisible(dev.off())}

if(Save_Plot){CairoPNG(filename ='Image/correctedbinder_reorder_app.png', 
                       width = 500, height = 500)}
heatmap(psm(sort(Binder_wine_50)), Rowv = NA, Colv = NA, scale='none',
        col = viridis(2)[2:1])
if(Save_Plot){invisible(dev.off())}

if(Save_Plot){CairoPNG(filename ='Image/VI_wine_reorder_app.png', 
                       width = 500, height = 500)}
heatmap(psm(sort(VI_wine)), Rowv = NA, Colv = NA, scale='none',
        col = viridis(2)[2:1])
if(Save_Plot){invisible(dev.off())}
if(Save_Plot){CairoPNG(filename ='Image/corrected_VI_wine_reorder_app.png', 
                       width = 500, height = 500)}
heatmap(psm(sort(VI_wine_50)), Rowv = NA, Colv = NA, scale='none',
        col = viridis(2)[2:1])
if(Save_Plot){invisible(dev.off())}

tab_VI     = table(wine[,1], VI_wine)    
sum(tab_VI) - sum(sort(tab_VI, decreasing=T)[1:3])

tab_VI_reg = table(wine[,1], VI_wine_50)
sum(tab_VI_reg) - sum(sort(tab_VI_reg, decreasing=T)[1:3])

tab_Binder    = table(wine[,1], Binder_wine) 
sum(tab_Binder) - sum(sort(tab_Binder, decreasing=T)[1:3])

tab_Binder_reg    = table(wine[,1], Binder_wine_50) 
sum(tab_Binder_reg) - sum(sort(tab_Binder_reg, decreasing=T)[1:3])
table_Binder = prop.table(table(Binder_wine))[order(table(Binder_wine), 
                                                    decreasing = TRUE)]
names(table_Binder) = sort(names(table_Binder))

if(Save_Plot){CairoPNG(filename ='Image/clustnocorr_app.png', 
                       width = 800, height = 600)}
par(mar = c(5,5,2,2))
barplot(table_Binder, xlab = "clusters", ylab = "frequency", col =viridis(4)[3],
        main ="",family = "serif", font = 1, font.lab = 2, cex.lab = 2.5, 
        cex.axis = 2.2)
if(Save_Plot){invisible(dev.off())}
table_Binder_50 = prop.table(table(Binder_wine_50))[
  order(table(table(Binder_wine_50)), decreasing = TRUE)]
names(table_Binder_50) = sort(names(table_Binder_50))
if(Save_Plot){CairoPNG(filename ='Image/cluster_corr_app.png', 
                       width = 800, height = 600)}
par(mar = c(5,5,2,2))
barplot(table_Binder_50, xlab = "clusters", ylab = "frequency", 
        col=viridis(4)[3], main ="", family = "serif", font = 1, font.lab = 2, 
        cex.lab = 2.5, cex.axis = 2.2)
if(Save_Plot){invisible(dev.off())}
table_VI = prop.table(table(VI_wine))[order(table(VI_wine), 
                                                    decreasing = TRUE)]
names(table_VI) = sort(names(table_VI))

if(Save_Plot){CairoPNG(filename ='Image/VI_clustnocorr_app.png', 
                       width = 800, height = 600)}
par(mar = c(5,5,2,2))
barplot(table_VI, xlab = "clusters", ylab = "frequency", col =viridis(4)[3],
        main ="",family = "serif", font = 1, font.lab = 2, cex.lab = 2.5, 
        cex.axis = 2.2)
if(Save_Plot){invisible(dev.off())}
table_VI_50 = prop.table(table(VI_wine_50))[order(table(VI_wine_50), 
                                                  decreasing = TRUE)]
names(table_VI_50) = sort(names(VI_wine_50))

if(Save_Plot){CairoPNG(filename ='Image/VI_cluster_corr_app.png', 
                       width = 800, height = 600)}
par(mar = c(5,5,2,2))
barplot(table_VI_50, xlab = "clusters", ylab = "frequency", col =viridis(4)[3],
        main ="",family = "serif", font = 1, font.lab = 2, cex.lab = 2.5, 
        cex.axis = 2.2)
if(Save_Plot){invisible(dev.off())}
