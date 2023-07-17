# The following code produces the following plots: 
# 1. Partition cost of DPM as function of entropy 
# 2. Partition cost of PYPM as function of entropy
# 3. Base cost with normal kernel with two different variances
################################################################################
# Load relevant libraries, functions and data ----------------------------------
rm(list=ls())
# Set the working directory to the current folder 
# Code to set the working directory to the current folder from RStudio
library(rstudioapi) # version 0.14
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(ggplot2)    # version 3.4.2
library(latex2exp)  # version 0.9.6
library(extrafont)  # version 1.0
loadfonts(device = "all")

# Load functions
source("ERC_fcts.R")

# compute all possible vectors of frequencies
n = 100 # set sample size

for (k in c(2, 3, 4)){
  # Compute all vectors of numerosity
  a = 1
  nn = matrix(nrow = n**k, ncol = k)
  if (k == 3){
    for (i in 1:n){
      for (j in 1:n){
        for (l in 1:n){
          nn[a,] = c(i,j,l)
          a = a+1
          }
        }
      }
  } else if (k == 2){
    for (j in 1:n){
      for (l in 1:n){
        nn[a,] = c(j,l)
        a = a+1
      }
    }
  } else if (k == 4){
    for(h in (1:n)){
      for (i in 1:n){
        for (j in 1:n){
          for (l in 1:n){
            nn[a,] = c(h,i,j,l)
            a = a+1
          }
        }
      }
    }
  }
  if (k == 3){
    nn3 = nn
  } else if (k == 2){
    nn2 = nn
  } else if (k == 4){
   nn4 = nn
  }
}

nn2_save = nn2 
nn3_save = nn3
nn4_save = nn4

# compute Dirichlet base cost 
for(k in c(2, 3, 4)){
  if (k == 4){
    nn = nn4_save
  } else if (k == 3){
    nn = nn3_save
  } else if (k == 2){
    nn = nn2_save
  }

  nn = cbind(nn, rowSums(nn))
  nn = nn[nn[,(k+1)]==n,1:k]
  nn_temp = cbind(nn, num_part(nn, k))
  entropy = 0 
  for (kk in 1:k){
    entropy = entropy - (nn[,kk] / n) * log(nn[,kk] / n, base = k) 
  }
  nn = cbind(nn_temp, entropy)
  nn = nn[order(entropy),]
  
  
  nn = cbind(nn, rep(1,dim(nn)[1]))
  for (ll in 1:(dim(nn)[1])){
    nn[ll,(k+3)] = DPeppf(nn[ll,1:k])
  }


  if (k == 4){
    nn4 = data.frame(nn)
    rw = nn4[,k+1] * nn4[,k+3]
    dat1 = aggregate(rw ~ entropy, data = nn4, sum)
    datw = aggregate(nn4[,k+1] ~ entropy, data = nn4, sum)
    nn4 = data.frame(cbind( aggregate(nn4, list(Name=nn4$entropy), mean)[, (k+3)] , dat1$rw / datw[,2]))
  } else if (k == 3){
    nn3 = data.frame(nn)
    rw = nn3[,k+1] * nn3[,k+3]
    dat1 = aggregate(rw ~ entropy, data = nn3, sum)
    datw = aggregate(nn3[,k+1] ~ entropy, data = nn3, sum)
    nn3 = data.frame( cbind( aggregate(nn3, list(Name=nn3$entropy), mean)[, (k+3)] , dat1$rw / datw[,2]))
  } else if (k == 2){
    nn2 = data.frame(nn)
    rw = nn2[,k+1] * nn2[,k+3]
    dat1 = aggregate(rw ~ entropy, data = nn2, sum)
    datw = aggregate(nn2[,k+1] ~ entropy, data = nn2, sum)
    nn2 = data.frame(cbind( aggregate(nn2, list(Name=nn2$entropy), mean)[, (k+3)] , dat1$rw / datw[,2]))
  }
}

# Figure 2 of the manuscript
p = ggplot() + theme_bw() + 
  geom_line(data = nn4, aes(x = X1, y = X2), color = 4, linewidth = 1.2) +
  geom_line(data = nn3, aes(x = X1, y = X2), color = 2, linewidth = 1.2) +
  geom_line(data = nn2, aes(x = X1, y = X2), color = 3, linewidth = 1.2) +
  xlab('Entropy') +
  ylab('Partition Cost') + 
  xlim(0,1) + 
  ylim(-360,-200) +
  theme(text = element_text(size=20, family="serif"))+ 
  annotate(geom="text", x=0.75, y=-245, label=TeX("$K_n = 4$"),
           color=4, size = 6) +
  annotate(geom="text", x=0.75, y=-268, label=TeX("$K_n = 3$"),
           color=2, size = 6) + 
  annotate(geom="text", x=0.75, y=-300, label=TeX("$K_n = 2$"),
           color=3, size = 6) 
p



# compute  Pitman-Yor base cost 
for(k in c(2, 3, 4)){
  if (k == 4){
    nn = nn4_save
  } else if (k == 3){
    nn = nn3_save
  } else if (k == 2){
    nn = nn2_save
  }
  
  nn = cbind(nn, rowSums(nn))
  nn = nn[nn[,(k+1)]==n,1:k]
  nn_temp = cbind(nn , num_part(nn, k))
  entropy = 0 
  for (kk in 1:k){
    entropy = entropy - (nn[,kk] / n) * log(nn[,kk] / n, base = k) 
  }
  nn = cbind(nn_temp, entropy)
  nn = nn[order(entropy),]
  
  
  nn = cbind(nn, rep(1,dim(nn)[1]))
  for (ll in 1:(dim(nn)[1])){
    nn[ll,(k+3)] = PYeppf(nn[ll,1:k])
  }
  
  
  if (k == 4){
    nn4 = data.frame(nn)
    rw = nn4[,k+1] * nn4[,k+3]
    dat1 = aggregate(rw ~ entropy, data = nn4, sum)
    datw = aggregate(nn4[,k+1] ~ entropy, data = nn4, sum)
    nn4 = data.frame(cbind( aggregate(nn4, list(Name=nn4$entropy), mean)[, (k+3)] , dat1$rw / datw[,2]))
  } else if (k == 3){
    nn3 = data.frame(nn)
    rw = nn3[,k+1] * nn3[,k+3]
    dat1 = aggregate(rw ~ entropy, data = nn3, sum)
    datw = aggregate(nn3[,k+1] ~ entropy, data = nn3, sum)
    nn3 = data.frame( cbind( aggregate(nn3, list(Name=nn3$entropy), mean)[, (k+3)] , dat1$rw / datw[,2]))
  } else if (k == 2){
    nn2 = data.frame(nn)
    rw = nn2[,k+1] * nn2[,k+3]
    dat1 = aggregate(rw ~ entropy, data = nn2, sum)
    datw = aggregate(nn2[,k+1] ~ entropy, data = nn2, sum)
    nn2 = data.frame(cbind( aggregate(nn2, list(Name=nn2$entropy), mean)[, (k+3)] , dat1$rw / datw[,2]))
  }
}

# Figure 3 of the manuscript
p = ggplot() + theme_bw() + 
  geom_line(data = nn4, aes(x = X1, y = X2), color = 4, linewidth = 1.2) +
  geom_line(data = nn3, aes(x = X1, y = X2), color = 2, linewidth = 1.2) +
  geom_line(data = nn2, aes(x = X1, y = X2), color = 3, linewidth = 1.2) +
  xlab('Entropy') +
  ylab('Partition Cost') + 
  xlim(0,1) + 
  ylim(-360,-200) +
  theme(text = element_text(size=20, family="serif"))+ 
  annotate(geom="text", x=0.75, y=-241.5, label=TeX("$K_n = 4$"),
           color=4, linewidth = 6) +
  annotate(geom="text", x=0.75, y=-264, label=TeX("$K_n = 3$"),
           color=2, linewidth = 6) + 
  annotate(geom="text", x=0.75, y=-297, label=TeX("$K_n = 2$"),
           color=3, linewidth = 6) 
p

# Heatmaps normal distribution for base cost
library(mnormt) # version 2.1.1
library(fields) # version 14.1
set.seed(0)
x1 <- seq(-3, 3, 0.1)
x2 <- seq(-3, 3, 0.1)
mean <- c(0, 0)

f <- function(x1, x2) -dmnorm(cbind(x1, x2), mean, cov, log = TRUE)

cov <- matrix(c(1, 0, 0, 1), nrow=2)
y1 <- outer(x1, x2, f)

cov <- matrix(c(3, 0, 0, 3), nrow=2)
y2 <- outer(x1, x2, f)

# create contour plot
contour(x1, x2, y1)

# Figure 1 (a) of the manuscript
hm_col_scale<-colorRampPalette(c("white","yellow","red"))(1000)
image.plot(y1,  
      col = hm_col_scale,
      zlim=c(min(y1), max(y1)),axes=F)
axis(1, at =  seq(0,1,1/6), label=seq(-3,3))
axis(2, at =  seq(0,1,1/6), label=seq(-3,3))

# Figure 1 (b) of the manuscript
image.plot(y2,  
      col = hm_col_scale,
      zlim=c(min(y1), max(y1)),axes=F)
axis(1, at =  seq(0,1,1/6), label=seq(-3,3))
axis(2, at =  seq(0,1,1/6), label=seq(-3,3))

# base cost as function of distance 
grid_theta1 = seq(-3,3,0.1)
grid_theta2 = seq(-3,3,0.1)

tab = expand.grid(grid_theta1, grid_theta2)
tab$dist = round(abs(tab$Var1 - tab$Var2),1)

sigma2 = 1
tab$cost = log(2*pi) + log(sigma2) + 
  1/2 * (tab$Var1)^2/sigma2 + 1/2 * (tab$Var2)^2/sigma2

# Figure 1 (c) of the manuscript
p = ggplot() + theme_bw() + 
  geom_point(data = tab, aes(x = dist, y = cost), color = 4, size = 1.2) +
  xlab('distance') +
  ylab('Base Cost') + 
  xlim(0,6) + 
  ylim(0,12) +
  theme(text = element_text(size=20, family="serif"))
p

sigma2is1 = aggregate(tab[,3:4], list(tab$dist), mean)


sigma2 = 3
tab$cost = log(2*pi) + log(sigma2) + 
  1/2 * (tab$Var1)^2/sigma2 + 1/2 * (tab$Var2)^2/sigma2


sigma2is2 = aggregate(tab[,3:4], list(tab$dist), mean)

# Figure 1 (d) of the manuscript
p = ggplot() + theme_bw() + 
  geom_point(data = tab, aes(x = dist, y = cost), color = 2, size = 1.2) +
  xlab('distance') +
  ylab('Base Cost') + 
  xlim(0,6) + 
  ylim(0,12) +
  theme(text = element_text(size=20, family="serif"))
p

