# Functions


# Simulate data from normal, t-stud or Bernoulli clusters
simdat <- function(n=200, Kn=3, mu = c(-4, 0, 4), sigma = rep(1, Kn), df = 10,
                   kernel="Gauss", p = 1){
  y = data.frame()
  c_truth = data.frame()
  nc = c(1, rep(floor(n/Kn), Kn - 1), n - (Kn-1) *  floor(n/Kn)) 
  
  if(kernel == "Gauss" & p == 1){
    for(c in 1:Kn){
      a = sum(c(nc[1:c])); b = sum(nc[1:(c+1)]) - 1
      y[a:b, 1] = rnorm(nc[c+1], mu[c], sigma[c])
      c_truth[a:b, 1] = rep(c, nc[c+1])
    }
  }
  
  if(kernel == "Gauss" & p > 1){
    for(c in 1:Kn){
      a = sum(c(nc[1:c])); b = sum(nc[1:(c+1)]) - 1
      y[a:b, 1:p] = matrix(rnorm(p*nc[c+1], mu[c,], sigma[c,]), ncol = p,
                         byrow = TRUE)
      c_truth[a:b, 1] = rep(c, nc[c+1])
    }
  }
  
  if(kernel =="tstud" & p == 1){
    for(c in 1:Kn){
      a = sum(c(nc[1:c])); b = sum(nc[1:(c+1)]) - 1
      y[a:b, 1] = rt(nc[c+1], df) + mu[c]
      c_truth[a:b, 1] = rep(c, nc[c+1])
    } 
  }
  
  return(list(y = y, c_truth = c_truth))
}

# Compute the marginal likelihood of a univariate norm-norm model 
# (marg lik of a cluster) 
margnn <- function(y, m = 0, s2 = 1, sigma2 = 1, log = TRUE){
  sigma = sqrt(sigma2)
  s = sqrt(s2)
  n = length(y)
  Y2 = sum(y^2)
  Y = sum(y)
  p = - n / 2 * log(2*pi) - (n-1) * log(sigma) - 1/2 * log(n*s2+sigma2) - 
    Y2/(2*sigma2) - m^2/(2*s2) + (s*Y/sigma + sigma*m/s)^2/(2*(n*s2 + sigma2))
  return(p)
} 

# Compute the marginal likelihood of a multivariate norm-norm model 
# (marg lik of a cluster) 
margnn_multi <- function(y, m = 0, s2 = 1, sigma2 = 1, log = TRUE){
  J = dim(y)[2]
  n = dim(y)[1]
  sigma = sqrt(sigma2)
  s = sqrt(s2)
  p = 0
  for (jj in (1:J)){
    Y2 = sum(y[,jj]^2); Y = sum(y[,jj])
    p = p - n/2 * log(2*pi) - (n-1) * log(sigma) - 1/2 * log(n*s2 + sigma2) - 
      Y2/(2*sigma2) - m^2/(2*s2) + (s*Y/sigma + sigma*m/s)^2/(2*(n*s2 + sigma2))
  }
  return(p)
} 


entropy <- function(labels){
  pc = prop.table(table(labels)); k = length(unique(labels))
  entropy = - sum( pc * log( pc, base = k) )
  return(list(entropy = entropy, k= k))
}

DPeppf <- function(nn, alpha = 1){
  k = length(nn)
  logp = - k * log(alpha) - sum(lgamma(nn))
  return(logp)
}

PYeppf <- function(nn, alpha = 1, sigma = 0.5){
  k = length(nn)
  logp = - sum( log(alpha + sigma * seq(1,k)) ) - 
    sum(lgamma(nn - sigma)) + k * lgamma(1- sigma)
  return(logp)
}

num_part <- function(nn, k){
  temp = factorial(apply(nn, 1, sum))
  for(kk in 1:k){
    temp = temp / factorial(nn[,kk])
  }
  return(temp)
}

# Compute all the marginal likelihood of a norm-norm model 
# (marg lik of a new cluster) univariate
marg_new_all <- function(y, m = 0, s2 = 1, sigma2 = 1, log = TRUE){
  sigma = sqrt(sigma2); s = sqrt(s2)
  p = - 1 / 2 * log(2 * pi) - 1/2 * log(s2 + sigma2) - y^2 / (2 * sigma2) - 
    m^2 / (2 * s2) + (s * y / sigma + sigma * m / s ) ^ 2 /(2*(s2+sigma2))
  return(p)
} 
# Compute all the marginal likelihood of a norm-norm model 
# (marg lik of a new cluster) multivariate
marg_new_all_multi <- function(y, m = 0, s2 = 1, sigma2 = 1, log = TRUE){
  J = dim(y)[2]; sigma = sqrt(sigma2); s = sqrt(s2)
  p = - J / 2 * log(2 * pi) - 
    J/2 * log(s2 + sigma2) - 
    rowSums(y^2) / (2 * sigma2) - 
    J*m^2 / (2 * s2) + 
    rowSums((s * y / sigma + sigma * m / s )^2) / (2 * (s2 + sigma2))
  return(p)
}

# Compute all univariate marginal lik in beta bern
marg_new_all_bern <- function(y, a, b, log = TRUE){
  # rowSums(lbeta(a + y, b + 1 - y) - lbeta(a, b))
  p=ncol(y)
  rowSums(ifelse(y==0, log(b), log(a)))-p*log(a+b)
}

mcmc_DPM_norm_norm <- function(y, 
                              hyper   = c(1, 0, 1, 1), 
                              clus    = FALSE, 
                              hprior  = FALSE,
                              totiter = Nsamples,
                              verbose = max(round(Nsamples/20),1)){
  alpha = hyper[1]
  m = hyper[2]
  s2 = hyper[3]
  sigma2 = hyper[4]
  sigma = sqrt(sigma2)
  s = sqrt(s2)
  c_saved = matrix(nrow=totiter,ncol=n)
  mnew    = marg_new_all(y)
  if(!clus[1]){clus = kmeans(y, centers = floor(n/10), nstart = 20)$cluster}
  print(paste("initialization to", length(unique(clus)), "clusters"))
  for (iter in (1:totiter)){
    for (i in (1:n)){
      clus[i] = NA
      c_unique = unique(clus[!is.na(clus)])
      prob = double(length(c_unique)+1) 
      j = 1
      for (cc in c_unique){
        which_cc = which(clus == cc & !is.na(clus))
        n_cc     = length(which_cc) 
        y_cc    = sum(y[which_cc])
        prob[j] = 
          # difference of marginal liks in cc
          - log(2*pi)/2 - 
          log(sigma) - 
          (log(((n_cc+1)*s2+sigma2)/(n_cc*s2+sigma2)))/2 - 
          y[i]^2/(2*sigma2) +
          (s*(y_cc+y[i])/sigma+sigma*m/s)^2/(2*((n_cc+1)*s2+sigma2)) -
          (s*y_cc/sigma+sigma*m/s)^2/(2*(n_cc*s2+sigma2)) +
          # prior part of marg lik
          log(n_cc)
        j = j + 1
      } 
      prob[j] = mnew[i] + log(alpha)
      prob    = exp(prob - max(prob))
      clus[i] = sample(c(c_unique, setdiff(1:n,c_unique)[1]), 1, prob = prob)
    }
    c_saved[iter,] = clus
    if(hprior[1]){
      k = length(unique(clus))
      eta  = rbeta(1, alpha+1, n)
      zeta = sample(c(0,1), 1, 
                    prob = c(hprior[1] + k + 1, n*(hprior[2] - log(eta))))
      alpha = rgamma(1, hprior[1] + k - zeta, hprior[2] - log(eta) )
    }
    if(iter%%verbose==0){
      print(iter)
      #print(table(clus))
      print(alpha)
    }
  }
  return(c_saved)
}

entropy_regularization = function(chain, lambda , plots=F){
  n_samp = dim(chain)[1]
  weights = double(n_samp)
  for (iter in (1:n_samp)){
    labels = chain[iter,]
    pc = prop.table(table(labels))
    entropy = - sum( pc * log( pc, base = length(pc)) )
    weights[iter] = exp(lambda*entropy)
  }
  which = sample.int(n_samp, replace = TRUE, prob = weights)
  if (plots){
    barplot(weights/sum(weights))
    print(weights/sum(weights))
  }
  return(chain[which,])
}

mcmc_DPM_norm_norm_multi <- function(y, 
                                     hyper=c(0.1, 0, 1, 1), 
                                     clus = FALSE,
                                     hprior = FALSE,
                                     totiter = Nsamples,
                                     verbose = max(round(Nsamples/20),1)){
  alpha  = hyper[1]
  m      = hyper[2]
  s2     = hyper[3]
  sigma2 = hyper[4]
  s      = sqrt(s2)
  sigma  = sqrt(sigma2)
  J      = ncol(y)
  if(!clus[1]){clus = kmeans(y, centers = floor(n/10), nstart = 20)$cluster}
  print(paste("initialization to", length(unique(clus)), "clusters"))
  c_saved = matrix(nrow=totiter,ncol=n)
  mnew = marg_new_all_multi(y, m, s2, sigma2, log = TRUE)
  for (iter in 1:totiter){
    for (i in 1:n){
      clus[i] = NA
      c_unique = unique(clus[!is.na(clus)])
      prob = double(length(c_unique)+1) 
      y_i  = y[i,]
      j    = 1
      temp = 0
      for (cc in c_unique){
        which_cc = which(clus == cc & !is.na(clus))
        n_cc     = length(which_cc)
        # difference of marginal liks in cc
        temp     = - J / 2 * log(2 * pi) -
          J * log(sigma) - 
          sum(y_i^2)/ (2 * sigma2) +
          log(n_cc)  +
          J/2 * log( (n_cc *s2 + sigma2)/((n_cc+1) * s2 + sigma2))
        for (jj in 1:J){
          Y_cc_jj = sum(y[which_cc,jj])
          temp = temp + 
            (s * (Y_cc_jj+y_i[jj])/sigma+sigma*m/s)^2 / (2*((n_cc+1)*s2+sigma2)) -
            (s * Y_cc_jj/sigma+sigma*m/s)^2 / (2*(n_cc*s2+sigma2))
        }
        prob[j] = temp
        j = j + 1
      } 
      # new cc
      prob[j] = mnew[i] + log(alpha)
      prob    = exp(prob - max(prob))
      clus[i] = sample(c(c_unique, setdiff(1:n,c_unique)[1]), 1, prob = prob)
    }
    c_saved[iter,] = clus
    if(hprior[1]){
      k = length(unique(clus))
      eta  = rbeta(1, alpha+1, n)
      zeta = sample(c(0,1), 1, 
                    prob = c(hprior[1] + k + 1, n*(hprior[2] - log(eta))))
      alpha = rgamma(1, hprior[1] + k - zeta, hprior[2] - log(eta) )
    }
    if(iter%%verbose==0){
      print(iter)
    }
  }
  return(c_saved)
}


mcmc_DPM_bern_beta <- function(y, 
                               hyper=c(0.1, 2, 2), 
                               clus = FALSE,
                               hprior = FALSE,
                               totiter = Nsamples,
                               verbose = max(round(Nsamples/20),1)){
  alpha = hyper[1]
  a     = hyper[2]
  b     = hyper[3]
  J     = ncol(y)
  if(!clus[1]){clus = kmeans(y, centers = min(floor(n/10), 50),
                             nstart = 20)$cluster}
  print(paste("initialization to", length(unique(clus)), "clusters"))
  c_saved = matrix(nrow=totiter,ncol=n)
  mnew = marg_new_all_bern(y, a, b, log = TRUE)
  for (iter in 1:totiter){
    for (i in 1:n){
      clus[i] = NA
      c_unique = unique(clus[!is.na(clus)])
      prob = double(length(c_unique)+1) 
      y_i  = y[i,]
      j    = 1
      temp = 0
      for (cc in c_unique){
        which_cc = which(clus == cc & !is.na(clus))
        n_cc     = length(which_cc)
        # difference of marginal liks in cc
        n1 = colSums(y[which_cc,,drop = FALSE])
        n0 = n_cc - n1
        temp  = sum(-lbeta(a + n1, b + n_cc - n1) + 
                      lbeta(a + n1 + y[i,], b + (n_cc + 1) - n1 - y[i,])) +
          log(n_cc)
        prob[j] = temp
        j = j + 1
      } 
      # new cc
      prob[j] = mnew[i] + log(alpha)
      prob    = exp(prob - max(prob))
      clus[i] = sample(c(c_unique, setdiff(1:n,c_unique)[1]), 1, prob = prob)
    }
    c_saved[iter,] = clus
    if(hprior[1]){
      k = length(unique(clus))
      eta  = rbeta(1, alpha+1, n)
      zeta = sample(c(0,1), 1, 
                    prob = c(hprior[1] + k + 1, n*(hprior[2] - log(eta))))
      alpha = rgamma(1, hprior[1] + k - zeta, hprior[2] - log(eta) )
    }
    if(iter%%verbose==0){
      print(iter)
      print(alpha)
      print(table(clus))
    }
  }
  return(c_saved)
}

mcmc_DPM_bern_beta_edited <- function(y, 
                                      hyper=c(0.1, 2, 2), 
                                      clus = FALSE,
                                      hprior = FALSE,
                                      totiter = Nsamples,
                                      verbose = max(round(Nsamples/20),1)){
  alpha = hyper[1]
  a     = hyper[2]
  b     = hyper[3]
  J     = ncol(y)
  n     = nrow(y)
  if(!clus[1]){clus = stats::kmeans(y, centers = min(floor(n/10), 50),
                             nstart = 20)$cluster}
  print(paste("initialization to", length(unique(clus)), "clusters"))
  c_saved = matrix(nrow=totiter,ncol=n)
  mnew = marg_new_all_bern(y, a, b, log = TRUE)
  for (iter in 1:totiter){
    for (i in 1:n){
      clus[i] = NA
      c_unique = unique(clus[!is.na(clus)])
      prob = double(length(c_unique)+1) 
      y_i  = y[i,]
      j    = 1
      for (cc in c_unique){
        which_cc = which(clus == cc & !is.na(clus))
        n_cc     = length(which_cc)
        # difference of marginal liks in cc
        n1 = colSums(y[which_cc,,drop=F])
        prob[j] = sum(
          ifelse(y_i==0, log(b + n_cc - n1), log(a + n1)) ) - 
          J*log(a + b + n_cc) + log(n_cc)
        j = j + 1
      } 
      # new cc
      prob[j] = mnew[i] + log(alpha)
      prob    = exp(prob - max(prob))
      clus[i] = sample(c(c_unique, setdiff(1:n,c_unique)[1]), 1, prob = prob)
    }
    c_saved[iter,] = clus
    if(hprior[1]){
      k = length(unique(clus))
      eta  = rbeta(1, alpha+1, n)
      zeta = sample(c(0,1), 1, 
                    prob = c(hprior[1] + k + 1, n*(hprior[2] - log(eta))))
      alpha = rgamma(1, hprior[1] + k - zeta, hprior[2] - log(eta) )
    }
    if(iter%%verbose==0){
      print(iter)
      print(table(clus))
    }
  }
  return(c_saved)
}

# Simulate data from spherical multivariate normal clusters
simdat_mult <- function(freq, dim=1, mu = c(-4, 0, 4), sigma = rep(1,3)){
  y  = matrix(nrow=n, ncol=J)
  Kn = length(freq)
  c_truth = double()
  for(c in 1:Kn){
    if(c>1){
      a = sum(c(freq[1:(c-1)]))+1
    } else {
      a=1
    }
    b = sum(freq[1:c])
    y[a:b, 1:dim] = rnorm(freq[c]*dim, mu[c], sigma[c])
    c_truth[a:b]  = c
  } 
  return(list(y = y, c_truth = c_truth))
}
