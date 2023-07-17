# Entopy-Regularized-Clustering

R codes for entropy regularized point estimates of clustering under random partition models.

**Authors**: [Beatrice Franzolini](https://beatricefranzolini.github.io) and [Giovanni Rebaudo](https://giovannirebaudo.github.io)

#### Overview 
This repository is associated with the article [Fanzolini and Rebaudo (2023) **Entropy regularization in probabilistic clustering.**]()
The key contribution of the paper is outlined below.
 
> [...] we propose a novel Bayesian estimator of the clustering configuration.
The proposed estimator is equivalent to a post-processing procedure that reduces the number of sparsely-populated clusters and enhances interpretability. 
The procedure takes the form of entropy-regularization of the Bayesian estimate. 
While being computationally convenient with respect to alternative strategies, it is also theoretically justified as a correction to the Bayesian loss function used for point estimation and, as such, can be applied to any posterior distribution of clusters, regardless of the specific model used.

This repository provides codes to replicate the results in Fanzolini and Rebaudo (2023) **Entropy regularization in probabilistic clustering.** *Statistical Methods & Applications*, in press.

More precisely, we provide the `R` code to implement **Algorithms 1 [Fanzolini and Rebaudo (2023+)]()** to obtain **entropy regularized point estimates** of the clustering in random partition models.

The repository contains the following:

1. `ERC_main.R` code to reproduce the main results in the article;
2. `ERC_fcts.R` functions needed to run the main code;
3. `ERC_cost.R` code to compute and plot the cost functions shown in the articles;
4. `Data-and-Results` folder with data and results of the analyses.

#### Questions or bugs
For bug reporting purposes, e-mail [Beatrice Franzolini](https://beatricefranzolini.github.io) (franzolini@pm.me).

#### Citation
Please cite the following publication if you use this repository in your research: Franzolini, B. and Rebaudo, G. (2023). [Entropy regularization in probabilistic clustering](). *Statistical Methods & Applications*, in press.




