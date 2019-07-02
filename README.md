# Het-SBM
This project is providing fitting procedures for the Het-SBM model described in the thesis ['Generalised Stochastic Blockmodels and their Applications in the Analysis of Brain Networks' ](https://core.ac.uk/download/pdf/42619639.pdf) Chapter 4, Pavlovic 2015 (see [the link to the paper](https://www.biorxiv.org/content/10.1101/672071v1.abstract)  based on this chapter). 
# Variables Specified by User
In this section, we give a list of variables supplied by a user.

X is a 3-dimensional array which encodes adjacency matrices for each subject. The dimensions of X are nodes (n) by nodes (n) by subjects (K), where n stands for a total number of nodes in a network and K stands for a total number of subjects in the study.  For example, X[ , ,1] is a square adjacency matrix corresponding to the first subject. 

D is a 2-dimensional array (matrix) which encodes covariates of interest (i.e., design matrix). Design matrix D is given as subjects (K) by covariates (P). It is worth mentioning that P is also counting intercept in a linear model as well as covariates of interest like age, gender, patient/control etc. 

Het-SBM fits data over a range of candidate models or candidate number of clusters. These are specified by parameters 'qmin' (indicating the smallest number of clusters in your data) and 'qmax' (indicating the largest number of clusters in your data). By setting qmin=2 and qmax=5, the function will fit 4 models, so that the first model has 2 clusters, second 3 clusters, third 4 clusters and fourth 5 clusters. As noted in the referenced work the model which maximises ICL criterion is the final. 


Next, we describe parameters related to iterations:  t_max is the maximum number of steps in the M-step (default value is 10), h_max  is the maximum number of steps in the E-step (default value is 10). The convergence criterion for the  update of taus is given by threshold_h (default value is 10^-10), while the convergence criterion for betas and alphas is given by threshold_psi (default value is set to 10^-10). 

To fit regression coefficients (betas) within M-step, the user will have two choices. 
- By setting method = "Firth", the user will obtain Firth type estimates for betas using the function firth_bs().
- By setting method = "MLE", the user will obtain Maximum Likelihood estimates (MLE) using the function no_firth_bs().
In both options, the convergence criterion for betas  is given by the threshold_lik (default value is 10^-10), the maximum number of steps in a logistic regression fitting is given by M_iter (default value is 10), the maximal absolute variation in beta value is given by the maxstep (default value is 5) while the step-halving parameter is half_iter (default value is 5).

Initialisation: There are 4 main initialisation strategies which can be chosen using the startType parameter. 

1. When startType = "KMeans", the starting point is based on the k-means algorithm from the package 'amap'. Linked to  this option, kmeanType can be used to pass the distance measure for centres which can be "euclidean", "maximum", "manhattan", "canberra", "binary", "pearson", "abspearson", "correlation", "abscorrelation", "spearman" or "kendall" (default value is set to "correlation"), and  kmeanMax which sets the maximum number of iterations (default is 30). 

2. When startType  is set to  "StartingPoint", then the algorithm expects initialisation values for each candidate models. For example with qmin=2 and qmax=3, a variable Z_initialisation is needed as an argument. Z_initialisation is expected to be a list such that Z_initialisation[[1]]= NULL, Z_initialisation[[2]] is a partition vector with 2 clusters and Z_initialisation[[3]] is a partition vector with 3 clusters.

3. When startType  is set to "Random", then the algorithm uses uniform sampling to generate some cluster labels based on the given cluster numbers (i.e. qmin and qmax). 

4. When startType  is set  to  "Hclust", then the algorithm will use hierarchical clustering from stats package with distance metric set to "manhattan", while method is set to "ward.D2". 

It is worth highlighting that for "KMeans" and "Hclust" options it is also possible to pick a specific subject whose data will be used for initialisation. This can be handled by iSubjStartingPoint parameter. For example, iSubjStartingPoint = 10 indicates that the network corresponding to the 10th subject will be used for initialisation. 


