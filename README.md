# Het-SBM
This project is providing fitting procedures for the Het-SBM model described in the thesis ['Generalised Stochastic Blockmodels and their Applications in the Analysis of Brain Networks' ](https://core.ac.uk/download/pdf/42619639.pdf) Chapter 4, Pavlovic 2015 (see [the link to the paper](https://www.biorxiv.org/content/10.1101/672071v1.abstract)  based on this chapter). 
# Variables Specified by User
In this section, we give a list of variables supplied by a user.

X is a 3-dimensional array which encodes adjacency matrices for each subject. The dimensions of X are nodes (n) by nodes (n) by subjects (K), where n stands for a total number of nodes in a network and K stands for a total number of subjects in the study.  For example, X[ , ,1] is a square adjacency matrix corresponding to the first subject. 

D is a 2-dimensional array (matrix) which encodes covariates of interest (i.e., design matrix). Design matrix D is given as subjects (K) by covariates (P). It is worth mentioning that P is also counting intercept in a linear model as well as covariates of interest like age, gender, patient/control etc. 

Het-SBM fits data over a range of candidate models or candidate number of clusters. These are specified by parameters 'qmin' (indicating the smallest number of clusters in your data) and 'qmax' (indicating the largest number of clusters in your data). By setting qmin=2 and qmax=5, the function will fit 4 models, so that the first model has 2 clusters, second 3 clusters, third 4 clusters and fourth 5 clusters. As noted in the referenced work the model which maximises ICL criterion is the final. 


Next, we describe parameters related to iterations:  t_max is the maximum number of steps in the M-step (default value is 10), h_max  is the maximum number of steps in the E-step (default value is 10). Convergence criteria for the  update of tau is given by threshold_h (default value is 10^-10), while convergence criteria for beta and alpha is given by threshold_psi (default value is set to 10^-10). 

The estimates of beta are handled internally by a function firth_bs(). Convergence criteria for beta is given by threshold_lik (default value is 10^-10), maximum number of steps in a logistic regression is given by M_iter (default value is 10), maxstep is the maximal absolute variation in beta value (default value is 5), step-halving parameter is half_iter (default value is 5). 
