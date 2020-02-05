Het-SBM <- function(X,D,qmin,qmax, t_max = 10, h_max = 10,threshold_h = 10^-10,
                     tau_min = 10^-10, threshold_lik = 10^-10, M_iter = 10,
                     threshold_psi = 10^-10, startType = "KMeans",
                     kmeanType = "correlation", kmeanMax = 30,startZ = NULL,
                     iSubjStartingPoint = NULL, method = "Firth"){
  n   = dim(X)[1]  # Number of nodes
  K   = dim(X)[3]  # Number of subjects
  P   = dim(D)[2]  # Number of covariates
  # h is concerened with update of tau
  # t is concerened with maximisation of all other parameters
  rezultat =  list()
  for(Q in qmin:qmax){
    convergence  = FALSE
    fail_rate    = 0
    varBound     = numeric(t_max)
    while (!convergence & fail_rate < 1 ){	
      # Initialise tau or update tau
      t         = 1
      delta_psi = 1
      print(paste("Model_Q = ",Q))
      while(t <= t_max & delta_psi > threshold_psi){
        print(paste("t_iteration=",t))
        if(t == 1){
          X_all = apply(X,c(1,2),function(x)sum(x))      # Sum data over subjects
          if(startType == "KMeans"){
            ##############################################
            # Initialise Z, using the k-means algorithm  #
            ##############################################
            attempts = 0
            while(attempts < 30){
                Z_all  = tryCatch(expr = {amap::Kmeans(X_all,Q,method=kmeanType,iter.max = kmeanMax + attempts)$cluster},
                                error = function(e) e,
                                warning = function(w) w)
              if(inherits(Z_all, "error") || inherits(Z_all, "warning")){
                attempts = attempts + 1
              }else{
                break()
              }
            }
          }
          if(startType == "StartingPoint"){
            Z_all =  startZ[[Q]]
          }
          if(startType == "Random"){
            #######################################################################
            #                   Initialise Z, using random point                  #
            #######################################################################
            attempts2 = 0
            while(attempts2<11){
              Z_all = sample(1:Q,n,replace=TRUE)
              Z_nqs = as.data.frame(table(Z_all))$Freq
              if(any(Z_nqs<3)){
                attempts2 = attempts2 + 1
                Z_all     = sample(1:Q,n,replace=TRUE)        	      
              }else{
                break()
              }
            }
          }
          # if(startType=="Hclust"){
          #   #######################################################################
          #   #      Initialise Z, using Hieararchical Clustering algorithm         #
          #   #######################################################################
          # X_all_d=dist(X_all,method="manhattan")
          # tmp=hclust(X_all_d,method="ward.D2")
          # Z_all=cutree(tmp,k=Q)
          # }
          #if(startType=="Hclust"){
            #######################################################################
            #      Initialise Z, using Hieararchical Clustering algorithm         #
            #######################################################################
          #   if (is.null(iSubjStartingPoint)){
          #     X_all_d = dist(X_all,method="manhattan")
          #   }else{
          #     X_all_d = dist(X[,,iSubjStartingPoint],method="manhattan")
          #   }
          #   tmp = hclust(X_all_d,method="ward.D2")
          #   Z_all = cutree(tmp,k=Q)
          # }
          if(startType == "Hclust"){
            #######################################################################
            #      Initialise Z, using Hieararchical Clustering algorithm         #
            #######################################################################
            if (is.null(iSubjStartingPoint)){
              X_all_d = stats::dist(X_all,method="manhattan")
            } else {
              if (length(iSubjStartingPoint) == 1) {
                X_all_d = stats::dist(X[,,iSubjStartingPoint],method="manhattan")
              } else {
                X_all_d = stats::dist(apply(X[,,iSubjStartingPoint],c(1,2),function(x)sum(x)),method="manhattan")
              }
            }
            tmp   = stats::hclust(X_all_d, method ="ward.D2")
            Z_all = stats::cutree(tmp, k = Q)
          }
          # Build tau using the first guess of Z
          tau = matrix(tau_min,n,Q)              #populate tau with minimum values we cannot deal with zero
          for(q in 1:Q){
            tau[Z_all==q,q] = 1-(Q-1)*tau_min    #normalise tau
          }
        }
        else{ # for all t that are >=2 we eneter to optimise tau
          h          = 1
          delta_h    = 1
          monoBlock  = which(sapply(1:Q, function(x) length(which(apply(tau, 1, which.max)==x))) == 1) # Check for a block with one node
          
          # do stuff not affected by upadtes of tau
          tmp1 = array(0, dim = c(K, Q, Q))
          for (l in 1:Q){									
            tmp1[,,l]  =  tcrossprod(D,beta[,l,])
          }
          # logOnePlusExpTmp1 <- log(1+exp(tmp1))  
          sumLogOnePlusExpTmp1 = apply(log(1+exp(tmp1)), 2:3, sum)
          Xtmp1                = array(0, dim = c(n,n*Q,Q))
          for(i in 1:n){
            for(q in 1:Q){
              Xtmp1[i,,q] =  as.numeric(X[i,,] %*% tmp1[,,q])
            }
          }          
          # Update Taus
          while(h <= h_max & delta_h > threshold_h){
            tau_previous =  tau	 					
            for(i in 1:n){
              posTerm                       = crossprod(as.numeric(tau), Xtmp1[i,,])
              negTerm                       = tcrossprod(colSums(tau[-i,]), sumLogOnePlusExpTmp1)
              tmp                           = log(alpha) + posTerm - negTerm
              tau[i,]                       =  1/colSums(exp(outer(as.numeric(tmp), as.numeric(tmp), FUN="-")))
              tau[i,tau[i,]< tau_min]       = tau_min                                 # smaller than smallest is set to our view of small :)
              tau[i,which.max(tau[i,])[1]]  = (1-sum(tau[i,-which.max(tau[i,])[1]]))  # to keep normalization in check!
            }
            delta_h      = max(abs(tau-tau_previous)/tau_previous)   
            # print(delta_h)
            h            = h + 1	
          }# Here we update tau					
        }# for all t that are >=2 we eneter to optimised tau
        monoBlock =  which(sapply(1:Q, function(x) length(which(apply(tau, 1, which.max)==x))) == 1) # Check for the monoblock
        if(t > 1) {
          alpha_previous =  alpha
        }  
        alpha     =  apply(tau,2,sum)/n
        # print(alpha)
        if (t == 1){
          Emp_Pis = matrix(0,Q,Q)  # Emp_Pis are initail guess for intercept. 
          for(q in 1:Q){
            for(l in q:Q){
              privr1 = as.matrix(X_all[Z_all==q,Z_all==l])
              if(q == l){
                if(dim(privr1)[1]==1){
                  Emp_Pis[q,l] = sum(X_all[Z_all==q,])/(K*(n-1))
                }
                else{
                  Emp_Pis[q,l] = sum(privr1)/(K*dim(privr1)[1]*(dim(privr1)[2]-1))	
                }
              }
              else{
                Emp_Pis[q,l]   =  sum(privr1)/(K*dim(privr1)[1]*dim(privr1)[2])
              }
            }
          }
          Emp_Pis[Emp_Pis == 1]                     =  0.9
          beta                                      =  array(0,dim=c(Q,Q,P))
          beta_init                                 =  matrix(0,Q,Q)
          beta_init[upper.tri(Emp_Pis,diag = TRUE)] = log(Emp_Pis[upper.tri(Emp_Pis,diag=TRUE)]/(1-Emp_Pis[upper.tri(Emp_Pis,diag=TRUE)]))
          beta_init[beta_init == -Inf]                = 0
          beta[,,1]                                 = beta_init
        }
        beta_previous <- beta
        if(method == "Firth"){
          tmp   = firth_bs(beta, tau, D,X,threshold_lik=10^-10,M_iter=10,maxstep=5,half_iter=5) 
        }
        if(method == "MLE"){
          tmp   = no_firth_bs(beta, tau, D,X,threshold_lik=10^-10,M_iter=10,maxstep=5,half_iter=5)
        }
        beta = tmp$beta
        FI   = tmp$FI
        if (t>1) {
          delta_psi  =  max(abs((alpha-alpha_previous)/alpha_previous),abs((beta-beta_previous)/beta_previous), na.rm = TRUE)
        }  
        varBound[t]  =  sum(tau%*%log(alpha))-sum(tau*log(tau))
        for (q in 1:Q){
          for (l in q:Q){
            if(q == l){
              gamma = tcrossprod(tau[,q])
            }else{
              gamma = tcrossprod(tau[,q],tau[,l])
              gamma = gamma + t(gamma)    
            }
            # if (q==l & any(l== monoBlock)){
            #   
            # } else {
              tmp1        =  D%*%beta[q,l,]			
              varBound[t] = varBound[t] + 0.5*(crossprod(apply(X,3, function(x) sum(gamma[x==1])), tmp1) - (sum(gamma)-sum(diag(gamma))) * sum(log(1+exp(tmp1))))
            # }					
          }
        }
        t =  t +1
      }#end of while loop for t
      Qlik = 2 * sum(tau%*%log(alpha))
      for (q in 1:Q){
        for (l in q:Q){
          if(q == l){
            gamma = tcrossprod(tau[,q])
          }else{
            gamma = tcrossprod(tau[,q],tau[,l])
            gamma = gamma + t(gamma)    
          }
            tmp1 =  D %*% beta[q,l,]			
            Qlik =  Qlik + crossprod(apply(X,3, function(x) sum(gamma[x==1])), tmp1) - (sum(gamma)-sum(diag(gamma))) * sum(log(1+exp(tmp1)))
        }
      }
      Qlik    = 0.5*Qlik
      ICL     = Qlik - P*Q*(Q+1)/4*log(K*n*(n-1)/2) - (Q-1)/2*log(n)
      part    = apply(tau,1,function(x)which.max(x))
      part_sz = sapply(1:Q,function(x)sum(part==x))
      if(any(part_sz < 1)){
        fail_rate = fail_rate + 1
      }else{
        convergence = TRUE
      }
    }
    if(!convergence){
      warning(paste("Het-SBM with ", Q, " classes: no convergence after ", fail_rate," trials",sep="" ))
    }
    if(length(monoBlock) >0) warning(paste("ERMM with ", Q, "classes: presence of mono-blocks",sep="" ))
    rezultat[[Q]]             = list()
    rezultat[[Q]]$tau         = tau
    rezultat[[Q]]$alpha       = alpha
    rezultat[[Q]]$beta        = beta
    rezultat[[Q]]$ICL         = ICL
    rezultat[[Q]]$VB          = varBound
    rezultat[[Q]]$FI          = FI
    rezultat[[Q]]$convergence = convergence
    rezultat[[Q]]$monoBlock   = monoBlock
  }
  return(rezultat)	
}