no_firth_bs<-function(beta, tau, D,X,threshold_lik=10^-10,M_iter=10,maxstep=5,half_iter=5){
	# use lik without firth correction
 	# maxstep is the difference in beta updates
	# halfiter the number of half step when update is worst than in the previous step
	P<-dim(D)[2]
	n<-dim(X)[1]
	K<-dim(X)[3]
	Q<-dim(tau)[2]
	
    # Defining objects:
   
    beta[rep(lower.tri(beta[,,1],diag=FALSE),P)]<-0
	beta[is.nan(beta)] <- 0
	
	#Constant terms depending on pi_{qlk}
	
	U_const<-array(1,c(Q,Q,K))	    # pi_{qlk}
	pr0<-array(1,c(Q,Q,K))             # (1-pi_{qlk})  Constant term in likelihood when x_ijk=0 or (1-pi_qlk)
    FI_sqr<-array(0,c(Q,Q,K))       	# pi_{qlk}(1-pi_qlk)
    # dFI_const<-array(0,c(Q,Q,K))	#pi_{qlk} (1-pi_{qlk}) (1-2*pi_{qlk})  Constant term in the derivative of the FI matrix
	
	#Objects to be computed
	
	# U_star<-array(0,c(Q,Q,P)) 	    # Modified score function
	FI<-array(0,c(Q,Q,P,P))		    	# Fisher Information (FI) matrix
	FI_inv<-array(0,c(Q,Q,P,P))	    # Inverse of Fisher Information matrix
	
	# Terms which depend on the current tau
	
# # 	gamma<-array(0,c(n,n,Q,Q))	
	# sum_gamma_without_diag <- matrix(0,Q,Q) # sum gamma but not its diagonal
	# sum_connection <- array(0,c(Q,Q,K))            # sum gamma when x_ijk=1, gamma_ijql*x_{ijk}=gamma_ijql otherwise zero
	# sum_no_connection <- array(0,c(Q,Q,K))      # sum gamma when x_ijk=0, gamma_ijql*x_{ijk}=gamma_ijql otherwise zero
	monoBlock <- which(sapply(1:Q, function(x) length(which(apply(tau, 1, which.max)==x))) == 1) # isolate single node block
		
	# Design matrix
	
	tcross_D <- array(0,dim=c(P,P,K))                       #outer products of column vectors in design matrix D
	for (k in 1:K) tcross_D[,,k] <- tcrossprod(D[k,])    #outer products of column vectors in design matrix D
	
	# Two likelihood parts, so like_beta_star = lik_beta + lik_beta_star
	
	# lik_beta_star <- 0 														# log of determinant of FI
	lik_beta <- 0																	# non-modified likelihood
	# block_lik_beta_star_tmp<-matrix(0,Q,Q) 					# blockwise likelihood contribution
	block_lik_beta_tmp<-matrix(0,Q,Q) 					# blockwise likelihood contribution
    
	# Compute terms which depend on tau (modified on 5/8/15)
	indDiag_ij <- as.numeric(outer(1:n + (0:(n-1)) * n, (0:(Q^2-1)) * n^2, FUN ="+" ))
	indDiag_ql <- as.numeric(outer(1:n^2, 0:(Q-1) * (Q+1) * n^2, FUN ="+" ))
	
	gamma <-  array(kronecker(tau, tau), dim=c(n, n, Q, Q)) # This is the term tau_{iq}tau_{jq} if q=l and tau_{iq}tau_{jl}+tau_{il}tau{jq} if q<l
	gamma <- gamma + aperm(gamma, c(1,2,4,3))
	gamma[indDiag_ij] <- 0
	gamma[indDiag_ql] <- gamma[indDiag_ql] / 2
	
	sum_gamma_without_diag <- apply(gamma, 3:4, sum)
	
	sum_connection <- array(crossprod(matrix(gamma, n^2, Q^2), matrix(X, n^2, K)), dim = c(Q, Q, K))

	sum_no_connection <-  array(sum_gamma_without_diag, dim = c(Q, Q, K)) - sum_connection
	
	
	# for (iM in monoBlock){
	# 	sum_connection[iM,iM,] = sum_no_connection[iM,iM,] = sum_gamma_without_diag[iM,iM] <-0
	# }
	
	# Initialisation (FIRST ROUND)

	delta_lik<-1 # difference between likelihoods
	m=1   
	# Compute FI and lik_beta star         
	for(q in 1:Q){
		for(l in q:Q){
			# if (q==l  & any(q== monoBlock)){
			# 	FI[q,l,,] = beta[q,l,] <- NaN               # Single nodes are not interpretable by the model
			# } else {	
				U_const[q,l,]=1/(1+exp(-D%*%beta[q,l,]))	
				FI_sqr[q,l,]= U_const[q,l,]*(1-U_const[q,l,])													        #
				# FI_sqr[q,l,]=exp(D%*%beta[q,l,])/(1+exp(D%*%beta[q,l,]))^2                          #
				pr0[q,l,]=1/(1+exp(D%*%beta[q,l,]))                                                                 #
				for(k in 1:K){
					FI[q,l,,] <- FI[q,l,,] + sum_gamma_without_diag[q,l] *FI_sqr[q,l,k] *tcross_D[,,k]  # sum on k
				}
				FI[q,l,,] <- 0.5*FI[q,l,,]
			   # Protection against the negative determinants as log(-x)=NaN and also log(0)=Inf
				# det_tmp=det(FI[q,l,,])
					# if(det_tmp>0){
						# lik_beta_star=lik_beta_star+log(det_tmp)
						# block_lik_beta_star_tmp[q,l] <- log(det_tmp)    						# log of determinant of FI
					# }else{
						# lik_beta_star=lik_beta_star+ log(exp(det_tmp)+10^-10)
						# block_lik_beta_star_tmp[q,l] <- log(exp(det_tmp)+10^-10) 		# log of determinant of FI
					# }                 
			# }
		}
	}
	# Compute full likelihood   like_beta_star = lik_beta + lik_beta_star
	# Numerical stabilisation:   -log(1+exp(-x)) is x for large negative x
	tmp <- log(U_const)   # tmp is array QQK
	if (min(tmp)==-Inf){
		tmp2 <- which(tmp==-Inf, arr.ind= TRUE)
		for (iInf in  1:dim(tmp2)[1]){
			tmp[tmp2[iInf,1], tmp2[iInf,2], tmp2[iInf,3]] <- D[tmp2[iInf,3],]%*%beta[tmp2[iInf,1], tmp2[iInf,2],]  #positive 
		}
	}
	#Numerical stabilisation:  -log(1+exp(x)) is -x for large x
	tmp3 <- log(pr0)
	if (min(tmp3)==-Inf){
		tmp2 <- which(tmp3==-Inf, arr.ind= TRUE)
		for (iInf in  1:dim(tmp2)[1]){
			tmp3[tmp2[iInf,1], tmp2[iInf,2], tmp2[iInf,3]] <- -D[tmp2[iInf,3],]%*%beta[tmp2[iInf,1], tmp2[iInf,2],]  #negative
		}
	}
	# Compute non-modified likelihood lik_beta:
	tmpex=sum_connection * tmp + sum_no_connection * tmp3
    lik_beta <- 0.5 * sum(tmpex)    
    # lik_beta_star=lik_beta+0.5*lik_beta_star   
    # block_lik_beta_star=0.5*(apply(tmpex,c(1,2),sum) +  block_lik_beta_star_tmp)
    block_lik_beta =0.5*(apply(tmpex,c(1,2),sum) +  block_lik_beta_tmp)
    
    
    # print(paste("lik_beta_star:", lik_beta_star)) 
	print(paste("lik_beta:", lik_beta))
   
    # Main cycle: 
   while(m <=M_iter & delta_lik>threshold_lik){
        # print(paste("iteration",m))
        # Previous values:
        
		# lik_beta_star_previous =lik_beta_star  					    		# lik_beta_star not known at the first iteration (m=1)
		lik_beta_previous=lik_beta			       					    				# lik_beta not known at the first iteration (m=1)
		# block_lik_beta_star_previous = block_lik_beta_star 		# blockwise contribution to lik_beta_star
		block_lik_beta_previous = block_lik_beta
		beta_prev<-beta                         # previous beta
		FI_prev<-FI                                # previous FI
	   
	   # Two likelihood parts, so like_beta_star = lik_beta + lik_beta_star

		# lik_beta_star <- 0                      
		lik_beta <- 0
		
		#Define objects:  Score and dFI derivative of FI
		U<-array(0,c(Q,Q,P))			    		# Score function
		# dFI<-array(0,c(Q,Q,P,P,P))	    		# Derivative of Fisher Information matrix
		delta_beta<-array(0,c(Q,Q,P)) 	#difference between new and previous beta
		for(q in 1:Q){
			for(l in q:Q){
				# if (q==l  & any(q== monoBlock)){
				# 	
				# } else {	
					# dFI_const[q,l,]=FI_sqr[q,l,]*(1-2*U_const[q,l,])
					# dFI_const[q,l,]=FI_sqr[q,l,]*(1-exp(D%*%beta[q,l,]))/(1+exp(D%*%beta[q,l,]))	    #
					for(k in 1:K){
						U[q,l,] <- U[q,l,] + (sum_connection[q,l,k] - sum_gamma_without_diag[q,l] * U_const[q,l,k])*t(D[k,])
						# for(p in 1:P){
							# dFI[q,l,p,,] <- dFI[q,l,p,,] + D[k,p] * dFI_const[q,l,k] *  sum_gamma_without_diag[q,l] * tcross_D[,,k]
						# }
					}
					# dFI[q,l,,,] <- 0.5 * dFI[q,l,,,] 
					U[q,l,] <- 0.5 * U[q,l,] 
					FI_inv[q,l,,]=solve(FI[q,l,,],tol=.Machine$double.eps^2)
					# for(p in 1:P){
						# U_star[q,l,p]=U[q,l,p] + 0.5*sum(diag(tcrossprod(FI_inv[q,l,,],dFI[q,l,p,,])))        # Modified Score
					# }
					# beta[q,l,]=beta[q,l,]	+FI_inv[q,l,,]%*%U_star[q,l,]                                                 # Update BETA
				    beta[q,l,]=beta[q,l,]	+FI_inv[q,l,,]%*%U[q,l,]                                                 # Update BETA
					delta_beta[q,l,]=beta[q,l,]- beta_prev[q,l,]
					if(any(abs(delta_beta[q,l,])> maxstep)){
						delta_beta[q,l,abs(delta_beta[q,l,])> maxstep] = 5 * sign(delta_beta[q,l,abs(delta_beta[q,l,])> maxstep])
						beta[q,l,] <- beta_prev[q,l,] + delta_beta[q,l,]
					}					
					U_const[q,l,]=1/(1+exp(-D%*%beta[q,l,]))														      # 
					pr0[q,l,]=1/(1+exp(D%*%beta[q,l,]))   
					 FI_sqr[q,l,]=  U_const[q,l,]*(1-U_const[q,l,])                                                          #
					# FI_sqr[q,l,]=exp(D%*%beta[q,l,])/(1+exp(D%*%beta[q,l,]))^2                             # 
					FI[q,l,,] <- 0
					for(k in 1:K){
						FI[q,l,,] <- FI[q,l,,] + sum_gamma_without_diag[q,l] *FI_sqr[q,l,k] * tcross_D[,,k] 			
					}
					FI[q,l,,]<-0.5*FI[q,l,,]
					# Protection against the negative determinants as log(-x)=NaN and also log(0)=Inf
					# det_tmp=det(FI[q,l,,])
					# if(det_tmp>0){
						# lik_beta_star=lik_beta_star+log(det_tmp)
						# block_lik_beta_star_tmp[q,l] <- log(det_tmp)
					# }else{
						# lik_beta_star=lik_beta_star+ log(exp(det_tmp)+10^-10)
						# block_lik_beta_star_tmp[q,l] <- log(exp(det_tmp)+10^-10)
					# }
					# print(paste("lik_beta_star:",lik_beta_star))
					# print(paste("q=",q,"l=",l))
				# }			
			}
		}
				
		#Numerical stabilisation:  -log(1+exp(x)) is -x,  for large x
		tmp <- log(U_const)
		if (min(tmp)==-Inf){
			tmp2 <- which(tmp==-Inf, arr.ind= TRUE)
			for (iInf in  1:dim(tmp2)[1]){
				tmp[tmp2[iInf,1], tmp2[iInf,2], tmp2[iInf,3]] <- D[tmp2[iInf,3],]%*%beta[tmp2[iInf,1], tmp2[iInf,2],]
			}
		}
		#Numerical stabilisation:  -log(1+exp(x)) is -x for large x
		tmp3 <- log(pr0)
		if (min(tmp3)==-Inf){
			tmp2 <- which(tmp3==-Inf, arr.ind= TRUE)
			for (iInf in  1:dim(tmp2)[1]){
				tmp3[tmp2[iInf,1], tmp2[iInf,2], tmp2[iInf,3]] <- -D[tmp2[iInf,3],]%*%beta[tmp2[iInf,1], tmp2[iInf,2],]
			}
		}
		# Compute likelihood
		tmpex=sum_connection * tmp + sum_no_connection * tmp3
    		# lik_beta <- 0.5 * sum(tmpex)    
    		# lik_beta_star=lik_beta+0.5*lik_beta_star   
   		# block_lik_beta_star =0.5*(apply(tmpex,c(1,2),sum) +  block_lik_beta_star_tmp)
   		block_lik_beta=0.5*(apply(tmpex,c(1,2),sum))
    		
    		# Check improvement:
   		# b_id=(block_lik_beta_star-block_lik_beta_star_previous)<0  # logical value block elements which gave a decrease in likelihood
   		b_id=(block_lik_beta - block_lik_beta_previous)<0  # logical value block elements which gave a decrease in likelihood
    	ihalf_iter <- 0
    	while (ihalf_iter < half_iter & any(b_id)){
    		ql_id <- which(b_id, arr.ind=TRUE)
	        delta_beta=beta- beta_prev
	        delta_beta[rep(b_id, P)] <- delta_beta[rep(b_id, P)] *0.5 # halving of delta_beta if needed
	        beta <- beta_prev + delta_beta # update beta
	        for(ql in 1:dim(ql_id)[1]){
	        	q <- ql_id[ql,1]		
	        	l <- ql_id[ql,2]		
				# if (q==l  & any(q== monoBlock)){
				# 	
				# } else {	
					U_const[q,l,]=1/(1+exp(-D%*%beta[q,l,]))														      # 
					pr0[q,l,]=1/(1+exp(D%*%beta[q,l,]))   
					FI_sqr[q,l,]=  U_const[q,l,]*(1-U_const[q,l,])                                                          #
					FI[q,l,,] <- 0
					for(k in 1:K){
						FI[q,l,,] <- FI[q,l,,] + sum_gamma_without_diag[q,l] *FI_sqr[q,l,k] * tcross_D[,,k] 			
					}
					FI[q,l,,]<-0.5*FI[q,l,,]
					# Protection against the negative determinants as log(-x)=NaN and also log(0)=Inf
					# det_tmp=det(FI[q,l,,])
					# if(det_tmp>0){
						# block_lik_beta_star_tmp[q,l] <- log(det_tmp)
					# }else{
						# block_lik_beta_star_tmp[q,l] <- log(exp(det_tmp)+10^-10)
					# }
					# print(paste("lik_beta_star:",lik_beta_star))
					# print(paste("q=",q,"l=",l))
					
					#Numerical stabilisation:  -log(1+exp(x)) is -x,  for large x
					tmp <- log(U_const[q,l,])
					if (min(tmp)==-Inf){
						tmp[tmp==-Inf] <- D[tmp==-Inf,]%*%beta[q, l,]
					}
					#Numerical stabilisation:  -log(1+exp(x)) is -x for large x
					tmp3 <- log(pr0[q,l,])
					if (min(tmp3)==-Inf){
						tmp[tmp2==-Inf] <- -D[tmp==-Inf,]%*%beta[q, l,]
					}
					# block_lik_beta_star[q,l] <- 0.5 * sum(sum_connection[q,l,] * tmp + sum_no_connection[q,l,]  * tmp3) + block_lik_beta_star_tmp[q,l]										
					block_lik_beta[q,l] <- 0.5 * sum(sum_connection[q,l,] * tmp + sum_no_connection[q,l,]  * tmp3) 									
				# }			
			}
			b_id <- (block_lik_beta - block_lik_beta_previous)<0 
			ihalf_iter <- ihalf_iter +1
		}
        
        # lik_beta_star <- sum(block_lik_beta_star)
        lik_beta <- sum(block_lik_beta)
       
		# delta_lik = lik_beta_star-lik_beta_star_previous
		delta_lik = lik_beta - lik_beta_previous
		# Condition: if no improvemnet return to previous values
		if (delta_lik < 0) {
			beta 	<- beta_prev
			FI 		<- FI_prev
			# lik_beta_star <- lik_beta_star_previous
			# block_lik_beta_star<-block_lik_beta_star_previous 
			lik_beta <- lik_beta_previous
			block_lik_beta <- block_lik_beta_previous 
		}
		# print(paste("lik_beta_star:", lik_beta_star))
		print(paste("lik_beta:", lik_beta))

		m=m+1
	}
	for (p in 1:P){
		beta[,,p]<-beta[,,p] + t(beta[,,p]) - diag(diag(beta[,,p]))
	}
	# return(list(beta=beta,FI=FI, lik_beta_star= lik_beta_star, lik_beta= lik_beta,block_lik_beta_star=block_lik_beta_star))
	return(list(beta=beta, FI=FI, lik_beta= lik_beta, block_lik_beta=block_lik_beta))
}