# Corresponds to LKMR-Main and LKMR-Supplement

# MCMC LOOPING
for (M in 1:num.reps)  {
    #Full Conditional Posteriors  
    
    # 1. Block diagonal matrix
    list.G = list()
	poly = polydot(degree=2, offset=1)
	
	for (g in 1:numtime) {
		list.G[[g]] = solve(kernelMatrix(poly, Z[,(grpsize*g-(grpsize-1)):(grpsize*g)]) + cofactor*diag(n)) / tau.sqf[g]
	}
    cov.bf1 = do.call(adiag, list.G)
        
    list.sig.sqf = c()
	list.sig.sqf[1] = 1/sig.sqf[1]
	list.sig.sqf[numtime] = 1/sig.sqf[(numtime-1)]
	for (g in 2:(numtime - 1)) {
		list.sig.sqf[g] = 1/sig.sqf[g] + 1/sig.sqf[g-1]
	}

	list.mat.sig.sqf = list()
	for (g in 1:numtime) {
		list.mat.sig.sqf[[g]] = (list.sig.sqf[g])*diag(n)
	}

	cov.bf2 = do.call(adiag, list.mat.sig.sqf)
    
    #3. Off-diagonals
	new.sig.sqf = rep(1/sig.sqf, each=n)
	cov.bf3 	= offDiagonal((-1)*new.sig.sqf)
    
	# 4. Sigma_h covariance matrix
	cov.bf = cov.bf1 + cov.bf2 + cov.bf3
   
   	# h - updated
    cov.bef      <- sig.sq0 * solve(XX+cov.bf+cofactor*diag(n*numtime))
    mean.bef     <- 1/sig.sq0 * cov.bef%*%t(X)%*%(Y - U%*%conf)
    beta.f[M+1,] <- rmvnorm(1,mean=mean.bef,sigma=cov.bef)

    # sig.sq0 - updated
    sh.sig      <- n*(numtime+1)/2 + sig.shape
    sc.sig      <- 1/2*t(Y-X%*%beta.f[M+1,] - U%*%conf)%*%(Y-X%*%beta.f[M+1,] - U%*%conf)+ 1/2*t(beta.f[M+1,])%*%cov.bf%*%beta.f[M+1,] + sig.scale
    sig.sq0     <- rinvgamma(1, shape=sh.sig, scale=sc.sig)
    sigsq0.post <- c(sigsq0.post, sig.sq0)
    
    # tausq - updated 
    gam <- c()
    for (j in 1:numtime){
    	term.beta = as.matrix(beta.f[M+1, ((j-1)*n+1):(j*n) ])
    	term.gam = t(term.beta) %*% solve(kernelMatrix(poly, Z[,(grpsize*j-(grpsize-1)):(grpsize*j)]) + cofactor*diag(n)) %*% term.beta 
    	gam[j]  = rinvGauss(1, nu=sqrt(lambda1.sq*sig.sq0/term.gam), lambda=lambda1.sq)            		
    	tau.sqf[j] = 1/gam[j] 
    }  	
    tausqf.post[M+1,] <- tau.sqf
    
    # omegasq - updated
    et <- c()
    for (k in 1:(numtime-1)){
    	nu.k     <- sqrt(lambda2.sq * sig.sq0/sum((beta.f[M+1,((k-1)*n+1):(k*n)]-beta.f[M+1,(k*n+1):((k+1)*n)])^2))
    	et[k]  <- rinvGauss(1, nu=nu.k, lambda=lambda2.sq)
    	sig.sqf[k] <- 1/et[k]
    	}
    sigsqf.post[M+1,] <- sig.sqf 	
	    
    # lambda1 - updated
    sh.lam1       <- numtime*(n+1)/2 + a.1 
    sc.lam1       <- 1/2*sum(tau.sqf) + b.1 
    lambda1.sq    <- rgamma(1, shape=sh.lam1, rate=sc.lam1)
    lambda1f.post <- c(lambda1f.post, lambda1.sq)
	
	# lambda2 - udpated
    sh.lam2       <- numtime - 1 + a.2 
    sc.lam2       <- 1/2*sum(sig.sqf) + b.2 
    lambda2.sq    <- rgamma(1, shape=sh.lam2, rate=sc.lam2)
    lambda2f.post <- c(lambda2f.post, lambda2.sq) 
    
    # Beta - updated
    mean.conf = solve(t(U) %*% U) %*% t(U) %*% (Y - X%*%beta.f[M+1,])
    sig.conf = sig.sq0 * solve(t(U) %*% U) 
    conf = rmvnorm(1,mean=mean.conf,sigma=sig.conf)
    conf.post[M+1,] <- conf  
    conf = matrix(conf, nrow=c)   
}