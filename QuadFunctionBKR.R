## function to fit Bayesian kernel machine regression
# quadratic kernel

library(kernlab)
beta.update <- function(X, Vinv, y, sigsq.eps) {
	XVinv <- crossprod(X, Vinv)
	Vbeta <- chol2inv(chol(XVinv %*% X))
	cholVbeta <- chol(Vbeta)
	betahat <- Vbeta %*% XVinv %*% y
	n01 <- rnorm(ncol(X))
	betahat + crossprod(sqrt(sigsq.eps)*cholVbeta, n01)
}

sigsq.eps.update <- function(a.eps, b.eps, y, X, beta, Vinv) {
	mu <- y - X%*%beta
	prec.y <- rgamma(1, shape=a.eps + nrow(X)/2, rate=b.eps + 1/2*crossprod(mu, Vinv)%*%mu)
	1/prec.y
}

rdelta.update <- function(r, delta, lambda, y, X, beta, sigsq.eps, cholV, In, Kpart, Vinv, n, nz, Z, ztest, control.params) {
	mu.r <- control.params$mu.r
	sigma.r <- control.params$sigma.r
	r.muprop <- control.params$r.muprop
	a.p0 <- control.params$a.p0
	b.p0 <- control.params$b.p0
	delta.star <- delta
	r.star <- r
	
	move.type <- ifelse(all(delta[ztest] == 0), 1, sample(c(1,2),1))
	if(move.type == 1) {
		comp <- ifelse(length(ztest) == 1, ztest, sample(ztest, 1))
		r.jump <- control.params$r.jump1
		delta.star[comp] <- 1 - delta[comp]
		r.star[comp] <- ifelse(delta.star[comp] == 0, 0, rgamma(1, shape=r.muprop^2/r.jump^2, rate=r.muprop/r.jump^2))
		diffpriors <- (lgamma(sum(delta.star[ztest]) + a.p0) + lgamma(length(ztest) - sum(delta.star[ztest]) + b.p0) - lgamma(sum(delta[ztest]) + a.p0) - lgamma(length(ztest) - sum(delta[ztest]) + b.p0)) + ifelse(delta[comp] == 1, -1, 1)*dgamma(ifelse(delta[comp] == 1, r[comp], r.star[comp]), shape=mu.r^2/sigma.r^2, rate=mu.r/sigma.r^2, log=TRUE)
		negdifflogproposal <- -ifelse(delta[comp] == 1, -1, 1)*dgamma(ifelse(delta[comp] == 1, r[comp], r.star[comp]), shape=r.muprop^2/r.jump^2, rate=r.muprop/r.jump^2, log=TRUE)
	} else if(move.type == 2) {
		comp <- ifelse(length(which(delta == 1)) == 1, which(delta == 1), sample(which(delta == 1), 1))
		r.jump <- control.params$r.jump2
		r.star[comp] <- rgamma(1, shape=(r[comp])^2/r.jump^2, rate=(r[comp])/r.jump^2)
		diffpriors <- dgamma(r.star[comp], shape=mu.r^2/sigma.r^2, rate=mu.r/sigma.r^2, log=TRUE) - dgamma(r[comp], shape=mu.r^2/sigma.r^2, rate=mu.r/sigma.r^2, log=TRUE)
		negdifflogproposal <- -dgamma(r.star[comp], shape=(r[comp])^2/r.jump^2, rate=(r[comp])/r.jump^2, log=TRUE) + dgamma(r[comp], shape=(r.star[comp])^2/r.jump^2, rate=(r.star[comp])/r.jump^2, log=TRUE)
	}

	lambda.star <- lambda
	poly = polydot(degree=2, offset=1)
	rZ = t(sqrt(r.star)*t(Z))
	Kpartstar = kernelMatrix(poly, rZ)
	
	## M-H step
	return(MHstep(In=In, r=r, lambda=lambda, lambda.star=lambda.star, r.star=r.star, delta=delta, delta.star=delta.star, Kpart=Kpart, Kpartstar=Kpartstar, y=y, X=X, beta=beta, sigsq.eps=sigsq.eps, diffpriors=diffpriors, negdifflogproposal=negdifflogproposal, cholV=cholV, Vinv=Vinv, move.type=move.type))
}

r.update <- function(r, whichcomp, delta, lambda, y, X, beta, sigsq.eps, cholV, In, Kpart, Vinv, n, nz, Z, control.params) {
	r.jump <- control.params$r.jump
	mu.r <- control.params$mu.r
	sigma.r <- control.params$sigma.r
	rcomp <- unique(r[whichcomp])
	if(length(rcomp) > 1) stop("rcomp should only be 1-dimensional")
	
	## generate a proposal
	rcomp.star <- rgamma(1, shape=rcomp^2/r.jump^2, rate=rcomp/r.jump^2)
	lambda.star <- lambda
	delta.star <- delta
	move.type <- NA
	
	## part of M-H ratio that depends on the proposal distribution
	negdifflogproposal <- -dgamma(rcomp.star, shape=rcomp^2/r.jump^2, rate=rcomp/r.jump^2, log=TRUE) + dgamma(rcomp, shape=rcomp.star^2/r.jump^2, rate=rcomp.star/r.jump^2, log=TRUE)
	
	## prior distribution
	diffpriors <- dgamma(rcomp.star, shape=mu.r^2/sigma.r^2, rate=mu.r/sigma.r^2, log=TRUE) - dgamma(rcomp, shape=mu.r^2/sigma.r^2, rate=mu.r/sigma.r^2, log=TRUE)
	
	r.star <- r
	r.star[whichcomp] <- rcomp.star
	poly = polydot(degree=2, offset=1)
	rZ = t(sqrt(r.star)*t(Z))
	Kpartstar = kernelMatrix(poly, rZ)
	
	## M-H step
	return(MHstep(In=In, r=r, lambda=lambda, lambda.star=lambda.star, r.star=r.star, delta=delta, delta.star=delta.star, Kpart=Kpart, Kpartstar=Kpartstar, y=y, X=X, beta=beta, sigsq.eps=sigsq.eps, diffpriors=diffpriors, negdifflogproposal=negdifflogproposal, cholV=cholV, Vinv=Vinv, move.type=move.type))
}

lambda.update <- function(r, delta, lambda, y, X, beta, sigsq.eps, cholV, In, Kpart, Vinv, control.params) {
	lambda.jump <- control.params$lambda.jump
	mu.lambda <- control.params$mu.lambda
	sigma.lambda <- control.params$sigma.lambda

	## generate a proposal
	lambda.star <- rgamma(1, shape=lambda^2/lambda.jump^2, rate=lambda/lambda.jump^2)
	r.star <- r
	delta.star <- delta
	Kpartstar <- Kpart
	move.type <- NA
	
	## part of M-H ratio that depends on the proposal distribution
	negdifflogproposal <- -dgamma(lambda.star, shape=lambda^2/lambda.jump^2, rate=lambda/lambda.jump^2, log=TRUE) + dgamma(lambda, shape=lambda.star^2/lambda.jump^2, rate=lambda.star/lambda.jump^2, log=TRUE)
	
	## prior distribution
	diffpriors <- dgamma(lambda.star, shape=mu.lambda^2/sigma.lambda^2, rate=mu.lambda/sigma.lambda^2, log=TRUE) - dgamma(lambda, shape=mu.lambda^2/sigma.lambda^2, rate=mu.lambda/sigma.lambda^2, log=TRUE)
	
	## M-H step
	return(MHstep(In=In, r=r, lambda=lambda, lambda.star=lambda.star, r.star=r.star, delta=delta, delta.star=delta.star, Kpart=Kpart, Kpartstar=Kpartstar, y=y, X=X, beta=beta, sigsq.eps=sigsq.eps, diffpriors=diffpriors, negdifflogproposal=negdifflogproposal, cholV=cholV, Vinv=Vinv, move.type=move.type))
}

MHstep <- function(In, r, lambda, lambda.star, r.star, delta, delta.star, Kpart, Kpartstar, y, X, beta, sigsq.eps, diffpriors, negdifflogproposal, cholV, Vinv, move.type) {
	## compute log M-H ratio
	Vstar <- In + lambda.star*Kpartstar # Altered
	cholVstar <- chol(Vstar)
	Vinvstar <- chol2inv(cholVstar)
	mu <- y - X%*%beta
	diffliks <- -sum(log(diag(cholVstar))) + sum(log(diag(cholV))) - 1/2/sigsq.eps*crossprod(mu, Vinvstar - Vinv)%*%mu
	logMHratio <- diffliks + diffpriors + negdifflogproposal
	logalpha <- min(0,logMHratio)

	## return value
	acc <- FALSE
	if( log(runif(1)) <= logalpha ) {
		r <- r.star
		delta <- delta.star
		Kpart <- Kpartstar
		lambda <- lambda.star
		Vinv <- Vinvstar
		cholV <- cholVstar
		acc <- TRUE
	}
	return(list(r=r, lambda=lambda, delta=delta, Kpart=Kpart, acc=acc, cholV=cholV, Vinv=Vinv, move.type=move.type))
}	

h.update <- function(lambda, Kpart, Vinv, sigsq.eps, y, X, beta) {
	K = Kpart
	lamKVinv <- lambda*K%*%Vinv
	h.postmean <- lamKVinv%*%(y-X%*%beta)
	h.postvar <- sigsq.eps*lamKVinv
	h.postvar.sqrt <- try(chol(h.postvar), silent=TRUE)
	if(class(h.postvar.sqrt) == "try-error") {
		sigsvd <- svd(h.postvar)
		h.postvar.sqrt <- t(sigsvd$v %*% (t(sigsvd$u) * sqrt(sigsvd$d)))
	}
	h.postmean + crossprod(h.postvar.sqrt, rnorm(length(h.postmean)))
}

newh.update <- function(Z, Znew, Vinv, lambda, sigsq.eps, r, y, X, beta, n0, n1, nz) {
	Z.star = rbind(Z, Znew)
	poly = polydot(degree=2, offset=1)
	rZ = t(sqrt(r)*t(Z.star))
	Kpartall = kernelMatrix(poly, rZ)

	nall <- n0 + n1
	Kmat <- Kpartall 
	Kmat0 <- Kmat[1:n0,1:n0 ,drop=FALSE]
	Kmat1 <- Kmat[(n0+1):nall,(n0+1):nall ,drop=FALSE]
	Kmat10 <- Kmat[(n0+1):nall,1:n0 ,drop=FALSE]

	if(is.null(Vinv)) Vinv <- chol2inv(chol(diag(1,n0,n0) + lambda*Kmat0))
	lamK10Vinv <- lambda*Kmat10 %*% Vinv
	Sigma.hnew <- lambda*sigsq.eps*(Kmat1 - lamK10Vinv %*% t(Kmat10))
	mu.hnew <- lamK10Vinv %*% (y - X%*%beta)
	root.Sigma.hnew <- try(chol(Sigma.hnew), silent=TRUE)
	if(class(root.Sigma.hnew) == "try-error") {
		sigsvd <- svd(Sigma.hnew)
		root.Sigma.hnew <- t(sigsvd$v %*% (t(sigsvd$u) * sqrt(sigsvd$d)))
	}
	mu.hnew + crossprod(root.Sigma.hnew, rnorm(n1))
}

kmbayes <- function(nsamp=1000, y, X, Z, id, quiet=TRUE, Znew, starting.values=list(), control.params=list(), modsel=FALSE, ztest) {
	n <- nrow(X)
	nbeta <- ncol(X)
	In <- diag(1,n,n)
	nz <- ncol(Z)

	## create empty matrices to store the posterior draws in
	chain <- list(h.hat=matrix(0,nsamp,n), beta=matrix(0,nsamp,nbeta), lambda=rep(NA,nsamp), sigsq.eps=rep(NA,nsamp), r=matrix(NA,nsamp,nz), acc.r=rep(0,nsamp), acc.lambda=rep(0,nsamp), delta=matrix(1,nsamp,nz))
	if(modsel) {
		chain$acc.rdelta <- rep(0,nsamp)
		chain$move.type <- rep(0,nsamp)
	}

	## components to predict h(Znew)
	if(!missing(Znew)) {
		if(is.null(dim(Znew))) Znew <- matrix(Znew, nrow=1)
		if(class(Znew) == "data.frame") Znew <- data.matrix(Znew)
		if(ncol(Z) != ncol(Znew)) {
			stop("Znew must have the same number of columns as Z")
		}
		chain$hnew <- matrix(0,nsamp,nrow(Znew))
		colnames(chain$hnew) <- rownames(Znew)
	}

	## components if model selection is being done
	if(modsel) {
		if(missing(ztest)) {
			ztest <- 1:nz
		}
	} else {
		ztest <- NULL
	}
	
	## control parameters
	control.params <- modifyList(list(r.jump=1, mu.r=5, sigma.r=5, lambda.jump=1, mu.lambda=1, sigma.lambda=10, a.p0=1, b.p0=1, r.muprop=1, a.eps=1e-3, b.eps=1e-3), control.params)
	
	## initial values
	starting.values <- modifyList(list(h.hat=1, beta=10, sigsq.eps=1, r=1, lambda=1, delta=1), starting.values)
	chain$h.hat[1,] <- starting.values$h.hat
	chain$beta[1,] <- starting.values$beta
	chain$lambda[1] <- starting.values$lambda
	chain$sigsq.eps[1] <- starting.values$sigsq.eps
	chain$r[1,] <- starting.values$r
	if(modsel) {
		chain$delta[1,ztest] <- starting.values$delta
	}
	
	poly = polydot(degree=2, offset=1)
	rZ = t(sqrt(chain$r[1,])*t(Z))
	Kpart = kernelMatrix(poly, rZ)
	
	cholV <- chol(In + chain$lambda[1]*Kpart)
	Vinv <- chol2inv(cholV)
	
	chain$time1 <- Sys.time()
	for(s in 2:nsamp) {
		
		###################################################
		## generate posterior samples from marginalized distribution P(beta, sigsq.eps, lambda, r | y)
		
		## beta
		chain$beta[s,] <- beta.update(X=X, Vinv=Vinv, y=y, sigsq.eps=chain$sigsq.eps[s-1])
		
		## \sigma_\epsilon^2
		chain$sigsq.eps[s] <- sigsq.eps.update(a.eps=control.params$a.eps, b.eps=control.params$b.eps, y=y, X=X, beta=chain$beta[s,], Vinv=Vinv)
		
		## lambda
		varcomps <- lambda.update(r=chain$r[s-1,], delta=chain$delta[s-1,], lambda=chain$lambda[s-1], y=y, X=X, beta=chain$beta[s,], sigsq.eps=chain$sigsq.eps[s], cholV=cholV, In=In, Kpart=Kpart, Vinv=Vinv, control.params=control.params)
		chain$lambda[s] <- varcomps$lambda
		if(varcomps$acc) {
			Kpart <- varcomps$Kpart
			Vinv <- varcomps$Vinv
			cholV <- varcomps$cholV
			chain$acc.lambda[s] <- varcomps$acc
		}

		## r
		rSim <- chain$r[s-1,]
		comp <- which(!1:nz %in% ztest)
		if(length(comp) != 0) {
			varcomps <- r.update(r=rSim, whichcomp=comp, delta=chain$delta[s-1,], lambda=chain$lambda[s], y=y, X=X, beta=chain$beta[s,], sigsq.eps=chain$sigsq.eps[s], cholV=cholV, In=In, Kpart=Kpart, Vinv=Vinv, n=n, nz=nz, Z=Z, control.params=control.params)
			rSim <- varcomps$r
			if(varcomps$acc) {
				Kpart <- varcomps$Kpart
				Vinv <- varcomps$Vinv
				cholV <- varcomps$cholV
				chain$acc.r[s] <- varcomps$acc
			}
		}
		## for those variables being selected: joint posterior of (r,delta)
		if(modsel) {
			varcomps <- rdelta.update(r=rSim, delta=chain$delta[s-1,], lambda=chain$lambda[s], y=y, X=X, beta=chain$beta[s,], sigsq.eps=chain$sigsq.eps[s], cholV=cholV, In=In, Kpart=Kpart, Vinv=Vinv, n=n, nz=nz, Z=Z, ztest=ztest, control.params=control.params)
			chain$delta[s,] <- varcomps$delta
			rSim <- varcomps$r
			chain$move.type[s] <- varcomps$move.type
			if(varcomps$acc) {
				Kpart <- varcomps$Kpart
				Vinv <- varcomps$Vinv
				cholV <- varcomps$cholV
				chain$acc.rdelta[s] <- varcomps$acc
			}
		}
		chain$r[s,] <- rSim

		###################################################
		## generate posterior sample of h(z) from its posterior P(h | beta, sigsq.eps, lambda, r, y)
		
		chain$h.hat[s,] <- h.update(lambda=chain$lambda[s], Kpart=Kpart, Vinv=Vinv, sigsq.eps=chain$sigsq.eps[s], y=y, X=X, beta=chain$beta[s,])
		
		###################################################
		## generate posterior samples of h(Znew) from its posterior P(hnew | beta, sigsq.eps, lambda, r, y)
		
		if(!missing(Znew)) {
			chain$hnew[s,] <- newh.update(Z=Z, Znew=Znew, Vinv=Vinv, lambda=chain$lambda[s], sigsq.eps=chain$sigsq.eps[s], r=chain$r[s,], y=y, X=X, beta=chain$beta[s,], n0=n, n1=nrow(Znew), nz=nz)
		}
		
		###################################################
		## print details of the model fit so far
		if(s%%(nsamp/10)==0 & !quiet) {
			print(s)
			cat(mean(chain$acc.lambda[2:s]), "   lam accept rate\n")
			cat(mean(chain$acc.r[2:s]), "   r nosel accept rate\n")
			if(modsel) {
				cat(mean(chain$acc.rdelta[2:s]), "   rdelt accept rate\n")
				cat(mean(chain$acc.rdelta[2:s][chain$move.type[2:s] == 1]), "   rdelt[move 1] accept rate\n")
				cat(mean(chain$acc.rdelta[2:s][chain$move.type[2:s] == 2]), "   rdelt[move 2] accept rate\n")
				cat(round(colMeans(chain$delta[1:s,ztest ,drop=FALSE]),4), "   post incl probs\n")
				cat(round(colMeans(chain$r[2:s,], na.rm=TRUE),4), "   post mean of r\n")
			}
			print(difftime(Sys.time(), chain$time1))
		}	
	}
	chain$time2 <- Sys.time()
	chain$nsamp <- nsamp
	chain$starting.values <- starting.values
	chain$control.params <- control.params	
	chain$X <- X
	chain$Z <- Z
	chain$y <- y
	chain$ztest <- ztest
	if(!missing(Znew)) chain$Znew <- Znew
	chain
}

## function to obtain posterior samples of h(znew) from fit of Bayesian kernel machine regression
predz.samps <- function(fit, Znew, quiet=FALSE) {
	if(is.null(dim(Znew))) Znew <- matrix(Znew, nrow=1)
	if(class(Znew) == "data.frame") Znew <- data.matrix(Znew)
	Z <- fit$Z
	nsamp <- fit$nsamp
	if(ncol(Z) != ncol(Znew)) {
		stop("Znew must have the same number of columns as Z")
	}
	
	hnew.samps <- sapply(1:fit$nsamp, function(s) {
		if(s%%(nsamp/10)==0 & !quiet) print(s)
		newh.update(Z=Z, Znew=Znew, Vinv=NULL, lambda=fit$lambda[s], sigsq.eps=fit$sigsq.eps[s], r=fit$r[s,], y=fit$y, X=fit$X, beta=fit$beta[s,], n0=nrow(Z), n1=nrow(Znew), nz=ncol(Z))
	})
	rownames(hnew.samps) <- rownames(Znew)
	t(hnew.samps)
}

## function to approximate the posterior mean and variance as a function of the estimated tau, lambda, beta, and sigsq.eps
newh.postmean <- function(fit, Znew, sel) {
	if(is.null(dim(Znew))) Znew <- matrix(Znew, nrow=1)
	if(class(Znew) == "data.frame") Znew <- data.matrix(Znew)

	Z <- fit$Z
	X <- fit$X
	y <- fit$y
	n0 <- nrow(Z)
	n1 <- nrow(Znew)
	nz <- ncol(Z)
	nall <- n0+n1
	lambda <- mean(fit$lambda[sel])
	sigsq.eps <- mean(fit$sigsq.eps[sel])
	r <- colMeans(fit$r[sel,])
	beta <- colMeans(fit$beta[sel, ,drop=FALSE])

	#Kpartall <- as.matrix(dist(sqrt(matrix(r, byrow=TRUE, n0+n1, nz))*rbind(Z,Znew)))^2
	#Kmat <- exp(-Kpartall)
	
	Z.star = rbind(Z, Znew)
	poly = polydot(degree=2, offset=1)
	rZ = t(sqrt(r)*t(Z.star))
	Kpartall = kernelMatrix(poly, rZ)
	Kmat = Kpartall
	
	Kmat0 <- Kmat[1:n0,1:n0 ,drop=FALSE]
	Kmat1 <- Kmat[(n0+1):nall,(n0+1):nall ,drop=FALSE]
	Kmat10 <- Kmat[(n0+1):nall,1:n0 ,drop=FALSE]

	Vinv <- chol2inv(chol(diag(1,n0,n0) + lambda*Kmat0))
	lamK10Vinv <- lambda*Kmat10 %*% Vinv
	Sigma.hnew <- lambda*sigsq.eps*(Kmat1 - lamK10Vinv %*% t(Kmat10))
	mu.hnew <- lamK10Vinv %*% (y - X%*%beta)
	list(postmean=drop(mu.hnew), postvar=Sigma.hnew)
}