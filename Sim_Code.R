library(lars); library(lasso2); library(mvtnorm); library(SuppDists); library(MCMCpack); library(grplasso); library(magic); library(kernlab); library(MASS); library(fields)

#jobid <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
#set.seed(jobid)

set.seed(1)

####################################################################
n      		= 100
grpsize 	= 5
numtime 	= 4
p			    = n*numtime
REP 		  = 1
cofactor	= 1e-7
num.reps	= 10000
sel 		  = seq(2000, num.reps, by=5)
a.1 		  = 60 
b.1 		  = 10
a.2 		  = 45
b.2 		  = 10
sig.shape 	= 5
sig.scale 	= 5
c 			  = 2 
ar 			  = 0.8 # modify this for varying levels of AR-1. 

autocorr.mat <- function(p, rho) {
  mat <- diag(p)
  return(rho^abs(row(mat)-col(mat)))
}
time_mat = autocorr.mat(numtime,ar)

metal_mat <- matrix( c(1, 0.19, 0.32, 0.15, 0.4, 0.19, 1, 0.11, 0.25, 0.4, 0.32, 0.11, 1, 0.17, 0.24, 0.15, 0.25, 0.17, 1, 0.35, 0.4, 0.4, 0.24, 0.35, 1), 5,5, byrow=TRUE) 

cmat = kronecker(time_mat, metal_mat)
obs <- matrix( mvrnorm(n, mu=rep(0, 5*4), Sigma=cmat), nrow = n)

Z = obs

Z.1 = Z[,1:5]
Z.2 = Z[,6:10]
Z.3 = Z[,11:15]
Z.4 = Z[,16:20]

conf = rep(1, c)
U = cbind(scale(matrix(rnorm(n, mean = 10, sd = 1))), scale(matrix(sample(1:2, n, replace = TRUE))))

X.Cov = U
	    	         
res.sd = 1 
Y = matrix(rnorm(n=n, mean=0, sd=res.sd)) + U%*%conf + 0.5*(Z.2[,1]^2 - Z.2[,2]^2 + 1/2*Z.2[,1]*Z.2[,2] + Z.2[,1] + Z.2[,2]) + 0.8*(Z.3[,1]^2 - Z.3[,2]^2 + 1/2*Z.3[,1]*Z.3[,2] + Z.3[,1] + Z.3[,2]) + 1*(Z.4[,1]^2 - Z.4[,2]^2 + 1/2*Z.4[,1]*Z.4[,2] + Z.4[,1] + Z.4[,2])

oX = diag(n)
counter = 1
while(counter < numtime) {
	oX = cbind(oX, diag(n))
	counter = counter+1
}
X      = oX 
XX     = t(X)%*%X

offDiagonal = function(x) {
		diagBelow = diag(x)
		i = 1
		while (i <= n) {
			diagBelow=rbind(rep(0,length(x)	+i),cbind(diagBelow,rep(0,length(x) + i - 1)))
			i = i + 1
		}
		mat <- diagBelow + t(diagBelow) - diag(diag(diagBelow))
		return(mat)
}

##############################################################################
# PRIORS 
lambda1.sq <- rgamma(1,shape=a.1, rate=b.1) 
lambda2.sq <- rgamma(1,shape=a.2, rate=b.2) 
sig.sq0  <- rinvgamma(1, shape = sig.shape, scale = sig.scale)
mean.bef <- rep(0,p)
tau.sqf  <- rgamma(numtime, shape = (n + 1)/2, rate=lambda1.sq/2)
sig.sqf  <- rexp(numtime-1, rate=lambda2.sq/2)


#Full Conditional Posteriors  
# beta
    
# 1. Block diagonal matrix 
list.G = list()
poly = polydot(degree=2, offset=1)
for (g in 1:numtime) {
	list.G[[g]] = solve(kernelMatrix(poly, Z[,(grpsize*g-(grpsize-1)):(grpsize*g)]) + cofactor*diag(n)) / tau.sqf[g]
	}  
cov.bf1 = do.call(adiag, list.G)

# 2. Diagonal

list.sig.sqf = c()
list.sig.sqf[1] = 1/sig.sqf[1]
list.sig.sqf[numtime] = 1/sig.sqf[(numtime-1)]
for (g in 2:(numtime - 1)) {
	list.sig.sqf[g] = 1/sig.sqf[g] + 1/sig.sqf[g-1]
}

list.mat.sig.sqf = list()
for (g in 1:numtime) {
	list.mat.sig.sqf[[g]] = list.sig.sqf[g] * diag(n)
}

cov.bf2 = do.call(adiag, list.mat.sig.sqf)

#3. Off-diagonals
new.sig.sqf = rep(1/sig.sqf, each=n)
cov.bf3 	= offDiagonal((-1)*new.sig.sqf)
    
# 4. Beta covariance matrix
cov.bf = cov.bf1 + cov.bf2 + cov.bf3

# A few ways to invert cov.bf
SIG       	= sig.sq0*solve(cov.bf + cofactor*diag(n*numtime))

beta.fp    <- rmvnorm(1,mean=mean.bef,sigma=SIG)     
    
# FOR POSTERIOR
sigsq0.post <- lambda1f.post <- lambda2f.post <- NULL 
beta.f      <- rbind( beta.fp,matrix(rep(NA,num.reps*p),ncol=p) )
tausqf.post <- rbind( tau.sqf,matrix(rep(NA,num.reps*numtime),ncol=numtime) )
sigsqf.post <- rbind( sig.sqf,matrix(rep(NA,num.reps*(numtime-1)),ncol=numtime-1) )
conf.post = rbind(conf, matrix(rep(NA,num.reps*c),nrow=num.reps))
conf = matrix(conf, nrow = c)

######################################################
######################################################    
# POSTERIOR - No omegas or taus - just linear kernel
cat(c("Job started at:",date()),fill=TRUE)
source("MCMC.R")
cat(c("Job finished at:",date()),fill=TRUE)

resultsGFLK = list(Sigsq = sigsq0.post, Lam1 = lambda1f.post, Lam2 = lambda2f.post, Beta = beta.f, Tau = tausqf.post, Omega = sigsqf.post, Z = Z, Y=Y, U = U, sel = sel, Conf = conf.post, cofactor = cofactor)

mat.res.gfl = matrix(NA, nrow=numtime, ncol=5)

Z.1 = Z[,1:5]
Z.2 = Z[,6:10]
Z.3 = Z[,11:15]
Z.4 = Z[,16:20]

ind = 1:100
hgrid.bapprox.T.1 <- apply(resultsGFLK$Beta[sel,], 2, mean)[ind]
true_hgrid.linear.1 = rep(0, n)

# Linear h
mod_grid.linear.1 = lm(hgrid.bapprox.T.1 ~ true_hgrid.linear.1)

# coverage
gf.cb.1 = apply(resultsGFLK$Beta[sel,ind], 2, quantile, c(0.025, 0.975), na.rm = TRUE)

mat.res.gfl[1,] = c(summary(mod_grid.linear.1)$coef[1,1], 0, 0, sqrt( sum( (hgrid.bapprox.T.1 - true_hgrid.linear.1)^2) / n), sum(true_hgrid.linear.1 >= gf.cb.1[1,] & true_hgrid.linear.1 <= gf.cb.1[2,] )/length(true_hgrid.linear.1))

# Time 2
ind = 101:200

hgrid.bapprox.T.2 <- apply(resultsGFLK$Beta[sel,], 2, mean)[ind]

true_hgrid.linear.2 = 0.5*(Z.2[,1]^2 - Z.2[,2]^2 + 1/2*Z.2[,1]*Z.2[,2] + Z.2[,1] + Z.2[,2])

# Linear h
mod_grid.linear.2 = lm(hgrid.bapprox.T.2 ~ true_hgrid.linear.2)

# coverage
gf.cb.2 = apply(resultsGFLK$Beta[sel,ind], 2, quantile, c(0.025,0.975), na.rm = TRUE)

mat.res.gfl[2,] = c(summary(mod_grid.linear.2)$coef[1,1], summary(mod_grid.linear.2)$coef[2,1], summary(mod_grid.linear.2)$r.sq, sqrt( sum( (hgrid.bapprox.T.2 - true_hgrid.linear.2)^2) / n), sum(true_hgrid.linear.2 >= gf.cb.2[1,] & true_hgrid.linear.2 <= gf.cb.2[2,] )/length(true_hgrid.linear.2))

###########################
###########################
###########################
# Estimation at time 3.

ind = 201:300
hgrid.bapprox.T.3 <- apply(resultsGFLK$Beta[sel,], 2, mean)[ind]

true_hgrid.linear.3 = 0.8*(Z.3[,1]^2 - Z.3[,2]^2 + 1/2*Z.3[,1]*Z.3[,2] + Z.3[,1] + Z.3[,2])

mod_grid.linear.3 = lm(hgrid.bapprox.T.3 ~ true_hgrid.linear.3)
gf.cb.3 = apply(resultsGFLK$Beta[sel,ind], 2, quantile, c(0.025,0.975), na.rm = TRUE)
mat.res.gfl[3,] = c(summary(mod_grid.linear.3)$coef[1,1], summary(mod_grid.linear.3)$coef[2,1], summary(mod_grid.linear.3)$r.sq, sqrt( sum( (hgrid.bapprox.T.3 - true_hgrid.linear.3)^2) / n), sum(true_hgrid.linear.3 >= gf.cb.3[1,] & true_hgrid.linear.3 <= gf.cb.3[2,] )/length(true_hgrid.linear.3))

###########################
###########################
###########################

# Estimation at time 4.

ind = 301:400
hgrid.bapprox.T.4 <- apply(resultsGFLK$Beta[sel,], 2, mean)[ind]

true_hgrid.linear.4 = 1*(Z.4[,1]^2 - Z.4[,2]^2 + 1/2*Z.4[,1]*Z.4[,2] + Z.4[,1] + Z.4[,2])

mod_grid.linear.4 = lm(hgrid.bapprox.T.4 ~ true_hgrid.linear.4)

gf.cb.4 = apply(resultsGFLK$Beta[sel,ind], 2, quantile, c(0.025,0.975), na.rm = TRUE)

mat.res.gfl[4,] = c(summary(mod_grid.linear.4)$coef[1,1], summary(mod_grid.linear.4)$coef[2,1], summary(mod_grid.linear.4)$r.sq, sqrt( sum( (hgrid.bapprox.T.4 - true_hgrid.linear.4)^2) / n), sum(true_hgrid.linear.4 >= gf.cb.4[1,] & true_hgrid.linear.4 <= gf.cb.4[2,] )/length(true_hgrid.linear.4))

#save(mat.res.gfl, file=paste0("Output/matresGFL", jobid, ".RData"))

#################################

mat.res.bkr = matrix(NA, nrow=numtime, ncol=5)

#################################
# Time 1
Z.Time = Z.1

source('QuadFunctionBKR.R')
source('postmeanbkmrcode.R')

preds_g.linear.1 = apply(fit_g.linear$h.hat[sel,], 2, mean)
mod_g.linear.1 = lm(preds_g.linear.1 ~ true_hgrid.linear.1)

# coverage
bkr.cb.1 = apply(fit_g.linear$h.hat[sel,], 2, quantile, c(0.025,0.975), na.rm = TRUE)
mat.res.bkr[1,] = c(summary(mod_g.linear.1)$coef[1,1], 0, 0, sqrt( sum( (preds_g.linear.1 - true_hgrid.linear.1)^2) / n), sum(true_hgrid.linear.1 >= bkr.cb.1[1,] & true_hgrid.linear.1 <= bkr.cb.1[2,] )/length(true_hgrid.linear.1))

######################################################
# Time 2
Z.Time = Z.2

source('postmeanbkmrcode.R')

preds_g.linear.2 = apply(fit_g.linear$h.hat[sel,], 2, mean)
mod_g.linear.2 = lm(preds_g.linear.2 ~ true_hgrid.linear.2)
# coverage
bkr.cb.2 = apply(fit_g.linear$h.hat[sel,], 2, quantile, c(0.025,0.975), na.rm = TRUE)
mat.res.bkr[2,] = c(summary(mod_g.linear.2)$coef[1,1], summary(mod_g.linear.2)$coef[2,1], summary(mod_g.linear.2)$r.sq, sqrt( sum( (preds_g.linear.2 - true_hgrid.linear.2)^2) / n), sum(true_hgrid.linear.2 >= bkr.cb.2[1,] & true_hgrid.linear.2 <= bkr.cb.2[2,] )/length(true_hgrid.linear.2))

######################################################
# Time 3
Z.Time = Z.3

source('postmeanbkmrcode.R')

preds_g.linear.3 = apply(fit_g.linear$h.hat[sel,], 2, mean)
mod_g.linear.3 = lm(preds_g.linear.3 ~ true_hgrid.linear.3)
# coverage
bkr.cb.3 = apply(fit_g.linear$h.hat[sel,], 2, quantile, c(0.025,0.975), na.rm = TRUE)
mat.res.bkr[3,] = c(summary(mod_g.linear.3)$coef[1,1], summary(mod_g.linear.3)$coef[2,1], summary(mod_g.linear.3)$r.sq, sqrt( sum( (preds_g.linear.3 - true_hgrid.linear.3)^2) / n), sum(true_hgrid.linear.3 >= bkr.cb.3[1,] & true_hgrid.linear.3 <= bkr.cb.3[2,] )/length(true_hgrid.linear.3))

######################################################
# Time 4
Z.Time = Z.4

source('postmeanbkmrcode.R')

preds_g.linear.4 = apply(fit_g.linear$h.hat[sel,], 2, mean)
mod_g.linear.4 = lm(preds_g.linear.4 ~ true_hgrid.linear.4)

# coverage
bkr.cb.4 = apply(fit_g.linear$h.hat[sel,], 2, quantile, c(0.025,0.975), na.rm = TRUE)
mat.res.bkr[4,] = c(summary(mod_g.linear.4)$coef[1,1], summary(mod_g.linear.4)$coef[2,1], summary(mod_g.linear.4)$r.sq, sqrt( sum( (preds_g.linear.4 - true_hgrid.linear.4)^2) / n), sum(true_hgrid.linear.4 >= bkr.cb.4[1,] & true_hgrid.linear.4 <= bkr.cb.4[2,] )/length(true_hgrid.linear.4))

#save(mat.res.bkr, file=paste0("Output/matresBKR", jobid, ".RData"))

######################################################
######################################################

Z.Time 					= cbind(Z.1, Z.2, Z.3, Z.4)
mat.res.bkr.1kernel = matrix(NA, nrow=1, ncol=5)

source('postmeanbkmrcode.R')

#################

mat.res.1kernel.all = matrix(NA, nrow=numtime, ncol=5)

# Time 1:

cross.sec = cbind(Z.1, z6=median(Z.Time[,6]), z7=median(Z.Time[,7]), z8=median(Z.Time[,8]), z9=median(Z.Time[,9]), z10=median(Z.Time[,10]), z11=median(Z.Time[,11]), z12=median(Z.Time[,12]), z13=median(Z.Time[,13]), z14=median(Z.Time[,14]), z15=median(Z.Time[,15]), z16=median(Z.Time[,16]), z17=median(Z.Time[,17]), z18=median(Z.Time[,18]), z19=median(Z.Time[,19]), z20=median(Z.Time[,20]))
preds_g.linear.1kernel.T1 <- newh.postmean(fit_g.linear, Znew = cross.sec, sel = sel)$postmean
preds_g.linear.1kernel.T1.sd 	= sqrt(diag(newh.postmean(fit_g.linear, Znew = cross.sec, sel = sel)$postvar))
mod_g.linear.1kernel.T1 = lm(preds_g.linear.1kernel.T1 ~ true_hgrid.linear.1)

bkr.cb.1kernel.low = preds_g.linear.1kernel.T1 - 1.96*preds_g.linear.1kernel.T1.sd
bkr.cb.1kernel.high = preds_g.linear.1kernel.T1 + 1.96*preds_g.linear.1kernel.T1.sd
mat.res.1kernel.all[1,] = c(summary(mod_g.linear.1kernel.T1)$coef[1,1], 0, 0, sqrt( sum( (preds_g.linear.1kernel.T1 -  true_hgrid.linear.1)^2) / n), sum( true_hgrid.linear.1 >= bkr.cb.1kernel.low &  true_hgrid.linear.1 <= bkr.cb.1kernel.high )/length( true_hgrid.linear.1))

# Time 2:

cross.sec = cbind(z1=median(Z.Time[,1]), z2=median(Z.Time[,2]), z3=median(Z.Time[,3]), z4=median(Z.Time[,4]), z5=median(Z.Time[,5]), Z.2, z11=median(Z.Time[,11]), z12=median(Z.Time[,12]), z13=median(Z.Time[,13]), z14=median(Z.Time[,14]), z15=median(Z.Time[,15]), z16=median(Z.Time[,16]), z17=median(Z.Time[,17]), z18=median(Z.Time[,18]), z19=median(Z.Time[,19]), z20=median(Z.Time[,20]))
preds_g.linear.1kernel.T2 <- newh.postmean(fit_g.linear, Znew = cross.sec, sel = sel)$postmean
preds_g.linear.1kernel.T2.sd 	= sqrt(diag(newh.postmean(fit_g.linear, Znew = cross.sec, sel = sel)$postvar))
mod_g.linear.1kernel.T2 = lm(preds_g.linear.1kernel.T2 ~ true_hgrid.linear.2)

bkr.cb.1kernel.low = preds_g.linear.1kernel.T2 - 1.96*preds_g.linear.1kernel.T2.sd
bkr.cb.1kernel.high = preds_g.linear.1kernel.T2 + 1.96*preds_g.linear.1kernel.T2.sd
mat.res.1kernel.all[2,] = c(summary(mod_g.linear.1kernel.T2)$coef[1,1], summary(mod_g.linear.1kernel.T2)$coef[2,1], summary(mod_g.linear.1kernel.T2)$r.sq, sqrt( sum( (preds_g.linear.1kernel.T2 -  true_hgrid.linear.2)^2) / n), sum( true_hgrid.linear.2 >= bkr.cb.1kernel.low &  true_hgrid.linear.2 <= bkr.cb.1kernel.high )/length( true_hgrid.linear.2))

# Time 3:

cross.sec = cbind(z1=median(Z.Time[,1]), z2=median(Z.Time[,2]), z3=median(Z.Time[,3]), z4=median(Z.Time[,4]), z5=median(Z.Time[,5]), z6=median(Z.Time[,6]), z7=median(Z.Time[,7]), z8=median(Z.Time[,8]), z9=median(Z.Time[,9]), z10=median(Z.Time[,10]), Z.3, z16=median(Z.Time[,16]), z17=median(Z.Time[,17]), z18=median(Z.Time[,18]), z19=median(Z.Time[,19]), z20=median(Z.Time[,20]))
preds_g.linear.1kernel.T3 <- newh.postmean(fit_g.linear, Znew = cross.sec, sel = sel)$postmean
preds_g.linear.1kernel.T3.sd 	= sqrt(diag(newh.postmean(fit_g.linear, Znew = cross.sec, sel = sel)$postvar))
mod_g.linear.1kernel.T3 = lm(preds_g.linear.1kernel.T3 ~ true_hgrid.linear.3)

bkr.cb.1kernel.low = preds_g.linear.1kernel.T3 - 1.96*preds_g.linear.1kernel.T3.sd
bkr.cb.1kernel.high = preds_g.linear.1kernel.T3 + 1.96*preds_g.linear.1kernel.T3.sd
mat.res.1kernel.all[3,] = c(summary(mod_g.linear.1kernel.T3)$coef[1,1], summary(mod_g.linear.1kernel.T3)$coef[2,1], summary(mod_g.linear.1kernel.T3)$r.sq, sqrt( sum( (preds_g.linear.1kernel.T3 -  true_hgrid.linear.3)^2) / n), sum( true_hgrid.linear.3 >= bkr.cb.1kernel.low &  true_hgrid.linear.3 <= bkr.cb.1kernel.high )/length( true_hgrid.linear.3))

# Time 4:

cross.sec = cbind(z1=median(Z.Time[,1]), z2=median(Z.Time[,2]), z3=median(Z.Time[,3]), z4=median(Z.Time[,4]), z5=median(Z.Time[,5]), z6=median(Z.Time[,6]), z7=median(Z.Time[,7]), z8=median(Z.Time[,8]), z9=median(Z.Time[,9]), z10=median(Z.Time[,10]), z11=median(Z.Time[,11]), z12=median(Z.Time[,12]), z13=median(Z.Time[,13]), z14=median(Z.Time[,14]), z15=median(Z.Time[,15]), Z.4)
preds_g.linear.1kernel.T4 <- newh.postmean(fit_g.linear, Znew = cross.sec, sel = sel)$postmean
preds_g.linear.1kernel.T4.sd 	= sqrt(diag(newh.postmean(fit_g.linear, Znew = cross.sec, sel = sel)$postvar))
mod_g.linear.1kernel.T4 = lm(preds_g.linear.1kernel.T4 ~ true_hgrid.linear.4)

bkr.cb.1kernel.low = preds_g.linear.1kernel.T4 - 1.96*preds_g.linear.1kernel.T4.sd
bkr.cb.1kernel.high = preds_g.linear.1kernel.T4 + 1.96*preds_g.linear.1kernel.T4.sd
mat.res.1kernel.all[4,] = c(summary(mod_g.linear.1kernel.T4)$coef[1,1], summary(mod_g.linear.1kernel.T4)$coef[2,1], summary(mod_g.linear.1kernel.T4)$r.sq, sqrt( sum( (preds_g.linear.1kernel.T4 -  true_hgrid.linear.4)^2) / n), sum( true_hgrid.linear.4 >= bkr.cb.1kernel.low &  true_hgrid.linear.4 <= bkr.cb.1kernel.high )/length( true_hgrid.linear.4))

#save(mat.res.1kernel.all, file=paste0("Output/matresBKR1kernelall", jobid, ".RData"))

colnames(mat.res.gfl) = colnames(mat.res.bkr) = colnames(mat.res.1kernel.all) = c("Intercept", "Slope", "Rsquared", "RMSE", "Coverage")

round(mat.res.gfl,2) #LKMR
round(mat.res.bkr,2) #BKMR
round(mat.res.1kernel.all,2) #JKBKMR