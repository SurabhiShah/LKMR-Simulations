# postmean bkmr code

set.seed(101)

lm.pdi0 <- lm(Y ~ X.Cov)
lm.pdi <- lm(Y ~ Z.Time + X.Cov)

zranges <- diff(apply(Z.Time, 2, range))
Drange <- max(zranges)
r.prior.params <- data.frame(alpha.r=0.1, DD=Drange, qbound1=c(1,1,1/2), qbound2=c(1/8,1/4,1/8))
r.prior.params$rbound1 <- with(r.prior.params, -log(1-0.50)/(qbound1*DD)^2)
r.prior.params$rbound2 <- with(r.prior.params, -log(1-0.50)/(qbound2*DD)^2)
r.prior.params$mu.r <- 0.15 
r.prior.params$sigma.r <- 0.1

fit_g.linear <- with(r.prior.params[1,], kmbayes(nsamp=max(sel), y=Y, X=X.Cov, Z=Z.Time, quiet=FALSE, starting.values=list(r=0.05, lambda=10, sigsq.eps=summary(lm.pdi)$sigma^2, beta=coef(lm.pdi)[paste0("X",colnames(X.Cov))]), control.params=list(mu.r=mu.r, sigma.r=sigma.r, mu.lambda=100, sigma.lambda=100, r.jump1=2, r.jump2=0.05, r.muprop=1, lambda.jump=15), modsel=FALSE))
