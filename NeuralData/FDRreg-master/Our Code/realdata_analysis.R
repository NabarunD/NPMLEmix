# gather results
# setwd('sim_output/')
# filelist = list.files(pattern = 'sim_realdata')
# reslist = list()
# loadenv = new.env()
# for (i in filelist){
#   load(i, envir = loadenv)
#   reslist = append(reslist, list(as.list(loadenv, all.names = TRUE)))
# }

# setwd('..')
require(FDRreg)
source('logistic_gaussian_updated.R')
source('altmax.R')
source('altmax_logistic_decreasing.R')
# results on real data
synchrony_smithkohn2008 = read.csv('../data/synchrony_smithkohn2008.csv', header=TRUE)

ddfull = cbind(synchrony_smithkohn2008$Dist, synchrony_smithkohn2008$TuningCor, synchrony_smithkohn2008$z)
ddfull = as.data.frame(ddfull)
names(ddfull) = c('Dist', 'TuningCor', 'z')

# the null distribution is not modeled as N(0,1) but rather as N(mu,sig)
# where mu and sig are computed using Efron's ML estimate in a two-groups model
efron_res = FDRreg::efron(ddfull$z, nulltype='empirical', df=10)
mu0hat = efron_res$mu0; sig0hat = efron_res$sig0;

# using Scott et al's Empirical Bayes FDR Regression
X = ddfull[,1:2]
# Set up spline basis functions (df=3)
df = 3
b1 = bs(synchrony_smithkohn2008$Dist, df=df)
b2 = bs(synchrony_smithkohn2008$TuningCor, df=df)
Xs = model.matrix( ~  b1 + b2 - 1)
# produces a matrix with 6 column, first three are splines corresponding to Dist, next three corr. to TuningCor
# FDRreg analysis uses a scaled and centered Xs
# lets use the same
Xs = scale(Xs)

# for our methods, need to rescale z so that nulls are N(0,1)
zscaled = (ddfull$z - mu0hat)/sig0hat
lines(scott_res$x_grid, (1-scott_res$pi0)*scott_res$f1_grid, col='blue', lwd=2, lty='dotted')
##########################################
# running all methods

# running scott's algorithm
scott_res = FDRreg(ddfull$z, Xs, nulltype='empirical', control=list(lambda = 1))
#scott_res = FDRreg(newdfull, Xs, nulltype='empirical', control=list(lambda = 1))
# running our methods
# for our methods, need to rescale z so that nulls are N(0,1)
zscaled = (ddfull$z - mu0hat)/sig0hat
# marginal method 1
gridy_res = gridy(zscaled, Xs, verbose = TRUE)
# fullmle starting from marginal1
am1_res = lgamd(zscaled, cbind(1, Xs), gridy_res$b, verbose = TRUE)

# running altmax with marginal method 2 estimate of beta
m2_res = bmlogist(zscaled, cbind(1, Xs), verbose = TRUE)
am2_res = lgamd(zscaled, cbind(1, Xs), m2_res$bt, verbose = TRUE)

# save all results
save(ddfull, scott_res, gridy_res, am1_res, m2_res, am2_res, file = 'save_new.Rdata')

##########################################################
# finished running and saving all results

# loading saved results for analysis
# load saved results
load('save_new.Rdata')
########################################################
# adding one step kwmle to marginal2 - dont need to run this again
m2p = 1/(1 + exp(- cbind(1, Xs) %*% m2_res$bt))
m2_onestep = kwmle(zscaled, 1 - m2p, maxit = 1e2)
m2_localfdr = (1 -m2p) * dnorm(zscaled)/((1 -m2p) * dnorm(zscaled) + m2p * m2_onestep$f1u)

m2_res = list(b = m2_res$bt, mu = m2_res$mu, p = m2p, 
              atoms = m2_onestep$atoms, probs = m2_onestep$probs, localfdr = m2_localfdr, f1u = m2_onestep$f1u)

########################################################

# plot to see how closely scott and our estimates of pi match up
# scott is actually using a ridge regression coefficient in the fit
pdf('pi_e_new.pdf', width = 10, height = 10)
par(mar = c(5.1, 5.1, 2.1, 2.1))

xtb = as.vector(cbind(1, Xs) %*% am1_res$b)
xtbs = sort(xtb, index.return = TRUE)
#xtbs = sort(ddfull$z, index.return = TRUE)
plot(xtbs$x, am1_res$p[xtbs$ix], type = 'l', ylim = c(0,1), col = 'darkcyan',
     ylab = expression(hat(pi)), xlab = expression(x^T~hat(beta)),cex.lab=1,cex.axis=1,cex.main=1,ps=12)
lines(xtbs$x, am2_res$p[xtbs$ix], col = 'blue')
lines(xtbs$x, scott_res$priorprob[xtbs$ix], col = 'red')
# write a legend for this plot
legend(x = 'topleft',
       legend = c('scott', 'marginal1+fullmle', 'marginal2+fullmle'),
       col = c('red', 'darkcyan', 'blue'), lty = 1, bty = 'n',cex=2,pt.cex=1)

legend(x = 'topright',
       legend = c('lFDRs (Scott)', 'lFDRs (fMLE)'),
       col = c('blue', 'black'), lty = 1, bty = 'n',cex=1.3,pt.cex=1)

dev.off()


####################################################################################
# fullmle gives the same answer when starting from either initial estimate
# that is encouraging

# functions for rejection
bh = function(q, alpha){
  n = length(q)
  qs = sort(q, index.return = TRUE)
  bestt = tail(which(qs$x <= alpha*(1:n)/n), 1)
  if(length(bestt)){
    rejects = qs$ix[1:bestt]
  } else {
    rejects = integer(0)
  }
  return(rejects)
}

lfdr_reject = function(lfdr, fdr_nominal = 0.1){
  sl = sort(lfdr)
  k = sum(cumsum(sl)/seq_along(sl) <= fdr_nominal)
  if(k) rejects = which(lfdr <= sl[k]) else rejects = numeric(0)
}

fdr_nominal=0.1
lfdr=fmlelfdr
sl = sort(lfdr)
k = sum(cumsum(sl)/seq_along(sl) <= fdr_nominal)
sl[k]
plot(ddfull$z,scottlfdr,pch=16,col="blue",cex=0.01,xlim=c(-2.5,7))
points(ddfull$z,fmlelfdr,pch=16,col="red",cex=0.01)
abline(h=sl[k],lty=2,col="red")
# we can verify this by checking differences between solutions
# compute f1 and lfdr for each solution
# marginal1+fullmle
am1_f1_zscaled = sapply(am1_res$atoms, function(ai) dnorm(zscaled - ai)) %*% am1_res$probs
am1_f11_zscaled = sapply(am1_res$atoms, function(ai) (1/sig0hat)*dnorm((zscaled - ai - mu0hat)/(sig0hat))) %*% am1_res$probs
# computing local fdr based on am1
am1_lfdr = (1 - am1_res$p) * dnorm(zscaled)/((1 - am1_res$p) * dnorm(zscaled) + am1_res$p * am1_f1_zscaled)
# marginal2+fullmle
am2_f1_zscaled = sapply(am2_res$atoms, function(ai) dnorm(zscaled - ai)) %*% am2_res$probs
# computing local fdr based on am1
am2_lfdr = (1 - am2_res$p) * dnorm(zscaled)/((1 - am2_res$p) * dnorm(zscaled) + am2_res$p * am2_f1_zscaled)

# difference between estimates of pi
mean(abs(am1_res$p - am2_res$p))
# difference between estimates of f1
mean(abs(am1_f1_zscaled - am2_f1_zscaled))
# difference between estimates of localfdr
mean(abs(am1_lfdr - am2_lfdr))

# getting rejection sets
am1_rejects = lfdr_reject(am1_lfdr)
am2_rejects = lfdr_reject(am2_lfdr)
# the same neuronal pairs are rejected
setdiff(am1_rejects, am2_rejects)
identical(am1_rejects, am2_rejects)

am1r2 = intersect(am1_rejects, which(ddfull$z > 0))
am2r2 = intersect(am2_rejects, which(ddfull$z > 0))

hist(zscaled[both_methods],col="grey",breaks=80,xlim=c(-0.5,8),xlab="Test statistic",cex.lab=1.3,cex.main=1.5,main="Histogram on rejection")
hist(zscaled[only_am1],col="red",add=T,xlim=c(-0.5,8))
hist(zscaled[only_scott],col="blue",add=T,xlim=c(-0.5,8))

lines((scott_res$x_grid-mu0hat)/(sig0hat), (1-scott_res$p0)*scott_res$f1_grid, col='green',lwd=5)
lines((scott_res$x_grid-mu0hat)/(sig0hat), scott_res$p0*scott_res$f0_grid, col='violet',lwd=5)
lines((scott_res$x_grid-mu0hat)/(sig0hat), scott_res$p0*scott_res$f0_grid+(1-scott_res$p0)*scott_res$f1_grid,col="yellow",lwd=5)

####################################################################################
# plot showing extra rejection under likelihood methods

#I am unable to find scott_rejects
#Can I get it in the following way?
scott_rejects = which(getFDR(scott_res$postprob)$FDR <= 0.1 & ddfull$z > 0) #lfdr_reject(1 - scott_res$postprob)

# let's make the plot showing extra rejections
only_scott = setdiff(scott_rejects, am1_rejects)
both_methods = intersect(scott_rejects, am1_rejects)
only_am1 = setdiff(am1_rejects, scott_rejects)

# need to bin covariate space to create matrices
# lets pick bl bins on each axes
bl = 100

x1breaks = seq(from = min(X[,1]) - 1e-2, to = max(X[,1]) + 1e-2, length.out = bl)
x2breaks = seq(from = min(X[,2]) - 1e-2, to = max(X[,2]) + 1e-2, length.out = bl)

mat_only_scott = matrix(0, bl, bl)
for (i in only_scott){
  p1 = sum(X[i,1] >= x1breaks)
  p2 = sum(X[i,2] >= x2breaks)
  mat_only_scott[p1, p2] = mat_only_scott[p1, p2] + 1
}

mat_both_methods = matrix(0, bl, bl)
for (i in both_methods){
  p1 = sum(X[i,1] >= x1breaks)
  p2 = sum(X[i,2] >= x2breaks)
  mat_both_methods[p1, p2] = mat_both_methods[p1, p2] + 1
}


mat_only_am1 = matrix(0, bl, bl)
for (i in only_am1){
  p1 = sum(X[i,1] >= x1breaks)
  p2 = sum(X[i,2] >= x2breaks)
  mat_only_am1[p1, p2] = mat_only_am1[p1, p2] + 1
}

colors1 = colorRampPalette(c('red','red4'))
colors2 = colorRampPalette(c('cyan','darkblue'))
colors3 = colorRampPalette(c('greenyellow','green4'))

temp_fun = function(mat) length(unique(as.vector(mat))) - 1

colors_list = c('white', colors1(temp_fun(mat_only_scott)),
                colors2(temp_fun(mat_both_methods)),
                colors3(temp_fun(mat_only_am1)))

# combine all matrices into one
mat_combined = mat_only_scott
mat_combined[which(mat_both_methods != 0)] = 10 + mat_both_methods[which(mat_both_methods != 0)]
mat_combined[which(mat_only_am1 != 0)] = 20 + mat_only_am1[which(mat_only_am1 != 0)]
mat_combined[which(mat_only_scott != 0)] = mat_only_am1[which(mat_only_scott != 0)]

pdf('extra_rejections_new.pdf', width = 10, height = 10)
image(x1breaks, x2breaks, mat_combined, xlab = 'Dist', ylab = 'TuningCor',
      col = colors_list,cex.lab=1.5,cex.axis=1.5)
legend(x = 2200, y=0, legend = c('both methods', 'only fullmle'), col = c('blue', 'green4'), lty = 1, bty = 'n',cex=2)
#legend(x = 'topright', legend = c('both methods', 'only fullmle','only scott'), col = c('grey','red','blue'), lty = 1, bty = 'n',cex=2)
dev.off()


# histogram of data with scaled versions of f0 and f1
###################################################################################

hist(zscaled, 200, prob=TRUE, col='lightgrey', border='grey', xlim=c(-3,6),
     main='Fitted marginals for Scott', xlab='Test statistic',cex.lab=1.3,cex.main=1.3)
rug(zscaled)
lines((scott_res$x_grid-mu0hat)/(sig0hat), (1-scott_res$p0)*scott_res$f1_grid, col='red')
lines((scott_res$x_grid-mu0hat)/(sig0hat), scott_res$p0*scott_res$f0_grid, col='blue')
lines((scott_res$x_grid-mu0hat)/(sig0hat), scott_res$p0*scott_res$f0_grid+(1-scott_res$p0)*scott_res$f1_grid)
legend(x = 2.5, y=0.3, bty = 'n',
       legend = c(expression(bar(pi)~f[1] + (1-bar(pi))~f[0]),
                  expression(bar(pi)~f[1]),
                  expression((1-bar(pi))~f[0])),
       col = c('black', 'red', 'blue'), lty= 1,cex=1.3)
# lets plot the histogram overlapped with scaled versions of f0 and f1

hist(ddfull$z, 200, prob=TRUE, col='lightgrey', border='grey', xlim=c(-2,6),
     main='Fitted marginals for Scott', xlab='Test statistic (no scaling)',cex.lab=1.3,cex.main=1.3)
rug(ddfull$z)
lines(scott_res$x_grid, (1-scott_res$p0)*scott_res$f1_grid, col='red')
lines(scott_res$x_grid, scott_res$p0*scott_res$f0_grid, col='blue')
lines(scott_res$x_grid, scott_res$p0*scott_res$f0_grid+(1-scott_res$p0)*scott_res$f1_grid)
legend(x = 2.5, y=0.3, bty = 'n',
       legend = c(expression(bar(pi)~f[1] + (1-bar(pi))~f[0]),
                  expression(bar(pi)~f[1]),
                  expression((1-bar(pi))~f[0])),
       col = c('black', 'red', 'blue'), lty= 1,cex=1.3)
# lets plot the histogram overlapped with scaled versions of f0 and f1
pdf('histogram_fitted_densities_new.pdf', width = 10, height = 10)

hist(ddfull$z, freq=F, breaks=200, col = 'lightgrey', border = 'grey',xlim=c(-2,6), xlab="Test statistic (no scaling)",
     pr=T, main="Fitted marginals for fMLE",cex.main=1.3,cex.lab=1.3)
rug(ddfull$z)
curve((1 - mean(am1_res$p)) * dnorm((x-mu0hat)/sig0hat) * (1/sig0hat), from = min(ddfull$z), to = max(ddfull$z),
      col = 'red', add = TRUE)
points(ddfull$z, mean(am1_res$p) * am1_f11_zscaled, col = 'blue', cex = 0.05)
# also showing the mixture of these two - but this is not necessarily the correct notion of fit
points(ddfull$z, (1 - mean(am1_res$p)) * dnorm((ddfull$z-mu0hat)/(sig0hat))*(1/sig0hat) + mean(am1_res$p) * am1_f11_zscaled,
       cex = 0.05)
legend(x = 2.5, y=0.3, bty = 'n',
       legend = c(expression(bar(pi)~f[1] + (1-bar(pi))~f[0]),
                  expression(bar(pi)~f[1]),
                  expression((1-bar(pi))~f[0])),
       col = c('black', 'blue', 'red'), lty= 1,cex=1.3)

hist(zscaled, breaks=200, freq = FALSE, col = 'lightgrey', border = 'grey',xlim=c(-3,6), xlab="Test statistic", pr=T, main="Fitted marginals for fMLE",cex.main=1.3,cex.lab=1.3)
rug(zscaled)
curve((1 - mean(am1_res$p)) * dnorm(x), from = min(zscaled), to = max(zscaled),
      col = 'red', add = TRUE)
points(zscaled, mean(am1_res$p) * am1_f1_zscaled, col = 'blue', cex = 0.05)
# also showing the mixture of these two - but this is not necessarily the correct notion of fit
points(zscaled, (1 - mean(am1_res$p)) * dnorm(zscaled) + mean(am1_res$p) * am1_f1_zscaled,
       cex = 0.05)
legend(x = 2.5, y=0.3, bty = 'n',
       legend = c(expression(bar(pi)~f[1] + (1-bar(pi))~f[0]),
                                             expression(bar(pi)~f[1]),
                                             expression((1-bar(pi))~f[0])),
      col = c('black', 'blue', 'red'), lty= 1,cex=1.3)

dev.off()

jdistfmle=(1-am1_res$p) * dnorm((ddfull$z-mu0hat)/(sig0hat))*(1/sig0hat) + am1_res$p * am1_f11_zscaled
fmlelik=-mean(jdistfmle)
fmlelfdr=((1-am1_res$p) * dnorm((ddfull$z-mu0hat)/(sig0hat))*(1/sig0hat))/jdistfmle
jdistscott=scott_res$priorprob*scott_res$M1+(1-scott_res$priorprob)*scott_res$M0
scottlik=-mean(jdistscott)
scottlfdr=scott_res$localfdr
plot(ddfull$z,fmlelfdr,pch=16,cex=0.01,col="blue",xlim=c(-2,7))
points(ddfull$z,scottlfdr,pch=16,cex=0.01,col="red")
##############################################################
# making the regression surface
regres = lm(ddfull$z ~ Xs)
zfitted = regres$fitted.values
if(.Platform$OS.type=="windows") {
  quartz<-function() windows()
}
quartz()
# marginal plots
plot(X[,1], zfitted)
plot(X[,2], zfitted)

# plotting the 3d regression surface
x1grid = seq(from = min(X[,1]), to = max(X[,1]), length.out = 1e2)
x2grid = seq(from = min(X[,2]), to = max(X[,2]), length.out = 1e2)
newdatamatrix = matrix(0, length(x1grid) * length(x2grid),2)
for (i in 1:length(x1grid)){
  newdatamatrix[(i-1)*length(x2grid) + 1:length(x2grid), ] = cbind(x1grid[i], x2grid)
}

ndms = cbind(1, bs(newdatamatrix[,1], df = 3), bs(newdatamatrix[,2], df = 3))
zndms = ndms %*% regres$coefficients

# shows clear dependence on x1 and x2
# pdf('reg_surface.pdf', width = 10, height = 10)
persp(x1grid, x2grid, matrix(zndms, ncol = length(x1grid), byrow = T),
      xlab = 'distance', ylab = 'correlation', zlab = 'z',
      theta = 40, phi = 20, col = 'blue')
# dev.off()
# finished the regression surface
###############################################################







# save(synchrony_smithkohn2008, Xs,
#      efron_res, scott_res, gridy_res, m2_res, am2_res,
#      file = 'realdata_analysis.Rdata')

# load('realdata_analysis.Rdata')
plotFDR(scott_res, breaks=150, mar=c(2,4,1,1))


# differences in esitmate of pi
mean((scott_res$priorprob - gridy_res$p)^2)
mean((scott_res$priorprob - am2_res$p)^2)

# differences in local fdr
mean((scott_res$localfdr - gridy_res$localfdr)^2)
mean((scott_res$localfdr - am2_res$localfdr)^2)

scott_rejects = lfdr_reject(scott_res$localfdr)
gridy_rejects = lfdr_reject(gridy_res$localfdr)
am2_rejects = lfdr_reject(am2_res$localfdr)

length(scott_rejects)
length(gridy_rejects)
length(am2_rejects)

length(intersect(gridy_rejects, scott_rejects))
length(intersect(gridy_rejects, am2_rejects))
length(intersect(scott_rejects, am2_rejects))

####### All Plots #######

### 1 ###

fdr_nominal=0.1
lfdr1=fmlelfdr
sl1 = sort(lfdr1)
k1 = sum(cumsum(sl1)/seq_along(sl1) <= fdr_nominal)
sl1[k1]
sl1vec=which(fmlelfdr<sl1[k1])
lfdr2=scott_res$localfdr
sl2 = sort(lfdr2)
k2 = sum(cumsum(sl2)/seq_along(sl2) <= fdr_nominal)
sl2[k2]
sl2vec=which(scott_res$localfdr<sl2[k2])
plot(ddfull$z,scottlfdr,pch=20,col="blue",cex=0.8,xlim=c(-2.5,7),xlab="Test statistic",ylab="lFDR",cex.lab=1.5)
points(ddfull$z,fmlelfdr,pch=20,col="red",cex=0.8)
abline(h=sl1[k1],lty=2,col="red")
abline(h=sl2[k2],lty=2,col="blue")
legend(x = "topright", legend = c('Scott lFDR', 'fMLE lFDR'), col = c('blue', 'red'), lty = 1, lwd = 3, bty = 'n',cex=2)

### 2 ###
ttx=seq(0,1,0.001)
tty=seq(0,1,0.001)
plot(fmlelfdr,scottlfdr,col="red",pch=20,cex=0.8,cex.lab=1.5)
lines(ttx,tty,col="blue",lty=2,lwd=4)
ind=which(fmlelfdr>scottlfdr)
points(fmlelfdr[ind],scottlfdr[ind],pch=20,cex=0.8)
legend(x = "topleft", legend = c('y=x line'), col = c('blue'), lty = 2, lwd=4,  bty = 'n',cex=2)

### 3 ###

hist(ddfull$z, 200, prob=TRUE, col='lightgrey', border='grey', xlim=c(-2,6),
     main='Fitted marginals for FDRreg', xlab='Test statistic (no scaling)',cex.lab=1.5,cex.main=1.5)
rug(ddfull$z)
lines(scott_res$x_grid, (1-scott_res$p0)*scott_res$f1_grid, col='blue')
lines(scott_res$x_grid, scott_res$p0*scott_res$f0_grid, col='red')
lines(scott_res$x_grid, scott_res$p0*scott_res$f0_grid+(1-scott_res$p0)*scott_res$f1_grid)
legend(x = 2.5, y=0.3, bty = 'n',
       legend = c(expression(bar(pi)~f[1] + (1-bar(pi))~f[0]),
                  expression(bar(pi)~f[1]),
                  expression((1-bar(pi))~f[0])),
       col = c('black', 'red', 'blue'), lty= 1,lwd=2,cex=1.3)
# lets plot the histogram overlapped with scaled versions of f0 and f1
#pdf('histogram_fitted_densities_new.pdf', width = 10, height = 10)

### 4 ###

hist(ddfull$z, freq=F, breaks=200, col = 'lightgrey', border = 'grey',xlim=c(-2,6), xlab="Test statistic (no scaling)",
     pr=T, main="Fitted marginals for fMLE",cex.main=1.5,cex.lab=1.5)
rug(ddfull$z)
curve((1 - mean(am1_res$p)) * dnorm((x-mu0hat)/sig0hat) * (1/sig0hat), from = min(ddfull$z), to = max(ddfull$z),
      col = 'red', add = TRUE)
points(ddfull$z, mean(am1_res$p) * am1_f11_zscaled, col = 'blue', cex = 0.05)
# also showing the mixture of these two - but this is not necessarily the correct notion of fit
points(ddfull$z, (1 - mean(am1_res$p)) * dnorm((ddfull$z-mu0hat)/(sig0hat))*(1/sig0hat) + mean(am1_res$p) * am1_f11_zscaled,
       cex = 0.05)
legend(x = 2.5, y=0.3, bty = 'n',
       legend = c(expression(bar(pi)~f[1] + (1-bar(pi))~f[0]),
                  expression(bar(pi)~f[1]),
                  expression((1-bar(pi))~f[0])),
       col = c('black', 'blue', 'red'), lty= 1, lwd=2,cex=1.3)

### 5 ###

plot(xtbs$x, am1_res$p[xtbs$ix], type = 'l', ylim = c(0,1), col = 'darkcyan',
     ylab = expression(hat(pi)), xlab = expression(x^T~hat(beta)),cex.lab=1.1,cex.axis=1,cex.main=1,ps=12)
lines(xtbs$x, am2_res$p[xtbs$ix], col = 'blue')
lines(xtbs$x, scott_res$priorprob[xtbs$ix], col = 'red')
# write a legend for this plot
legend(x = 'topleft',
       legend = c('scott', 'marginal1+fullmle', 'marginal2+fullmle'),
       col = c('red', 'darkcyan', 'blue'), lty = 1, bty = 'n',cex=2,pt.cex=1)

### 6 ###

image(x1breaks, x2breaks, mat_combined, xlab = 'Dist', ylab = 'TuningCor',
      col = colors_list,cex.lab=1.5,cex.axis=1.5)
legend(x = 2200, y=0.7, legend = c('both methods', 'only fullmle'), col = c('blue', 'green4'), lty = 1, bty = 'n',lwd=3,cex=2)

### 7 ###

os = setdiff(sl2vec,sl1vec)
bm = intersect(sl2vec,sl1vec)
oa = setdiff(sl1vec,sl2vec)
hist(ddfull$z[bm],col="lightgrey",breaks=50,xlim=c(1,6),xlab="Test statistic",cex.lab=1.3,cex.main=1.5,
     ylab="Frequency of rejection",main="Histogram on rejection")
hist(ddfull$z[oa],col="red",breaks=50,add=T,xlim=c(1,6))
hist(ddfull$z[os],col="blue",breaks=50,add=T,xlim=c(1,6))
legend(x="topright", legend = c('both methods', 'only fullmle','only FDRreg'),col=c('lightgrey', 'red','blue'), lty = 1, bty = 'n',lwd=4,cex=2)
