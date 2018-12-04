# solving carina
# need a lookup for f0(v)
hobj = hist(milkyway, breaks = 100, plot = FALSE)
# make a function by linearly interpolating between this fit
f0v = approx(hobj$mids, hobj$density, v)$y

# hist(milkyway, breaks = 50, freq =FALSE)
# points(v[order(v)], f0v[order(v)], cex = 0.2, col = 'blue', type = 'l')

# EM solution for carina
# using only R
emcarina = function(y, x, f0y, tol = 1e-6, maxit = 1e3, verbose = FALSE){
  require(Iso)
  require(Matrix)
  
  # sorting by x
  ox = order(x)
  yy = y[ox]
  f0yy = f0y[ox]
  iox = invPerm(ox)
  
  runFlag = 1
  itcount = 0
  w = rep(1/2, length(yy))
  if (verbose) pb = txtProgressBar(max = maxit, style = 3)
  while(runFlag){
    
    # updating f1
    mu = sum(w*yy)/sum(w)
    sig = sqrt(sum(w*(yy - mu)^2)/sum(w))
    # updating pi
    pivec = pava(w, decreasing = TRUE)
    # updating w
    wnew = pivec * dnorm(yy, mu, sig)/(pivec * dnorm(yy, mu, sig) + (1-pivec)*f0yy)
    
    itererr = mean(abs(wnew - w))
    runFlag = min(runFlag, !(itererr <= tol))
    runFlag = min(runFlag, !(itcount == maxit))
    
    w = wnew
    itcount = itcount + 1
    
    if (verbose) setTxtProgressBar(pb, itcount)
  }
  
  return(list(p = pivec[iox], mu = mu, sig = sig, w = w[iox]))
}

eo = emcarina(v, R, f0v, verbose = TRUE)
labels_cov = 1 + (eo$w <= 0.5)

# mean(abs(eo$w - postp))



plot(R[order(R)], eo$p[order(R)], type = 'l', ylim = c(0,1),xlab="Distance from center (X)", ylab="Estimated frequency of Carina Stars",cex.lab=1.5,cex.main=1.5)

# carina no covariates
emcarina_twog = function(y, f0y, tol = 1e-6, maxit = 1e3, verbose = FALSE){
  require(Iso)
  require(Matrix)
  
  runFlag = 1
  itcount = 0
  
  # initial values
  w = rep(1/2, length(y))
  if (verbose) pb = txtProgressBar(max = maxit, style = 3)
  while(runFlag){
    
    # updating p
    p = mean(w)
    
    # updating f1
    mu = sum(w*y)/sum(w)
    sig = sqrt(sum(w*(y - mu)^2)/sum(w))
    # updating w
    wnew = p * dnorm(y, mu, sig)/(p * dnorm(y, mu, sig) + (1-p)*f0y)
    
    itererr = mean(abs(wnew - w))
    runFlag = min(runFlag, !(itererr <= tol))
    runFlag = min(runFlag, !(itcount == maxit))
    
    w = wnew
    itcount = itcount + 1
    
    if (verbose) setTxtProgressBar(pb, itcount)
  }
  
  return(list(p = p, mu = mu, sig = sig, w = w))
}

eo2g = emcarina_twog(v, f0v)
labels_2g = 1 + (eo2g$w <= 0.5)
abline(h=eo2g$p,col="red")
legend(x="topright", legend = c('Covariates model', 'Two-groups'),col=c('blue', 'red'), lty = 1, bty = 'n',lwd=4,cex=2)


# difference in the two solutions
sum(labels_cov != labels_2g)
labels_cov[which(labels_cov != labels_2g)]
labels_2g[which(labels_cov != labels_2g)]

cols = c('red', 'blue')

pdf('carina_class_cov.pdf', width = 7, height = 7)
par(bty = 'l')
plot(R, v, cex = 0.7, col = cols[labels_cov],
     xlab = 'Distance from Center (X)',
     ylab = 'Line-of-Sigth Velocity (Y)',
     main = 'Two-groups with Covariates',pch=20,cex.lab=1.5,cex.main=1.5)
legend(x = 'topleft', legend = c('Carina', 'Milky Way'), col = cols, lwd = 4, bty = 'n',cex=2)
dev.off()

pdf('carina_class_2g.pdf', width = 7, height = 7)
par(bty = 'l')
plot(R, v, cex = 0.7, col = cols[labels_2g],
     xlab = 'Distance from Center (X)',
     ylab = 'Line-of-Sigth Velocity (Y)',
     main = 'Two-groups',pch=20,cex.lab=1.5,cex.main=1.5)
legend(x = 'topleft', legend = c('Carina', 'Milky Way'), col = cols, lwd = 4, bty = 'n', cex=2)
dev.off()

pdf('carina_pi_est.pdf', width = 7, height = 7)
par(bty = 'l')
plot(R[order(R)], eo$p[order(R)], type = 'l',
     xlab = 'Distance from Center (X)',
     ylab = 'Estimated Frequency of Carina')
abline(h = eo2g$p, col = 'red')
legend(x = 'topright', legend = c('Covariates', 'Two Groups'), col = c('black', 'red'), lwd = 2, bty = 'n')
dev.off()

####### All Plots #######
fmlelfdr=eo$w
scottlfdr=eo2g$w
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
plot(v,1-eo$w,pch=20,col="blue",cex=0.8,xlab="Test statistic",ylab="lFDR",cex.lab=1.5)
points(v,1-eo2g$w,pch=20,col="red",cex=0.8)
abline(h=sl1[k1],lty=2,col="red")
abline(h=sl2[k2],lty=2,col="blue")
legend(x = "topright", legend = c('2g lFDR', 'fMLE lFDR'), col = c('blue', 'red'), lty = 1, lwd = 3, bty = 'n',cex=2)

### 2 ###
ttx=seq(0,1,0.001)
tty=seq(0,1,0.001)
plot(fmlelfdr,scottlfdr,col="red",pch=20,cex=0.8,cex.lab=1.5)
lines(ttx,tty,col="blue",lty=2,lwd=4)
ind=which(fmlelfdr>scottlfdr)
points(fmlelfdr[ind],scottlfdr[ind],pch=20,cex=0.8)
legend(x = "topleft", legend = c('y=x line'), col = c('blue'), lty = 2, lwd=4,  bty = 'n',cex=2)

### 3 ###

hist(v, breaks = 100, prob=TRUE, col='lightgrey', border='grey',
     main='Fitted marginals for 2 groups', xlab='line of sight velocity',xlim=c(-100,350),ylim=c(0,0.03),cex.lab=1.5,cex.main=1.5)
rug(v)
m1=eo2g$p*sapply(((v-eo2g$mu)/(eo2g$sig)),dnorm)
m2=(1-eo2g$p)*f0v
m3=m1+m2
points(v[order(v)], m1[order(v)], col='red',ty='l')
points(v[order(v)], m2[order(v)], col='blue',ty='l',lty=2)
points(v[order(v)], m3[order(v)],ty='l',lwd=2)
legend(x = 2.5, y=0.3, bty = 'n',
       legend = c(expression(bar(pi)~f[1] + (1-bar(pi))~f[0]),
                  expression(bar(pi)~f[1]),
                  expression((1-bar(pi))~f[0])),
       col = c('black', 'red', 'blue'), lty= 1,lwd=2,cex=1.3)
# lets plot the histogram overlapped with scaled versions of f0 and f1
#pdf('histogram_fitted_densities_new.pdf', width = 10, height = 10)

### 4 ###

hist(v, breaks = 100, prob=TRUE, col='lightgrey', border='grey',
     main='Fitted marginals for covariate model', xlab='line of sight velocity',xlim=c(-100,350),ylim=c(0,0.03),cex.lab=1.5,cex.main=1.5)
rug(v)
m1=eo$p*sapply(((v-eo$mu)/(eo$sig)),dnorm)
m2=(1-eo$p)*f0v
m3=m1+m2
points(v[order(v)], m1[order(v)], col='red',ty='l')
points(v[order(v)], m2[order(v)], col='blue',ty='l',lty=2)
points(v[order(v)], m3[order(v)],ty='l',lwd=2)
legend(x = 2.5, y=0.3, bty = 'n',
       legend = c(expression(bar(pi)~f[1] + (1-bar(pi))~f[0]),
                  expression(bar(pi)~f[1]),
                  expression((1-bar(pi))~f[0])),
       col = c('black', 'red', 'blue'), lty= 1,lwd=2,cex=1.3)
# lets plot the histogram overlapped with scaled versions of f0 and f1
#pdf('histogram_fitted_densities_new.pdf', width = 10, height = 10)

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


### new code ###

dim(x)
y1=sort(x[,1])
y2=x[,4][order(x[,1])]
ourdat=cbind(y1,y2)
#ourdat=x[,c(1,4)]
ubd=5
lbd=1
l=length(ourdat[,1])
w=rep(0.5,l)
niter=50
f0v=f0v[order(x[,1])]
pivec=NULL
denvec=NULL
for(i in 1:niter)
{
 pivec=pava(w,dec=T)
 f=function(p)
 {
   sst=max(log(p[1]*sapply((ourdat[,2]-p[2])/p[3], dnorm)+(1-p[1])*sapply((ourdat[,2]-p[2])/p[4],dnorm)),-200)
   #return(-mean(log(sst)))
   return(-mean(sst))
 }
 est=optim(c(0.3,210,2,2),f,method="L-BFGS-B",lower=c(0.02,200,lbd,lbd),upper=c(0.98,250,ubd,ubd),hessian=FALSE)$par
 denvec=est[1]*sapply((ourdat[,2]-est[2])/est[3], dnorm)+(1-est[1])*sapply((ourdat[,2]-est[2])/est[4],dnorm)
 w=pivec*denvec+(1-pivec)*f0v
 print(i)
}
est
plot(ourdat[,1],pivec)
sst=ii[1]*sapply((ourdat[,2]-ii[2])/ii[3], dnorm)+(1-ii[1])*sapply((ourdat[,2]-ii[2])/ii[4],dnorm)
