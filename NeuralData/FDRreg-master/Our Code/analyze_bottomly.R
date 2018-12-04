# require(pbapply)
# 
# reslist = pbreplicate(1e4, {
#   madd = arima.sim(list(ma = c(1.5, 0.75)), n = 100)
#   mleres = arima(madd, order = c(0,0,2), include.mean = FALSE)
#   mleres$coef
# })
# 
# reslist = t(reslist)
# 
# colMeans(reslist)
# cov(reslist)
# 
# 
# hist(reslist[,2], breaks = 100, freq = FALSE)
# curve(dnorm(x, mean(reslist[,2]), sd(reslist[,2])), col = 'blue', add = TRUE)

source('../NeuralData/FDRreg-master/Our Code/logistic_gaussian_updated.R')
source('../NeuralData/FDRreg-master/Our Code/altmax.R')
source('../NeuralData/FDRreg-master/Our Code/altmax_logistic_decreasing.R')

ds1 = read.delim('datasetS1.txt', header = TRUE)

# let's remove gene's with low read count
yy = ds1$p.value[which(!ds1$Low_read_count)]
xx = ds1$logConc[which(!ds1$Low_read_count)]

# 13 breaks so that there are equal number of hypotheses in each bin
breakhere = (1:13)/13
xbreaks = quantile(xx, breakhere)

# scatter plot
plot(xx, yy, col = rgb(0,0,0, alpha = 0.1), cex = 0.2)
npr = ksmooth(xx, yy, bandwidth = 2)
lines(npr, col = 'blue')
abline(v = xbreaks, col = 'red')

# let's try to do logistic with dummy variables for each dummy
# remember to set intercept to FALSE
cs = c(0, cumsum(table(sapply(xx, function (xi) sum(xi <= xbreaks)))))
dummyx = matrix(0, nrow = length(xx), ncol = 13)
for (i in 1:13){dummyx[(cs[i]+1):cs[i+1],i] = 1}

gridy_logistic_decreasing = function(u, x, intercept = TRUE, verbose = FALSE){
  
  # intercept is added to x by default
  n = length(u)
  if(!is.matrix(x)) x = matrix(x, ncol = 1)
  #x = matrix(ns(x, df = df), ncol = ncol(x) * df)
  if(intercept) x = cbind(rep(1, n), x)
  
  # construct a small grid
  ptcutoff = cutoff_estimate(u, punif)$pt
  
  print('got cutoff')
  
  spread = n ^ -0.25
  resolution = n ^ -0.5
  ptlen = 2 * spread/resolution + 1
  smallgrid = seq(from = ptcutoff - spread,
                  to = ptcutoff + spread,
                  length.out = ptlen)
  # if smallgrid contains points outside of 0.01 and 0.99 they are cleansed
  smallgrid = smallgrid[which((smallgrid >= 0.01) & (smallgrid <= 0.99))]
  ptlen = length(smallgrid)
  
  ll = numeric(ptlen)
  wg_list = NULL
  lp_list = NULL
  
  f0u = dunif(u)
  
  # evaluating likelihood at grid points
  for (gi in 1:ptlen){
    pt = smallgrid[gi]
    ef = em_f1(u, rep(pt, n))
    
    print('got f1')
    
    lp = logistic_pi(f0u, ef$f1y, x)
    
    print('got pi')
    
    ll[gi] = mean(log((1 - lp$p) * f0u + lp$p * ef$f1y))
    
    wg_list[[gi]] = ef
    lp_list[[gi]] = lp
    
    if (verbose) print(paste0(gi , '/', ptlen))
  }
  
  best_index = which.max(ll)
  pt = smallgrid[best_index]
  ef = wg_list[[best_index]]
  lp = lp_list[[best_index]]
  
  lfdr = ((1 - lp$p) * f0u)/((1 - lp$p) * f0u + lp$p * ef$f1y)
  
  return(list(b = lp$b, p = lp$p, f1 = ef$f1, f1u = ef$f1y,
              localfdr = lfdr, pi0 = pt))
}

em_f1 = function(y, p, tol = 1e-6, maxit = 1e3){
  # y are observations coming from mixture of uniform and decreasing density
  # p are probability that u come from non-null (decreasing density)
  
  n = length(y)
  # we are going to work with sorted y's
  sy = sort(y, index.return = TRUE)
  y = sy$x
  p = p[sy$ix]
  
  # initial estimate
  d = diff(c(0, y))
  f1eval = pwg(rep(1/n, n), d)
  # pn here stores the prior probability of non-null
  pn = p
  
  flag = 1; cnt = 0; loglik = numeric(0);
  while(flag){
    
    # update pn
    # pn now stores the posterior probability of non-null
    pnnew = p*f1eval/((1-p) + p*f1eval)
    
    # update f1
    w = pnnew/sum(pnnew)
    f1evalnew = pwg(w, d)
    
    relerr = sum(abs(pn - pnnew))/sum(pn) + 
      sum(abs(f1evalnew - f1eval))/sum(f1eval)
    flag = min(flag, relerr > tol)
    flag = min(flag, cnt < maxit)
    cnt = cnt + 1
    
    pn = pnnew
    f1eval = f1evalnew
    loglik = append(loglik, mean(log((1-p) + p*f1eval)))
  }
  
  f1 = lc(y, f1eval)
  
  result = list(localfdr = 1 - pn[order(sy$ix)],
                f1sorted = f1eval, f1y = f1eval[order(sy$ix)],
                f1 = f1, loglik = loglik)
  return(result)
}



# now let's run
m1 = gridy_logistic_decreasing(yy, dummyx, intercept = FALSE, verbose = TRUE)







