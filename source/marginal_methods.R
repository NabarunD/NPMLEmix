m1 = function(y, x, pi0grid,
              mt = abs(y) - mean(abs(rnorm(1e4))), 
              grid_len = max(100, round(sqrt(length(y)))), blambda = 1e-6/length(y), verbose = FALSE){
  
  n = length(y); f0y = dnorm(y);
  kwm = kwprimal_weights(y, num_atoms = grid_len);
  fmy = kwm$f1y;
  
  m1step = function (pt) {
    weights = pmax(pmin(1 - (1-pt)*f0y/fmy, 0.99), 0.01)
    muinit = mean(mt)/pt
    binit = lm(mosaic::logit(pmax(pmin(mt/muinit, 0.9), 0.1)) ~ 0 + x)$coefficients;
    robj = lgstep(y, x, weights, binit, grid_len, blambda)
    robj$pi0 = pt
    return(robj)
  }
  
  if (verbose) {res = pblapply(pi0grid, m1step)}
  else {res = lapply(pi0grid, m1step)}
  
  ll_list = sapply(res, function (robj) robj$ll)
  bi = which.max(ll_list)
  robj = res[[bi]]
  robj$ll_list = ll_list
  
  return(robj)
}


lgstep = function(y, x, w, b, num_atoms, blambda){
  f0y = dnorm(y)
  # take one step of EM on data
  kwo = kwprimal_weights(y, weights = w, num_atoms)
  lp = lregoptim(f0y, kwo$f1y, x, b, lambda = blambda)
  # compute results
  den = (1 - lp$p) * f0y + lp$p * kwo$f1y; ll = mean(log(den)); localfdr = (1 - lp$p) * f0y/(den);
  return(list(p = lp$p, b = lp$b, f1y = kwo$f1y, kwo=kwo,
              localfdr = localfdr, den = den, ll = ll))
}

# solve logistic problem given alternate density (AM version)
lregoptim = function(f0y, f1y, x, binit, lambda = 1e-2/length(f0y)){
  # by default this provides a l2-regularization
  # equivalent to gaussian prior N(0,100) on all parameters
  # lambda = c(1e-2/length(f0y), rep(lambda, length(binit)-1))
  
  # defining objective function
  llobj = function(bb, f0y, f1y, x, lambda){
    pivec = as.vector(1/(1 + exp(-x %*% bb)))
    - mean(log(pivec*f1y + (1-pivec)*f0y)) + sum(lambda*bb^2)
  }
  
  # defining first gradient
  llfirst = function(bb, f0y, f1y, x, lambda){
    pivec = as.vector(1/(1 + exp(-x %*% bb)))
    wts = (f1y - f0y)/(pivec*f1y + (1-pivec)*f0y) * pivec * (1-pivec)
    - apply(x, 2, function(vec) mean(vec * wts)) + 2*lambda*bb
  }
  
  optimres = optim(par = binit, fn = llobj, gr = llfirst,
                   f0y = f0y, f1y = f1y, x = x, lambda = lambda,
                   method = 'BFGS')
  
  b = optimres$par
  p = as.vector(1/(1 + exp(-x %*% b)))
  
  return(list(b = b, p = p))
}



# marginal method 2, using a grid of potential values of pi0
# remember to use a good choice of mt (marginal transform)
m2 = function(y, x, pi0grid,
              mt = abs(y) - mean(abs(rnorm(1e4))),
              nlslambda = 1e-6/length(y), verbose = FALSE, solvef1 = FALSE){
  # for an user given sequence pi0 (overall non-null proportion)
  # solve non-linear least squares on marginal_transform to get initial value of beta
  
  n = length(y); f0y = dnorm(y); 
  mugrid = mean(mt)/pi0grid;
  
  computefn = function (mu) {
    # this function uses non-linear least squares to solve for beta for a fixed value of mu
    binit = lm(mosaic::logit(pmax(pmin(mt/mu, 0.9), 0.1)) ~ 0 + x)$coefficients
    nlso = m2boptim(mt/mu, x, binit, nlslambda)
    er = mean((mt - mu * nlso$p)^2) 
    return(list(nlso = nlso, er = er, mu = mu, pi0 = mean(mt)/mu))
  }
  
  if (verbose) {
    res = pblapply(mugrid, computefn)
  } else {
    res = lapply(mugrid, computefn)
  }
  
  ers = sapply(res, function (robj) robj$er)
  bi = which.min(ers)
  regres = res[[bi]]
  
  
  if(solvef1){
    robj = m1(y, x, pi0grid = pi0grid[bi], mt = mt)
  } else {
    robj = list(pi0grid = pi0grid, mugrid = mugrid,
                pi0 = pi0grid[bi], mu = mugrid[bi],
                regres = regres, res = res, ers = ers)
  }
  
  robj$pi0grid = pi0grid; robj$ers = ers;
  
  return(robj)
}

# solve marginal regression (logistic model) assuming mean under alternate is 1
m2boptim = function(y, x, binit, lambda = 1e-2/length(y)){
  # by default this provides a l2-regularization
  # equivalent to gaussian prior N(0,100) on all parameters
  
  # defining objective function
  lsobj = function(bb, y, x, lambda){
    pivec = as.vector(1/(1 + exp(-x %*% bb)))
    mean((y - pivec)^2) + lambda*sum(bb^2)
  }
  
  # defining first gradient
  lsfirst = function(bb, y, x, lambda){
    pivec = as.vector(1/(1 + exp(-x %*% bb)))
    wts = -2 * (y - pivec) * pivec * (1-pivec)
    apply(x, 2, function(vec) mean(vec * wts)) + 2*lambda*bb
  }
  
  optimres = optim(par = binit, fn = lsobj, gr = lsfirst,
                   y = y, x = x, lambda = lambda,
                   method = 'BFGS')
  
  b = optimres$par
  p = as.vector(1/(1 + exp(-x %*% b)))
  
  return(list(b = b, p = p))
}



