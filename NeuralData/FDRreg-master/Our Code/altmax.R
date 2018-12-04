###########################################################
# logistic solver when x has one column, no intercept term
###########################################################

lpi = function(f1y, f0y, x){
  # for one-dimensional x
  
  tol = 1e-6
  maxit = 1e2
  
  b = 0 # initialization
  
  conv_flag = 1; it_count = 0;
  while(conv_flag){
    exb = exp(x * b)
    
    t1 = f1y * exb + f0y
    t2 = exb + 1
    t3 = f1y * exb^2 - f0y
    
    first = mean((x * exb * (f1y - f0y))/(t1 * t2))
    second = mean((f1y - f0y) * x^2 * exb * t3/(t1^2 * t2^2))
    
    bnew = b + first/second
    
    it_count = it_count + 1
    conv_flag = ifelse((sum(abs(bnew - b)) <= tol) || (it_count == maxit), 0, 1)
    
    b = bnew
  }
  
  p = 1/(1 + exp(-x * b))
  return(list(b = b, p = p))
}


#############################################################################
# alternating maximization of likelihood
# when pi is modeled as logistic and f1 is modeled as a gaussian convolution
#############################################################################

# assume both initializations are given
# assume one-dimensional x

lgam = function(y, x, b, f1y, tol = 1e-5, maxit = 1e2, verbose = FALSE){
  
  ct = proc.time()
  n = length(y)

  f0y = dnorm(y)
  p = 1/(1 + exp(-x * b))

  conv_flag = 1
  it_count = 0
  if (verbose) print(it_count)
  loglik = mean(log((1-p) * f0y + p * f1y))
  
  while(conv_flag){
    
    #print(loglik)
    
    # alternating optimization
    kwo = kwmle_rcpp(y, p)
    f1y = kwo$f1u
    bnew = lpi(f1y, f0y, x)$b
    pnew = 1/(1 + exp(-x * bnew))
    
    #print(bnew)
    
    it_count = it_count + 1
    if (verbose) print(it_count)
    conv_flag = ifelse((mean(abs(pnew - p)) <= tol) || it_count == maxit, 0, 1)
    
    p = pnew
    loglik = append(loglik, mean(log((1-p) * f0y + p * f1y)))
  }
  
  return(list(atoms = kwo$atoms, probs = kwo$probs,
              b = bnew, p = p, loglik = loglik,
              runtime = proc.time() - ct))
}

##############################################
# one-dimensional version of gridy
##############################################

gridy_one = function(u, x){
  
  # intercept is added to x by default
  n = length(u)
  if(!is.matrix(x)) x = matrix(x, ncol = 1)
  
  # construct a small grid
  ptcutoff = cutoff_estimate(u, pnorm)$pt
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
  kw_list = NULL
  lp_list = NULL
  
  f0u = dnorm(u)
  
  # evaluating likelihood at grid points
  for (gi in 1:ptlen){
    pt = smallgrid[gi]
    
    kwo = kwmle_rcpp(u, pt)
    lp = lpi(kwo$f1u, f0u, x)
    
    ll[gi] = mean(log((1 - lp$p) * f0u + lp$p * kwo$f1u))
    
    kw_list[[gi]] = kwo
    lp_list[[gi]] = lp
  }
  
  best_index = which.max(ll)
  pt = smallgrid[best_index]
  kwo = kw_list[[best_index]]
  lp = lp_list[[best_index]]
  
  lfdr = (lp$p * f0u)/((1 - lp$p) * f0u + lp$p * kwo$f1u)
  
  return(list(b = lp$b, atoms = kwo$atoms, probs = kwo$probs,
              f1u = kwo$f1u, p = lp$p, localfdr = lfdr, pi0 = pt))
}


#############################################################################
# alternating maximization (general dimension)
#############################################################################

lgamd = function(y, x, b, tol = 1e-5, maxit = 1e2, verbose = FALSE){
  y = matrix(y, ncol = 1)
  ct = proc.time()
  n = length(y)
  if (is.vector(x)) x = matrix(x, ncol = 1)
  
  f0y = dnorm(y)
  p = 1/(1 + exp(-x %*% b))
  
  conv_flag = 1
  it_count = 0
  if (verbose) print(it_count)
  loglik = NULL
  
  while(conv_flag){
    # alternating optimization
    kwo = kwmle_rcpp(y, p)
    f1y = kwo$f1u
    bnew = logistic_pi(f0y, f1y, x)$b
    pnew = 1/(1 + exp(-x %*% bnew))
    
    #print(bnew)
    
    it_count = it_count + 1
    if (verbose) print(it_count)
    conv_flag = ifelse((mean(abs(pnew - p)) <= tol) || it_count == maxit, 0, 1)
    
    p = pnew
    loglik = append(loglik, mean(log((1-p) * f0y + p * f1y)))
  }
  
  return(list(atoms = kwo$atoms, probs = kwo$probs,
              b = bnew, p = p, loglik = loglik,
              runtime = proc.time() - ct))
}

# implementing second marginal method
blogist2 = function(y, x, tol = 1e-10, maxit = 1e3, bt = NULL, mu = NULL){
  if(is.null(mu)) mu = 1
  if(is.null(bt)) bt = unname(lm(y ~ 0 + x)$coefficients)/mu
  blist = bt
  mulist = mu
  objlist = mean((y - mu/(1 + exp(- x * bt)))^2)
  
  flag = 1
  itcount = 0 
  while(flag){
    pibx = 1/(1 + exp(-x * bt))
    pi2 = mean(pibx^2)
    ypi = mean(y * pibx)
    
    mulist = append(mulist, ypi/pi2)
    objlist = append(objlist, mean((y - tail(mulist,1) * pibx)^2))
    
    # want to maximize (ypi^2)/pi2
    # needed for first gradient
    dypi = mean(x * y * pibx * (1-pibx))
    dpi2 = 2 * mean(x * pibx^2 * (1-pibx))
    
    # first gradient
    first_num = 2 * pi2 * ypi * dypi - ypi^2 * dpi2
    first = first_num/pi2^2
    
    # needed for second gradient
    ddypi = mean(x^2 * y * pibx * (1 - pibx) * (1 - 2 * pibx))
    ddpi2 = 2 * mean(x^2 * pibx^2 * (1-pibx) * (2 - 3 * pibx))
    
    # second gradient
    second_num11 = pi2^2 * 2 * (dpi2 * ypi * dypi + pi2 * dypi^2 + pi2 * ypi * ddypi)
    second_num12 = pi2^2 * (ypi^2 * ddpi2 + 2 * ypi * dypi * dpi2)
    second_num1 = second_num11 - second_num12
    second_num2 = first_num * 2 * pi2 * dpi2
    second = (second_num1 - second_num2)/pi2^4
    
    update = -first/second
    # seems like we need a line search at this point
    # check overshooting etc
    # starting from truth seems to produce the right answer
    flagsearch = 1
    a = 1
    while(flagsearch){
      btnew = bt - a * first/second
      pibxnew = 1/(1 + exp(-x * btnew))
      pi2new = mean(pibxnew^2)
      ypinew = mean(y * pibxnew)
      munew = ypinew/pi2new
      
      objnew = mean((y - munew * pibxnew)^2)
      if(objnew <= tail(objlist, 1)) {flagsearch = 0} else {a = a/2}
      if(a < 1e-20) {flagserach = 0; flag = 0}
    }
    
    btnew = bt - a * first/second
    if (flag) flag = (!(abs(btnew - bt)/abs(bt) <= tol)) + (itcount == maxit)
    bt = btnew
    blist = append(blist, bt)
    itcount = itcount + 1
  }
  pibx = 1/(1 + exp(-x * bt))
  mulist = append(mulist, mean(y * pibx)/mean(pibx^2))
  objlist = append(objlist, mean((y - tail(mulist,1) * pibx)^2))
  return(list(bt = bt, mu = tail(mulist, 1), blist = blist, mulist = mulist, objlist = objlist))
}

# run this, this functions call the other two as required
bmlogist = function(y, x, tol = 1e-10, maxit = 1e3, nstart = 50, verbose = FALSE){
  reslist = list()
  objlist = NULL
  mulist = mean(y)/seq(from = 0, to = 1, length.out = nstart + 2)[1 + 1:nstart]
  for (i in 1:nstart){
    try({
      mu = mulist[i] #sign(mean(y)) * abs(rnorm(1))
      if(is.vector(x) || (ncol(x) == 1)){
        res = blogist2(y, x, tol, maxit, mu = mu)
      } else {
        res = blogistd(y, x, tol, maxit, mu = mu)
      }
      reslist = append(reslist, list(res))
      objlist = append(objlist, tail(res$objlist,1))
    }, silent = TRUE)
    if(verbose) print(i)
  }
  bestrun = which.min(objlist)
  return(reslist[[bestrun]])
}

blogistd = function(y, x, tol = 1e-10, maxit = 1e3, bt = NULL, mu = NULL){
  y = matrix(y, ncol = 1)
  if(is.null(mu)) mu = 1
  if(is.null(bt)) bt = unname(lm(y ~ 0 + x)$coefficients)/mu
  blist = bt
  mulist = mu
  objlist = mean((y - mu/(1 + exp(- x %*% bt)))^2)
  
  d = ncol(x)
  flag = 1
  itcount = 0 
  while(flag){
    pibx = 1/(1 + exp(-x %*% bt))
    pi2 = mean(pibx^2)
    ypi = mean(y * pibx)
    
    mulist = append(mulist, ypi/pi2)
    objlist = append(objlist, mean((y - tail(mulist,1) * pibx)^2))
    
    # want to maximize (ypi^2)/pi2
    # needed for first gradient
    dypi_temp = y * pibx * (1-pibx)
    dypi = apply(x, 2, function (vec) mean(vec * dypi_temp))
    dpi2_temp = 2 * pibx^2 * (1-pibx)
    dpi2 = apply(x, 2, function (vec) mean(vec * dpi2_temp))
    
    # first gradient
    first_num = 2*pi2*ypi*dypi - ypi^2*dpi2
    first = first_num/pi2^2
    
    # needed for second gradient
    ddypi_temp = as.vector(y * pibx * (1 - pibx) * (1 - 2 * pibx))
    ddypi_f = Vectorize(function(ii,jj) mean(x[,ii] * x[,jj] * ddypi_temp))
    ddypi = outer(1:d, 1:d, ddypi_f)
    
    ddpi2_temp = 2 * pibx^2 * (1-pibx) * (2 - 3 * pibx)
    ddpi2_f = Vectorize(function(ii,jj) mean(x[,ii] * x[,jj] * ddpi2_temp))
    ddpi2 = outer(1:d, 1:d, ddpi2_f)
    
    dypi_2 = outer(dypi, dypi)
    dpi2_2 = outer(dpi2, dpi2)
    dypi_dpi2 = outer(dypi, dpi2)
    
    second_num1 = pi2^2 * (2*ypi*pi2*ddypi + 2*ypi*dypi_dpi2 + 2*pi2*dypi_2 -
                             2*ypi*t(dypi_dpi2) - ypi^2*ddpi2)
    second_num2 = 2*pi2 * (2*ypi*pi2*dypi_dpi2 - ypi^2*dpi2_2)
    second = (second_num1 - second_num2)/pi2^4
    
    # second gradient
    second1 = (2*pi2*outer(dypi,dypi) + 2*ypi*pi2*ddypi + 2*ypi*outer(dypi,dpi2) - 2*ypi*outer(dpi2,dypi) - ypi^2*ddpi2)/pi2^2
    second2 = 2 * (2*ypi*pi2*outer(dpi2,dypi) - ypi^2*outer(dpi2,dpi2))/pi2^3
    second = second1 - second2
    
    update = - solve(second, first)
    # line search
    flagsearch = 1
    a = 1
    while(flagsearch){
      btnew = bt + a * update
      pibxnew = 1/(1 + exp(-x %*% btnew))
      pi2new = mean(pibxnew^2)
      ypinew = mean(y * pibxnew)
      munew = ypinew/pi2new
      
      objnew = mean((y - munew * pibxnew)^2)
      if(objnew <= tail(objlist, 1)) {flagsearch = 0} else {a = a/2}
      if(a < 1e-20) {flagserach = 0; flag = 0}
    }
    
    btnew = bt + a * update
    if (flag) flag = (!(mean(abs(btnew - bt))/mean(abs(bt)) <= tol)) + (itcount == maxit)
    bt = btnew
    blist = rbind(blist, bt)
    itcount = itcount + 1
  }
  pibx = 1/(1 + exp(-x %*% bt))
  mulist = append(mulist, mean(y * pibx)/mean(pibx^2))
  objlist = append(objlist, mean((y - tail(mulist,1) * pibx)^2))
  return(list(bt = bt, mu = tail(mulist, 1), blist = blist, mulist = mulist, objlist = objlist))
}



