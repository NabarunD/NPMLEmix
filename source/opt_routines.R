# this function is used by REBayes to construct the mids
histbin <- function(x, m = histm, eps = 1e-06, weights = weights) {
  u <- seq(min(x) - eps, max(x) + eps, length = m)
  xu <- findInterval(x, u)
  txu <- tabulate(xu)
  midu <- (u[-1] + u[-m])/2
  wnz <- (txu > 0)
  if (length(weights)) {
    if (length(weights) == length(x)) 
      w <- as.numeric(tapply(weights, xu, sum))
    else stop("length(weights) not equal length(x)")
  }
  else w <- txu[wnz]/sum(txu[wnz])
  list(midu = midu[wnz], w = w)
}

# let's try implementing primal formulation of kw in mosek
# this can be used to solve both the weighted and pi-constrained problems

kwprimal = function(y, pivec = rep(1,length(y)), weights = rep(1,length(y)), grid_len = 100){
  LIMIT_THREAD = TRUE
  m = grid_len
  n = length(y)
  
  # if (length(pivec) == 1) pivec = rep(pivec, n)
  # if (length(weights) == 1) pivec = rep(weights, n)
  
  atoms = seq(from = min(y), to = max(y), length.out = m)
  fmat = dnorm(outer(y, atoms, "-"))
  
  # setting up problem
  # want to maximize so sense = max
  # the variables which need to be optimized are (p_j, v_i)
  moseko = list(sense = "max")
  # the problem is stated by mosek in the form:
  # sum of non-linear functions separable in arguments + linear part (c^T\pi)
  # for us there is no linear part so c = 0
  moseko$c = c(rep(0, m + n),weights)
  # monotone constraints in pi
  A = rbind(cbind(fmat, -diag(1/pivec), matrix(0,n,n)), c(rep(c(1,0), times = c(m,2*n))))
  # sparsify to increase speed
  moseko$A = as(A, "CsparseMatrix")
  # moseko$A = Matrix(A, sparse = TRUE)
  # matrix of two rows
  # first row (blc) is lower bound on each linear relationship defined by A
  # second row (blu) is upper bound on each linear relationship defined by A
  moseko$bc = rbind(c((((pivec-1) * dnorm(y))/(pivec)), 1), c((((pivec-1) * dnorm(y))/(pivec)), 1))
  # box constraints: individual constraints on each pi_i
  # it seems to be easier to specify these separately
  # similar to above
  # first row (bxl) is lower bound on each pi_i
  # second row (bxu) is upper bound on each pi_i
  moseko$bx = rbind(c(rep(0, m + n), rep(-Inf, n)), c(rep(1, m), rep(Inf, 2*n)))
  # this is the interesting part, specifies the non-linear function in the objective
  # the function needs to separable in pi_i,
  # \sum_k f_k phi_k( g_k \pi_{j_k} + h_k)
  # where k is the k-th column of the following matrix
  # the number of columns is the number of non-linear functions the objective is separated into
  # the k-th column is a function of the j_k variable
  # the non-linearity is determined by 'type' (written as phi above)
  # you can multiply phi_k by some constant: f_k
  # also phi_k can be evaluated at a linear function of the variable in question: (g_k, h_k)
  #opro = matrix(list(), nrow = 5, ncol = n)
  #rownames(opro) = c("type", "j", "f", "g", "h")
  #opro[1, ] = rep("LOG", n)
  #opro[2, ] = m + 1:n
  #opro[3, ] = weights # coefficients outside log, in this case all 1
  #opro[4, ] = pivec # coefficient of v_i inside log
  #opro[5, ] = (1-pivec) * dnorm(y) #constants inside log
  #moseko$scopt = list(opro = opro)
  #if(LIMIT_THREAD) {moseko$iparam$num_threads <- 1}
  FE <- sparseMatrix(i=c(seq(1,3*n,3),seq(3,3*n,3)),j=c(((m+1):(m+n)),((m+n+1):(m+2*n))),x=rep(1,2*n),dims=c(3*n,m+2*n))
  gE <- rep(c(0, 1, 0), n)
  
  # Assemble input data
  moseko$F <- FE
  moseko$g <- gE
  moseko$cones <- cbind(matrix(list("PEXP", 3), nrow=2, ncol=n))
  rownames(moseko$cones) <- c("type","dim")
  
  ans <- mosek(moseko)
  
  probs = ans$sol$itr$xx[1:m]
  f1y = (ans$sol$itr$xx[(m + 1):n]-(1-pivec)*dnorm(y))/(pivec)
   
  return(list(atoms = atoms, probs = probs, f1y = f1y, ll = mean(log(f1y))))
}
# kwprimal using weights and histogram
kwprimal_weights = function(y, weights = 1, num_atoms = 100, hist_flag = TRUE, num_hist = num_atoms){
  
  n = length(y);
  if (length(weights) == 1) {weights = rep(weights, n)}
  # constructing histogram object
  if(hist_flag){
    histo = histbin(y, m = num_hist, weights = weights)
    # the compressed dataset is midu and w
    # print(names(histo))
    yy = histo$midu; ww = histo$w;
  } else {
    yy = y; ww = weights;
  }
  
  kwo = kwprimal(yy, weights = ww, grid_len = num_atoms);
  
  if (hist_flag){
    fmat = dnorm(outer(y, kwo$atoms, "-"))
    kwo$f1y = as.vector(fmat %*% kwo$probs)
    kwo$ll = mean(log(kwo$f1y))
  }
  
  return(kwo)
}

# solve logistic problem given alternate density (EM version)
lregem = function(weights, x, binit, lambda = 1e-2/length(weights)){
  # by default this provides a l2-regularization
  # equivalent to gaussian prior N(0,100) on all parameters
  # lambda = c(1e-2/length(f0y), rep(lambda, length(binit)-1))
  
  # defining objective function
  llobj = function(bb, weights, x, lambda){
    pivec = as.vector(1/(1 + exp(-x %*% bb)))
    - mean(weights * log(pivec) + (1-weights) * log(1-pivec)) + sum(lambda*bb^2)
  }
  
  # defining first gradient
  llfirst = function(bb, weights, x, lambda){
    pivec = as.vector(1/(1 + exp(-x %*% bb)))
    tt = pivec * (1-pivec) * (weights/pivec - (1-weights)/(1-pivec))
    - apply(x, 2, function(vec) mean(vec * tt)) + 2*lambda*bb
  }
  
  optimres = optim(par = binit, fn = llobj, gr = llfirst,
                   weights = weights, x = x, lambda = lambda,
                   method = 'BFGS')
  
  b = optimres$par
  p = as.vector(1/(1 + exp(-x %*% b)))
  
  return(list(b = b, p = p))
}

# EM
lgem = function(y, x,
                  weights, binit,
                  grid_len = 3*max(100, round(sqrt(length(y)))), histFlag = TRUE, timed = 60,
                  maxit = 600, tol = 1e-6, 
                  blambda = 1e-6/length(y)){
  
  st = proc.time()['elapsed']
  
  n = length(y); f0y = dnorm(y);
  
  # initial values, may need to revisit this
  lp = NULL; lp$b = binit; lp$p = as.vector(1/(1 + exp(-x %*% lp$b)));
  
  # starting iterations
  convFlag = 1; ll = NULL; lp_list = NULL; kw_list = NULL; time_list = NULL; itcount = 1; err_list = NULL;
  
  # if(verbose) pb = txtProgressBar(max = maxit+1, style = 3)
  while(convFlag){
    kwo = kwprimal_weights(y, weights = weights, num_atoms = grid_len, hist_flag = histFlag)
    lp = lregem(weights, x, lp$b, lambda = blambda)
    
    weights_new = lp$p*kwo$f1y/((1-lp$p)*f0y + lp$p*kwo$f1y)
    iter_error = rmse(weights_new, weights); err_list = append(err_list ,iter_error);
    
    lp_list = append(lp_list, list(lp)); kw_list = append(kw_list, list(kwo));
    ll = append(ll, mean(log((1-lp$p)*f0y + lp$p*kwo$f1y)))
    
    ct = proc.time()['elapsed'] - st
    time_list = append(time_list, ct)
    convFlag = (ct <= timed) & (itcount < maxit) & (iter_error > tol)
    weights = weights_new
    itcount = itcount + 1
    
  }
  
  localfdr = 1 - weights
  den = (1-lp$p)*f0y + lp$p*kwo$f1y
  
  return(list(atoms = kwo$atoms, probs = kwo$probs, f1y = kwo$f1y,
              b = lp$b, p = lp$p,
              f0y = f0y, den = den,
              localfdr = localfdr,
              ll = ll,  lp_list = lp_list, kw_list = kw_list,
              err_list = err_list, time_list = time_list,
              runtime = proc.time()['elapsed'] - st))
}


lla = function(a, pix, y, x, fmy){
  # given a scaling a, pi_hat pix, data y,x compute the solution at a * pix
  # this uses the approximate weighted scheme
  # fmy is required for this function
  # fmy = (1-pix) * f0y + pix * f1y
  # this can be estimated by running KWMLE on y
  f0y = dnorm(y)
  eps = 1e-10
  paraw = pmax(pmin(pix * a, 1-eps), eps)
  lmo = lm(mosaic::logit(paraw)~ 0 + x)
  ba = lmo$coefficients; pa = 1/(1 + exp(-lmo$fitted.values));
  weights = pmax(pmin(1 - (1-pa)*f0y/fmy, 0.99), 0.01)
  kwa = kwprimal_weights(y, weights = weights, num_atoms = 100)
  f1ya = kwa$f1y
  
  obj = list(b = ba, p = pa, atoms = kwa$atoms, probs = kwa$probs, f1y = f1ya, f0y = f0y)
  obj$den = with(obj, p * f1y + (1-p) * f0y)
  obj$ll = with(obj, mean(log(den)))
  obj$localfdr = with(obj, (1-p)*f0y/den)
  return(obj)
}

moveleft = function(pix, y, x, fmy, threshold_ll, atol = 1e-2){
  # this version uses the weighted lla for speed
  # function to move to the left from starting point until likelihood dips below threshold_ll
  # it is ensured that the final solution has higher likelihood than threshold
  ahi = 1; alo = 0;
  
  while((ahi - alo) > atol){
    acurrent = (alo + ahi)/2
    rr = lla(acurrent, pix, y, x, fmy); rr$alpha = acurrent;
    if(rr$ll < threshold_ll){alo = acurrent}
    else {ahi = acurrent}
  }
  rr = lla(ahi, pix, y, x, fmy); rr$alpha = ahi;
  return(rr)
}

moveright = function(pix, y, x, fmy, threshold_ll, atol = 1e-2){
  # this version uses the weighted lla for speed
  # function to move to the right from starting point until likelihood dips below threshold_ll
  # it is ensured that the final solution has higher likelihood than threshold
  ahi = 1/quantile(pix, 0.95); alo = 1;
  # in most cases, we are dealing with scenarios where 
  # the solution is at ahi
  # so we shall begin with a check there
  # if that check fails we shall do a standard bisection
  rrhi = lla(ahi, pix, y, x, fmy); rrhi$alpha = ahi;
  if(rrhi$ll >= threshold_ll){rr = rrhi} else {
    while((ahi - alo) > atol){
      acurrent = (alo + ahi)/2
      rr = lla(acurrent, pix, y, x, fmy); rr$alpha = acurrent;
      if(rr$ll < threshold_ll){ahi = acurrent}
      else {alo = acurrent}
    }
    rr = lla(alo, pix, y, x, fmy); rr$alpha = alo;
  }
  return(rr)
}

compute_emc = function(dd, em_time = 30){
  
  # computing EM, based on FDRreg
  # fdrreg solution
  # with splines
  f_ = with(dd, modf(FDRreg(y, xs, control = list(center = FALSE, scale = FALSE))))
  # we then compute the EM solution starting from PR
  # b should be padded with 2 zeros if splines are not used
  # should be padded with 6 zeros if splines are used
  pr_init_ = with(f_, list(b = c(mosaic::logit(1 - p0), rep(0, 6)),
                           w = (1 - p0) * f1y/((1 - p0) * f1y + p0 * f0y))) 
  # print('Predictive Recursion Initialization')
  # EM init by prfdr
  em_ = lgem(dd$y, cbind(1, dd$xs),
               weights = pr_init_$w, binit = pr_init_$b, grid_len = 300,
               timed = em_time, maxit = Inf, tol = 0)
  em_ = extract_index(em_, length(em_$ll))

  # parameters corresponding to EM solution
  sx_hat = function(xx){
    xs = spline_expand(xx)
    as.vector(cbind(1, xs) %*% em_$b)
  }
  tdparams_hat = clean_tdp(em_)
  
  # bootstrapped data
  dd_new = makedata(n, sx_hat, tdparams_hat)
  
  # conservative solution
  emc_new_ = with(dd_new, moveleft(pix, y, cbind(1, xs), den, ll))
  emc_ = extrapolate_dd(emc_new_, dd)

  return(list(f_ = f_, pr_init_ = pr_init_, em_ = em_, emc_ = emc_))
}



