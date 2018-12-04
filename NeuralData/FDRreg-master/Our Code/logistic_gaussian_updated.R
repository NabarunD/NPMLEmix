######################################################
# to do isotonic regression
require(Iso) # use function pava from this
######################################################

# truncate a non-cdf object to the interval [0,1]
truncate_01 = function(xx) pmax(pmin(xx, 1), 0)

###########################################################
# calculate Cramer Von Mises statistic for given pi-tilde
###########################################################

cvm = function(u = NULL, F0 = NULL, xx = NULL, pt){
  # to the sample u, fit a cdf of the form (1-pt) * F0 + pt * F
  # using l2-error
  # F0 is assumed to be a known cdf, should be supplied as a function
  # the solution depends only on F0(sort(u))
  # this object can directly be passed to the function, which will increase efficiency
  # otherwise this is computed based on u and F0
  
  require(Iso)
  n = ifelse(is.null(xx), length(u), length(xx))
  yy = (1:n)/n
  if(is.null(xx)) xx = F0(sort(u))
  not_cdf = (yy - (1-pt)*xx)/pt
  d = truncate_01(not_cdf)
  f_hat = pava(d)
  dpt = sqrt(mean((not_cdf - f_hat)^2)) * pt
  return(list(cvm = dpt, F_hat = f_hat))
}

############################################################
# bisection implementation of cutoff estimate
############################################################

cutoff_estimate = function(u, F0, cutoff = NULL, pad = 1e-2, tol = 1e-5, verbose = FALSE){
  
  # u are observations
  # F0 is the known cdf of the null component
  # the default cutoff is 0.1 * log(log(n))/sqrt(n)
  # but this can be overridden to any use defined choice
  # to do Cramer von-Mises limiting distribution, use cutoff = 0.6972/sqrt(n)
  
  n = length(u); su = sort(u);
  if(is.null(cutoff)) cutoff = 0.1 * log(log(n))/sqrt(n);
  ptlo = 0 + pad; pthi = 1 - pad; flag = 1;
  yy = (1:n)/n; xx = F0(su)
  
  while(flag){
    pt = (ptlo + pthi)/2
    cvmo = cvm(xx = xx, pt = pt)
    dpt = cvmo$cvm
    
    if(dpt >= cutoff) ptlo = pt else pthi = pt
    flag = ifelse(pthi - ptlo <= tol, 0, 1)
    
    if(verbose) print(c(ptlo, pthi))
  }
  return(list(pt = pt, f_hat =cvmo$F_hat))
}

########################################################
# solve for pi using a logistic model given f0 and f1
########################################################

logistic_pi = function(f0u, f1u, x, intercept = FALSE, tol = 1e-6, maxit = 1e2){
  
  # Newton-Rhapson method
  # x is a n times p matrix
  # update = - (1/n * sum(w_i x_i x_i^T) )^inv (1/n * sum(v_i x_i))
  # where x_i are i-th row of x
  
  n = length(f0u)
  if(intercept) x = cbind(rep(1,length(n)), x)
  if(!is.matrix(x)) x = matrix(x, ncol = 1)
  
  # f1u may contain exact zero's when a grenander type estimate is used to f1
  if (min(f1u) == 0) f1u = pmax(f1u, 1e-6)
  c = f0u/f1u
  
  # try until
  error_flag = 1; has_failed_once = 0;
  while(error_flag){
    tryCatch({
      
      # if inverting runs into problems, perturb x a little bit
      if(has_failed_once) x = x + array(1e-4 * rnorm(prod(dim(x))), dim = dim(x))
      
      b = solve(crossprod(x), apply(x, 2, mean)) # initialization
      
      conv_flag = 1; it_count = 0;
      while(conv_flag){
        exb = as.vector(exp(- x %*% b))
        pp = 1/(1 + c * exb)
        qq = 1/(1 + exb)
        v = pp - qq
        w = pp * (1-pp) - qq * (1-qq)
        
        wx = apply(x, 2, function(zz) w * zz)
        second = crossprod(x, wx)
        first = apply(x, 2, function(zz) sum(v * zz))
        
        bnew = b - solve(second, first)
        
        it_count = it_count + 1
        conv_flag = ifelse( (sum(abs(bnew - b)) <= tol) || (it_count == maxit), 0, 1)
        
        b = bnew
      }
      
      
      error_flag = 0
    }, error = function(e) has_failed_once <<- 1)
  }
  
  p = as.vector(1/(1 + exp(-(x %*% b))))
  return(list(p = p, b = b))
}

############################################################################################
# approximate implementation of Kiefer-Wolfowitz MLE
# for one-dimensional data
# can solve mixture model given a vector of probability i-th obs came from standard normal
############################################################################################
require(inline)

code = '
// parsing into Rcpp objects
Rcpp::NumericVector rpf0(spf0);
Rcpp::NumericMatrix rpf1(spf1);
Rcpp::NumericVector rq(sq);
int n = rpf1.nrow(), grid_len = rpf1.ncol();
double tol = as<double>(stol);
int maxit = as<int>(smaxit);

// parsing into arma objects
arma::colvec pf0(rpf0.begin(), rpf0.size(), false);
arma::mat pf1(rpf1.begin(), n, grid_len, false);
arma::colvec q(rq.begin(), rq.size(), false);

// initializing more arma objects required during iterations
arma::colvec qnew(grid_len, arma::fill::zeros);
arma::colvec weights(grid_len, arma::fill::zeros);
arma::colvec m(n, arma::fill::zeros);

// setting up convergence
bool conv_flag = TRUE;
int it_count = 0;
while(conv_flag){
  m = pf0 + pf1 * q;
  weights = mean(pf1.each_col()/m).t();
  qnew = q % weights;
  qnew = qnew/sum(qnew);
  ++it_count;
  if ((mean(abs(qnew - q)) <= tol) || (it_count == maxit)) conv_flag = FALSE;
  q = qnew;
}

// let\'s get out of here
return Rcpp::wrap(q);
'

kwr = cxxfunction(signature(spf0 = 'numeric', spf1 = 'numeric',
                            sq = 'numeric', stol = 'numeric',
                            smaxit = 'numeric'), 
 code, plugin = 'RcppArmadillo')

kwmle_rcpp = function(x, pp = rep(1, length(x)),
                      grid_type = 'static',
                      grid_len = 10 * round(sqrt(length(x))),
                      init_atoms = NULL, init_probs = NULL,
                      tol = 1e-5/grid_len, maxit = 1e3){
    
  # approximate solution to the following problem:
  # given observations coming from (1-pp[i])dnorm + pp[i]f
  # where pp is a known vector, optimize likelihood over the unknown function f
  # which is modeled as an arbitrary mixture of unit-variance gaussians
  
  # the default value of tolerance should be 1e-3/length(atoms)

    ct = proc.time()[3]

    n = length(x)
    atoms = switch(grid_type,
                  static = quantile(x, seq(from = 1/grid_len, to = 1 - 1/grid_len, length.out = grid_len)),
                  exemplar = x,
                  specified = init_atoms)

    grid_len = length(atoms)
    stopifnot(grid_len > 0)
    if((grid_type == 'specified') & !is.null(init_probs)) qinit = init_probs else qinit = rep(1/length(atoms), length(atoms))

    fmat = sapply(atoms, function(tt) dnorm(x - tt))
    f0vec = dnorm(x)

    pf0 = (1-pp) * f0vec
    pf1 = sapply(1:length(atoms), function(j) pp * fmat[, j])

    qconverged = kwr(pf0, pf1, qinit, tol, maxit)

    f1u = as.vector(fmat %*% qconverged)

    rt = proc.time()[3] - ct 
    return(list(atoms = atoms, probs = as.vector(qconverged), f1u = f1u, runtime = rt))
}

# not to be used, kept for systems where rcpp is not installed

kwmle = function(x, pp = rep(1, length(x)),
				 atoms = NULL, probs = NULL,
                 grid_type = 'static', grid_len = 10 * round(sqrt(length(x))), 
                 tol = 1e-3/length(atoms), maxit = 1e3){

  # approximate solution to the following problem:
  # given observations coming from (1-pp[i])dnorm + pp[i]f
  # where pp is a known vector, optimize likelihood over the unknown function f
  # which is modeled as an arbitrary mixture of unit-variance gaussians
  
  # the default value of tolerance should be 1e-3/length(atoms)

  # setting up a grid of atoms
  tgrid =switch(grid_type,
                static = quantile(x, seq(from = 1/grid_len, to = 1 - 1/grid_len, length.out = grid_len)),
                exemplar = x,
                specified = atoms)

  if((grid_type == 'specified') & !is.null(probs)) q = probs else q = rep(1/length(tgrid), length(tgrid))
  
  # iterations
  # q's determine mixture sizes in f1
  # pp is a vector, every observations comes from dnorm with different probs - given by 1-pp
  flag = 1
  fmat = sapply(tgrid, function(tt) dnorm(x - tt))
  f0vec = dnorm(x)
  it_count = 0
  while(flag){
    # updating
    marginal_vec = (1-pp) * f0vec + pp * as.vector(fmat %*% q)
    weight = sapply(1:length(tgrid), function(j) mean(pp * fmat[, j]/marginal_vec))
    qnew = q * weight
    qnew = qnew/sum(qnew)
    
    # checking for covergence
    it_count = it_count + 1
    flag = ifelse((mean(abs(qnew - q)) <= tol) || (it_count == maxit), 0, 1)
        
    # accepting new solution
    q = qnew
  }
  
  # the estimate f1 can always be written in terms of atoms and probs
  # it is also efficient to return this estimate evaluated at all sample points
  # since it is often used for other methods we are working in conjunction with  
  f1u = as.vector(fmat %*% q)
  
  return(list(atoms = tgrid, probs = q, f1u = f1u))
}

############################################################
# utility to make a density given atoms, probs and vars
############################################################

make_density = function(atoms, probs,
                          vars = rep(1, length(atoms))){
  
  f = function(xx) sum(sapply(1:length(atoms), function(j) 
    probs[j] * dnorm(xx, atoms[j], sqrt(vars[j]))))
  return(Vectorize(f))
}

########################################################################
# function to produce a suitable grid to evaluate f1 and its estimates
########################################################################

make_grid = function(tdp, sample_size = 1e4, len = 3000){
  n = sample_size;
  seq(from = min(tdp[, 2] - sqrt(1 + tdp[, 3]) * sqrt(2 * log(n))),
      to = max(tdp[, 2] + sqrt(1 + tdp[, 3]) * sqrt(2 * log(n))),
      length.out = len)
}

##########################################################################
# a solution of the joint optimization problem based on
# a grid expansion around the cutoff estimate
# a) use the cutoff estimate coded above
# b) construct a small grid around this
# c) take a step starting from each grid-point and compute full likelihood
# d) choose the solution with highest likelihood
##########################################################################

require(splines)
gridy = function(u, x, intercept = TRUE, verbose = FALSE,
                 kwmle_maxit = 1e3){
  
  # intercept is added to x by default
  n = length(u)
  if(!is.matrix(x)) x = matrix(x, ncol = 1)
  #x = matrix(ns(x, df = df), ncol = ncol(x) * df)
  if(intercept) x = cbind(rep(1, n), x)
  
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
    
    kwo = kwmle_rcpp(u, pt, maxit = kwmle_maxit)
    lp = logistic_pi(f0u, kwo$f1u, x)
    
    ll[gi] = mean(log((1 - lp$p) * f0u + lp$p * kwo$f1u))
    
    kw_list[[gi]] = kwo
    lp_list[[gi]] = lp
    
    if (verbose) print(paste0(gi , '/', ptlen))
  }
  
  best_index = which.max(ll)
  pt = smallgrid[best_index]
  kwo = kw_list[[best_index]]
  lp = lp_list[[best_index]]
  
  lfdr = ((1 - lp$p) * f0u)/((1 - lp$p) * f0u + lp$p * kwo$f1u)
  
  return(list(b = lp$b, atoms = kwo$atoms, probs = kwo$probs,
         f1u = kwo$f1u, p = lp$p, localfdr = lfdr, pi0 = pt))
}

######################################################################
# a faster method based on debiasing of the cutoff estimate
#
# a) get pi-cutoff from data and then f1-hat corresponding to this
# b) resample from the two groups mixture (pi-cutoff, f1-hat ; f0 is always known) a few times
# c) compute pi-cutoff-star for these samples and compute bias: mean(pi-cutoff-star) - pi-cutoff
# d) remove bias from original pi-cutoff
# e) take a step from the bias-corrected pi-cutoff

######################################################################

debiased = function(u, x, df = 1, B = 10, intercept = TRUE){
  
  # intercept is added to x by default
  n = length(u)
  if(!is.matrix(x)) x = matrix(x, ncol = 1)
  x = matrix(ns(x, df = df), ncol = ncol(x) * df)
  if(intercept) x = cbind(rep(1, n), x)
  
  # initial estimate of pi0
  pt = cutoff_estimate(u, pnorm)$pt
  kwo = kwmle_rcpp(u, pt)

  # then we draw samples from such a kwo
  ptstar = numeric(B)
  for (i in 1:B){
    ustar = rnorm(n)
    nnstar = runif(n) >= pt
    ustar[nnstar] = ustar[nnstar] + sample(kwo$atoms, sum(nnstar), prob = kwo$probs, replace = TRUE)

    ptstar[i] = cutoff_estimate(ustar, pnorm)$pt
  }

  bias_estimate = mean(ptstar) - pt
  pt = pt - bias_estimate

  # final kwmle_rcpp
  kwo = kwmle_rcpp(u, pt)
  # final logistic
  lp = logistic_pi(dnorm(u), kwo$f1u, x)

  # compute local fdr corresponding to this  
  localfdr = (lp$p * dnorm(u))/((1 - lp$p) * dnorm(u) + lp$p * kwo$f1u)
  
  return(list(b = lp$b, atoms = kwo$atoms, probs = kwo$probs,
         f1u = kwo$f1u, p = lp$p, localfdr = localfdr, pi0 = pt))
}

##########################################################################################
# a cutoff-less cross-validated method using bisection
# with flipped folds, take ss to be only one fold but holdout to be all other folds
##########################################################################################

heldout_gridy = function(u, x, df = 1, intercept = TRUE){

  # intercept is added to x by default
  n = length(u)
  if(!is.matrix(x)) x = matrix(x, ncol = 1)
  x = matrix(ns(x, df = df), ncol = ncol(x) * df)
  if(intercept) x = cbind(rep(1, n), x)

  ptlo = 1/sqrt(n); pthi = 1 - 1/sqrt(n); tol = 1/sqrt(n);
  conv_flag = 1;

  # construct some folds
  k = 5
  fold_id = rep(1:k, each = ceiling(n/5))[sample(n)]

  while(conv_flag){

  	pt = (ptlo + pthi)/2
    ptleft = pt - tol/2
    ptright = pt + tol/2

    heldout_ll_left = 0; heldout_ll_right = 0;

    for(i in 1:k){
    	ss = which(fold_id == i)
    	holdout = setdiff(1:n, ss)

	    kwo_left = kwmle_rcpp(u[ss], ptleft)
	    lp_left = logistic_pi(dnorm(u[ss]), kwo_left$f1u, x[ss, ])

	    fh_left = as.vector(sapply(kwo_left$atoms, function(tt) dnorm(u[holdout] - tt)) %*% kwo_left$probs)
	    ph_left = as.vector(1/(1 + exp(x[holdout, ] %*% lp_left$b)))
  	  	heldout_ll_left = heldout_ll_left + mean(log((1 - ph_left) * dnorm(u[holdout]) + (ph_left) * fh_left))

	    kwo_right = kwmle_rcpp(u[ss], ptright)
	    lp_right = logistic_pi(dnorm(u[ss]), kwo_right$f1u, x[ss, ])

	    fh_right = as.vector(sapply(kwo_right$atoms, function(tt) dnorm(u[holdout] - tt)) %*% kwo_right$probs)
	    ph_right = as.vector(1/(1 + exp(x[holdout, ] %*% lp_right$b)))
  	  	heldout_ll_right = heldout_ll_right + mean(log((1 - ph_right) * dnorm(u[holdout]) + (ph_right) * fh_right))
    }

    if(heldout_ll_left < heldout_ll_right) ptlo = pt else pthi = pt
    if(pthi - ptlo <= tol) conv_flag = 0
  }

  # final kwmle_rcpp
  kwo = kwmle_rcpp(u, pt)
  # final logistic
  lp = logistic_pi(dnorm(u), kwo$f1u, x)

  # compute local fdr corresponding to this  
  localfdr = (lp$p * dnorm(u))/((1 - lp$p) * dnorm(u) + (lp$p) * kwo$f1u)
  
  return(list(b = lp$b, atoms = kwo$atoms, probs = kwo$probs,
         f1u = kwo$f1u, p = lp$p, localfdr = localfdr, pi0 = pt))

}


