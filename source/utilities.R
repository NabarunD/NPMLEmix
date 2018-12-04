fdp_list = seq(from = 0.05, to = 0.30, by = 0.05)
cols = c('orange', 'red', 'blue', 'green4', 'cyan3','violetred1', 'slateblue1')

##########################################################
# expanding each column of a matrix via splines
##########################################################

spline_expand = function(xx, df = 3){
  require(splines)
  b1 = bs(xx[,1], df = df)
  b2 = bs(xx[,2], df = df)
  model.matrix( ~  b1 + b2 - 1)
}

##############################################################
# main function for generating data
# from logistic pix and gaussian mixture f1
# sx denotes the link within logistic
# tdparams determines the mixing density underlying f1
# covariates are drawn from Unif([0,1]^2)
# third order spline expansion of covariates is also returned
##############################################################

makedata = function(n, sx, tdparams){
  
  # constructing f1
  atoms = tdparams$atoms; probs = tdparams$probs; variances = tdparams$variances;
  # checking validity of tdparams
  stopifnot(length(unique(c(length(atoms), length(probs), length(variances)))) == 1)
  # number of components
  nc = length(atoms)
  # f1 function
  f1 = make_density(tdparams)
  
  # here we generate data
  ####################################################################
  # starting data generation process
  # covariate matrix
  x = cbind(runif(n), runif(n))
  # spline expansion
  xs = spline_expand(x)
  
  # generating pi(.)
  pix = 1/(1 + exp(-sx(x)))
  
  # generating observations
  nulls = (runif(n) >= pix);  nind = which(nulls); nnind = setdiff(1:n, nind);
  # populating data
  theta = array(0, dim = n)
  wc = sample(nc, length(nnind), replace = TRUE, prob = probs)
  theta[nnind] = atoms[wc] + rnorm(length(nnind)) * sqrt(variances[wc])
  y = rnorm(n) + theta
  # completed data generation process
  #####################################################################
  
  f0y = dnorm(y)
  f1y = f1(y)
  den = (1-pix)*f0y + pix*f1y
  localfdr = (1-pix)*f0y/den
  ll = mean(log(den))
  
  return(list(y = y, x = x, xs = xs,
              sx = sx, pix = pix,
              f0y = f0y, f1y = f1y, localfdr = localfdr, ll = ll, f1 = f1,
              nnind = nnind, den = den,
              sx = sx, tdparams = tdparams))
}

##########################################################
# construct (approximate) distribution function
# from density f
##########################################################

makecdf = function(f){
  xx = seq(from = 0, to = 1, length.out = 1e5)
  fx = f(xx)
  Fx = cumsum(fx)/length(xx)
  Fxx = approxfun(xx, Fx, yleft = 0, yright = 1)
  return(Fxx)
}

##########################################################
# functions to compute the True Positive Proportion
# and False Discovery Proportion 
# between predictions (preds) and truth (truth)
##########################################################

calculateTPR <- function(preds, truth) {
  # Calculate true positive rate for given model predictions at specified
  # threshold
  # Args:
  #  preds: the predicted set s from knockoff
  #  truth: true set s   
  if (is.null(preds) | length(preds) == 0) return(0)  
  as.numeric(sum(truth %in% preds) / length(truth))
}

calculateFDR <- function(preds, truth) {
  # Calculate false discovery rate for given model predictions at specified
  # threshold
  # Args:
  #  preds: the predicted set s from knockoff
  #  truth: true set s     
  if (is.null(preds) | length(preds) == 0) return(0)
  as.numeric(sum(!preds %in% truth) / length(preds))
}

##########################################################
# functions to compute the conditional TPR and FDR
# based on true localfdr's and preds (set of rejected hypotheses)
##########################################################

calculateCTPR <- function(preds, lfdr) {
  if (is.null(preds) | length(preds) == 0) return(0)  
  sum(1-lfdr[preds])/sum(1-lfdr)
}

calculateCFDR <- function(preds, lfdr) {
  if (is.null(preds) | length(preds) == 0) return(0)
  mean(lfdr[preds])
}

############################################################
# utility to make a density given atoms, probs and vars
############################################################

make_density = function(tdp){
  f <- function(xx) sum(sapply(1:length(tdp$atoms), function(j) 
    tdp$probs[j] * dnorm(xx, tdp$atoms[j], sqrt(1 + tdp$variances[j]))))
  return(Vectorize(f))
}

########################################################################
# function to produce a suitable grid to evaluate f1 and its estimates
########################################################################

make_grid = function(tdp, sample_size = 1e4, len = 3000){
  n = sample_size; inflate = max(sqrt(1+tdp$variances) * sqrt(2*log(n)));
  seq(from = min(tdp$atoms) - inflate, to = max(tdp$atoms) + inflate, length.out = len)
}

########################################################################
# functions to calculate MSE and deciles
########################################################################

rmse = function(vec, truth) sqrt(mean((vec - truth)^2))
deciles = function(vec) quantile(vec, seq(from = 0.1, to = 0.9, by = 0.1))

####################################
# Benjamini-Hochberg
####################################

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

####################################
# lfdr based rejection
####################################

lfdr_reject = function(lfdr, fdr_nominal = 0.1){
  sl = sort(lfdr)
  k = tail(which(cumsum(sl)/seq_along(sl) <= fdr_nominal), 1)
  if(length(k)){if (k) rejects = which(lfdr <= sl[k])} else {rejects = numeric(0)}
}

# from a full mle solution compute the mse of estimated lfdr
# compared to true lfdr at each iteration

rmsepath = function(obj, dd){
  rmses = numeric(length(obj$ll))
  with(obj , {
    for (i in 1:length(rmses)){
      p = obj$lp_list[[i]]$p
      f1 = obj$kw_list[[i]]$f1y
      lfdr = 1 - (p * f1)/((1-p) * f0y + p * f1)
      rmses[i] <<- rmse(dd$localfdr, lfdr)
    }
  })
  return(rmses)
}


# from a full mle solution with lots of iterations, extract an earlier iteration for comparison
extract_index = function(obj, ii){
  with(obj, {
    rr <<- list(atoms = kw_list[[ii]]$atoms, probs = kw_list[[ii]]$probs,
                f1y = kw_list[[ii]]$f1y, f0y = f0y,
                b = lp_list[[ii]]$b, p = lp_list[[ii]]$p,
                ll = ll[ii],
                runtime = time_list[ii])
  })
  rr$den = with(rr, (1-p)*f0y + p * f1y)
  rr$localfdr = with(rr, (1-p)*f0y/((1-p)*f0y + p * f1y))
  return(rr)
}

add_alpha = function(cols, alpha = 0.7) 
  rgb(red = t(col2rgb(cols))/255, alpha = alpha)

rmlast = function(ss) substr(ss, 1, nchar(ss) - 1)

# find good indices
fgi = function(vec, alpha = 0.5) which(vec >= max(vec) - alpha * sd(vec))

# find index upto which vec increases
fc = function(vec, ivec = 1:length(vec)){
  first_decrease = which(sign(c(diff(vec), -1)) == -1)[1]
  ivec[first_decrease]
}

# relocate vec so that minimum is zero
rmmin = function(vec) {vec - min(vec)}

# mod fdrreg
modf = function(obj){
  obj$b = obj$model$coef
  obj$p = obj$priorprob
  obj$f0y = obj$M0
  obj$f1y = obj$M1
  obj$den = (1-obj$p) * obj$f0y + obj$p * obj$f1y
  obj$ll = mean(log(obj$den))
  return(obj)
}

# extract initialization
extract_init = function(obj){
  return(list(p = obj$p, b = obj$b, w = 1 - obj$localfdr))
}

# construct sumamries to be used with errbar plots
require(Hmisc)
geteo = function(fmat){
  x = fdp_list
  y = apply(fmat, 2, median)
  yplus = apply(fmat, 2, quantile, 3/4)
  yminus = apply(fmat, 2, quantile, 1/4)
  return(list(x = x, y = y, yplus = yplus, yminus = yminus))
}

# list of colors
cl = c('orangered', 'dodgerblue', 'darkorchid', 'green4', 'violetred1', 'red', 'blue')

extrapolate_dd = function(res, dd){
  
  # res is a solution from some optim method
  # containing b, atoms and probs
  # need to calculate p, f1y, ... etc on the dataset dd
  # based on these estimates
  
  out = list(b = res$b, atoms = res$atoms, probs = res$probs)
  
  if (1 + ncol(dd$x) == length(res$b)){
    out$p = as.vector(1/(1 + exp(- cbind(1, dd$x) %*% res$b)))
  } else {
    out$p = as.vector(1/(1 + exp(- cbind(1, dd$xs) %*% res$b)))
  }
  
  fmat = dnorm(outer(dd$y, res$atoms, "-"))
  out$f1y = as.vector(fmat %*% res$probs)
  out$f0y = dnorm(dd$y)
  
  out$den = with(out, (1-p)*f0y + p*f1y)
  out$ll = mean(log(out$den))
  out$localfdr = with(out, (1-p)*f0y/den)
  
  return(out)
}


# clean tdparams to remove very small probabilities

clean_tdp = function(tdp, eps = 1e-6){
  bi = which(tdp$probs >= eps); out = NULL;
  out$atoms = tdp$atoms[bi]; out$probs = tdp$probs[bi]/sum(tdp$probs[bi]);
  if (length(tdp$variances)) {out$variances = tdp$variances[bi]} else {out$variances = rep(0, length(bi))}
  return(out)
}




# need to compute squared Hellinger based on a new sample
compute_f1 = function(yy, kwo){
  tdp = list(atoms = kwo$atoms, probs = kwo$probs, variances = rep(0, length(kwo$atoms)))
  tdp = clean_tdp(tdp)
  f1 = make_density(tdp)
  f1(yy)
}

compute_hel = function(dd, kwo){
  ygrid = make_grid(dd$tdparams, length(dd$y))
  dd_f1yy = dd$f1(ygrid)
  kwo_f1yy = compute_f1(ygrid, kwo)
  return(mean((sqrt(dd_f1yy) - sqrt(kwo_f1yy))^2))
}

