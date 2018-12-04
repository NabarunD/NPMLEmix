require(Rcpp)
sourceCpp('C:/Users/Nabarun Deb/Dropbox/CovariateTesting/Real Data Analysis/NeuralData/FDRreg-master/Our Code/pava_weighted_grenander.cpp')

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

# lets write em for logistic-decreasing
ldam = function(y, x, b, tol = 1e-7, maxit = 1e2, verbose = FALSE){
  
  ct = proc.time()['elapsed']
  n = length(y)
  
  p = 1/(1 + exp(-x * b))
  
  conv_flag = 1
  it_count = 0
  loglik = NULL
  
  while(conv_flag){
    
    # alternating optimization
    ef = em_f1(y, p)
    f1y = ef$f1y
    bnew = lpi(f1y, rep(1,n), x)$b
    pnew = 1/(1 + exp(-x * bnew))
    
    it_count = it_count + 1
    if (verbose) print(it_count)
    conv_flag = ifelse((mean(abs(pnew - p)) <= tol) || it_count == maxit, 0, 1)
    
    p = pnew
    loglik = append(loglik, mean(log((1-p) + p * f1y)))
  }
  
  return(list(f1 = ef$f1, b = bnew, p = p, loglik = loglik,
              runtime = proc.time()['elapsed'] - ct))
}

# left continuous f1
lc = function(x, fx)
  Vectorize(approxfun(c(0,x), c(fx[1],fx), method = 'constant', yleft = 0, yright = 0, f = 1))

##########################################
# EM for f1 when pi is known
#########################################

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




# # pava to solve weighted grenander
# pava_weighted_grenander = function(w, d){
#   
#   n = length(d)
#   fxx = 1/n * w/d
#   fxx = fxx/(sum(fxx * d))
#   block_counts = rep(1, n)
#   
#   findblock = function(f){
#     return(head(which(diff(f) > 0), 1))
#   }
#   
#   merge_blocks = function(a, k){
#     a[k] = a[k] + a[k+1]
#     l = length(a)
#     if (l > k+1){
#       a[(k+1):(l-1)] = a[(k+2):(l)]
#     }
#     a = a[1:(l-1)]
#     return(a)
#   }
#   
#   update_block = function(a, k, v){
#     a[k] = v
#     l = length(a)
#     if (l > k+1){
#       a[(k+1):(l-1)] = a[(k+2):(l)]
#     }
#     a = a[1:(l-1)]
#     return(a)
#   }
#   
#   flag = 1
#   cnt = 0
#   while(flag){
#     b = findblock(fxx)
#     if(length(b)){
#       block_counts = merge_blocks(block_counts, b)
#       block_upper = sum(block_counts[1:b])
#       block_lower = block_upper - block_counts[b] + 1
#       block_sumw = sum(w[block_lower:block_upper])
#       block_sumd = sum(d[block_lower:block_upper])
#       
#       newfb = block_sumw/block_sumd
#       
#       fxx = update_block(fxx, b, newfb)
#     } else {
#       flag = 0
#     }
#     cnt = cnt + 1
#     if (cnt > 2e5){
#       flag = 0
#     }
#   }
#   
#   #print(cnt)
#   fxx = unlist(lapply(1:length(block_counts), function(ii) rep(fxx[ii], block_counts[ii])))
#   return(fxx)
# }
# 


