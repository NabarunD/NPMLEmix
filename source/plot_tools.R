do_everything = function(file_list,
                         out_prefix,
                         splines_flag = FALSE,
                         method_ids = c('ce_', 'f_', 'emc_'),
                         method_names = c('conservative', 'FDRreg', 'EM (conservative)')){
  
  aaa = pbsapply(file_list, function (fname)
  {le = new.env(); load(fname, envir = le); as.list(le);})
  
  # want to plot the estimates of pi
  # first of all, let us draw some new x's
  m = 5e2; x_ = cbind(runif(m), runif(m));
  sxx_ = aaa[[1]]$dd$sx(x_); osxx_ = order(sxx_);
  pix_formal = as.vector(1/(1 + exp(-sxx_))) # the formal pi function at fresh points
  # spline expansion of the new x's
  xs_ = spline_expand(x_)
  pix_formal_obj = list(sxx_ = sxx_,
                        osxx_ = osxx_,
                        pix_formal = pix_formal)
  
  # estimates of pi evaluated at x_
  if (splines_flag){
    pix_est = lapply(aaa, function (ro) sapply(method_ids, function (mi)
      1/(1 + exp(-cbind(1, xs_) %*% get(mi, ro)$b))))
  } else {
    pix_est = lapply(aaa, function (ro) sapply(method_ids, function (mi)
      1/(1 + exp(-cbind(1, x_) %*% get(mi, ro)$b))))
  }
  
  pix_est = array(unlist(pix_est),
                  dim = c(nrow(pix_est[[1]]), ncol(pix_est[[1]]), length(pix_est)))
  dimnames(pix_est)[[2]] = method_names
  
  tdparams = aaa[[1]]$dd$tdparams
  f1 = make_density(tdparams)
  y_ = seq(from = -10, to = 10, by = 0.05)
  f1y_formal = f1(y_)
  f1y_formal_obj = list(y_ = y_,
                        f1y_formal = f1y_formal)
  
  f1y_est = lapply(method_ids, function (mn) extract_f1(aaa, mn, y_))
  f1y_est = array(unlist(f1y_est),
                  dim = c(nrow(f1y_est[[1]]), ncol(f1y_est[[1]]), length(f1y_est)))
  dimnames(f1y_est)[[3]] = method_names
  
  # list of colors
  cl = c('orange', 'dodgerblue', 'darkorchid', 'red', 'blue', 'violetred1')
  
  # removing folder specification from out_prefix
  plot_title_prefix = substring(out_prefix, tail(gregexpr("\\/", out_prefix)[[1]], 1) + 1)
  
  pdf(paste0(out_prefix, '_estimates.pdf'), width = 15, height = 7)
  
  par(mfrow = c(1,2))
  plot_pis(pix_est, pix_formal_obj, cl)
  title(paste0(plot_title_prefix, '_pix'))
  
  plot_f1s(f1y_est, f1y_formal_obj, cl)
  title(paste0(plot_title_prefix, '_f1y'))
  dev.off()
  
}

# this is an utility function for linearly interpolating an f1y type output to fresh data points
extract_f1 = function(aaa, mn, ynew){
  sapply(1:length(aaa), function (ii) {
    f1y = with(aaa[[ii]], get(mn))$f1y
    yy = aaa[[ii]]$dd$y
    approx(yy, f1y, ynew, rule = 2)$y
  })
}

plot_pis = function(pmat, po, cl){
  with(po, {
    par(bty = 'l')
    plot(sxx_[osxx_], pix_formal[osxx_], type = 'l', ylim = c(0,1),
         ylab = 'formal pi(x)', xlab = 's(x)', lwd = 2)
    for (j in 1:dim(pmat)[2]){
      for (i in 1:dim(pmat)[3]){
        lines(sxx_[osxx_], pmat[osxx_,j,i],  col = add_alpha(cl[j], 0.1))
      }
    }
    lines(sxx_[osxx_], pix_formal[osxx_], lwd = 2)
    
    legend(x = 'topleft',
           legend = unlist(c('formal', dimnames(pmat)[2])),
           lwd = 4, col = c('black', cl), bty = 'n')
  })
}

plot_f1s = function(fmat, fo, cl){
  with(fo, {
    par(bty = 'l')
    plot(y_, f1y_formal,
         type ='l', ylim = c(0, 1.1 * max(fmat, f1y_formal)),
         xlab = 'y', ylab = 'density')
    for (j in 1:dim(fmat)[3]){
      for (i in 1:dim(fmat)[2]){
        lines(y_, fmat[,i,j], cex = 0.3, col = add_alpha(cl[j], 0.2))
      }
    }
    points(y_, f1y_formal, type ='l', lwd = 2)
    legend(x = 'topleft',
           legend = unlist(c('formal', dimnames(fmat)[3])),
           lwd = 4, col = c('black', cl), bty = 'n')
  })
}


process_result = function(robj, dd,
                          metrics = c('p', 'localfdr', 'f1y', 'lfdr_rej')){
  # this can also handle errors in terms of denominator (den)
  # and weighted BH rejection (wtbh_rej)
  # this are not returned by default
  # add these options to metrics if you wish
  
  ro = NULL
  
  if('p' %in% metrics) ro$pierr = rmse(robj$p, dd$pix)
  if('localfdr' %in% metrics) ro$lfdrerr = rmse(robj$localfdr, dd$localfdr)
  if('f1y' %in% metrics) ro$f1err = rmse(robj$f1y, dd$f1y)
  if('den' %in% metrics) ro$denerr = rmse(robj$den, dd$den)
  
  if('lfdr_rej' %in% metrics) {
    # computing rejections based on lfdr
    lrej = sapply(fdp_list, function (ff) lfdr_reject(robj$localfdr, ff))
    # computing fdp and tpp at each nominal level
    ro$lfdp = sapply(lrej, function(rr) calculateFDR(rr, dd$nnind))
    ro$ltpp = sapply(lrej, function(rr) calculateTPR(rr, dd$nnind))
  }
  if('wtbh_rej' %in% metrics) {
    # computing rejections based on weighted BH
    pvals = 2 * (1 - pnorm(abs(dd$y)))
    wrej = sapply(fdp_list, function (ff) bh((1-robj$p)*pvals, ff))
    # computing fdp and tpp at each nominal level
    ro$wfdp = sapply(wrej, function(rr) calculateFDR(rr, dd$nnind))
    ro$wtpp = sapply(wrej, function(rr) calculateTPR(rr, dd$nnind))
  }
  
  return(ro)
}

# analyze all methods in a simulation run
pr = function(rri, ignore_strings = c('dd', 'gs')){
  # this object contains original data in dd
  # alongwith results from several methods
  # extracting method names
  mn = setdiff(names(rri), ignore_strings)
  ro = sapply(mn, function (mni) process_result(get(mni, rri), rri$dd))
  # colnames(ro) = mn
  colnames(ro) = sapply(mn, rmlast)
  return(ro)
}

# this is for averaging FDP to get FDR (but can be used for other tasks following the same pattern)
extract_avg = function(pmat, aggfunc = mean){
  # pmat is num_methods times niter, but each entry is a list of size length(fdp_list)
  # out = apply(pmat, 1, function(prow) colMeans(t(sapply(prow, unlist))))
  
  # there seems to be a strange bug somewhere
  tempFunc = function(prow){
    bad_indices = which(sapply(prow, length) != length(fdp_list))
    if(length(bad_indices)){
      for (bi in bad_indices){
        prow[[bi]] = rep(0, length(fdp_list))
      }
    }
    apply(t(sapply(prow, unlist)), 2, aggfunc)
  }
  
  out = apply(pmat, 1, tempFunc)
  # out = apply(pmat, 1, function(prow) apply(t(sapply(prow, unlist)), 2, aggfunc))
  return(out)
} 

# plots performance curves of all methods
# makes two separate plots par(mfrow = ...) should be set manually
perfcurves = function(pmat, cols, fdp_list, aggfunc = mean){
  
  # pmat is 2 times num_method times niter
  fdrmat = extract_avg(pmat[1,,], aggfunc)
  tprmat = extract_avg(pmat[2,,], aggfunc)
  
  # FDR plot
  plot(fdp_list, fdrmat[,1], lwd = 4, type = 'l', col = cols[1], ylim = c(0,0.6),
       main = 'False Discovery Rate', ylab = '', xlab = 'Nominal level')
  for(i in 2:ncol(fdrmat)) lines(fdp_list, fdrmat[,i], lwd = 4,  col = cols[i])
  lines(fdp_list, fdp_list, lty = 2, col = 'grey', lwd = 3)
  legend(x = 'topleft', legend = colnames(fdrmat), col = cols, lwd = 4, bty = 'n')
  
  # TPR plot
  plot(fdp_list, tprmat[,1], lwd = 4, type = 'l', col = cols[1], ylim = c(0,1),
       main = 'True Positive Rate', ylab = '', xlab = 'Nominal level')
  for(i in 2:ncol(tprmat)) lines(fdp_list, tprmat[,i], lwd = 4,  col = cols[i])
  legend(x = 'topleft', legend = colnames(tprmat), col = cols, lwd = 4, bty = 'n')
  
}


plot_aggregate = function(amat, aggfunc = mean){
  
  a = matrix(as.numeric(amat[1,,]), byrow = T, ncol = dim(amat)[2])
  boxplot(a, names = dimnames(amat)[[2]], col = cols, las = 2, main = 'RMSE pix')
  
  a = matrix(as.numeric(amat[2,,]), byrow = T, ncol = dim(amat)[2])
  boxplot(a, names = dimnames(amat)[[2]], col = cols, las = 2, main = 'RMSE f1')
  
  a = matrix(as.numeric(amat[3,,]), byrow = T, ncol = dim(amat)[2])
  boxplot(a, names = dimnames(amat)[[2]], col = cols, las = 2, main = 'RMSE localfdr')
  
  lfdp_index = which(sapply(dimnames(amat)[[1]] ,function(ss) grepl('lfdp', ss)))
  ltpp_index = which(sapply(dimnames(amat)[[1]] ,function(ss) grepl('ltpp', ss)))
  
  perfcurves(amat[c(lfdp_index, ltpp_index),,], cols, fdp_list = fdp_list, aggfunc)
  
}

# for plotting results of one replicate
# it is assumed that fdrreg is coerced into carrying our variable names: p, f1y
plotfn = function(robj,
                  cols = c('orange', 'red', 'blue', 'green4',
                           'cyan3','violetred1', 'slateblue1'),
                  ignore_strings = c('dd', 'gs', 'pr_init_')){
  with(robj, {
    
    mn = setdiff(names(robj), ignore_strings)
    
    f1err = round(sapply(mn, function (mni) rmse(get(mni)$f1y, dd$f1y)), 4)
    pierr = round(sapply(mn, function (mni) rmse(get(mni)$p, dd$pix)), 4)
    lfdrerr = round(sapply(mn, function (mni) rmse(get(mni)$localfdr, dd$localfdr)), 4)
    denerr = round(sapply(mn, function (mni) rmse(get(mni)$den, dd$den)), 4)
    
    f1main = paste(c(rbind(rmlast(mn), f1err)), collapse =', ')
    pimain = paste(c(rbind(rmlast(mn), pierr)), collapse =', ')
    lfdrmain = paste(c(rbind(rmlast(mn), lfdrerr)), collapse =', ')
    denmain = paste(c(rbind(rmlast(mn), denerr)), collapse =', ')
    
    par(mfrow = c(2,2), bty = 'l')
    
    # estimates of f1
    plot(dd$y, dd$f1y, cex = 0.3,
         ylim = c(0, 1.3 * max(dd$f1y)),
         ylab = 'f1', xlab = 'y', main = f1main)
    for (i in 1:length(mn)){
      points(dd$y, get(mn[i])$f1y, col = cols[(i-1) %% length(cols) + 1], cex = 0.3)
    }
    
    # estimates of pi
    plot(sx(dd$x), dd$pix, cex = 0.3,
         ylim = c(0,1), ylab = 'pi', xlab = 'sx(x)', type ='n', main = pimain)
    for (i in 1:length(mn)){
      points(sx(dd$x), get(mn[i])$p, col = cols[(i-1) %% length(cols) + 1], cex = 0.3)
    }
    points(sx(dd$x), dd$pix, cex = 0.3)
    legend(x = 'topleft', bty = 'n', legend = sapply(mn, rmlast), col = cols, lwd = 4)
    
    # estimates of localfdr
    plot(dd$localfdr, get(mn[1])$localfdr,
         cex = 0.3, col = cols[1], xlim = c(0,1), ylim = c(0,1),
         ylab = 'estimated lfdr', xlab = 'true lfdr', main = lfdrmain)
    for (i in 2:length(mn)){
      points(dd$localfdr, get(mn[i])$localfdr, col = cols[(i-1) %% length(cols) + 1], cex = 0.3)
    }
    abline(c(0,1))
    
    # this is only for setting 76
    # plot(dd$x[,1], get(mn[1])$localfdr,
    #      cex = 0.3, col = cols[1], xlim = c(0,1), ylim = c(0,1),
    #      ylab = 'localfdr', xlab = 'x1', main = lfdrmain)
    # for (i in 2:length(mn)){
    #   points(dd$x[,1], get(mn[i])$localfdr, col = cols[(i-1) %% length(cols) + 1], cex = 0.3)
    # }
    # points(dd$x[,1], dd$localfdr, cex = 0.3)
    
    # estimates of denominator
    plot(dd$den, get(mn[1])$den,
         cex = 0.3, col = cols[1],
         ylab = 'estimated denominator', xlab = 'true denominator', main = denmain)
    for (i in 2:length(mn)){
      points(dd$den, get(mn[i])$den, col = cols[(i-1) %% length(cols) + 1], cex = 0.3)
    }
    abline(c(0,1))
    
    
    
    par(mfrow = c(1,1))
  })
}

