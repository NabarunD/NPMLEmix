library(foreach)
library(doParallel)
library(R.utils)

nCores <- Sys.getenv("SLURM_CPUS_PER_TASK")
registerDoParallel(cores = nCores)

require(FDRreg, quietly = TRUE)
require(R.utils, quietly = TRUE)
require(Rcpp, quietly = TRUE)
require(Rmosek, quietly = TRUE)
require(REBayes, quietly = TRUE)
require(CAMAN, quietly = TRUE)
require(progress, quietly = TRUE)
require(pbapply, quietly = TRUE)
require(Hmisc, quietly = TRUE)

source('C:/Users/Nabarun Deb/Dropbox/CovariateTesting/Code/Code_Dump_May_2018/source/opt_routines.R')
source('C:/Users/Nabarun Deb/Dropbox/CovariateTesting/Code/Code_Dump_May_2018/source/utilities.R')
source('C:/Users/Nabarun Deb/Dropbox/CovariateTesting/Code/Code_Dump_May_2018/source/my_list.R')
source('C:/Users/Nabarun Deb/Dropbox/CovariateTesting/Code/Code_Dump_May_2018/source/marginal_methods.R')
LIMIT_THREAD = TRUE

niter = 1
n = 1e3

runSim = function(n, sx, tdparams, file_prefix, iter){

  # constructing data
  dd = makedata(n, sx, tdparams)

  # FDRreg solution
  ff=FDRreg(dd$y,dd$xs)
  f_ = modf(ff)

  # marginal methods
  pi0grid = seq(from = 0.01, to = 0.49, by = 0.01);
  m1_ = with(dd, m1(y, cbind(1, xs), pi0grid))
  m2_ = with(dd, m2(y, cbind(1, xs), pi0grid, solvef1 = TRUE))

  # we have three possible initializations
  # we should also keep track of how many times which initialization achieved a higher likelihood
  #init_list = list(f_ = f_)
  init_list = list(f_=f_, m1_ = m1_, m2_ = m2_)
  init_bi = which.max(sapply(init_list, function (ro) ro$ll))
  init_best = extract_init(init_list[[init_bi]])
  init_best_name = names(init_bi)

  start_time <- Sys.time()
  # EM starting from here, run at most 500 iterations
  em_ = with(dd, lgem(y, cbind(1, xs), weights = init_best$w, binit = init_best$b, timed = Inf, maxit = 100))
  em_ = extract_index(em_, length(em_$time_list))
  end_time <- Sys.time()
  print(end_time-start_time)

  robj = list(dd = dd, init_best_name = init_best_name,
              f_ = f_, m1_ = m1_, m2_ = m2_, em_ = em_)

  base::save(robj, file = paste0(file_prefix, '_iter', iter, '.Rdata'))
}

for (ii in 1:4){
  for (jj in 1:4){

    #indices = foreach(i = 1:niter, .combine = 'rbind') %dopar% {

      # running
      sx = slist[[ii]]
      tdparams = tdlist[[jj]]
      #file_prefix = paste0('./sim_output/n', n, '/finalsim', ii, jj)

      tt = system.time({
        runFlag = 1;
        while(runFlag){
          tryCatch({
            runSim(n, sx, tdparams, file_prefix, i)
            runFlag = 0
          }, error = function (e) runFlag <<- 1)
        }
      })['elapsed'];
      print(c(ii, jj, i))
      print(tt)

    }

  } # end loop over jj
} # end loop over ii


###############################################

marg1=NULL
marg2=NULL
truell=NULL
c1good=0
c2good=0
for(ii in 1:4)
{
  for(jj in 1:4)
  {
    print(ii)
    print(jj)
    sx = slist[[ii]]
    tdparams = tdlist[[jj]]
    dd = makedata(n, sx, tdparams)
    pi0grid = seq(from = 0.01, to = 0.49, by = 0.01);
    m1_ = with(dd, m1(y, cbind(1, xs), pi0grid))
    m2_ = with(dd, m2(y, cbind(1, xs), pi0grid, solvef1 = TRUE))
    marg1=c(marg1,m1_$ll)
    marg2=c(marg2,m2_$ll)
    truell=c(truell,dd$ll)
    if(m1_$ll>dd$ll)
      c1good=c1good+1
    if(m2_$ll>dd$ll)
      c2good=c2good+1
  }
}
plot(marg1,marg2,ty='l')
ii=1
jj=1
dd = makedata(n, sx, tdparams)
sx = slist[[ii]]
tdparams = tdlist[[jj]]
pi0grid = seq(from = 0.01, to = 0.49, by = 0.01);
m1_ = with(dd, m1(y, cbind(1, xs), pi0grid))
m2_ = with(dd, m2(y, cbind(1, xs), pi0grid, solvef1 = TRUE))
m1_$ll
m2_$ll

n=1e02
ii=1
jj=1
sx = slist[[ii]]
tdparams = tdlist[[jj]]
dd = makedata(n, sx, tdparams)
pi0grid = seq(from = 0.01, to = 0.49, by = 0.01);
start_time <- proc.time()
ff=FDRreg(dd$y,dd$xs)
f_ = modf(ff)
st=proc.time() - start_time
print(end_time-start_time)

start_time <- proc.time()
# marginal methods
pi0grid = seq(from = 0.01, to = 0.49, by = 0.01);
m1_ = with(dd, m1(y, cbind(1, xs), pi0grid))
st = proc.time() - start_time
as.numeric(st[3])
print(end_time-start_time)

start_time <- Sys.time()
m2_ = with(dd, m2(y, cbind(1, xs), pi0grid, solvef1 = TRUE))
end_time <- Sys.time()
print(end_time-start_time)

init_list = list(f_=f_, m1_ = m1_, m2_ = m2_)
init_bi = which.max(sapply(init_list, function (ro) ro$ll))
init_best = extract_init(init_list[[init_bi]])
init_best_name = names(init_bi)

start_time <- Sys.time()
# EM starting from here, run at most 500 iterations
em_ = with(dd, lgem(y, cbind(1, xs), weights = init_best$w, binit = init_best$b, timed = Inf, maxit = 100))
em_ = extract_index(em_, length(em_$time_list))
end_time <- Sys.time()
print(end_time-start_time)

init_list = list(f_=f_)
init_bi = which.max(sapply(init_list, function (ro) ro$ll))
init_best = extract_init(init_list[[init_bi]])
init_best_name = names(init_bi)

# EM starting from here, run at most 500 iterations
em1_ = with(dd, lgem(y, cbind(1, xs), weights = init_best$w, binit = init_best$b, timed = Inf, maxit = 100))
em1_ = extract_index(em1_, length(em1_$time_list))

init_list = list(m1_ = m1_)
init_bi = which.max(sapply(init_list, function (ro) ro$ll))
init_best = extract_init(init_list[[init_bi]])
init_best_name = names(init_bi)

# EM starting from here, run at most 500 iterations
em2_ = with(dd, lgem(y, cbind(1, xs), weights = init_best$w, binit = init_best$b, timed = Inf, maxit = 100))
em2_ = extract_index(em2_, length(em2_$time_list))

init_list = list(m2_ = m2_)
init_bi = which.max(sapply(init_list, function (ro) ro$ll))
init_best = extract_init(init_list[[init_bi]])
init_best_name = names(init_bi)

# EM starting from here, run at most 500 iterations
em3_ = with(dd, lgem(y, cbind(1, xs), weights = init_best$w, binit = init_best$b, timed = Inf, maxit = 100))
em3_ = extract_index(em3_, length(em3_$time_list))

sam=1000
t1=rnorMix(sam,norMix(mu=em1_$atoms,sigma=rep(1,length(em1_$atoms)),w=em1_$probs))
t2=rnorMix(sam,norMix(mu=em2_$atoms,sigma=rep(1,length(em2_$atoms)),w=em2_$probs))
t3=rnorMix(sam,norMix(mu=em3_$atoms,sigma=rep(1,length(em3_$atoms)),w=em3_$probs))

x <- data.frame(v1=t1,v2=t2,v3=t3)
library(ggplot2);library(reshape2)
data<- melt(x)
ggplot(data,aes(x=value, fill=variable)) + geom_density(alpha=0.25)
ggplot(data,aes(x=value, fill=variable)) + geom_histogram(alpha=0.25)
ggplot(data,aes(x=variable, y=value, fill=variable)) + geom_boxplot()


n=1e02
ii=1
jj=1
sx = slist[[ii]]
tdparams = tdlist[[jj]]
dd = makedata(n, sx, tdparams)
pi0grid = seq(from = 0.01, to = 0.49, by = 0.01);
ff=FDRreg(dd$y,dd$xs)
f_ = modf(ff)

# marginal methods
pi0grid = seq(from = 0.01, to = 0.49, by = 0.01);
m1_ = with(dd, m1(y, cbind(1, xs), pi0grid))
m2_ = with(dd, m2(y, cbind(1, xs), pi0grid, solvef1 = TRUE))

init_list = list(f_=f_, m1_ = m1_, m2_ = m2_)
init_bi = which.max(sapply(init_list, function (ro) ro$ll))
init_best = extract_init(init_list[[init_bi]])
init_best_name = names(init_bi)

# EM starting from here, run at most 500 iterations
em_ = with(dd, lgem(y, cbind(1, xs), weights = init_best$w, binit = init_best$b, timed = Inf, maxit = 100))
em_ = extract_index(em_, length(em_$time_list))

init_list = list(f_=f_)
init_bi = which.max(sapply(init_list, function (ro) ro$ll))
init_best = extract_init(init_list[[init_bi]])
init_best_name = names(init_bi)

# EM starting from here, run at most 500 iterations
em12_ = with(dd, lgem(y, cbind(1, xs), weights = init_best$w, binit = init_best$b, timed = Inf, maxit = 100))
em12_ = extract_index(em12_, length(em12_$time_list))

init_list = list(m1_ = m1_)
init_bi = which.max(sapply(init_list, function (ro) ro$ll))
init_best = extract_init(init_list[[init_bi]])
init_best_name = names(init_bi)

# EM starting from here, run at most 500 iterations
em22_ = with(dd, lgem(y, cbind(1, xs), weights = init_best$w, binit = init_best$b, timed = Inf, maxit = 100))
em22_ = extract_index(em22_, length(em22_$time_list))

init_list = list(m2_ = m2_)
init_bi = which.max(sapply(init_list, function (ro) ro$ll))
init_best = extract_init(init_list[[init_bi]])
init_best_name = names(init_bi)

# EM starting from here, run at most 500 iterations
em32_ = with(dd, lgem(y, cbind(1, xs), weights = init_best$w, binit = init_best$b, timed = Inf, maxit = 100))
em32_ = extract_index(em32_, length(em32_$time_list))

sam=1000
t1=rnorMix(sam,norMix(mu=em12_$atoms,sigma=rep(1,length(em12_$atoms)),w=em12_$probs))
t2=rnorMix(sam,norMix(mu=em22_$atoms,sigma=rep(1,length(em22_$atoms)),w=em22_$probs))
t3=rnorMix(sam,norMix(mu=em32_$atoms,sigma=rep(1,length(em32_$atoms)),w=em32_$probs))

x <- data.frame(v1=st[,1],v2=st[,2],v3=st[,3])
library(ggplot2);library(reshape2)
data<- melt(x)
#ggplot(data,aes(x=value, fill=variable)) + geom_density(alpha=0.25)
#ggplot(data,aes(x=value, fill=variable)) + geom_histogram(alpha=0.25)
ggplot(data,aes(x=variable, y=value, fill=variable)) + geom_boxplot()

cutoff=function(vec1,alph)
{
  vecc=sort(vec1)
  mvec1=cumsum(vecc)/(1:length(vecc))
  cc=which(mvec1<=alph)
  rejset=as.numeric(vec1<=vecc[max(cc)])
  truerej=numeric(n)
  truerej=replace(truerej,dd$nnind,1)
  st=table(rejset,truerej)
  return(c(st[2]/sum(rejset),st[4]/sum(truerej)))
}
tt=c(0.05,0.10,0.15,0.20,0.25,0.30)
cutoff(f_$localfdr,tt[4])
cutoff(m1_$localfdr,tt[4])
cutoff(m2_$localfdr,tt[4])
cutoff(em_$localfdr,tt[4])
mean(f_$localfdr[(which(f_$localfdr<=0.1))])
mean(m1_$localfdr[(which(m1_$localfdr<=0.1))])
sum(1-f_$localfdr)
sum(1-m1_$localfdr)
