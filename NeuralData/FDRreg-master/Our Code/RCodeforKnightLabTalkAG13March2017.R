#Open the data
setwd("/Users/aditya/Dropbox/CovariateTesting/NeuralData/FDRreg-master/data")
synchrony_smithkohn2008 = read.csv('synchrony_smithkohn2008.csv', header=TRUE)
sys = synchrony_smithkohn2008
#Only the last three columns seem relevant to us. It is not clear though as to how the Dist variable is calculated.
z = sys$z
#Histogram of the z scores
#The total number of data points seems to 7004 which is smaller than 128*127/2. So not all pairs of neurons are recorded?
#Histogram of the z values
hist(z, 250)

ddfull = cbind(synchrony_smithkohn2008$Dist, synchrony_smithkohn2008$TuningCor, synchrony_smithkohn2008$z)
ddfull = as.data.frame(ddfull)
names(ddfull) = c('Dist', 'TuningCor', 'z')

#Load the FDR Regression package:
require(FDRreg)

# the null distribution is not modeled as N(0,1) but rather as N(mu,sig)
# where mu and sig are computed using Efron's ML estimate in a two-groups model
efron_res = FDRreg::efron(ddfull$z, nulltype='empirical', df=10)
mu0hat = efron_res$mu0; sig0hat = efron_res$sig0;

# using Scott et al's Empirical Bayes FDR Regression
X = ddfull[,1:2]
# Set up spline basis functions (df=3)
df = 3
b1 = bs(synchrony_smithkohn2008$Dist, df=df)
b2 = bs(synchrony_smithkohn2008$TuningCor, df=df)
Xs = model.matrix( ~  b1 + b2 - 1)
# produces a matrix with 6 column, first three are splines corresponding to Dist, next three corr. to TuningCor
# FDRreg analysis uses a scaled and centered Xs
# lets use the same
Xs = scale(Xs)

# for our methods, need to rescale z so that nulls are N(0,1)
zscaled = (ddfull$z - mu0hat)/sig0hat

#Let us just work with these rescaled values throughout.
#Histogram of these values

hist(zscaled, breaks = 500, freq = F, main = "Histogram of the test statistics (z)", xlab = "Test Statistics (z)")
#Superimpose the standard normal density
points(sort(zscaled), dnorm(sort(zscaled)), type = "l", col = "red")
#points(sort(zscaled), dnorm(sort(zscaled), mu0, sigma0), type = "l", col = "red")

legend("topright", c("N(0, 1) density"), lty = c(1), lwd = c(1), col = c("red"))


#Run two methods (Benjamini-Hochberg and local FDR with predictive recursion) which do not use covariates. 
#Benjamini-Hochberg
help(BenjaminiHochberg)
BH1 = BenjaminiHochberg(zscaled, 0.1)
#Basically here the z-values are converted to p-values by Prob that |Z| is at least |z| where Z is standard normal
#and then one applies the usual BH procedure to the p-values. 
sum(BH1 == 1 & z > 0)
#There seem to be 329 rejections.
min(zscaled[BH1 == 1])
#This minimum is 2.82922. So every z-value larger than 2.82922 is rejected. 
summary(zscaled[BH1 == 1])
#Looks like all the rejected p-values are positive. There is no significant p-value (according to BH) that is negative.
hist(zscaled, 500, main = "Rejection Regions (both methods use level = 0.1)", xlab = "Test Statistics (z)")
min(zscaled[BH1 == 1 & z > 0]) #This is 2.836328
abline(v = min(zscaled[BH1 == 1 & z > 0]), col = "blue")
sum(BH1 == 1 & z > 0) #299 rejections

#Efron's local FDR with (1) the empirical null estimated according to Efron's method, (2) alternative density estimated via predictive recursion, and (3) pi_0 given also by the predictive recursion algorithm. Can the Patra and Sen method be used for this? 
#e1 = efron(zscaled, nulltype = 'empirical', df = 10)
#There seems to be no help file for this efron function. The code for this function can be got by typing efron in R.
#names(e1)
#[1] "mids"    "breaks"  "zcounts" "zdens"   "z"       "fz"      "mu0"    
# [8] "sig0"    "p0"      "fdr"    
#This procedure for estimating the empirical null comes from Section 4 in the paper "Size, Power and False Discovery Rates" 
#by Bradley Efron (published in the Annals of Statistics). 
#Is the empirical null Gaussian (with mean mu0 and standard deviation sig0)? Here e1$mu0 = 0.6081196 and e1$sig0 = 0.8140629. 
#With this f0, the next step seems to be to estimate the alternative density f1. For this, they are using predictive recursion.
#mu0 = e1$mu0
#sigma0 = e1$sig0
fhat = prfdr(zscaled, mu0 = 0, sig0 = 1)
#The most important output of fhat is postprob which gives the probability of being a signal for each observation i.e., the probability of the null hypothesis being false. 
postprob = fhat$postprob
BFDR = getFDR(fhat$postprob)
#This getFDR function outputs two quantities: localfdr and FDR. localfdr is simply 1 - postprob. FDR(z) is defined as the average of lfdr conditional on Z \leq z. It is essentially computed by ordering the postprob, taking cumulative sums and then reordering them. Scott et al seem to put cutoffs on FDR this as their ultimate procedure to reject p-values (instead of putting cutoffs on local fdr). 
sum(BFDR$FDR <= 0.1 & z > 0) #This gives the number of rejections by the local FDR method. 
#Every time I run this, I seem to be getting different answers. The predictive recursion algorithm seems to give different answers for each run. We can check this easily. 
#f1 =  prfdr(z, mu0 = mu0, sig0 = sigma0, control = list(npasses = 20, gridsize = 500))
#f2 =  prfdr(z, mu0 = mu0, sig0 = sigma0, control = list(npasses = 20, gridsize = 500))
#head(cbind(f1$postprob, f2$postprob))

min(zscaled[BFDR$FDR <= 0.1 & z > 0]) #2.406
abline(v = min(zscaled[BFDR$FDR <= 0.1 & z > 0]), col = "red")
sum(BFDR$FDR <= 0.1 & z > 0) #497 rejections. 

legend("topright", c("B-H (299 rejections, z > 2.836)", "Local FDR (497 rejections, z > 2.406)"), lty = c(1, 1), lwd = c(1, 1), col = c("blue", "red"))

#Plotting the test statistics against the covariates 

plot(X[,1], zscaled, main = "z vs Distance", xlab = "Inter-neuron distances (X1)", ylab = "Test Statistic for Synchrony", cex = 0.2)
abline(lm(zscaled ~ X[,1]), col = "blue")

plot(X[,2], zscaled, main = "z vs Tuning correlations", xlab = "Tuning Curve Correlations (X2)", ylab = "Test Statistic for Synchrony", cex = 0.2)
m1 = lm(zscaled ~ X[,2] + I(X[,2]^2))
summary(m1)
points(X[order(X[,2]), 2], m1$fitted.values[order(X[,2])], type = "l", col = "red")




#Run Scott's method:
scott_res = FDRreg(ddfull$z, Xs, nulltype='empirical', control=list(lambda = 1))
names(scott_res)
plot(c(scott_res$x_grid, scott_res$x_grid), c(scott_res$f1_grid, scott_res$f0_grid), type = "n")
#This does seem like a density:
sum(diff(scott_res$x_grid)*scott_res$f1_grid[-1])
sum(diff(scott_res$x_grid)*scott_res$f0_grid[-1])
points(scott_res$x_grid, scott_res$f0_grid, type = "l", col = "blue")
points(scott_res$x_grid, scott_res$f1_grid, type = "l", col = "red")
legend("topright", c("Null density (f0)", "Alternative density (f1)"), lty = c(1, 1), lwd = c(1, 1), col = c("blue", "red"))
#Note that we do not work with the original data z. But we rather work with zscaled. 
#Let us re-plot f0 and f1 after rescaling. 
scal_grid = (scott_res$x_grid - mu0hat)/sig0hat
scott_f1_scal = scott_res$f1_grid*sig0hat
scott_f0_scal = scott_res$f0_grid*sig0hat
plot(c(scal_grid, scal_grid), c(scott_f0_scal, scott_f1_scal), type = "n")
points(scal_grid, scott_f0_scal, type = "l", col = "blue")
points(scal_grid, scott_f1_scal, type = "l", col = "red")
legend("topright", c("Null density (f0)", "Alternative density (f1)"), lty = c(1, 1), lwd = c(1, 1), col = c("blue", "red"))
#It is easy to check that f0 now corresponds to the standard normal density:
#points(scal_grid, dnorm(scal_grid), col = "green")

#The results for our methods are all saved and can be retrieved via load('save_for_now.Rdata'). 
setwd("/Users/aditya/Dropbox/CovariateTesting/NeuralData/FDRreg-master/Our Code")
load('save_for_now.Rdata')

#The results for our methods are saved in four objects: gridy_res, am1_res, m2_res, am2_res:
#1) Marginal Method One: gridy_res
#2) Marginal Method Two: m2_res (this apparently does not work well). 
#3) Full MLE with Marginal One as Initialization: am1_res
#4) Full MLE with Marginal Two as Initialization: am2_res

names(gridy_res)
names(m2_res)
names(am1_res)
names(am2_res)

#Plotting the estimates for f1 from all the methods:
plot(c(scal_grid, scal_grid, scal_grid, scal_grid), c(scott_f1_scal, f1_m1, f1_fullm1, f1_fullm2), type = "n", xlab = "z", ylab = "Alternative Density (f1)", main = "Estimates of f1")

#1) Scott's method
points(scal_grid, scott_f1_scal, type = "l")

#2) gridy_res
f1_m1 = rep(0, length(scal_grid))
for(i in 1:length(scal_grid))
{
    f1_m1[i] = sum(gridy_res$probs*dnorm(scal_grid[i] - gridy_res$atoms))
}
points(scal_grid, f1_m1, type = "l", col = "blue")

#3) m2_res
f1_m2 = rep(0, length(scal_grid))
for(i in 1:length(scal_grid))
{
    f1_m2[i] = sum(m2_res$probs*dnorm(scal_grid[i] - m2_res$atoms))
}
#points(scal_grid, f1_m2, type = "l", col = "purple")


#3) am1_res
f1_fullm1 =  rep(0, length(scal_grid))
for(i in 1:length(scal_grid))
{
    f1_fullm1[i] = sum(am1_res$probs*dnorm(scal_grid[i] - am1_res$atoms))
}
points(scal_grid, f1_fullm1, type = "l", col = "green")

#4) am2_res
f1_fullm2 =  rep(0, length(scal_grid))
for(i in 1:length(scal_grid))
{
    f1_fullm2[i] = sum(am2_res$probs*dnorm(scal_grid[i] - am2_res$atoms))
}
points(scal_grid, f1_fullm2, type = "l", col = "red")

legend("topright", c("Scott's Method", "Marginal One", "Full MLE - M1", "Full MLE - M2"), lty = c(1, 1, 1, 1), lwd = c(1, 1, 1, 1), col = c("black", "blue", "green", "red"))

#Plotting the estimates of the estimated function pi from all the methods: 
#Because this is an additive model, I will plot pi in two functions (as a function of the first variable and as a function of the second variable)

plot(c(X[or1,1], X[or1,1], X[or1,1], X[or1, 1]), c(scott_res$model$coef[2]*Xs[or1,1] + scott_res$model$coef[3]*Xs[or1,2] + scott_res$model$coef[4]*Xs[or1,3], gridy_res$b[2]*Xs[or1,1] + gridy_res$b[3]*Xs[or1,2] + gridy_res$b[4]*Xs[or1,3], am1_res$b[2]*Xs[or1,1] + am1_res$b[3]*Xs[or1,2] + am1_res$b[4]*Xs[or1,3], am2_res$b[2]*Xs[or1,1] + am2_res$b[3]*Xs[or1,2] + am2_res$b[4]*Xs[or1,3]), type = "n", xlab = "Inter-neuron Distance (X1)", ylab = "The (logit of) Pi function (as a function of X1)", main = "Plot of the (logit of) Pi Function in terms of X1") 

#For Scott's method:
or1 = order(X[,1])
points(X[or1,1], scott_res$model$coef[2]*Xs[or1,1] + scott_res$model$coef[3]*Xs[or1,2] + scott_res$model$coef[4]*Xs[or1,3], type = "o", cex = 0.5)

#For gridy_res
or1 = order(X[,1])
points(X[or1,1], gridy_res$b[2]*Xs[or1,1] + gridy_res$b[3]*Xs[or1,2] + gridy_res$b[4]*Xs[or1,3], type = "o", cex = 0.5, col = "blue")

#For m2_res
or1 = order(X[,1])
#points(X[or1,1], m2_res$b[2]*Xs[or1,1] + m2_res$b[3]*Xs[or1,2] + m2_res$b[4]*Xs[or1,3], type = "o", cex = 0.5, col = "purple")


#For am1_res
or1 = order(X[,1])
points(X[or1,1], am1_res$b[2]*Xs[or1,1] + am1_res$b[3]*Xs[or1,2] + am1_res$b[4]*Xs[or1,3], type = "o", cex = 0.5, col = "green")

#For am2_res
or1 = order(X[,1])
points(X[or1,1], am2_res$b[2]*Xs[or1,1] + am2_res$b[3]*Xs[or1,2] + am2_res$b[4]*Xs[or1,3], type = "o", cex = 0.5, col = "red")

legend("topright", c("Scott's Method", "Marginal One", "Full MLE - M1", "Full MLE - M2"), lty = c(1, 1, 1, 1), lwd = c(1, 1, 1, 1), col = c("black", "blue", "green", "red"))


#Plotting pi as a function of the tuning curve correlation

plot(c(X[or2,2], X[or2,2], X[or2,2], X[or2, 2]), c(scott_res$model$coef[5]*Xs[or2,4] + scott_res$model$coef[6]*Xs[or2,5] + scott_res$model$coef[7]*Xs[or2,6], gridy_res$b[5]*Xs[or2,4] + gridy_res$b[6]*Xs[or2,5] + gridy_res$b[7]*Xs[or2,6], am1_res$b[5]*Xs[or2,4] + am1_res$b[6]*Xs[or2,5] + am1_res$b[7]*Xs[or2,6], am2_res$b[5]*Xs[or2,4] + am2_res$b[6]*Xs[or2,5] + am2_res$b[7]*Xs[or2,6]), type = "n", xlab = "Tuning Curve Correlation (X2)", ylab = "The (logit of) Pi function (as a function of X2)", main = "Plot of the (logit of) Pi function in terms of X2") 

#For Scott's method:
or2 = order(X[,2])
points(X[or2,2], scott_res$model$coef[5]*Xs[or2,4] + scott_res$model$coef[6]*Xs[or2,5] + scott_res$model$coef[7]*Xs[or1,6], type = "o", cex = 0.5)

#For gridy_res
points(X[or2,2], gridy_res$b[5]*Xs[or2,4] + gridy_res$b[6]*Xs[or2,5] + gridy_res$b[7]*Xs[or2,6], type = "o", cex = 0.5, col = "blue")

#For m2_res
#points(X[or2,2], m2_res$b[5]*Xs[or2,4] + m2_res$b[6]*Xs[or2,5] + m2_res$b[7]*Xs[or2,6], type = "o", cex = 0.5, col = "purple")


#For am1_res
points(X[or2,2], am1_res$b[5]*Xs[or2,4] + am1_res$b[6]*Xs[or2,5] + am1_res$b[7]*Xs[or2,6], type = "o", cex = 0.5, col = "green")

#For am2_res
points(X[or2,2], am2_res$b[5]*Xs[or2,4] + am2_res$b[6]*Xs[or2,5] + am2_res$b[7]*Xs[or2,6], type = "o", cex = 0.5, col = "red")

legend("top", c("Scott's Method", "Marginal One", "Full MLE - M1", "Full MLE - M2"), lty = c(1, 1, 1, 1), lwd = c(1, 1, 1, 1), col = c("black", "blue", "green", "red"))


#Plotting the pi(xi) values as a function of xi^Tbeta for a fixed beta. 
Xs1 = cbind(rep(1, length(Xs[,1])), Xs)
xaxs = Xs1%*%am1_res$b
or = order(xaxs)
plot(c(xaxs[or],xaxs[or],xaxs[or],xaxs[or]), c(scott_res$priorprob[or], gridy_res$p[or], am1_res$p[or], am2_res$p[or]), type = "n", xlab = "X * beta", ylab = "Prior Probabilities", main = "Plot of Estimates of the Pi function")

#Scott's method
points(xaxs[or], scott_res$priorprob[or], type = "l")

#gridy_res
points(xaxs[or], gridy_res$p[or], type = "l", col = "blue")

#m2_res
#points(xaxs[or], m2_res$p[or], type = "l", col = "purple")

#am1_res
points(xaxs[or], am1_res$p[or], type = "l", col = "green")

#am2_res
points(xaxs[or], am2_res$p[or], type = "l", col = "red")

legend("topleft", c("Scott's Method", "Marginal One", "Full MLE - M1", "Full MLE - M2"), lty = c(1, 1, 1, 1), lwd = c(1, 1, 1, 1), col = c("black", "blue", "green", "red"))

#Understanding the rejections:
#Finding the rejection sets for each of the methods (nominal level is 0.1)
#Scott's method: 
scott_rejects = which(getFDR(scott_res$postprob)$FDR <= 0.1 & ddfull$z > 0)
length(scott_rejects)

#Sujayam's function
lfdr_reject = function(lfdr, fdr_nominal = 0.1){
  sl = sort(lfdr)
  k = sum(cumsum(sl)/seq_along(sl) <= fdr_nominal)
  if(k) rejects = which(lfdr <= sl[k]) else rejects = numeric(0)
}

#gridy_res
gridy_f1_zscaled = sapply(gridy_res$atoms, function(ai) dnorm(zscaled - ai)) %*% gridy_res$probs
gridy_lfdr = (1 - gridy_res$p) * dnorm(zscaled)/((1 - gridy_res$p) * dnorm(zscaled) + gridy_res$p * gridy_f1_zscaled)
gridy_rejects = lfdr_reject(gridy_lfdr)
length(gridy_rejects)
#Slightly fewer rejections compared to Scott's method (makes sense from the plot of the pi's)

#m2_res
m2_f1_zscaled = sapply(m2_res$atoms, function(ai) dnorm(zscaled - ai)) %*% m2_res$probs
m2_lfdr = (1 - m2_res$p) * dnorm(zscaled)/((1 - m2_res$p) * dnorm(zscaled) + m2_res$p * m2_f1_zscaled)
m2_rejects = lfdr_reject(m2_lfdr)
length(m2_rejects)
#This does not work at all. 

#am1_res
am1_f1_zscaled = sapply(am1_res$atoms, function(ai) dnorm(zscaled - ai)) %*% am1_res$probs
am1_lfdr = (1 - am1_res$p) * dnorm(zscaled)/((1 - am1_res$p) * dnorm(zscaled) + am1_res$p * am1_f1_zscaled)
am1_rejects = lfdr_reject(am1_lfdr)
length(am1_rejects)
#This (and the next method) unsurprisingly give the most rejections. 

#am2_res
am2_f1_zscaled = sapply(am2_res$atoms, function(ai) dnorm(zscaled - ai)) %*% am2_res$probs
am2_lfdr = (1 - am2_res$p) * dnorm(zscaled)/((1 - am2_res$p) * dnorm(zscaled) + am2_res$p * am2_f1_zscaled)
am2_rejects = lfdr_reject(am2_lfdr)
length(am2_rejects)
#This is identical to am1. 

c(length(scott_rejects), length(gridy_rejects), length(m2_rejects), length(am1_rejects), length(am2_rejects))

#Additional Rejects (of Full MLE compared to Local FDR without Covariates)
noco_rejects = which(BFDR$FDR <= 0.1 & z > 0)
par(mfrow = c(2, 1))
sid = setdiff(am1_rejects, noco_rejects)
plot(X[,1], X[,2], xlab = "Inter-neuron Distances (X1)", ylab = "Tuning Curve Correlation (X2)", main = "Extra rejections compared to No Covariate method (575)", type = "n" )
points(X[sid, 1], X[sid, 2])
length(sid)

sid = setdiff(noco_rejects, am1_rejects)
plot(X[,1], X[,2], xlab = "Inter-neuron Distances (X1)", ylab = "Tuning Curve Correlation (X2)", main = "Rejected by the No-Covariate method but not by our method (52)", type = "n" )
points(X[sid, 1], X[sid, 2])
length(sid)

#Histogram
par(mfrow = c(1, 1))
sid = setdiff(am1_rejects, noco_rejects)
hist(zscaled[sid], breaks = 150)
plot(X[,1], X[,2], xlab = "Inter-neuron Distances (X1)", ylab = "Tuning Curve Correlation (X2)", main = "Extra rejections compared to No Covariate method", type = "n" )
points(X[sid, 1], X[sid, 2])

#Additional rejects (of full MLE compared to Scott's method)
par(mfrow = c(2, 1))
sid = setdiff(am1_rejects, scott_rejects)
plot(X[,1], X[,2], xlab = "Inter-neuron Distances (X1)", ylab = "Tuning Curve Correlation (X2)", main = "Extra rejections compared to Scott's method (276)", type = "n" )
points(X[sid, 1], X[sid, 2])
length(sid)


sid = setdiff(scott_rejects, am1_rejects)
plot(X[,1], X[,2], xlab = "Inter-neuron Distances (X1)", ylab = "Tuning Curve Correlation (X2)", main = "Rejected by Scott's method but not by our method (8)", type = "n" )
points(X[sid, 1], X[sid, 2])
length(sid)

#Rejection Regions of the four methods (Scott et al, no covariate lfdr, marginal one and am1)
par(mfrow = c(2, 2))
si = noco_rejects
length(si)
plot(X[,1], X[,2], xlab = "Inter-neuron Distances (X1)", ylab = "Tuning Curve Correlation (X2)", main = "Rejections of the No-Covariate method (497)", type = "n")
points(X[si, 1], X[si, 2])
si = scott_rejects
length(si)
plot(X[,1], X[,2], xlab = "Inter-neuron Distances (X1)", ylab = "Tuning Curve Correlation (X2)", main = "Rejections of the Scott et al. (2015) method (752)", type = "n")
points(X[si, 1], X[si, 2])
si = gridy_rejects
length(si)
plot(X[,1], X[,2], xlab = "Inter-neuron Distances (X1)", ylab = "Tuning Curve Correlation (X2)", main = "Rejections of the First Marginal Method (722)", type = "n")
points(X[si, 1], X[si, 2])
si = am1_rejects
length(si)
plot(X[,1], X[,2], xlab = "Inter-neuron Distances (X1)", ylab = "Tuning Curve Correlation (X2)", main = "Rejections of the Full MLE initialized by M1 (1020)", type = "n")
points(X[si, 1], X[si, 2])
par(mfrow = c(1, 1))



#Box Plot of the estimated pi values: 
boxplot(cbind(scott_res$priorprob, gridy_res$p, am1_res$p, am2_res$p), names = c("Scott", "Marg 1", "Full MLE - 1", "Full MLE - 2"), main = "Box Plots of the Estimated Pi values")
abline(h = 1 - fhat$pi0, col = "red") #This is the value of the pi for the no-covariate local FDR method. 














