#carina = read.table('car_published.dat')
cold = read.table('carina.dat')

cold = cold[order(cold[,3]), ]
x = cold[cold[,8] + cold[,9] > 0, ]
x = x[x[,6] < 3,]

R = x[,1]
Th = x[,2] * pi/180

# v is y
# x1, x2 are 'coordinates'
x1 = R * cos(Th)
x2 = R * sin(Th)
v = x[,4]

# this is the histogram of f0
milkyway = read.table('carina.besancon')
milkyway = drop(as.matrix(milkyway))
hist(milkyway, breaks = 100)


# this is the preprocessing
# x[,10] contains EM probs from previous run

plot(x1, x2, cex = 0.3)
plot(R, v, pch = 19, cex = 0.3)

pdf('v_vs_R.pdf', width = 7, height = 7)
par(bty = 'l')
plot(R, v, pch = 19, cex = 0.3,
     xlab = 'Distance from Center',
     ylab = 'Velocity')
dev.off()

pdf('densities.pdf', width = 7, height = 7)
par(bty = 'l')
plot(density(v, adjust = 0.9), main = ' ', xlab = 'Line-of-sight Velocity', ylab = 'Estimated Density')
hmw = hist(milkyway, breaks = 100, plot = FALSE)
hmw$density = (1-mean(eo$p)) * hmw$density
plot(hmw, add = TRUE, freq = FALSE, col = rgb(0,0,1,alpha = 0.7), border = rgb(0,0,1,alpha = 0))
dev.off()




# x[,6] seems interesting
plot(x[,6], v, pch = 19, cex = 0.3)
# this is magnesium strength

# one of the components
mg1 = x[which(x[,10] >= 0.5),6]
hist(mg1, freq = FALSE, breaks = 30)
points(mg1, dnorm(mg1, mean(mg1), sd(mg1)), cex = 0.2)

mg2 = x[which(x[,10] < 0.5),6]
hist(mg2, freq = FALSE, breaks = 30)
points(mg2, dnorm(mg2, mean(mg2), sd(mg2)), cex = 0.2)

# running KW on magnesium
kwo = kw(x[,6], sd = 0.1)
hist(x[,6], breaks = 50, freq = FALSE, ylim = c(0, 2.3))
points(x[,6], kwo$f1y, col = 'blue', cex = 0.2)

plot(kwo$atoms, kwo$probs, type = 'h')

# running mclust
mo = Mclust(x[,6], G = 2, modelNames = 'V')
# looking at fitted density
mohat = with(mo$parameters,
             pro[1] * dnorm(x[,6], mean[1], sqrt(variance$sigmasq[1])) + 
               pro[2] * dnorm(x[,6], mean[2], sqrt(variance$sigmasq[2]))) 
hist(x[,6], breaks = 50, freq = FALSE, ylim = c(0, 2.3))
points(x[,6], mohat, col = 'red', cex = 0.2)

plot(x[,10], v, pch = 19, cex = 0.3)


hist(v, breaks = 50)

milkyway = read.table('carina.besancon')
hist(milkyway[,1], breaks = 100)

par(mfrow = c(1,2))
plot(density(fb[,1]))
plot(density(v))
par(mfrow = c(1,1))

# this is the old data
# we shall need to rerun this with full data
# v is velocity
# R is distance
plot(R, v, cex = 0.3)

hist(R, breaks = 50)
plot(density(R))

ub = max(R) * 0.9

nc1 = 2/(ub)^2


hist(R, breaks = 50, freq = FALSE)
rgrid = seq(from = 0, to = ub, length.out = 100)
f1 = 0.65 * rgrid * nc1
lamb = 1/10
f2 = 0.35 * lamb * exp(-lamb * rgrid)
lines(rgrid, f1+f2, col = 'blue')

hist(R, breaks = 50, freq = FALSE)
lines(rgrid, f1, col = 'blue')


# look at modelling slide in bodhi-da's talk


cnew = read.table('car_published.dat')
# this is line of sight velocity
vnew = cnew[,11]

hist(vnew, breaks = 50)

# gathering coordinates
# this was R
hist(x[,1], breaks = 50)
summary(x[,1])

hist(cnew[,1], breaks = 50)

for (i in 2:18) {print(i); print(summary(cnew[,i]))}
# 5,7,8
hist(cnew[,9], breaks = 50)



