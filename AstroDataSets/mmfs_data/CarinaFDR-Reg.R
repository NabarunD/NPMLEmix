
data = read.table('carina.dat')


x = data[data$V8 + data$V9> 0,]
# x = x[x$V6 < 3,]
rv = x$V4; # represents the Radial velocity data of stars in the  Carina galaxy

rr = x$V1;
th = x$V2;
mg = x$V6;
x1 = rr*cos(th);  x2 = rr*sin(th);  
r = sqrt(x1^2 +x2^2);
plot(x1,x2)
plot(r,rv)
plot(mg, rv)

#plot(r,rv)



dens = density(rv);
plot(dens)

