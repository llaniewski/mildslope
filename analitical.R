
library(Bessel)

k = 5
cn = 256
theta = seq(0,2*pi,len=cn+1)[-(cn+1)]

# wave = exp(1i*k*x)
dwave_dn = 1i*k*exp(1i*k*cos(theta))*cos(theta)

plot_cpx = function(x,z,...) matplot(x,cbind(Re(z),Im(z)),pch=c("R","I"),...)
plot_cpx(theta, dwave_dn)

dwave_dn_coef = fft(dwave_dn)/length(dwave_dn)
dwave_dn_coef[Mod(dwave_dn_coef) < 1e-13] = 0
plot_cpx(theta, dwave_dn_coef)

N=300
x = seq(-4,4,len=N)
p = expand.grid(x=x,y=x)
p$r = sqrt(p$x*p$x+p$y*p$y)
p$theta = atan2(p$y,p$x)

image_cpx = function(f,zlim=c(-1,1)*max(abs(Re(f)),abs(Im(f)),na.rm=TRUE)) {
    par(mfrow=c(2,1))
    image(x,x,matrix(Re(f),N,N),zlim=zlim,asp=1)
    image(x,x,matrix(Im(f),N,N),zlim=zlim,asp=1)
}

n = 1
f = exp(p$theta*1i*n) * BesselJ(k*p$r,n)
image_cpx(f)

seq_circ = function(n) 0:(n-1)-ifelse(1:n>n/2,n,0)

i = seq_circ(cn)
dJ_dR = sapply(i,function(n) (BesselY(k,n-1)) - BesselY(k,n+1))/2

plot_cpx(i,dwave_dn_coef / dJ_dR)

coef = dwave_dn_coef / dJ_dR

sel = coef != 0
bigmat = sapply(i[sel], function(n) exp(p$theta*1i*n) * BesselY(k*p$r,n))

f = bigmat %*% coef[sel]
f[p$r<1] = NA
image_cpx(f)
range(Re(f),na.rm=TRUE)
range(Im(f),na.rm=TRUE)

range(Re(exp(1i*k*p$y) - f/2))
image_cpx(exp(1i*k*p$y) - f/2)

