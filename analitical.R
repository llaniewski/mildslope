
library(Bessel)

k = 1
cn = 256
theta = seq(0,2*pi,len=cn+1)[-(cn+1)]

seq_circ = function(n) 0:(n-1)-ifelse(1:n>n/2,n,0)
i = seq_circ(cn)


# wave = exp(1i*k*x)
dwave_dn = 1i*k*exp(1i*k*cos(theta))*cos(theta)

plot_cpx = function(x,z,...) matplot(x,cbind(Re(z),Im(z)),pch=c("R","I"),...)
plot_cpx(theta, dwave_dn)

dwave_dn_coef = fft(dwave_dn)/length(dwave_dn)
dwave_dn_coef[Mod(dwave_dn_coef) < 1e-13] = 0
plot_cpx(i, dwave_dn_coef)
plot(i, Mod(dwave_dn_coef),log="y")

N=300
x = seq(-10,10,len=N)
y = seq(-4,4,len=N)
p = expand.grid(x=x,y=y)
p$r = sqrt(p$x*p$x+p$y*p$y)
p$theta = atan2(p$y,p$x)

image_cpx = function(f,zlim=c(-1,1)*max(abs(Re(f)),abs(Im(f)),na.rm=TRUE)) {
    par(mfrow=c(2,1))
    image(x,y,matrix(Re(f),N,N),zlim=zlim,asp=1)
    image(x,y,matrix(Im(f),N,N),zlim=zlim,asp=1)
}

n = 3
f = exp(p$theta*1i*n) * BesselY(k*p$r,n)
f[p$r<1] = NA
image_cpx(f)

selected_bessel = BesselY


dJ_dR = sapply(i,function(n) k*(selected_bessel(k,n-1) - selected_bessel(k,n+1))/2)

plot(i,Mod(dJ_dR),log="y",ylim=c(1e-13,1e6),xlim=c(-25,25))
points(i,Mod(dwave_dn_coef),pch=16,col=2)

points(i,Mod(dwave_dn_coef/dJ_dR),pch=16,col=3)
range(Mod(dwave_dn_coef))

plot_cpx(i,dwave_dn_coef / dJ_dR)

coef = dwave_dn_coef / dJ_dR

sel = coef != 0
bigmat = sapply(i[sel], function(n) exp(p$theta*1i*n) * selected_bessel(k*p$r,n))

f = bigmat %*% coef[sel]
f[p$r<1] = NA
image_cpx(f)
range(Re(f),na.rm=TRUE)
range(Im(f),na.rm=TRUE)

range(Im(exp(1i*k*p$x) - f),na.rm=TRUE)


g = exp(1i*k*p$x) - f
image_cpx(g)
range(Re(g),na.rm=TRUE)
range(Im(g),na.rm=TRUE)

