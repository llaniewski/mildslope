r = 0.6/2
npar = 20
xpar = 2
ncirc = 30

k = 0

a = seq(0,pi,len=ncirc)
points = data.frame(
    x=c(-5,-5,seq(-xpar,xpar,len=npar+2),5,5,r*cos(a)),
    y=c(0,rep(0.5,npar+4),0,r*sin(a)),
    par=c(rep(FALSE,3),rep(TRUE,npar),rep(FALSE,3+ncirc))
)


n = nrow(points)
segments = data.frame(
    i1 = 1:n,
    i2 = c(2:n,1),
    tag = c(1,rep(3,npar+3),2,3,rep(4,ncirc-1),3)
)

plot(points$x,points$y,asp=1)
segments(
    points$x[segments$i1],points$y[segments$i1],
    points$x[segments$i2],points$y[segments$i2],
    col=segments$tag
)

f = file("circ_20.poly",open="w")
np = sum(points$par)
p = cbind(points$x,points$y,matrix(0,nrow(points),np))
p[points$par,-(1:2)] = diag(np)
p = cbind(1:nrow(p),p)
write.table(t(c(nrow(p),2,np,0)), file=f,row.names=FALSE,col.names=FALSE)
write.table(p, file=f,row.names=FALSE,col.names=FALSE)
s = cbind(segments$i1,segments$i2,segments$tag)
s = cbind(1:nrow(s),s)
write.table(t(c(nrow(s),2,1)), file=f,row.names=FALSE,col.names=FALSE)
write.table(s, file=f,row.names=FALSE,col.names=FALSE)
h = data.frame()
write.table(t(c(nrow(h))), file=f,row.names=FALSE,col.names=FALSE)
close(f)
