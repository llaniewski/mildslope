write.poly = function(filename, points, segments,limits = c(-0.2,1)) {
    f = file(filename,open="w")
    p = cbind(points$x,points$y,ifelse(points$par,limits[2],0),ifelse(points$par,limits[1],0))
    p = cbind(1:nrow(p),p)
    write.table(t(c(nrow(p),2,ncol(p)-3,0)), file=f,row.names=FALSE,col.names=FALSE)
    write.table(p, file=f,row.names=FALSE,col.names=FALSE)
    s = cbind(segments$i1,segments$i2,segments$tag)
    s = cbind(1:nrow(s),s)
    write.table(t(c(nrow(s),1)), file=f,row.names=FALSE,col.names=FALSE)
    write.table(s, file=f,row.names=FALSE,col.names=FALSE)
    h = data.frame()
    write.table(t(c(nrow(h))), file=f,row.names=FALSE,col.names=FALSE)
    close(f)
}


r = 0.6/2
npar = 60
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
    tag = c(2,rep(1,npar+3),3,1,rep(1,ncirc-1),1)
)

plot(points$x,points$y,asp=1)
segments(
    points$x[segments$i1],points$y[segments$i1],
    points$x[segments$i2],points$y[segments$i2],
    col=segments$tag
)

write.poly("circ_20.poly",points, segments)


k = 0

npar_half = floor(npar/2)
a = seq(0,pi,len=ncirc)
points = data.frame(
    x=c(-5,-5,seq(-0.5-xpar,-0.5,len=npar_half+1),-0.5,0.5,seq(0.5,0.5+xpar,len=npar_half+1),5,5),
    y=c(0,rep(0.5,npar_half+2),5,5,rep(0.5,npar_half+2),0),
    par=c(rep(FALSE,3),rep(TRUE,npar_half),rep(FALSE,2),rep(TRUE,npar_half),rep(FALSE,3))
)

n = nrow(points)
segments = data.frame(
    i1 = 1:n,
    i2 = c(2:n,1),
    tag = c(1,rep(3,npar_half+2),3,rep(3,npar_half+2),2,4)
)

plot(points$x,points$y,asp=1)
segments(
    points$x[segments$i1],points$y[segments$i1],
    points$x[segments$i2],points$y[segments$i2],
    col=segments$tag
)

write.poly("side_channel.poly",points, segments)



k = 0
ncirc=5
a = seq(0,pi,len=ncirc)
eps = 0.02
d = xpar
points = data.frame(
    x=c(
        -5,seq(-xpar,xpar,len=npar+2),5,
        5,d-eps*sin(a),5,
        5,seq(xpar,-xpar,len=npar+2),-5
    ),
    y=c(
        rep(0.5,npar+4),
        eps,eps*cos(a),-eps,
        rep(-0.5,npar+4)
    ),
    par=c(
        rep(FALSE,2),rep(TRUE,npar),rep(FALSE,2),
        rep(FALSE,2+ncirc),
        rep(FALSE,2),rep(TRUE,npar),rep(FALSE,2)
    )
)

n = nrow(points)
segments = data.frame(
    i1 = 1:n,
    i2 = c(2:n,1),
    tag = c(
        rep(1,npar+3),3,rep(1,ncirc+1),4,rep(1,npar+3),2)
)

plot(points$x,points$y,asp=1,pch=16,cex=ifelse(points$par,1,0.3))
segments(
    points$x[segments$i1],points$y[segments$i1],
    points$x[segments$i2],points$y[segments$i2],
    col=segments$tag
)

write.poly("split.poly",points, segments, limits=c(-0.1,0.1))


k = 0

a = seq(0,pi,len=ncirc)
x = seq(-xpar,xpar,len=npar+2)
points = data.frame(
    x=c(-5.0,-5.0,ifelse(x<0,x,0)+0.5,+0.5,-0.5,-0.5),
    y=c(-0.5, 0.5,-ifelse(x<0,0,x)+0.5,-5,-5,-0.5),
    par=c(rep(FALSE,3),rep(TRUE,npar),rep(FALSE,4))
)

n = nrow(points)
segments = data.frame(
    i1 = 1:n,
    i2 = c(2:n,1),
    tag = c(1,rep(3,npar+3),2,rep(3,2))
)

a = -pi/4
points[,1:2] = as.matrix(points[,1:2]) %*% matrix(c(cos(a),sin(a),-sin(a),cos(a)),2,2)

plot(points$x,points$y,asp=1)
segments(
    points$x[segments$i1],points$y[segments$i1],
    points$x[segments$i2],points$y[segments$i2],
    col=segments$tag
)

write.poly("turn.poly",points, segments)
