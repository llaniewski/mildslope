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



npar = 60
xpar = 2
for (d in seq(0.1,0.9,0.1)) {
    r = d/2
    ncirc = ceiling(r*pi/0.03)+1
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

    write.poly(sprintf("circ_%02.0f.poly",d*100),points, segments)
}



make_wall = function(w,h) {
    npar = 30
    xpar = 2
    k = 0
    a = seq(0,pi,len=ncirc)
    x = seq(-1,1,len=npar+2)
    y = (cos(x*pi)+1)/2
    points = data.frame(
        x=c(-5, -5,    xpar*x,  5, 5,w/2,w/2,-w/2,-w/2),
        y=c( 0,0.5,0.5+h/4*y,0.5, 0,  0,h/2, h/2,   0),
        par=c(rep(FALSE,3),rep(TRUE,npar),rep(FALSE,7))
    )

    n = nrow(points)
    segments = data.frame(
        i1 = 1:n,
        i2 = c(2:n,1),
        tag = c(2,rep(1,npar+3),3,rep(1,5))
    )

    plot(points$x,points$y,asp=1)
    segments(
        points$x[segments$i1],points$y[segments$i1],
        points$x[segments$i2],points$y[segments$i2],
        col=segments$tag
    )

    write.poly(sprintf("wall_%03.0f_%03.0f.poly",w*100,h*100), points, segments)
}

for (h in seq(0.1, 1,0.1)) make_wall(0.01,h)
for (h in seq(0.1, 1,0.1)) make_wall(0.6,h)
for (h in seq(0.1, 1,0.1)) make_wall(1,h)
for (w in seq(0.1, 1,0.1)) make_wall(w,0.6)

npar = 100
xpar = 2
eps = 0.005
ncirc = ceiling(eps*pi/0.03)+1
for (d in seq(0.1,0.9,0.1)) {
    r = d/2
    
    k = 0
    a = seq(0,pi,len=ncirc)
    x = seq(-xpar,xpar,len=npar+2)
    alpha = 8.7
    y = 0.5 + r*sin(x*alpha)/(x*alpha)
    points = data.frame(
        x=c(-5,-5,x,5,5,eps,eps*cos(a),-eps),
        y=c(0,0.5,y,0.5,0,0,r-max(eps*sin(a))+eps*sin(a),0),
        par=c(rep(FALSE,3),rep(TRUE,npar),rep(FALSE,3+ncirc+2))
    )

    n = nrow(points)
    segments = data.frame(
        i1 = 1:n,
        i2 = c(2:n,1),
        tag = c(2,rep(1,npar+3),3,1,rep(1,ncirc-1+2),1)
    )

    plot(points$x,points$y,asp=1)
    segments(
        points$x[segments$i1],points$y[segments$i1],
        points$x[segments$i2],points$y[segments$i2],
        col=segments$tag
    )

    write.poly(sprintf("sinc_%02.0f.poly",d*100),points, segments)
}


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

write.poly("split.poly",points, segments, limits=c(-0.03,0.03))


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
    tag = c(2,rep(1,npar+3),3,rep(1,2))
)

a = -pi/4
points[,1:2] = as.matrix(points[,1:2]) %*% matrix(c(cos(a),sin(a),-sin(a),cos(a)),2,2)

plot(points$x,points$y,asp=1)
segments(
    points$x[segments$i1],points$y[segments$i1],
    points$x[segments$i2],points$y[segments$i2],
    col=segments$tag
)

write.poly("turn.poly",points, segments,limits=c(-1,1)*0.3)
