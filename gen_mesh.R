#install.packages("devtools")
#devtools::install_github("davidcsterratt/RTriangle", subdir="pkg")

k = 0
# points = matrix(c(
#     0, 0, 2, 2, 1.5, 1.5, 1, 1,
#     0, 1, 1, 0,   0, 0.5, 0.5,   0),ncol=2)
points = matrix(c(
    0, 0, 10, 10,
    0, 1,  1,  0),ncol=2)
n = nrow(points)
segments = matrix(c(1:n,2:n,1),ncol=2)
k = k + n


p = RTriangle::pslg(points, S = segments)

ret = RTriangle::triangulate(p,a=0.01,q=30)
write.table(ret$P, "mesh/empty1_points.txt", row.names=FALSE,col.names=FALSE)
write.table(ret$T-1, "mesh/empty1_triangles.txt", row.names=FALSE,col.names=FALSE)
ret = RTriangle::triangulate(p,a=0.001,q=30)
write.table(ret$P, "mesh/empty2_points.txt", row.names=FALSE,col.names=FALSE)
write.table(ret$T-1, "mesh/empty2_triangles.txt", row.names=FALSE,col.names=FALSE)



n = 20
a = seq(0,2*pi,len=n+1)[-1]
points = rbind(points, matrix(c(cos(a)*0.3+5,sin(a)*0.3+0.5),ncol=2))
segments = rbind(segments, matrix(c(1:n,2:n,1),ncol=2)+k)
k = k + n

p = RTriangle::pslg(points,S = segments, H=matrix(c(5,0.5),ncol=2))

ret = RTriangle::triangulate(p,a=0.01,q=30)
write.table(ret$P, "mesh/mesh1_points.txt", row.names=FALSE,col.names=FALSE)
write.table(ret$T-1, "mesh/mesh1_triangles.txt", row.names=FALSE,col.names=FALSE)
ret = RTriangle::triangulate(p,a=0.001,q=30)
write.table(ret$P, "mesh/mesh2_points.txt", row.names=FALSE,col.names=FALSE)
write.table(ret$T-1, "mesh/mesh2_triangles.txt", row.names=FALSE,col.names=FALSE)

