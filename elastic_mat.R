library(polyAlgebra)

midx = c("00","10","01","11")

X0 = PV("P0[i0*2+",1:2-1,"]")
X1 = PV("P0[i1*2+",1:2-1,"]")
X2 = PV("P0[i2*2+",1:2-1,"]")
x0 = PV("P1[i0*2+",1:2-1,"]")
x1 = PV("P1[i1*2+",1:2-1,"]")
x2 = PV("P1[i2*2+",1:2-1,"]")

zeta_X_val = c(X1-X0,X2-X0);
zeta_x_val = c(x1-x0,x2-x0);
zeta_X = PV("zX",midx); dim(zeta_X) = c(2,2)
zeta_x = PV("zx",midx); dim(zeta_x) = c(2,2)
C(zeta_X, zeta_X_val)
C(zeta_x, zeta_x_val)
# zeta_X %*% zeta = X
# zeta_x %*% solve(zeta_X) %*% X = x
det_val = zeta_X[1]*zeta_X[4] - zeta_X[2]*zeta_X[3]
det = PV("det")
C(det,det_val)

zeta_X_inv = zeta_X[c(4,2,3,1)] * c(1,-1,-1,1) * det^{-1}
dim(zeta_X_inv) = c(2,2)

X_x_val = zeta_x %*% zeta_X_inv
X_x = PV("Xx",midx); dim(X_x) = c(2,2)
C(X_x,X_x_val)

# X_x %*% X = x
# t(x) %*% x = t(X_x %*% X) %*% X_x %*% X = t(X) %*% t(X_x) %*% X_x %*% X
S_X_X_val = t(X_x) %*% X_x
S_X_X = PV("SXX",midx); dim(S_X_X) = c(2,2)
C(S_X_X,S_X_X_val)

F = PV("F",midx); dim(F) = c(2,2)
C(F,S_X_X_val - PV(c(1,0,0,1)))

