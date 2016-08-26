###########################################################
######### Created by Spencer Woody on 23 Aug 2016 #########
###########################################################

library(Matrix)
library(microbenchmark)

### No. 1 pt B

# Set N, P, X, W, and y

N <- 2000
P <-  500

X <- matrix(rnorm(N * P), nrow = N)
W <- diag(rep(1, N))
y <- matrix(rnorm(N), nrow = N)

# Inversion method

bhat.Inv <- solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% y


# LU decomp method

D <- t(X) %*% W %*% y
C <- t(X) %*% W %*% X

LUdecomp <- lu(C)
L <- expand(LUdecomp)$L
U <- expand(LUdecomp)$U

A <- solve(L, D)
bhat.LU <- solve(U, A)

# Cholesky decomp method

Cho <- chol(C)

L.Cho <- t(Cho)
U.Cho <- Cho

A.Cho <- solve(L.Cho, D)
bhat.Cho <- solve(U.Cho, A.Cho)



### No. 1 pt C

### No. 1 pt D

# N <- 2000
# P <-  500
#
# X <- matrix(rnorm(N * P), nrow = N)
# mask <- matrix(rbinom(N * P, 1, 0.05), nrow = N)
# X <- mask * X

# END