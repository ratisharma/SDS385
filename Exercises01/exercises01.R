###########################################################
######### Created by Spencer Woody on 24 Aug 2016 #########
###########################################################

library(Matrix)
library(microbenchmark)

### No. 1 pt C

# Set N, P, X, W, and y

N <- 4000
P <- 1000

X <- matrix(rnorm(N * P), nrow = N)
y <- matrix(rnorm(N), nrow = N)
W <- diag(rep(1, N))

# Inversion method

Inv.method <- function(X.Inv, W.Inv, y.Inv) {
	XtWX <- (t(X.Inv)*diag(W.Inv)) %*% X.Inv
	XtWY <- (t(X.Inv)*diag(W.Inv)) %*% y.Inv
	bhat.Inv <- solve(XtWX) %*% XtWY
	return(bhat.Inv)
}

Cho.decomp <- function(X.Cho, W.Cho, y.Cho) {
	D.Cho <- (t(X.Cho)*diag(W.Cho)) %*% y.Cho
	C.Cho <- (t(X.Cho)*diag(W.Cho)) %*% X.Cho
	
	U.Cho <- chol(C.Cho)
	L.Cho <- t(U.Cho)
	
	u <- forwardsolve(L.Cho, D.Cho)
	bhat.Cho <- backsolve(U.Cho, u)
	
	return(bhat.Cho)
}

microbenchmark(
	Inv.method(X, W, y),
	Cho.decomp(X, W, y),
	times = 5, unit = "ms") # N = 4000, P = 1000


### No. 1 pt D

N <- 2000
P <- 1000

# Sparsity measure # 0.01, 0.05, 0.25
alpha <- 0.25

X <- matrix(rnorm(N * P), nrow = N)
mask <- matrix(rbinom(N * P, 1, alpha), nrow = N)
X <- mask * X
W <- diag(rep(1, N))


Cho.decompSPARSE <- function(X.Cho, W.Cho, y.Cho) {
	X.Cho <- Matrix(X.Cho, sparse = T)
	D.Cho <- (t(X.Cho)*diag(W.Cho)) %*% y.Cho
	C.Cho <- (t(X.Cho)*diag(W.Cho)) %*% X.Cho
	
	U.Cho <- chol(C.Cho)
	L.Cho <- t(U.Cho)
	
	u <- forwardsolve(L.Cho, D.Cho)
	bhat.Cho <- backsolve(U.Cho, u)
	
	return(bhat.Cho)
}



microbenchmark(
	Inv.method(X, W, y),
	Cho.decomp(X, W, y),
	Cho.decompSPARSE(X, W, y), 
	times=5, unit = "ms") # N = 2000, P = 1000, alpha = 0.25

# END
