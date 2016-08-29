###########################################################
######### Created by Spencer Woody on 24 Aug 2016 #########
###########################################################

library(Matrix)
library(microbenchmark)

### No. 1 pt C

# Set N, P, X, W, and y

N <- 2000
P <-  500

X <- matrix(rnorm(N * P), nrow = N)
y <- matrix(rnorm(N), nrow = N)

# Inversion method

Inv.method <- function(X.Inv, W.Inv, y.Inv) {
	bhat.Inv <- solve(t(X.Inv) %*% W.Inv %*% X.Inv) %*% t(X.Inv) %*% W.Inv %*% y.Inv
	
	return(bhat.Inv)
}

LU.decomp <- function(X.LU, W.LU, y.LU) {
	D.LU <- t(X.LU) %*% W.LU %*% y.LU
	C.LU <- t(X.LU) %*% W.LU %*% X.LU

	LUdecomp <- lu(C.LU)
	L.LU <- expand(LUdecomp)$L
	U.LU <- expand(LUdecomp)$U

	A.LU <- solve(L.LU, D.LU)
	bhat.LU <- solve(U.LU, A.LU)
	
	return(bhat.LU)
}

microbenchmark(
	inversion_solver(y,X),
	chol_solver(y,X),
	times=5)

microbenchmark(t(X)%*%W%*%y, crossprod(X, W, y), times = 5)

microbenchmark(
mymat1 <- t(X) %*% W %*% y, 
mymat2 <- crossprod(X, W) %*% y, 
times = 2)


microbenchmark(
t(X) %*% W %*% y, 
crossprod(X, W) %*% y, 
times = 2)

Cho.decomp <- function(X.Cho, W.Cho, y.Cho) {
	D.Cho <- t(X.Cho) %*% W.Cho %*% y.Cho
	C.Cho <- t(X.Cho) %*% W.Cho %*% X.Cho
	
	Cho <- chol(C.Cho)
	
#	L.Cho <- t(Cho)
#	U.Cho <- Cho

#	A.Cho <- solve(L.Cho, D.Cho)
#	bhat.Cho <- solve(U.Cho, A.Cho)
	bhat.Cho <- forwardsolve(t(Cho), D.Cho)
	return(bhat.Cho)
}



microbenchmark(
bhat.Inv <- Inv.method(X, W, y), 
bhat.LU  <- LU.decomp(X, W, y), 
bhat.Cho <- Cho.decomp(X, W, y), 
times = 2
)

microbenchmark(
bhat.Inv <- Inv.method(X, W, y), 
bhat.Cho <- Cho.decomp(X, W, y), 
times = 2
)

### No. 1 pt D

N <- 2000
P <-  500

X <- matrix(rnorm(N * P), nrow = N)
mask <- matrix(rbinom(N * P, 1, 0.05), nrow = N)
X <- mask * X

# END

myfun <- function(x) {
	dim(combn(5,x))[2] * dim(combn(25, 4-x))[2] / dim(combn(30,4))[2]
}

