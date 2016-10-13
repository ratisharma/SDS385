###########################################################
######### Created by Spencer Woody on 22 Sep 2016 #########
###########################################################

library(Matrix)
library(Rcpp)
library(RcppEigen)

# sourceCpp("fastsgd.cpp")

sourceCpp("fastsgdridge.cpp")

# Note that file "url_Xt.rds" is already in column-oriented mode
tX <- readRDS('url_Xt.rds')
y <- readRDS('url_y.rds')

n <- length(y)
p <- nrow(tX)

M <- rep(1, n)
beta0 <- rep(0, p)

#sgd1 <- sgdcppr(tX, y, M, stepsize0 = 2, numpasses = 1, beta0, lambda = 0.0, edecay = 0.001)

sgd2 <- sgdcppr(tX, y, M, stepsize0 = 2, numpasses = 1, beta0, lambda = 0.01, edecay = 0.001)

png("nlltrace1.png", height = 580, width = 580)
plot(sgd2$logliktrace, 
	type = "l",
	xlab = "iteration",
	ylab = "Negative Log Likelihood",
	main = "Traceplot of Negative Log Likelihood with ell-2 penalization, lamba = 0.01")
dev.off()


# Goodness of fit

dotprod <- crossprod(tX, sgd2$beta)
yhats <- 1 / (1 + exp(-dotprod))

dotprod2 <- X %*% sgd1$beta

sens <- sum((y == 1) * (yhats > 0.5)) / sum(y==1) # true positive
spec <- sum((y == 0) * (yhats < 0.5)) / sum(y==0) # true negative

# Cross validation
numbins <- 5
jumble <- sample(1:n, n, replace = F)
bin.indices <- split(jumble, cut(1:n, numbins))

# sens.vec <- rep(NA, numbins)
# spec.vec <- rep(NA, numbins)

# Create lists for holding bins of data
tX.list <- list()
y.list  <- list()
M.list  <- list()


for (k in 1:numbins) {
	indices.k <- bin.indices[[k]]
	
	tX.list[[k]] <- tX[, indices.k]
	y.list[[k]]  <- y[indices.k]
	M.list[[k]]  <- M[indices.k]
	print(k)
}

lambdas <- c(0, 0.000001, 0.0001, 0.01)
sens.mat <- matrix(rep(NA, numbins * length(lambdas)), nrow = length(lambdas))
spec.mat <- matrix(rep(NA, numbins * length(lambdas)), nrow = length(lambdas))
accu.mat <- matrix(rep(NA, numbins * length(lambdas)), nrow = length(lambdas))

for (i in 1:length(lambdas)) {
	lambda.i <- lambdas[i]
	for (j in 1:numbins) { 
		indices.j <- c(1:numbins)[-j]
	
		tX.tr <- do.call(cbind, tX.list[ indices.j ])
		y.tr  <- do.call(cbind,  y.list[ indices.j ])
		M.tr  <- do.call(cbind,  M.list[ indices.j ])
	
		sgd.j <- sgdcppr(tX.tr, y.tr, M.tr, 	
				 stepsize0 = 2, numpasses = 1, beta0, lambda = lambda.i, edecay = 0.001)
	
		beta.j <- sgd.j$beta 
	
		tX.te <- tX.list[[ j ]]
		y.te  <-  y.list[[ j ]]
		M.te  <-  M.list[[ j ]]
	
		tXbeta.j <- crossprod(tX.te, beta.j)
		yhats.j <- 1 / (1 + exp(-tXbeta.j  - sgd.j$alpha))
	
		sens.mat[i, j] <- sum((y.te == 1) * (yhats.j > 0.5)) / sum(y.te==1)
		spec.mat[i, j] <- sum((y.te == 0) * (yhats.j < 0.5)) / sum(y.te==0)
		accu.mat[i, j] <- sum((y.te == 1) * (yhats.j > 0.5) + (y.te == 0) * (yhats.j < 0.5)) / length(y.te)
		sprintf("%i out of %i bins cycled through", j, numbins)
	}
	sprintf("i = %i out of %i lambdas cross-validated.", i, length(lambdas))	
}

sens.vec <- apply(sens.mat, 1, mean)
spec.vec <- apply(spec.mat, 1, mean)
accu.vec <- apply(accu.mat, 1, mean)

print(lambdas)
print(sens.vec)
print(spec.vec)
print(accu.vec)


# > print(lambdas)
# [1] 0e+00 1e-06 1e-04 1e-02
# > print(sens.vec)
# [1] 0.9901595 0.9897778 0.9888078 0.9726841
# > print(spec.vec)
# [1] 0.9934549 0.9936025 0.9930379 0.9922197
# > print(accu.vec)
# [1] 0.9923656 0.9923385 0.9916390 0.9857616

