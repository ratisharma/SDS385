######### Created by Spencer Woody on 23 Oct 2016 #########

library(MASS)
library(glmnet)
library(ggplot2)

source("myfuns.R")

# Read in data, scale everything
X <- as.matrix(read.csv("diabetesX.csv", header = TRUE))
Y <- as.numeric(unlist(read.csv("diabetesY.csv", header = FALSE)))

X <- scale(X)
Y <- scale(Y)

nll <- function(X, Y, beta, lambda) {
	f <- 0.5 * crossprod(Y - X %*% beta)
	g <- lambda * sum(abs(beta))
	
	mynll <- f + g
	
	return(as.numeric(mynll))
}

prox.ell1 <- function(u, lambda) {
	uprime <- abs(u) - lambda
	prox <- sign(u) * uprime * (uprime > 0)
	
	return(prox)
}

admmLasso <- function(X, Y, lambda = 0.01, rho = 5, maxiter = 1e4, tol = 1e-10, beta0 = NA) {
	N <- nrow(X)
	P <- ncol(X)
	
	lambda <- lambda * N
	
	if (is.na(beta0)) {
		beta0 <- rep(0, P)
	}
	
	XtX <- crossprod(X)
	XtY <- crossprod(X, Y)
	
	betas <- matrix(rep(NA, (maxiter + 1) * P), nrow = P)
	betas[, 1] <- beta0

	nlltrace <- rep(NA, maxiter + 1)
	nlltrace[1] <- t(beta0) %*% (0.5 * XtX %*% beta0 - XtY) + lambda * sum(abs(beta0))

	invmat <- solve(crossprod(X) + diag(rho, P))
	
	z <- rep(0, P)
	v <- rep(0, P)

	message <- "Convergence not reached..."

	for (i in 1:maxiter) {
		#betas[, i + 1] <- invmat %*% (crossprod(X, Y) + rho * (z - v))
		x <- invmat %*% (crossprod(X, Y) + rho * (z - v))
	
		z <- prox.ell1(x + v, lambda / rho)
	
		betas[, i + 1] <- z
	
		v <- v + x - z
	
		nlltrace[i + 1] <- t(z) %*% (0.5 * XtX %*% z - XtY) + lambda * sum(abs(z))
	
		nllchange <- nlltrace[i + 1] / nlltrace[i]  - 1
	
		if (abs(nllchange) < tol & i > 2) {
			betas <- betas[, -((i + 2):(maxiter+1))]
			nlltrace <- nlltrace[-((i + 2):(maxiter+1))]
			message <- sprintf("Convergence reached after %i iterations", i)
			break
		}
	}
	

	mylist <- list("finalbeta" = betas[, ncol(betas)], "beta" = betas, "nll" = nlltrace, "conv" = message)
	return(mylist)
}

fit1 <- glmnet(X, Y)
lambdas <- fit1$lambda

N = nrow(X)

system.time(myadmm <- admmLasso(X, Y, lambda = 0.00613679 * N, rho = 10, maxiter = 1e4))

which(myadmm$finalbeta == 0)
myadmm$conv
plot(myadmm$nll, type = "l", log = "x")

which(fit1$beta[, 50] == 0)



betamat <- matrix(rep(NA, length(lambdas)*ncol(X)), ncol = length(lambdas))

for (i in 1:length(lambdas)) {
	mylasso <- admmLasso(X, Y, lambda = lambdas[i] * N, rho = 5, maxiter = 1e4)
	betamat[, i]  <- mylasso$finalbeta
}

ybound <- c(min(betamat), max(betamat))

#pdf("myplot.pdf")
plot(log(lambdas), betamat[1, ], 
	ylim = ybound, 
	#log = "x",
	main = "Coefficients of Lasso Regression using my own lasso function",
	xlab = "",
	ylab = "Coefficient",
	type = "l", 
	col = 2)
for (j in 2:nrow(betamat)) {
	lines(log(lambdas), betamat[j, ], col = (j + 4))
}
#dev.off()






plot(fit1, xvar = "lambda")

nlltrace <- myadmm$nll
betas<-myadmm$beta

dif  <- rep(NA, length(nlltrace - 1))
bdif <- rep(NA, length(nlltrace - 1))

for (j in 2:length(nlltrace)) {
	dif[j-1] <- nlltrace[j] / nlltrace[j - 1] - 1.0
	bdif[j-1] <- sum(abs(betas[, j] - betas[, j - 1]))
}

plot(abs(dif), type = "l", log = "xy")
plot(bdif, type = "l", log = "xy")


