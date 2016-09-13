###########################################################
######### Created by Spencer Woody on 13 Sep 2016 #########
###########################################################

library(TTR)

# Read in data file, scale X (y = 1 represents a malignant tumor)

wdbc <- read.csv("wdbc.csv", header = FALSE)

X <- as.matrix(wdbc[, 3:12])
X <- scale(X)
X <- cbind(rep(1, nrow(X)), X)

y <- wdbc[, 2]
y <- y == "M"
m.i <- 1

# # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # #

# Function for computing w.i (logit transform of Xtbeta)

comp.wi <- function (X, beta) {
	wi <- 1 / (1 + exp(-X %*% beta))
	return(wi)
}

# Function for computing likelihood

loglik <- function(beta, y, X, m.i) {
	loglik <- apply((m.i - y) * (X %*% beta)+ m.i*log(1 + exp(-X %*% beta)), 2, sum)
	return(loglik)
}

# Function for computing gradient of likelihood

grad.loglik <- function(beta, y, X, mi){
  grad <- array(NA, dim = length(beta))
  wi   <- comp.wi(X, beta)
  grad <- apply(X*as.numeric(mi * wi - y), 2, sum)
  return(grad)
}

# # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # #

### NEWTON'S METHOD (results are referred to as beta.N)

beta.N <- as.matrix(rep(0, ncol(X)))

newton.steps <- 10

for (step in 1:newton.steps) {
	Newton.wi <- as.numeric(comp.wi(X, beta.N))
	Hessian   <- t(X) %*% diag(Newton.wi*(1-Newton.wi)) %*% X
	beta.N    <- beta.N - solve(Hessian, grad.loglik(beta.N, y, X, m.i))
}

# Function for making traceplots of all betas

graph.betatrace <- function(beta, beta.N, graphname) {
	num.vars <- nrow(beta)
	n.graphrows <- floor(sqrt(num.vars))
	n.graphcols <- ceiling(num.vars / n.graphrows)
	pdf(graphname, width = n.graphcols * 2, height = n.graphrows * 2.5)
	par(mfrow = c(n.graphrows, n.graphcols), oma=c(0,0,2,0))
	for (j in 1:num.vars) {
		plot(beta[j, ], 
			xlab = "iteration",
			ylab = paste("beta", sprintf("%i", j)),
			type = "l",
			ylim = c(min(beta.N[j], min(beta[j, ])) - 0.5,
			max(beta.N[j], max(beta[j, ])) + 0.5),
			log = "x"
			)
		abline(h = beta.N[j], col = "red")
	}
	title("Trace plot for all betas", outer = TRUE)
	dev.off()
}

# Stochastic gradient

SGD <- function(n.iterSGD,  beta.init, step.size, X, y, m.i) {
	#' Perform stochastic gradient descent for a binomial logistic regression
	#' 
	#' @param n.iterSGD  Number of iterations to loop through
	#' @param beta.init  Initial guess for coefficients (betas)
	#' @param step.size  Constant step size
	#' @param X  N by P matrix of covariate data
	#' @param y  P-vector of responses
	#' @param m.i  n-parameter of binomial (1 for case of binary logistic)
	#' 
	#' @return A matrix, each column is an iteration of computed betas.
	#
	#
		# Create zero matrix to store iterative values of beta; initialize beta
		betaSGD <- matrix(rep(0, ncol(X) * (n.iterSGD+ 1)), nrow = ncol(X))
		betaSGD[, 1] <- beta.init
		#
		for (step in 1:n.iterSGD) {
			
			# Draw random sample of single row of data (with replacement) 
			i <- sample(nrow(X), 1)
			
			# Compute w.i, fitted value of p-parameter of binomial
			w.i <- 1 / (1 + exp(-crossprod(X[i, ], betaSGD[, step])))
			
			# Compute gradient, the descent direction
			grad <- nrow(X) * (m.i * w.i - y[i]) * X[i, ]
			
			# Next set of betas
			betaSGD[, step + 1] <- betaSGD[, step] - step.size * grad
	}
		return(betaSGD)
}

# Create initial values of beta which are some distance away from "true" values
# of beta computed from Newton's method. At least 5 units away plus some 
# exponential noise, in either positive or negative direction
beta1 <- beta.N + (-1) ^ rbinom(ncol(X), 1, 0.5) * (5 + rexp(ncol(X), rate = 1))

# Set number of iterations and stepsize
my.iter <- 4e5
my.stepsize <- 0.0003

# Perform SGD, graph trace plot of betas
sgd1 <- SGD(my.iter, beta1, my.stepsize, X, y, m.i)
graph.betatrace(sgd1, beta.N, "SGD.pdf")

# Compute loglikelihood from results of SGD
trace1 <- loglik(sgd1, y, X, m.i)

# Set time period for moving average
my.time <- 500

# Create plot of EMA of log-likelihood
pdf("EMAloglik.pdf")
plot(EMA(trace1, n = my.time), 
	type = "l", 
	main = "Exponential moving average of log-likelihood",
	xlab = "iteration",
	ylab = "Log-likelihood",
	xlim = c(my.time, length(trace1)),
	log = "x")
dev.off()

# Decaying steps


SGD.decay <- function(n.iterSGD,  beta.init, C, t.0, alpha, X, y, m.i) {
	#' Perform stochastic gradient descent for a binomial logistic regression
	#' 
	#' @param n.iterSGD  Number of iterations to loop through
	#' @param beta.init  Initial guess for coefficients 
	#' @param C  Constant C in Robbins-Monro rule
	#' @param t.0  Constant t.0 in Robbins-Monro rule
	#' @param alpha Constant alpha in Robbins-Monro rule
	#' @param X  N by P matrix of covariate data
	#' @param y  P-vector of responses
	#' @param m.i  n-parameter of binomial (1 for case of binary logistic)
	#' 
	#' @return A matrix, each column is an iteration of computed betas.
	#
	#
		# Create NULL matrix to store iterative values of beta; initialize beta
		betaSGD.decay <- matrix(rep(0, ncol(X) * (n.iterSGD+ 1)), nrow = ncol(X))
		betaSGD.decay[, 1] <- beta.init
		#
		for (step.decay in 1:n.iterSGD) {
			
			# Draw random sample of single row of data (with replacement) 
			i <- sample(nrow(X), 1)
			
			# Compute w.i, fitted value of p-parameter of binomial
			w.i <- 1 / (1 + exp(-crossprod(X[i, ], betaSGD.decay[, step.decay])))
			
			# Compute gradient, the descent direction
			grad <- nrow(X) * (m.i * w.i - y[i]) * X[i, ]
			
			# Compute next step size
			stepsize.decay <- C * (step.decay + t.0) ^ -alpha
			
			# Next set of betas
			betaSGD.decay[, step.decay + 1] <- betaSGD.decay[, step] - stepsize.decay * grad
	}
		return(betaSGD.decay)
}

beta1 <- beta.N + (-1) ^ rbinom(ncol(X), 1, 0.5) * (3 + rexp(ncol(X), rate = 1))

decay1 <- SGD.decay(1e5, beta1, C = 10, t.0 = 1, alpha = 0.75, X, y, m.i)
graph.betatrace(decay1[, 10000:1e5], beta.N, "decay1.pdf")

# Running average of beta with decaying steps

burn.in <- 3e5
beta.burn <- sgd1[, burn.in:4e5]
beta.burnbar <- apply(beta.burn, 1, mean)

while (TRUE) {
	if (A & A | j > 1000) {
		break
	}
}
