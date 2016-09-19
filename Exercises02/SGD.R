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

# Function for computing full likelihood

loglik <- function(beta, y, X, m.i) {
	loglik <- apply((m.i - y) * (X %*% beta) + m.i * log(1 + exp(-X %*% beta)), 2, sum)
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

# # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # #

# Stochastic gradient

SGD <- function(X, y, m.i, beta.init, 
				n.iterSGD, step.size, a, burn.in) {
	#' Perform stochastic gradient descent for a binomial logistic regression
	#' 
	#' @param X  N by P matrix of covariate data
	#' @param y  P-vector of responses
	#' @param m.i  n-parameter of binomial (1 for case of binary logistic)
	#' @param beta.init  Initial guess for coefficients (betas)
	#' @param n.iterSGD  Number of iterations to loop through
	#' @param step.size  Constant step size
	#' @param a  Coefficient for weighting decrease for calculating loglik EMA 
	#              (larger a discounts older values faster)
	#' @param burn.in  Number of betas to burn before starting PR average
	#
	#' @return A list of four items
	#				1. A matrix, each column is an iteration of computed betas.
	#				2. A vector of EMA of loglikelihood of one data point for
	#					 each iteration
	#				3. A matrix of PR-averaged betas after burn in
	#				4. A vector of the full likelihood for each beta in No. 3
	#
		# Create zero matrix to store iterative values of beta; initialize beta
		betaSGD <- matrix(rep(NA, ncol(X) * (n.iterSGD+ 1)), nrow = ncol(X))
		betaSGD[, 1] <- beta.init
		
		EMAloglik <- rep(NA, n.iterSGD)
		
		PRbeta <- matrix(rep(NA, ncol(X) * (n.iterSGD - burn.in)), nrow = ncol(X))
		EMA.PR <- rep(NA, n.iterSGD - burn.in)
		#
		for (step in 1:n.iterSGD) {
			# Draw random sample of single row of data (with replacement) 
			i <- sample(nrow(X), 1)
			
			
			# Compute w.i, fitted value of p-parameter of binomial
			w.i <- 1 / (1 + exp(-crossprod(X[i, ], betaSGD[, step])))
			
			# Compute gradient, the descent direction
			grad <- nrow(X) * (m.i * w.i - y[i]) * X[i, ]
			
			# Next set of betas
			newbeta <- betaSGD[, step] - step.size * grad
			betaSGD[, step + 1] <- newbeta
			# EMA of loglik
			step.loglik <- loglik(newbeta, y[i], X[i, ], m.i)
			#
			if (step == 1) {
				EMAloglik[step] <- step.loglik
			}
			else {
				EMAloglik[step] <- a * step.loglik + (1 - a) * EMAloglik[step-1]
			}		
			
			if (step > burn.in) {
				j <- step - burn.in
				if (j == 1) {
					PRbeta[, j] <- newbeta
				}
				else {
					PRbeta[, j] <- (PRbeta[, j - 1] * (j - 1) + newbeta ) / j
				}
				EMA.PR[j] <- loglik(PRbeta[, j], y, X, m.i)
			}
	}
		mylist <- list(betaSGD, EMAloglik, PRbeta, EMA.PR)
		return(mylist)
}






# Create initial values of beta which are some distance away from "true" values
# of beta computed from Newton's method. At least 5 units away plus some 
# exponential noise, in either positive or negative direction
beta0 <- beta.N + (-1) ^ rbinom(ncol(X), 1, 0.5) * (5 + rexp(ncol(X), rate = 1))

# Set number of iterations and stepsize
numiter <- 1e5
my.stepsize <- 0.0001
a <- 0.5 / nrow(X)
burn <- floor(numiter / 1.5)

# Perform SGD, graph trace plot of betas
sgd.1 <- SGD(X, y, m.i, beta0, numiter, my.stepsize, a, burn)

beta.SGD  <- sgd.1[[1]]
EMAloglik <- sgd.1[[2]]
PR.beta   <- sgd.1[[3]]
PR.loglik <- sgd.1[[4]]

graph.betatrace(beta.SGD, beta.N, "SGDbetatrace.pdf")

pdf("SGDloglik.pdf")
plot(nrow(X)*EMAloglik, 
	log = "y", 
	main = "Exponential moving average of log-likelihood from SGD",
	ylab = "log-likelihood", 
	xlab = "iteration", type = "l")
dev.off()

graph.betatrace(PR.beta, beta.N, "SGDbetatracePR.pdf")

pdf("SGDloglikPR.pdf")
plot(PR.loglik, 
	main = "Full Log-likelihood from PR averaging of SGD",
	ylab = "log-likelihood", 
	xlab = "iteration", type = "l")
dev.off()

print(PR.beta[, ncol(PR.beta)])

# > print(PR.beta[, ncol(PR.beta)])
#  [1]  0.4253684 -6.8026923  1.6917163 -1.5052068 13.5295289  1.1387835
#  [7] -0.1433198  0.8158190  2.4698268  0.4751768 -0.5073092

# # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # #

# Decaying steps


SGD.decay <- function(X, y, m.i, beta.init, 
				n.iterSGD, C, t.0, alpha, 
				a, burn.in) {
	#' Perform stochastic gradient descent for a binomial logistic regression
	#' 
	#' @param X  N by P matrix of covariate data
	#' @param y  P-vector of responses
	#' @param m.i  n-parameter of binomial (1 for case of binary logistic)
	#' @param beta.init  Initial guess for coefficients (betas)
	#' @param n.iterSGD  Number of iterations to loop through
	#' @param C  Constant C in Robbins-Monro rule
	#' @param t.0  Constant t.0 in Robbins-Monro rule
	#' @param alpha Constant alpha in Robbins-Monro rule
	#' @param a  Coefficient for weighting decrease for calculating loglik EMA 
	#              (larger a discounts older values faster)
	#' @param burn.in  Number of betas to burn before starting PR average
	#
	#' @return A list of four items
	#				1. A matrix, each column is an iteration of computed betas.
	#				2. A vector of EMA of loglikelihood of one data point for
	#					 each iteration
	#				3. A matrix of PR-averaged betas after burn in
	#				4. A vector of the full likelihood for each beta in No. 3
	#
		# Create zero matrix to store iterative values of beta; initialize beta
		betaSGD <- matrix(rep(NA, ncol(X) * (n.iterSGD+ 1)), nrow = ncol(X))
		betaSGD[, 1] <- beta.init
		
		EMAloglik <- rep(NA, n.iterSGD)
		
		PRbeta <- matrix(rep(NA, ncol(X) * (n.iterSGD - burn.in)), nrow = ncol(X))
		EMA.PR <- rep(NA, n.iterSGD - burn.in)
		#
		for (step in 1:n.iterSGD) {
			# Draw random sample of single row of data (with replacement) 
			i <- sample(nrow(X), 1)
			
			
			# Compute w.i, fitted value of p-parameter of binomial
			w.i <- 1 / (1 + exp(-crossprod(X[i, ], betaSGD[, step])))
			
			# Compute gradient, the descent direction
			grad <- nrow(X) * (m.i * w.i - y[i]) * X[i, ]
			
			# Calculate decaying step size
			step.decay <- C * (step + t.0) ^ -alpha
			
			# Next set of betas
			newbeta <- betaSGD[, step] - step.decay * grad
			betaSGD[, step + 1] <- newbeta
			# EMA of loglik
			step.loglik <- loglik(newbeta, y[i], X[i, ], m.i)
			#
			if (step == 1) {
				EMAloglik[step] <- step.loglik
			}
			else {
				EMAloglik[step] <- a * step.loglik + (1 - a) * EMAloglik[step-1]
			}		
			
			if (step > burn.in) {
				j <- step - burn.in
				if (j == 1) {
					PRbeta[, j] <- newbeta
				}
				else {
					PRbeta[, j] <- (PRbeta[, j - 1] * (j - 1) + newbeta ) / j
				}
				EMA.PR[j] <- loglik(PRbeta[, j], y, X, m.i)
			}
	}
		mylist <- list(betaSGD, EMAloglik, PRbeta, EMA.PR)
		return(mylist)
}

# SGD.decay <- function(X, y, m.i, beta.init,
# 				n.iterSGD, C, t.0, alpha,
# 				a, burn.in)

beta1 <- beta.N + (-1) ^ rbinom(ncol(X), 1, 0.5) * (5 + rexp(ncol(X), rate = 1))

# Define number of iterations & parameters of Robbins-Monro


numiter2 <- 2e5
a2 <- 0.5 / nrow(X)
burn2 <- floor(numiter2 / 1.5)

C     <-  0.1
t.0   <-  1
alpha <-  0.6

decay.1 <- SGD.decay(X, y, m.i, beta1, numiter2, C, t.0, alpha, a2, burn2)



dbeta.SGD  <- decay.1[[1]]
dEMAloglik <- decay.1[[2]]
dPR.beta   <- decay.1[[3]]
dPR.loglik <- decay.1[[4]]

graph.betatrace(dbeta.SGD, beta.N, "dSGDbetatrace.pdf")

pdf("dSGDloglik.pdf")
plot(nrow(X)*dEMAloglik[10:length(dEMAloglik)], 
	log = "y", 
	main = "Exponential moving average of log-likelihood from decay SGD",
	ylab = "log-likelihood", 
	xlab = "iteration", type = "l")
dev.off()

graph.betatrace(dPR.beta, beta.N, "dSGDbetatracePR.pdf")

pdf("dSGDloglikPR.pdf")
plot(dPR.loglik, 
	main = "Full log-likelihood from PR averaging of decay SGD",
	ylab = "log-likelihood", 
	xlab = "iteration", type = "l")
dev.off()



# Running average of beta with decaying steps

# END