###########################################################
######### Created by Spencer Woody on 07 Sep 2016 #########
###########################################################

# Read in data file, standardize X

wdbc <- read.csv("wdbc.csv", header = FALSE)

X <- as.matrix(wdbc[, 3:12])
X <- scale(X)
X <- cbind(rep(1, nrow(X)), X)

y <- wdbc[, 2]
y <- y == "M"
beta <- as.matrix(rep(0, ncol(X)))
mi <- 1

# Function for computing w.i

comp.wi <- function (X, beta) {
	wi <- 1 / (1 + exp(-X %*% beta))
	return(wi)
}

# Function for computing likelihood, which handles the case that -X^T * beta is huge

loglik <- function(beta, y, X, mi) {
	XtBeta <- -X %*% beta
	if (max(XtBeta) > 700) {
		loglik <- apply((mi - y) * (X %*% beta)+ mi*XtBeta, 2, sum)
	}
	else {
		loglik <- apply((mi - y) * (X %*% beta)+ mi*log(1 + exp(XtBeta)), 2, sum)
		 }
	return(loglik)
}

# Function for computing gradient for likelihood

grad.loglik <- function(beta, y, X, mi){
  grad <- array(NA, dim = length(beta))
  wi <- comp.wi(X, beta)
  grad <- apply(X*as.numeric(mi * wi - y), 2, sum)
  return(grad)
}

### GRADIENT DESCENT

stepfactor <- 0.025
n.steps <- 50000
log.lik <- rep(NULL, n.steps + 1)

log.lik[1] <- loglik(beta, y, X, mi)

for (step in 1:n.steps) {
	beta <- beta - stepfactor * grad.loglik(beta, y, X, mi)
	log.lik[step + 1] <- loglik(beta, y, X, mi)
}

# Create trace plot of likelihood, check for convergence

png("beta_trace1.png")
plot(log.lik,
	 main = "Trace Plot for log likelihood(beta) using gradient descent", 
	 xlab = "Step (log scale)",
	 ylab = "log likelihood(beta) (log scale)",
	 log = "xy",
	 pch = 20)
dev.off()

# Compare results to R's glm function

mymodel <- glm(y ~ X[, c(-1)], family = "binomial")
summary(mymodel)

print(beta)

### NEWTON'S METHOD

beta.N <- as.matrix(rep(0, ncol(X)))

n.steps2 <- 10
log.lik2 <- rep(NULL, n.steps + 1)

log.lik2[1] <- loglik(beta.N, y, X, mi)

for (step in 1:n.steps2) {
	w.i <- as.numeric(comp.wi(X, beta.N))
	Hessian <- t(X) %*% diag(w.i*(1-w.i)) %*% X
	beta.N <- beta.N - solve(Hessian) %*% grad.loglik(beta.N, y, X, mi)
	log.lik2[step + 1] <- loglik(beta.N, y, X, mi)
}


png("beta_trace2.png")
plot(log.lik2,
	 main = "Trace Plot for log likelihood(beta) using Newton's method", 
	 xlab = "Step",
	 ylab = "log likelihood(beta)",
	 pch = 20)
dev.off()

# Show estimates from Newton's method

print(beta.N)