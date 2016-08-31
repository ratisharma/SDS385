###########################################################
######### Created by Spencer Woody on 24 Aug 2016 #########
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

# Function for computing likelihood

loglik <- function(beta, y, X, mi) {
	loglik <- apply((mi - y) * (X %*% beta)+ mi*log(1 + exp(-X %*% beta)), 2, sum)
	return(loglik)
}
# Function for computing gradient for likelihood

grad.loglik <- function(beta, y, X, mi){
  grad <- array(NA, dim = length(beta))
  wi <- comp.wi(X, beta)
  grad <- apply(X*as.numeric(mi * wi - y), 2, sum)
  return(grad)
}

###
### Gradient descent
###

stepfactor <- 0.025
n.steps <- 50000
log.lik <- NULL

for (step in 1:n.steps) {
	log.lik[step] <- loglik(beta, y, X, mi)
	beta <- beta - stepfactor * grad.loglik(beta, y, X, mi)
}

# Create trace plot of likelihood, check for convergence

png("beta_trace1.png")
plot(log.lik,
	 main = "Trace Plot for log likelihood(beta)", 
	 xlab = "Step",
	 ylab = "log likelihood(beta)")
dev.off()

# Compare results to R's glm function

mymodel <- glm(y ~ X[, c(-1)], family = "binomial")
summary(mymodel)

print(beta)

###
### Newton's method
###

beta.N <- as.matrix(rep(0, ncol(X)))

n.steps <- 10
log.lik2 <- NULL

for (step in 1:n.steps) {
	log.lik2[step] <- loglik(beta, y, X, mi)
	w.i <- as.numeric(comp.wi(X, beta.N))
	W <- diag(w.i*(1-w.i))
	Hessian <- t(X) %*% W %*% X
	beta.N <- beta.N - solve(Hessian) %*% grad.loglik(beta.N, y, X, mi)
}

# Show estimates from Newton's method

round(as.matrix(coef(mymodel)) - beta.N, 8)

