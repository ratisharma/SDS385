###########################################################
######### Created by Spencer Woody on 24 Aug 2016 #########
###########################################################

wdbc <- read.csv("wdbc.csv", header = FALSE)

X <- as.matrix(wdbc[, 3:12])
X <- scale(X)
X <- cbind(rep(1, nrow(X)), X)

y <- wdbc[, 2]
y <- y == "M"
beta <- as.matrix(rep(0, ncol(X)))
mi <- 1

comp.wi <- function (X, beta) {
	wi <- 1 / (1 + exp(-X %*% beta))
	return(wi)
}

loglik <- function(beta, y, X, mi) {
	loglik <- apply((mi - y) * (X %*% beta)+ mi*log(1 + exp(-X %*% beta)), 2, sum)
	return(loglik)
}

grad.loglik <- function(beta, y, X, mi){
  grad <- array(NA, dim = length(beta))
  wi <- comp.wi(X, beta)
  grad <- apply(X*as.numeric(mi * wi - y), 2, sum)
  return(grad)
}

stepfactor <- 0.01
n.steps <- 30000
log.lik <- NULL

for (step in 1:n.steps) {
	log.lik[step] <- loglik(beta, y, X, mi)
	grad <- grad.loglik(beta, y, X, mi)
	beta <- beta - stepfactor * grad
}

mymodel <- glm(y ~ X[, c(-1)], family = "binomial")
summary(mymodel)

print(beta)

