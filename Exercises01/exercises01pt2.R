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

stepfactor <- 0.025
n.steps <- 50000
log.lik <- NULL

for (step in 1:n.steps) {
	log.lik[step] <- loglik(beta, y, X, mi)
	beta <- beta - stepfactor * grad.loglik(beta, y, X, mi)
}

png("beta_trace1.png")
plot(log.lik,
	 main = "Trace Plot for log likelihood(beta)", 
	 xlab = "Step",
	 ylab = "log likelihood(beta)")
dev.off()

mymodel <- glm(y ~ X[, c(-1)], family = "binomial")
summary(mymodel)

print(beta)

# Newton's method

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

round(as.matrix(coef(mymodel)) - beta.N, 8)

# RESULTS
#
# FROM glm
# Coefficients:
#               Estimate Std. Error z value Pr(>|z|)
# (Intercept)    0.48702    0.56432   0.863   0.3881
# X[, c(-1)]V3  -7.22185   13.09494  -0.551   0.5813
# X[, c(-1)]V4   1.65476    0.27758   5.961  2.5e-09 ***
# X[, c(-1)]V5  -1.73763   12.27499  -0.142   0.8874
# X[, c(-1)]V6  14.00485    5.89090   2.377   0.0174 *
# X[, c(-1)]V7   1.07495    0.44942   2.392   0.0168 *
# X[, c(-1)]V8  -0.07723    1.07434  -0.072   0.9427
# X[, c(-1)]V9   0.67512    0.64733   1.043   0.2970
# X[, c(-1)]V10  2.59287    1.10701   2.342   0.0192 *
# X[, c(-1)]V11  0.44626    0.29143   1.531   0.1257
# X[, c(-1)]V12 -0.48248    0.60406  -0.799   0.4244
# ---
# 
# FROM GRADIENT DESCENT
# > print(beta)
#              [,1]
#  [1,]  0.48553491
#  [2,] -7.14617798
#  [3,]  1.65480926
#  [4,] -1.80712669
#  [5,] 13.99289624
#  [6,]  1.07426088
#  [7,] -0.07318783
#  [8,]  0.67573353
#  [9,]  2.59382838
# [10,]  0.44615349
# [11,] -0.48275720

# FROM NEWTON'S METHOD
# > round(beta.N, 5)
#         [,1]
#      0.48702
# V3  -7.22185
# V4   1.65476
# V5  -1.73763
# V6  14.00485
# V7   1.07495
# V8  -0.07723
# V9   0.67512
# V10  2.59287
# V11  0.44626
# V12 -0.48248


#END