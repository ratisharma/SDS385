###########################################################
######### Created by Spencer Woody on 18 Sep 2016 #########
###########################################################

# Import functions from linesearch.R

source("linesearch.R")

# Read in data file, scale X (y = 1 represents a malignant tumor)

wdbc <- read.csv("wdbc.csv", header = FALSE)

X <- as.matrix(wdbc[, 3:12])
X <- scale(X)
X <- cbind(rep(1, nrow(X)), X)

y <- wdbc[, 2]
y <- y == "M"
m.i <- 1

# Obtain estimates of beta from Newton's method

beta.N <- as.matrix(rep(0, ncol(X)))

newton.steps <- 10

for (step in 1:newton.steps) {
	Newton.wi <- as.numeric(comp.wi(X, beta.N))
	Hessian   <- t(X) %*% diag(Newton.wi*(1-Newton.wi)) %*% X
	beta.N    <- beta.N - solve(Hessian, grad(beta.N, y, X, m.i))
}

# Produce initial guesses for betas, based on values from Newton's method.

beta0 <- beta.N + (-1) ^ rbinom(ncol(X), 1, 0.5) * (5 + rexp(ncol(X), rate = 1))

# Perform minimization techniques

gd  <- grad.desc(X, y, m.i, beta0, 5e4, 0.025)
gd.l <- gd.line(X, y, m.i, beta0, 5e4)
qn.l <- qn.line(X, y, m.i, beta0, 5e4)


gd.lik   <- gd$Lik.trace
gd.l.lik <- gd.l$Lik.trace
qn.l.lik <- qn.l$Lik.trace

pdf("complikplot.pdf")
plot(gd.lik,
	log = "xy",
	type = "l",
	col = "blue",
	xlab = "iteration",
	ylab = "Log-likelihood",
	main = "Traceplot of convergence to log-likelihood")
lines(gd.l.lik, col = "red")
lines(qn.l.lik, col = "green")
points(length(gd.lik), bb[length(gd.lik)], col = "blue", pch = 19)
points(length(gd.l.lik), gd.l.lik[length(gd.l.lik)], col = "red", pch = 19)
points(length(qn.l.lik), qn.l.lik[length(qn.l.lik)], col = "green", pch = 19)
legend("topright",
	c("Gradient descent",
	"GD with line search",
	"Quasi-Newton with line search"),
	lty = c(1, 1, 1),
	col = c("blue", "red", "green"))
dev.off()	
