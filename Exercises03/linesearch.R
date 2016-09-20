###########################################################
######### Created by Spencer Woody on 18 Sep 2016 #########
###########################################################

# Function for computing w.i (logit transform of Xtbeta)

comp.wi <- function (X, beta) {
	wi <- 1 / (1 + exp(-X %*% beta))
	return(wi)
}

# Function for full likelihood

lik <- function(beta, y, X, m.i) {
	loglik <- apply((m.i - y) * (X %*% beta) + m.i * log(1 + exp(-X %*% beta)), 2, sum)
	return(loglik)
}

# Function for computing gradient of likelihood

grad <- function(beta, y, X, mi){
  grad <- array(NA, dim = length(beta))
  wi   <- comp.wi(X, beta)
  grad <- apply(X*as.numeric(mi * wi - y), 2, sum)
  return(grad)
}

# Hessian function

Hess <- function(beta, y, X, mi) {
	w.i <- as.numeric(comp.wi(X, beta))
	Hessian <- t(X) %*% diag(w.i*(1-w.i)) %*% X
	return(Hessian)
}


# Line search algorithm

linesearch <- function(beta.k, y, X, m.i, 
	direct, lik.k, grad.k,
	c = 0.001, max.alpha = 1, rho = 0.75) {
	best.alpha <- max.alpha
	# Initial Boolean value for armijo condition
	armijo.cond <-   ( lik(beta.k + best.alpha * direct, y, X, m.i) 
				     <= lik.k + c * best.alpha * crossprod(grad.k, direct) )
	while (!armijo.cond) {
		# Reduce alpha, recompute Boolean value for Armijo condition
		best.alpha <- rho * best.alpha
		armijo.cond <-   ( lik(beta.k + best.alpha * direct, y, X, m.i) 
					     <= lik.k + c * best.alpha * crossprod(grad.k, direct) )
	}
	return(best.alpha)
}


# Gradient descent function

grad.desc <- function(X, y, m.i, beta.init, n.iter, stepfactor) {
	p <- ncol(X)
	beta <- matrix(rep(NA, p * (n.iter + 1)), nrow = p)
	beta[, 1] <- beta.init
	
	lik.trace <- rep(NA, n.iter + 1)
	lik.trace[1] <- lik(beta[, 1], y, X, m.i)
	
	message <- "Convergence not reached!"
	for (iter in 1:n.iter) {
		beta.iter <- beta[, iter]
		
		# Calculate gradient and new beta, update beta and likelihood
		grad.iter <- grad(beta.iter, y, X, m.i)
		
		newbeta <- beta.iter - stepfactor * grad.iter
		beta[, iter + 1] <- newbeta
		lik.trace[iter + 1] <- lik(newbeta, y, X, m.i)
		
		# Convergence check
		lik.diff <- lik.trace[iter] - lik.trace[iter + 1]
		
		if (lik.diff < 1e-7) {
			beta <- beta[, -(iter:ncol(beta))] # remove excess cols of beta
			lik.trace <- lik.trace[-(iter:length(lik.trace))] # same for lik
			message <- sprintf("Stopped after %i iterations", iter) 
			break
		}
	}
	mylist <- list(Beta = beta, Lik.trace = lik.trace, Conv = message)
	return(mylist)
}

# Gradient descent function with line search

gd.line <- function(X, y, m.i, beta.init, n.iter) {
	# Initialize beta matrix and likelihood trace
	p <- ncol(X)
	beta <- matrix(rep(NA, p * (n.iter + 1)), nrow = p)
	beta[, 1] <- beta.init
	
	lik.trace <- rep(NA, n.iter + 1)
	lik.trace[1] <- lik(beta[, 1], y, X, m.i)
	
	message <- "Convergence not reached!"
	for (iter in 1:n.iter) {
		beta.iter <- beta[, iter]
		
		# Calculate gradient and direction
		grad.iter <- grad(beta.iter, y, X, m.i)
		direct.iter <- -(grad.iter)	
		
		# Store current likelihood
		lik.iter <- lik.trace[iter]
		
		# Perform linesearch
		stepsize.iter <- linesearch(beta.iter, y, X, m.i, 
								    direct.iter, lik.iter, grad.iter)	
		
		# Update beta using calculated stepsize and direction
		newbeta <- beta.iter + stepsize.iter * direct.iter
		
		# Store new values of beta and likelihood
		beta[, iter + 1] <- newbeta
		lik.trace[iter + 1] <- lik(newbeta, y, X, m.i)
		
		# Convergence check
		lik.diff <- lik.trace[iter] - lik.trace[iter + 1]
		if (lik.diff < 1e-7) {
			beta <- beta[, -(iter:ncol(beta))] # remove excess cols of beta
			lik.trace <- lik.trace[-(iter:length(lik.trace))] # same for lik
			message <- sprintf("Stopped after %i iterations", iter) 
			break
		}
	}
	mylist <- list(Beta = beta, Lik.trace = lik.trace, Conv = message)
	return(mylist)
}

# Quasi-Newton method

qn.line <- function(X, y, m.i, beta.init, n.iter) {
	# Initialize beta matrix, likelihood trace, and gradient
	p <- ncol(X)
	beta <- matrix(rep(NA, p * (n.iter + 1)), nrow = p)
	beta[, 1] <- beta.init
	
	lik.trace <- rep(NA, n.iter + 1)
	lik.trace[1] <- lik(beta[, 1], y, X, m.i)
	
	grad.new <- grad(beta.init, y, X, m.i)
	
	# Initialize inverse Hessian approximation with identity matrix
	id.mat <- diag(1, p)
	Hk <- id.mat
	
	message <- "Convergence not reached!"
	for (iter in 1:n.iter) {
		# "New" values for grad and beta from last iteration become "old" values
		grad.old <- grad.new
		beta.old <- beta[, iter]
		lik.iter <- lik.trace[iter]
		
		# Calculate direction from previous gradient and inverse Hessian
		direct.iter <- - Hk %*% grad.old
		
		# Perform linesearch
		stepsize.iter <- linesearch(beta.old, y, X, m.i, 
								    direct.iter, lik.iter, grad.old)	
		
		# Create new estimate of beta from calculated stepsize and direction
		beta.new <- beta.old + stepsize.iter * direct.iter
		
		# Update gradient
		grad.new <- grad(beta.new, y, X, m.i)
		
		# Update betas and likelihood trace
		beta[, iter + 1] <- beta.new
		lik.trace[iter + 1] <- lik(beta.new, y, X, m.i)		
		
		# Update Hk for next iteration
		y.k <- grad.new - grad.old
		s.k <- beta.new - beta.old
		
		ys.k <- tcrossprod(y.k, s.k)
		rho.k <- as.numeric(1 / (crossprod(y.k, s.k)))
		Hk <- (id.mat - rho.k * ys.k) %*% Hk %*% (id.mat - rho.k * t(ys.k)) + 
			   rho.k * s.k %*% t(s.k) # Formula from N&W p. 25
		
		# Convergence check
		lik.diff <- lik.trace[iter] - lik.trace[iter + 1]
		if (lik.diff < 1e-7) {
			beta <- beta[, -(iter:ncol(beta))] # remove excess cols of beta
			lik.trace <- lik.trace[-(iter:length(lik.trace))] # same for lik
			message <- sprintf("Stopped after %i iterations", iter) 
			break
		}
	}
	mylist <- list(Beta = beta, Lik.trace = lik.trace, Conv = message)
	return(mylist)
}


