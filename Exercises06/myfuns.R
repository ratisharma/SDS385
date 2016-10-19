######### Created by Spencer Woody on 17 Oct 2016 #########

nll <- function(X, Y, beta) {
	mynll <- t(beta) %*% (0.5 * crossprod(X) %*% beta - crossprod(X, Y))
	
	return(mynll)
}

grad <- function(X, Y, beta) {
	myprod <- Y - X%*%beta
	mygrad <- - crossprod(X, myprod)
	
	return(mygrad)
}

prox.ell1 <- function(u, lambda) {
	uprime <- abs(u) - lambda
	uprime.zeros <- cbind(rep(0, length(u)), uprime)
	prox <- sign(u) * apply(uprime.zeros, 1, max)
	
	return(prox)
}


proxlasso <- function(X, Y, lambda = 1e-4, gamma = 1e-4, beta0 = NA, numsteps = 50000, tol = 1e-10) {
	# If there is no beta0, make it a series of zeros	
	if (is.na(beta0)) {
		beta0 <- rep(0, ncol(X))
	}

	# Standardize lambda
	lambda = lambda * nrow(X)
	
	# Initialize matrix for 
	beta <- matrix(rep(NA, ncol(X) * (numsteps+1)), nrow = ncol(X))
	beta[, 1] <- beta0

	# Begin trace of nll
	nlltrace <- rep(NA, (numsteps + 1))
	nlltrace[1] <- nll(X, Y, beta0) + lambda * sum(abs(beta0))

	# Initialize convergence message in case convergence not reached
	message <- "Convergence not reached..."

	for (t in 1:numsteps) {
		u <- beta[, t] - gamma * grad(X, Y, beta[, t])
	
		beta[, t + 1] <- prox.ell1(u, gamma * lambda)
	
		nlltrace[t + 1] <- nll(X, Y, beta[, t + 1]) + lambda * sum(abs(beta[, t + 1]))
	
		# Convergence check
		nllchange <- nlltrace[t + 1] / nlltrace[t] - 1
		if (abs(nllchange) < tol) {
			# Remove excess betas and nll
			beta <- beta[, -((t+2):ncol(beta))]
			nlltrace <- nlltrace[-((t+2):length(nlltrace))]
			# Update convergence message
			message <- sprintf("Convergence reached after %i iterations", (t+1))
			break
		}
	}
	
	mylist <- list("finalbeta" = beta[, ncol(beta)], "beta" = beta, "nll" = nlltrace, "conv" = message)
	return(mylist)
}



proxlassofast <- function(X, Y, lambda = 1e-4, gamma = 1e-4, beta0 = NA, numsteps = 50000, tol = 1e-10) {
	# If there is no beta0, make it a series of zeros	
	if (is.na(beta0)) {
		beta0 <- rep(0, ncol(X))
	}

	# Standardize lambda
	lambda = lambda * nrow(X)

	# Create s vector
	s <- rep(NA, numsteps + 1)
	s[1] <- 10

	for (j in 2:length(s)) {
		s[j] <- (1 + (1 + 4 * s[j - 1] ^ 2) ^ 0.5) / 2
	}

	# Initialize z matrix
	z <- matrix(0, nrow = ncol(X), ncol = (numsteps + 1))
	
	# Initialize matrix for 
	beta <- matrix(rep(NA, ncol(X) * (numsteps+1)), nrow = ncol(X))
	beta[, 1] <- beta0

	# Begin trace of nll
	nlltrace <- rep(NA, (numsteps + 1))
	nlltrace[1] <- nll(X, Y, beta0) + lambda * sum(abs(beta0))

	# Initialize convergence message in case convergence not reached
	message <- "Convergence not reached..."

	for (t in 1:numsteps) {
		u <- z[, t] - gamma * grad(X, Y, z[, t])
	
		beta[, t + 1] <- prox.ell1(u, gamma * lambda)

		z[, t + 1] <- beta[, t + 1] + (s[t] - 1) / s[t + 1] * (beta[, t + 1] - beta[, t])

		nlltrace[t + 1] <- nll(X, Y, beta[, t + 1]) + lambda * sum(abs(beta[, t + 1]))
	
		# Convergence check
		nllchange <- nlltrace[t + 1] / nlltrace[t] - 1
		if (abs(nllchange) < tol) {
			# Remove excess betas and nll
			beta <- beta[, -((t+2):ncol(beta))]
			nlltrace <- nlltrace[-((t+2):length(nlltrace))]
			# Update convergence message
			message <- sprintf("Convergence reached after %i iterations", (t+1))
			break
		}
	}
	
	mylist <- list("finalbeta" = beta[, ncol(beta)], "beta" = beta, "nll" = nlltrace, "conv" = message)
	return(mylist)
}
