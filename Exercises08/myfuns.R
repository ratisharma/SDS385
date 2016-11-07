# James Scott
makeD2_sparse = function (dim1, dim2)  {
	require(Matrix)
    D1 = bandSparse(dim1 * dim2, m = dim1 * dim2, k = c(0, 1), 
        diagonals = list(rep(-1, dim1 * dim2), rep(1, dim1 * 
            dim2 - 1)))
    D1 = D1[(seq(1, dim1 * dim2)%%dim1) != 0, ]
    D2 = bandSparse(dim1 * dim2 - dim1, m = dim1 * dim2, k = c(0, 
        dim1), diagonals = list(rep(-1, dim1 * dim2), rep(1, 
        dim1 * dim2 - 1)))
    return(rBind(D1, D2))
}

# Hereafter all code written by Spencer Woody

prox.ell1 <- function(u, lambda) {
	u.prime <- abs(u) - lambda
	prox <- sign(u) * u.prime * (uprime > 0)
	
	return(prox)
}

defaultLaplace <- function(y, lambda = 0.01) {
	# y is a sparse image matrix 
	# lambda is the ell-2 penalty
	require(Matrix)
	
	Nx <- nrow(y)
	Ny <- ncol(y)
	N <- Nx * Ny
	
	y <- as.numeric(unlist(y))
	
	zerolist <- which(y == 0)
	
	D <- makeD2_sparse(Nx, Ny) 
	C <- Matrix::Diagonal(N) + lambda * crossprod(D)
	
	x <- solve(C, y)
	x <- as.numeric(x)
	x[zerolist] <- 0
	x <- Matrix::Matrix(x, nrow = Nx)
	
	return(x)
}


GSLaplace <- function(y, lambda = 0.01) {
	# y is a sparse image matrix 
	# lambda is the ell-2 penalty
	require(Matrix)
	
	Nx <- nrow(y)
	Ny <- ncol(y)
	N <- Nx * Ny
	
	y <- as.numeric(unlist(y))
	
	zerolist <- which(y == 0)
	
	D <- makeD2_sparse(Nx, Ny) 
	C <- Matrix::Diagonal(N) + lambda * crossprod(D)
	
	x <- solve(C, y)
	x <- as.numeric(x)
	x[zerolist] <- 0
	x <- Matrix::Matrix(x, nrow = Nx)
	
	return(x)
}

LOOCV <- function(y, yhat, S) {
	res.sq <- ((y - yhat) / diag(S))^2
	
	return(mean(res.sq))
}


admmFusedLasso <- function(y, maxiter = 1e3, rho = 0.5, lambda = 0.5, 
	e.abs = 1e-3, e.rel = 1e-3) {
	require(Matrix)
	
	Nx <- nrow(y)
	Ny <- ncol(y)
	y <- as.vector(y)
	
	zerolist <- which(y == 0)
	
	D <- makeD2_sparse(Nx, Ny) 
	
	N <- dim(D)[2]
	M <- dim(D)[1]

	C <- Diagonal(N) + rho * crossprod(D)
	
	x <- rep(0, N)
	r <- rep(0, M)
	u <- rep(0, M)
	
	for (i in 1:maxiter) {
		
		# Update and store x; update(but do not store) r)
		x <- solve(C, (y + rho * crossprod(D, r - u)), sparse = T)
		
		# Precache Dx = D
		Dx <- D %*% x
		r.new <- prox.ell1(Dx + u, lambda / rho)
		
		# Update residuals
		q <- Dx - r.new
		s <- - rho * t(D) %*% (r.new - r)
		
		# Update r and u
		r <- r.new
		u <- u + q
		
		# Convergence check
		e.prim <- sqrt(M) * e.abs + e.rel * max(sqrt(sum((Dx)^2)), sqrt(sum(r^2)))
		e.dual <- sqrt(N) * e.abs + e.rel * sqrt(sum((rho * t(D) %*% u)^2))

		cond.prim <- (sum(q^2) <= e.prim)
		cond.dual <- (sum(s^2) <= e.dual)
		
		if(cond.prim & cond.dual) {
			break
		}
	}
	x[zerolist] <- 0
	x <- Matrix(as.numeric(x), nrow = Nx)
	mylist <- list("x" = x, "iter" = i)
	return(mylist)
}

findneighbors <- function(xy, Nx = 128, Ny = 128) {
	# xy is from coordinate list
	
	neighborlist <- "something went wrong..."
	
	if (xy == 1) { # top-left corner
		
		neighborlist <- c(2, Nx + 1)
		
	} else if (xy == Nx) { # bottom-left corner
		
		neighborlist <- c(Nx - 1, 2 * Nx)
		
	} else if (xy == (Ny - 1) * Nx + 1) { # top-right corner
		
		neighborlist <- c((Ny - 2) * Nx + 1, (Ny - 1) * Nx + 2)
		
	} else if (xy == Nx * Ny) { # bottom-right corner
		
		#print("square")
		neighborlist <- c( (Ny - 1) * Nx, Nx * Ny - 1)
		
	} else if (xy %in% (2:(Nx-1)) ) { # left-hand edge (noncorner)
		
		neighborlist <- c(xy - 1, xy + 1, xy + Nx)
		
	} else if (xy %in% ((Ny - 1) * Nx + 2:(Nx - 1)) ) { # right-hand edge (noncorner)
		
		neighborlist <- c(xy - 1, xy + 1, xy - Nx)
		
	} else if (xy %in% (1:(Ny - 2) * Nx + 1)) { # upper edge (noncorner)
		
		neighborlist <- c(xy - Nx, xy + 1, xy + Nx)
		
	} else if (xy %in% (2:(Ny - 1) * Nx)) { # lower edge (noncorner)
		
		neighborlist <- c(xy - Nx, xy - 1, xy + Nx)
		
	} else { # all non-edge cases
		
		neighborlist <- c(xy - Nx, xy - 1, xy + 1, xy + Nx)
	
	}
	return(neighborlist)
}

