######### Created by Spencer Woody on 17 Oct 2016 #########

library(Matrix)
library(ggplot2)

source("makeD2_sparse.R")

brain <- read.csv("fmri_z.csv", header = T)

brain <- as.matrix(brain)
brain <- Matrix::Matrix(brain, sparse = T)
brain <- t(brain)
image(brain, cuts = 80)

sum(brain == 0) / length(brain)
# [1] 0.6452026

#brain.m <- melt(brain)
#ggplot(brain.)

D <- makeD2_sparse(128, 128) 


n <- length(brain)
lambda <- 0.01

y <- 1:n

C <- Matrix::Diagonal(n) + lambda * crossprod(D)
	
defaultsolve <- solve(C, y)	

#Gauss-Seidel

lower.tri(A, diag = F) * A + Diagonal(nrow(A)) * diag(A)
upper.tri(A, diag = F) * A

mychol <- chol(C)
decomp <- as(mychol, "sparseMatrix")
u <- forwardsolve(t(decomp), y)
x <- backsolve(decomp, u)

mychol <- chol(C)


# Jacobi 
DIAG <- Diagonal(nrow(A)) * diag(A)
ARRR <- A - DIAG

# library(reshape2)
# library(ggplot2)
# dat <- matrix(rnorm(100, 3, 1), ncol=10)
# names(dat) <- paste("X", 1:10)
# dat2 <- melt(dat, id.var = "X1")
#
# ggplot(brain.m, aes(as.factor(Var1), Var2, group=Var2)) +
#     geom_tile(aes(fill = value)) +
#     geom_text(aes(fill = brain.m$value, label = round(brain.m$value, 1))) +
#     scale_fill_gradient(low = "white", high = "red")
#
#
#
# 	library(reshape2)
# 	library(ggplot2)
# 	dat <- matrix(rnorm(100, 3, 1), ncol=10)
# 	names(dat) <- paste("X", 1:10)
# 	dat2 <- melt(dat, id.var = "X1")
# 	ggplot(dat2, aes(as.factor(Var1), Var2, group=Var2)) +
# 	    geom_tile(aes(fill = value)) +
# 	    geom_text(aes(fill = dat2$value, label = round(dat2$value, 1))) +
# 	    scale_fill_gradient(low = "white", high = "red")
