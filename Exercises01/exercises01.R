###########################################################
######### Created by Spencer Woody on 23 Aug 2016 #########
###########################################################

library(Matrix)
library(microbenchmark)

# No. 1 pt C

# No. 1 pt D

N <- 2000
P <-  500

X <- matrix(rnorm(N * P), nrow = N)
mask <- matrix(rbinom(N * P, 1, 0.05), nrow = N)
X <- mask * X

# END