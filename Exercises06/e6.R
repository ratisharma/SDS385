######### Created by Spencer Woody on 17 Oct 2016 #########

library(MASS)
library(glmnet)
library(ggplot2)

source("myfuns.R")

# Read in data, scale everything
X <- as.matrix(read.csv("diabetesX.csv", header = TRUE))
Y <- as.numeric(unlist(read.csv("diabetesY.csv", header = FALSE)))

X <- scale(X)
Y <- scale(Y)


lambdas <- seq(5e-5, 0.5, length.out = 50)

# Plot coefficients from lasso regression

betamat <- matrix(rep(NA, length(lambdas)*ncol(X)), ncol = length(lambdas))

for (i in 1:length(lambdas)) {
	mylasso <- proxlassofast(X, Y, lambda = lambdas[i])
	betamat[, i]  <- mylasso$finalbeta
}

ybound <- c(min(betamat), max(betamat))

pdf("myplot.pdf")
plot(lambdas, betamat[1, ], 
	ylim = ybound, 
	log = "x",
	main = "Coefficients of Lasso Regression using my own lasso function",
	xlab = "",
	ylab = "Coefficient",
	type = "l", 
	col = 2)
for (j in 2:nrow(betamat)) {
	lines(lambdas, betamat[j, ], col = (j + 4))
}
dev.off()


########################################################################

# Compare ordinary and accelerated versions of proximal gradient

m <- -20

mylasso <- proxlasso(X, Y, gamma = 0.0002, lambda = 0.001, tol = 10^m)
mylasso2 <- proxlassofast(X, Y, gamma = 0.0002, lambda = 0.001, tol = 10^m)

pdf("slowfast.pdf")
plot(mylasso$nll, type = "l", col = "black", 
log = "x",
xlab = "iteration",
ylab = "Negative log likelihood",
main = sprintf("NLL for Prox Method and Accelerated Prox method, tol = 1e%i", m))
lines(mylasso2$nll, col = "red")
points(length(mylasso$nll), mylasso$nll[length(mylasso$nll)], pch = 19, col = "black")
points(length(mylasso2$nll), mylasso2$nll[length(mylasso2$nll)], pch = 19, col = "red")
legend("topright", legend=c("Prox method", "Accelerated prox method"),
       col=c("black", "red"), lty=c(1,1), cex=0.8)
dev.off()

print(length(mylasso$nll))
print(length(mylasso2$nll))

sum(sign(mylasso$finalbeta) == sign(mylasso2$finalbeta))
ncol(X)

pdf("compcoefs.pdf")
plot(mylasso$finalbeta, mylasso2$finalbeta,
	main = "Comparison of Coefficients from two methods",
	xlab = "Coefficients from Prox Method",
	ylab = "Coefficients from Accel Prox Method")
abline(0, 1)
dev.off()

mylasso$finalbeta
mylasso2$finalbeta

# > print(length(mylasso$nll))
# [1] 19953
# > print(length(mylasso2$nll))
# [1] 616

# > mylasso$finalbeta
#  [1]  0.023560578 -0.153789540  0.282694838  0.208274783 -0.103033550
#  [6]  0.000000000 -0.113948268  0.053230187  0.407621844  0.041754067
# [11]  0.046650203  0.033109963 -0.008095386  0.057726330  0.000000000
# [16] -0.019305929  0.080070994  0.219445533  0.067286417  0.103984948
# [21] -0.004789084  0.011970890  0.000000000 -0.108030864  0.082766137
# [26]  0.086187888  0.034960689  0.028160431  0.045641503  0.050263648
# [31]  0.108286794 -0.101133491  0.000000000 -0.044642852 -0.037578683
# [36]  0.024963870  0.100191687 -0.113246424  0.075705929  0.038941220
# [41] -0.009111873  0.007940555  0.020034844  0.047921927  0.001962857
# [46]  0.000000000 -0.019820825 -0.003312034 -0.082569062  0.042338316
# [51]  0.000000000 -0.197988211 -0.401744208 -0.036250001 -0.066970573
# [56] -0.015327546  0.374411995  0.000000000  0.002665353  0.123333658
# [61]  0.096873190 -0.035289003  0.133837012  0.019053591
# > mylasso2$finalbeta
#  [1]  0.023538900 -0.153795619  0.282633738  0.208328306 -0.103014142
#  [6]  0.000000000 -0.113982017  0.053249091  0.407576274  0.041765184
# [11]  0.046624516  0.033198478 -0.008081280  0.060115879  0.000000000
# [16] -0.020599252  0.078493833  0.220584671  0.067223390  0.103999771
# [21] -0.004744999  0.011962682  0.000000000 -0.107915299  0.082692339
# [26]  0.085975450  0.035031585  0.028175283  0.045726185  0.050280271
# [31]  0.108269248 -0.101125076  0.000000000 -0.044626911 -0.037644834
# [36]  0.024955140  0.100180730 -0.110735176  0.073923786  0.037456934
# [41] -0.009942927  0.007065391  0.020099087  0.048478022  0.001293112
# [46]  0.000000000 -0.019492307 -0.003648320 -0.082558916  0.039746085
# [51]  0.000000000 -0.195383539 -0.404309136 -0.036472332 -0.066619737
# [56] -0.016430649  0.376181451  0.000000000  0.000000000  0.123918165
# [61]  0.097184059 -0.036062764  0.134341872  0.019007764

########################################################################
# > which(beta[, ncol(beta)] == 0)
#  [1] 13 14 15 16 17 21 23 26 31 35 36 38 39 40 41 42 43 45 47 48 50 54 55 56 60
# [26] 61 64


# Make plot tracking s 

numsteps = 1000

s <- rep(NA, numsteps + 1)
s[1] <- 10

for (j in 2:length(s)) {
	s[j] <- (1 + (1 + 4 * s[j - 1] ^ 2) ^ 0.5) / 2
}

weight <- 0

for (k in 2:length(s)) {
	weight[k - 1] <- (s[k] - 1) / (s[k-1])
}

pdf("splot.pdf")
plot(weight, type = "l", main = "Convergence of s", xlab = "k", ylab = "(s[k] - 1) / s[k-1]")
dev.off()


