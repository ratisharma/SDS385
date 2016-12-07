# Written by Spencer Woody on 6 Nov 2016

library(ggplot2)          # plotting
library(gridExtra)        # make a grid with ggplot
library(reshape2)         # melt function
library(caTools)          # trapz function
library(foreach)          # faster for-loops
library(doParallel)       # works with foreach
library(RDAVIDWebService) # run queries to DAVID API

# Use this to install RDAVIDWebService
# source("https://bioconductor.org/biocLite.R")
# biocLite("RDAVIDWebService") # requires java >=1.8, openssl >= 1.06


################################################################################
################################ RNA processing ################################
################################################################################


# Step 1.1 Normalize RNA to read depth
rna <- read.csv("rna.csv", header = TRUE, as.is = TRUE)
head(rna)

# Store names for later
rna_names <- rna[, 1]

# Read in scaled size factors
ssf <- read.csv("ssf.csv", header = T)
ssf
ssf <- unlist(ssf[1, 2:ncol(ssf)])
ssf

rna_norm <- sweep(as.matrix(rna[, -1]), 2, ssf, `/`)
rna_norm[1, ]

# find average and sd of each RNA for each 3 replicates
rna_norm_av <- matrix(rep(0, nrow(rna) * 9), nrow = nrow(rna))
rna_norm_sd <- matrix(rep(0, nrow(rna) * 9), nrow = nrow(rna))

for (t in 1:9) {
	rna_norm_t <- rna_norm[, c(t, t + 9, t + 18)]
	rna_norm_av[, t] <- apply(rna_norm_t, 1, mean)
	rna_norm_sd[, t] <- apply(rna_norm_t, 1, sd)
}

# Normalize to maximum value
rna_norm_av_max <- apply(rna_norm_av, 1, max)

rna_norm_av <- sweep(rna_norm_av, 1, rna_norm_av_max, '/') 
rna_norm_sd <- sweep(rna_norm_sd, 1, rna_norm_av_max, '/') 

# Create data frame with names of RNAs
rna_norm_av <- data.frame(rna_names, rna_norm_av)
colnames(rna_norm_av) <- c("gene_id", "t3", "t4", "t5", "t6", 
						   "t8", "t24", "t48", "t168", "t336")
						   
rna_norm_sd <- data.frame(rna_names, rna_norm_sd)
colnames(rna_norm_sd) <- colnames(rna_norm_av)

rna_norm_av <- na.omit(rna_norm_av)
rna_norm_sd <- na.omit(rna_norm_sd)

# Save them to csv's
write.csv(rna_norm_av, file = "Data/rna_norm_av.csv", row.names = FALSE)
write.csv(rna_norm_sd, file = "Data/rna_norm_sd.csv", row.names = FALSE)


# Note that we will discard insignificant RNA's later, but we will keep 
# all of them right now so that we can look at correlations b/t pro's and RNA's

################################################################################
############################## RNA sig #########################################
################################################################################


# Now, pick only significant RNAs

rna_norm_av <- read.csv("Data/rna_norm_av.csv", header = TRUE)
rna_norm_sd <- read.csv("Data/rna_norm_sd.csv", header = TRUE)

all_rna_names <- as.character(rna_norm_av[, 1])

sig_rna_names <- as.character(read.csv("sig_rnas.csv", header = TRUE)[, 1])

rna_keep <- which(all_rna_names %in% sig_rna_names)

rna_sig_norm_av <- rna_norm_av[rna_keep, ]
rna_sig_norm_sd <- rna_norm_sd[rna_keep, ]

write.csv(rna_sig_norm_av, file = "Data/rna_sig_norm_av.csv", row.names = FALSE)
write.csv(rna_sig_norm_sd, file = "Data/rna_sig_norm_sd.csv", row.names = FALSE)

################################################################################
########################### Protein processing #################################
################################################################################

# Step 1.2 Normalize proteins to read depth
pro <- read.csv("protein.csv", header = TRUE, as.is = TRUE)
head(pro)

# Look for proteins with no counts above 10 to discard later
pro_maxes <- apply(pro[, -1], 1, max)
pro_keep <- which(pro_maxes >= 10)

# Normalize each protein to total read depth at that time and replicate
pro_time_maxes <- colSums(pro[, -1])
pro_norm <- sweep(as.matrix(pro[, -1]), 2, pro_time_maxes, `/`)

# Now throw out insignificant proteins
pro_norm <- pro_norm[pro_keep, ]
pro_names <- pro[pro_keep, 1]

# find average and sd of each protein for each 3 replicates
pro_norm_av <- matrix(rep(0, nrow(pro_norm) * 9), nrow = nrow(pro_norm))
pro_norm_sd <- matrix(rep(0, nrow(pro_norm) * 9), nrow = nrow(pro_norm))

for (t in 1:9) {
	pro_norm_t <- pro_norm[, c(t, t + 9, t + 18)]
	pro_norm_av[, t] <- apply(pro_norm_t, 1, mean)
	pro_norm_sd[, t] <- apply(pro_norm_t, 1, sd)
}

# Normalize it so it's in (0, 1)
pro_norm_av_max <- apply(pro_norm_av, 1, max)

pro_norm_av <- sweep(pro_norm_av, 1, pro_norm_av_max, '/') 
pro_norm_sd <- sweep(pro_norm_sd, 1, pro_norm_av_max, '/') 

# Create a dataframe, and then save it
pro_norm_av <- data.frame(pro_names, pro_norm_av)
colnames(pro_norm_av) <- c("gene_id", "t3", "t4", "t5", "t6", 
						   "t8", "t24", "t48", "t168", "t336")
						   
pro_norm_sd <- data.frame(pro_names, pro_norm_sd)
colnames(pro_norm_sd) <- colnames(pro_norm_av)


write.csv(pro_norm_av, file = "Data/pro_norm_av.csv", row.names = FALSE)
write.csv(pro_norm_sd, file = "Data/pro_norm_sd.csv", row.names = FALSE)

############################################################################
############################## Clustering ##################################
############################################################################

## RNA k-means center plots

rna_sig_norm_av <- read.csv("Data/rna_sig_norm_av.csv", header = T)
pro_norm_av <- read.csv("Data/pro_norm_av.csv", header = T)

# Number of centers
rna_centers <- 15

# Run clustering
rna_k <- kmeans(rna_sig_norm_av[, -1], centers = rna_centers)
rna_k_centers <- rna_k$centers
rna_k_centers_m <- melt(rna_k_centers)

p <- ggplot(rna_k_centers_m, aes(x = Var2, y = Var1, fill = value)) + 
geom_tile() + 
scale_fill_gradient(low = "snow", high = "dodgerblue3") + 
xlab("Time point") + 
ylab("Cluster number") +
ggtitle("RNA Cluster Centroids")
p

## Protein k-means center plots

# Number of centers
pro_centers <- 25

# Run clustering
pro_k <- kmeans(pro_norm_av[, -1], centers = pro_centers)
pro_k_centers <- pro_k$centers
pro_k_centers_m <- melt(pro_k_centers)

q <- ggplot(pro_k_centers_m, aes(x = Var2, y = Var1, fill = value)) + 
geom_tile() + 
scale_fill_gradient(low = "snow", high = "dodgerblue3") + 
xlab("Time point") + 
ylab("Cluster number") +
ggtitle("Protein Cluster Centroids")
q

# Compare how many clusters we should choose for both RNAs and proteins

numclusters <- 1:30
rna_perc_exp <- c()
pro_perc_exp <- c()

for (i in numclusters) {
	rna_k_i <- kmeans(rna_sig_norm_av[, -1], centers = i)
	pro_k_i <- kmeans(pro_norm_av[, -1], centers = i)
	
	rna_perc_exp <- c(rna_perc_exp, rna_k_i$between / rna_k_i$totss * 100)
	pro_perc_exp <- c(pro_perc_exp, pro_k_i$between / pro_k_i$totss * 100)
}

h <- qplot(numclusters, geom = "blank") +
xlab("Number of Clusters") + 
ylab("Percentage of Variation Explained by Clusters") + 
labs(title = "Elbow Method for k-means clustering for RNAs vs. Proteins") +
geom_line(aes(x = numclusters, y = rna_perc_exp, colour = "RNA")) + 
geom_point(aes(x = numclusters, y =  rna_perc_exp, colour = "RNA")) + 
geom_line(aes(x = numclusters, y = pro_perc_exp, colour = "Proteins")) + 
geom_point(aes(x = numclusters, y = pro_perc_exp, colour = "Proteins")) +
scale_colour_manual(name = "", values = c("RNA" = "dodgerblue3", "Proteins" = "firebrick3"))  +
 geom_point(aes(x = numclusters[5], y = pro_perc_exp[5]), colour = "black", shape = 1, size = 4)
h

############################################################################
############################## Pro & RNAs ##################################
############################################################################

# Function for evaluating cumulative trapezoidal integral sums
ctrapz <- function(x, y) {
	require(caTools)
	n <- length(x)
	
	ctrapz <- c()
	
	for (i in 1:(n-1)) {
		ctrapz <- c(ctrapz, trapz(x[1: (i+1) ], y[1:(i+1)]))
	}
	return(ctrapz)
}

# rna_norm_av <- read.csv("Data/rna_norm_av.csv", header = T)
# pro_norm_av <- read.csv("Data/pro_norm_av.csv", header = T)
#
dict <- read.csv("Data/dictionary.csv", header = T)
head(dict)

# Dictionary names for RNAs and proteins
dict_rna_names <- as.character(dict[, 3])
dict_pro_names <- as.character(dict[, 1])

rna_names <- as.character(rna_norm_av[, 1])
pro_names <- as.character(pro_norm_av[, 1])

# Number of proteins
N <- nrow(pro_norm_av)

cor_mat <- matrix(rep(NA, N * 3), nrow = N)

times <- c(3, 4, 5, 6, 8, 24, 48, 24*7, 24*7*2)

# for-loop going over all proteins in file
for (i in 1:N) {
	# Match protein to its transcript
	pro_name_i <- pro_names[i]
	index_p2r <- which(dict_pro_names == pro_name_i) # index for protein -> rna 
	rna_name_i <- dict_rna_names[index_p2r]
	
	index_r2v <- which(rna_names == rna_name_i) # index for rna in rna datafile
	
	# Normalized RNA and protein counts for protein and its rna
	rna_vals <- as.numeric(rna_norm_av[index_r2v, -1])
	pro_vals <- as.numeric(pro_norm_av[i, -1])
	
	# Cumulative integral for RNA values and times
	rna_vals_int <- c(rna_vals[1], ctrapz(times, rna_vals))
	
	cor_mat[i, 1] <- pro_name_i
	cor_mat[i, 2] <- as.numeric(cor(rna_vals, pro_vals, method = "spearman"))
	cor_mat[i, 3] <- as.numeric(cor(rna_vals_int, pro_vals, method = "spearman"))
}

prop_reg <- as.numeric(cor_mat[, 2])
int_reg <- as.numeric(cor_mat[, 3])

# Histogram for Proportional Regulation
u <- qplot(prop_reg,
	xlab = "Spearman Correlation Coefficient",
	main = "Protein vs. RNA (Prop. Regulation)") + 
	geom_histogram(
	binwidth = diff(range(prop_reg)) / 8,
	col = "dodgerblue4",
	fill = "dodgerblue3")
u

# Histogram for Integral Regulation
v <- qplot(int_reg,
	xlab = "Spearman Correlation Coefficient",
	main = "Protein vs. Time Integral of RNA (Int. Regulation)") + 
	geom_histogram(
	binwidth = diff(range(int_reg)) / 8,
	col = "dodgerblue4",
	fill = "dodgerblue3")
v

# 2D Histogram
w <- qplot(prop_reg, int_reg, geom = "blank",
	xlab = "Correlation of Protein vs. RNA",
	ylab = "Correlation of Protein vs. Time Integral of RNA",
	main = "2D Histogram of Different Protein & RNA Regulation") + 
	stat_bin2d(bins = 10) +
	scale_fill_gradientn(colours = c("snow", "dodgerblue3"))
w

############################################################################
################################ Operons ##################################
############################################################################

rna_sig_norm_av <- read.csv("Data/rna_sig_norm_av.csv", header = T)
pro_norm_av <- read.csv("Data/pro_norm_av.csv", header = T)

dict <- read.csv("Data/dictionary.csv", header = T)
head(dict)

dict_rna_names <- as.character(dict[, 3])
dict_pro_names <- as.character(dict[, 1])

operonlist <- as.character(dict$operon.name)
unique_operons <- unique(operonlist)

rna_names <- as.character(rna_sig_norm_av[, 1])
pro_names <- as.character(pro_norm_av[, 1])

op_rna_cor <- c()
op_pro_cor <- c()

for (i in 1:length(unique_operons)) {
	operon_i <- unique_operons[i]
	operon_indices <- which(operonlist == operon_i)
	
	# Skip operons with fewer than 2 proteins
	if (length(operon_indices) < 2) {
		next
	}
	
	operon_rna_names <- dict_rna_names[operon_indices]
	operon_pro_names <- dict_pro_names[operon_indices]
	
	# Get values for RNAs and proteins in this operon
	rna_operon <- rna_sig_norm_av[which(rna_names %in% operon_rna_names), ]
	pro_operon <- pro_norm_av[which(pro_names %in% operon_pro_names), ]

	# Make correlation matrix
	pair_cor_rna <- cor(t(rna_operon[, -1]), method = "spearman")
	pair_cor_pro <- cor(t(pro_operon[, -1]), method = "spearman")

	# Add median pairwise correlation coefficient by selecting lower triangle
	op_rna_cor <- c(op_rna_cor, median(pair_cor_rna[lower.tri(pair_cor_rna)]))
	op_pro_cor <- c(op_pro_cor, median(pair_cor_pro[lower.tri(pair_cor_pro)]))
}

x <- qplot(na.omit(op_rna_cor),
	xlab = "Spearman Correlation Coefficient",
	main = "RNA Intra-operon Correlation") + 
	geom_histogram(
	binwidth = diff(range(na.omit(op_rna_cor))) / 9,
	col = "dodgerblue4",
	fill = "dodgerblue3")
x

y <- qplot(na.omit(op_pro_cor),
	xlab = "Spearman Correlation Coefficient",
	main = "Protein Intra-operon Correlation") + 
	geom_histogram(
	binwidth = diff(range(na.omit(op_pro_cor))) / 9,
	col = "dodgerblue4",
	fill = "dodgerblue3")
y

###############################################################
##################### FITTING A TIME COURSE ###################
###############################################################

detectCores()
cl <- makeCluster(4) # however many processors you have
registerDoParallel(cl)

rna.av <- read.csv("Data/rna_sig_norm_av.csv", header = T, as.is = T)
rna.sd <- read.csv("Data/rna_sig_norm_sd.csv", header = T, as.is = T)

pro.av <- read.csv("Data/pro_norm_av.csv", header = T, as.is = T)
pro.sd <- read.csv("Data/pro_norm_sd.csv", header = T, as.is = T)

dict <- read.csv("Data/dictionary.csv", header = T, as.is = T)

costfun <- function(vec) {
	# INPUTS
	# vec is a vector of guesses for the amplitude parameters (vec[1:3])
	# 	  and time parameters (vec[4:7])
	# AVG and SD must be PRE-DEFINED as averages of presences and 
	# 	  standard deviations of presences, respectively, for each time
	#
	# OUTPUT: Result of cost function defined in Houser et al. p. 21
	
	times <- c(3, 4, 5, 6, 8, 24, 48, 24*7, 24*7*2)
	
	A1 <- vec[1]; A2 <- vec[2]; A3 <- vec[3]
	tau1 <- vec[4]; tau2 <- vec[5]; tau3 <- vec[6]; tau4 <- vec[7]
	
	# Function for the fit
	fit <- (  A1 * (times <= tau1) + 

	( A1 + (A2 - A1) / (tau2 - tau1) * (times - tau1) ) * 
	(times > tau1 & times <= tau2) + 

	A2 * (times > tau2 & times <= tau3) + 

	( A2 + (A3 - A2) / (tau4 - tau3) * (times - tau3) ) * 
	(times > tau3 & times <= tau4) + 

	A3 * (times > tau4) )
	
	# Are the taus in the wrong order? is tau4 > final time? Is tau1 negative?
	check1 <- (tau1 < 0 || tau4 < tau3 || tau3 < tau2 || tau2 < tau1 || tau4 > times[9])
	# Are all the amplitudes within (0, 1)?
	check2 <- (A1 > 1 || A1 < 0 || A2 > 1 || A2 < 0 || A1 > 1 || A1 < 0)
	
	# Compute cost function, plus penalty if conditions are not all met
	# 0.001 added to denominator for numerical stability
	cost <- sum((fit - AVG)^2 / (SD + 0.01)) + (check1 + check2) * 1000

	return(cost)
}


fit <- function(vec, times) {
	##### *** Function for graphing estimated time course *** #####
	#
	# INPUTS
	# vec is a vector of estimates for the amplitude parameters (vec[1:3])
	# 	  and time parameters (vec[4:7])
	# times is a vector of continuos time points
	A1 <- vec[1]
	A2 <- vec[2]
	A3 <- vec[3]
	tau1 <- vec[4]
	tau2 <- vec[5]
	tau3 <- vec[6]
	tau4 <- vec[7]
	
	fit <- (  A1 * (times <= tau1) + 

	( A1 + (A2 - A1) / (tau2 - tau1) * (times - tau1) ) * 
	(times > tau1 & times <= tau2) + 

	A2 * (times > tau2 & times <= tau3) + 

	( A2 + (A3 - A2) / (tau4 - tau3) * (times - tau3) ) * 
	(times > tau3 & times <= tau4) + 

	A3 * (times > tau4) )
	
	return(fit)
}

plotfun <- function(est, AVG, SD, name, type) {
	# est is the estimated parameters
	# AVG is the average of each time point
	# SD is the standard deviation of each time point
	# name is the name of the gene's transcripts as a string
	# type is a string, either "RNA"
	times.cont <- seq(3, 24*7*2, length.out = 1000)
	if (name == "") {name = "NOT FOUND"}
	v <- qplot(times.cont, geom = "blank") +
	xlab("Time (h)") +
	ylab("Relative presence") + 
	labs(title = sprintf("Estimated time course function for %s of %s", type, name)) +
	geom_line(aes(x = times.cont, y = fit(est, times.cont)), col = "red") +
	geom_point(aes(x = times, y = as.numeric(AVG))) +
	geom_errorbar(aes(x=times, 
		ymin = as.numeric(AVG-SD), ymax = as.numeric(AVG+SD)), 
		width=0.05) + 
	coord_cartesian(ylim = c(0, 1)) +
	scale_x_log10()
	
	return(v)
}


###############################
## Take a random observation ##
###############################

set.seed(350)
n <- sample(1:nrow(rna.av), 1)

AVG <- rna.av[n, -1]
SD <- rna.sd[n, -1]

name <- dict[which(dict[, 3] == rna.av[n, 1]), 2]

# Randomize initial amplitude parameters, but keep taus constant for all fitting
times <- c(3, 4, 5, 6, 8, 24, 48, 24*7, 24*7*2)
A1 <- runif(1); A2 <- runif(1); A3 <- runif(1)
tau1 <- times[2]; tau2 <- times[4]; tau3 <- times[6]; tau4 <- times[8]

system.time(myoptim <- optim(c(A1, A2, A3, tau1, tau2, tau3, tau4), 
			costfun, 
			method = "Nelder-Mead", 
			control = list(maxit = 1000)))

# Display results of optim
myoptim$value
myoptim$counts
myoptim$convergence

# Store estimates of parameters and display them
est <- myoptim$par
est

plotfun(est, AVG, SD, name, type = "transcript")


# Make multiple plots of RNAs 
AVGmat <- matrix(rep(0, 6 * length(times)), nrow = 6)
 SDmat <- matrix(rep(0, 6 * length(times)), nrow = 6)
estmat <- matrix(rep(0, 6 * 7), nrow = 6, ncol = 7)
namevec <- vector(length = 6)
convvec <- vector(length = 6)

for (i in 1:6) {
	n <- sample(1:nrow(rna.av), 1)
	
	namevec[i] <- dict[which(dict[, 3] == rna.av[n, 1]), 2]
	
	AVG <- as.numeric(rna.av[n, -1])
	 SD <- as.numeric(rna.sd[n, -1])
	
	AVGmat[i, ] <- AVG 
	 SDmat[i, ] <- SD 
	
	A1 <- runif(1); A2 <- runif(1); A3 <- runif(1)
	tau1 <- times[2]; tau2 <- times[4]; tau3 <- times[6]; tau4 <- times[8]

	optim.i <- optim(c(A1, A2, A3, tau1, tau2, tau3, tau4), 
				     costfun, 
				     method = "Nelder-Mead", 
				     control = list(maxit = 1000))

	convvec[i] <- optim.i$convergence
	estmat[i, ]<- optim.i$par
}

convvec

a  <- plotfun(estmat[1, ], AVGmat[1, ], SDmat[1, ], namevec[1], type = "transcript")
b  <- plotfun(estmat[2, ], AVGmat[2, ], SDmat[2, ], namevec[2], type = "transcript")
c  <- plotfun(estmat[3, ], AVGmat[3, ], SDmat[3, ], namevec[3], type = "transcript")
d  <- plotfun(estmat[4, ], AVGmat[4, ], SDmat[4, ], namevec[4], type = "transcript")
ee <- plotfun(estmat[5, ], AVGmat[5, ], SDmat[5, ], namevec[5], type = "transcript")
ff <- plotfun(estmat[6, ], AVGmat[6, ], SDmat[6, ], namevec[6], type = "transcript")

grid.arrange(a, b, c, d, ee, ff, ncol=2, widths = c(3, 3))


# Make multiple plots of proteins now
AVGmat <- matrix(rep(0, 6 * length(times)), nrow = 6)
 SDmat <- matrix(rep(0, 6 * length(times)), nrow = 6)
estmat <- matrix(rep(0, 6 * 7), nrow = 6, ncol = 7)
namevec <- vector(length = 6)
convvec <- vector(length = 6)

for (i in 1:6) {
	n <- sample(1:nrow(pro.av), 1)
	
	namevec[i] <- dict[which(dict[, 1] == pro.av[n, 1]), 2]
	
	AVG <- as.numeric(pro.av[n, -1])
	 SD <- as.numeric(pro.sd[n, -1])
	
	AVGmat[i, ] <- AVG 
	 SDmat[i, ] <- SD 
	
	A1 <- runif(1); A2 <- runif(1); A3 <- runif(1)
	tau1 <- times[2]; tau2 <- times[4]; tau3 <- times[6]; tau4 <- times[8]

	optim.i <- optim(c(A1, A2, A3, tau1, tau2, tau3, tau4), 
				     costfun, 
				     method = "Nelder-Mead", 
				     control = list(maxit = 1000))

	convvec[i] <- optim.i$convergence
	estmat[i, ]<- optim.i$par
}

convvec

g   <- plotfun(estmat[1, ], AVGmat[1, ], SDmat[1, ], namevec[1], type = "protein")
h   <- plotfun(estmat[2, ], AVGmat[2, ], SDmat[2, ], namevec[2], type = "protein")
eye <- plotfun(estmat[3, ], AVGmat[3, ], SDmat[3, ], namevec[3], type = "protein")
jay <- plotfun(estmat[4, ], AVGmat[4, ], SDmat[4, ], namevec[4], type = "protein")
kay <- plotfun(estmat[5, ], AVGmat[5, ], SDmat[5, ], namevec[5], type = "protein")
ell <- plotfun(estmat[6, ], AVGmat[6, ], SDmat[6, ], namevec[6], type = "protein")

grid.arrange(g, h, eye, jay, kay, ell, ncol=2)


# Critique this approach too... small standard deviations overrule everything

############################
## Fit for entire RNA set ##
############################

# Fit for just 10; show advantages of foreach package
x1 <- matrix(nrow = 10, ncol = 7)

# regular for-loop
system.time(for (i in 1:10) {
	AVG <- rna.av[i, -1]
	SD  <- rna.sd[i, -1]
	
	A1 <- runif(1); A2 <- runif(1); A3 <- runif(1)
	
	optim.i <- optim(c(A1, A2, A3, tau1, tau2, tau3, tau4), 
	costfun, 
	method = "Nelder-Mead", 
	control = list(maxit = 1000))
	
	x1[i, ] <- optim.i$par
} )
#   user  system elapsed
# 24.737   0.223  25.119


# foreach
system.time (x2 <- 
	foreach(i=1:10, .combine='rbind', .inorder=TRUE) %dopar% {
	AVG <- rna.av[i, -1]
	SD  <- rna.sd[i, -1]
	
	A1 <- runif(1); A2 <- runif(1); A3 <- runif(1)
	
	optim.i <- optim(c(A1, A2, A3, tau1, tau2, tau3, tau4), 
	costfun, 
	method = "Nelder-Mead", 
	control = list(maxit = 1000))
	
	# parameter.matrix.example[i, ] <- optim.i$par
	optim.i$par
}
)
#  user  system elapsed
# 0.068   0.020  12.015
# ~ 40 minutes for all RNAs

#########
### Now do it for all RNAs & proteins. This will take a while...
#########


parameters.rna <- 
	foreach(i=1:nrow(rna.av), .combine='rbind', .inorder=TRUE) %dopar% {
	AVG <- rna.av[i, -1]
	SD  <- rna.sd[i, -1]
	
	A1 <- runif(1); A2 <- runif(1); A3 <- runif(1)
	
	optim.i <- optim(c(A1, A2, A3, tau1, tau2, tau3, tau4), 
	costfun, 
	method = "Nelder-Mead", 
	control = list(maxit = 1000))
	
	# parameter.matrix.example[i, ] <- optim.i$par
	optim.i$par
}

parameters.rna <- matrix(parameters.rna, nrow = nrow(parameters.rna))

rna.fit <- data.frame(cbind(rna.av[, 1], round(parameters.rna, 6)))
names(rna.fit) <- c("gene_id", "A1", "A2", "A3", "tau1", "tau2", "tau3", "tau4")
write.csv(file = "Results/rna_parameters.csv", rna.fit, row.names = FALSE)



##
# Proteins
##
parameters.pro <- 
	foreach(i=1:nrow(pro.av), .combine='rbind', .inorder=TRUE) %dopar% {
	AVG <- pro.av[i, -1]
	SD  <- pro.sd[i, -1]
	
	A1 <- runif(1); A2 <- runif(1); A3 <- runif(1)
	
	optim.i <- optim(c(A1, A2, A3, tau1, tau2, tau3, tau4), 
	costfun, 
	method = "Nelder-Mead", 
	control = list(maxit = 1000))
	
	# parameter.matrix.example[i, ] <- optim.i$par
	optim.i$par
}

# Take only numerics
parameters.pro <- matrix(parameters.pro, nrow = nrow(parameters.pro))

pro.fit <- data.frame(cbind(pro.av[, 1], round(parameters.pro, 6)))
names(pro.fit) <- c("gene_id", "A1", "A2", "A3", "tau1", "tau2", "tau3", "tau4")
write.csv(file = "Results/pro_parameters.csv", pro.fit, row.names = FALSE)

##
# RNA
##
 
###########################
#### Classify the fits ####
###########################

entrez.dict <- read.csv("Data/entrez_dictionary.csv", header = TRUE, as.is = TRUE)
dict <- read.csv("Data/dictionary.csv", header = TRUE, as.is = TRUE)
head(entrez.dict)
head(dict)

rna.fit <- read.csv("Results/rna_parameters.csv", header = TRUE, as.is = TRUE)

A1.r <- rna.fit$A1
A2.r <- rna.fit$A2
A3.r <- rna.fit$A3

pro.fit <- read.csv("Results/pro_parameters.csv", header = TRUE, as.is = TRUE)

A1.p <- pro.fit$A1
A2.p <- pro.fit$A2
A3.p <- pro.fit$A3

# Tolerence level
err <- 0.10

cond1.r  <- {(A3.r - err > A2.r) & (A2.r + err > A1.r)}
cond2.r  <- {(A2.r - err > A1.r) & (A3.r + err > A2.r)}
cond2B.r <- {(A3.r - err > A1.r) & (A2.r + err > A1.r) & (A2.r - err < A1.r)}

cond3.r  <- {(A3.r + err < A2.r) & (A2.r < A1.r + err)}
cond4.r  <- {(A2.r + err < A1.r) & (A3.r < A2.r + err)}
cond4B.r <- {(A3.r + err < A1.r) & (A2.r + err < A1.r) & (A2.r - err > A1.r)}

rna.up <- rna.fit[cond1.r | cond2.r | cond2B.r , 1]
rna.do <- rna.fit[cond3.r | cond4.r | cond4B.r , 1]

sum(rna.up %in% rna.do)

#
rna.up.entrez <- c()
for (name in rna.up) {
	short <- dict$short.name[which(dict$mRNA_ID == name)]
	new.entrez <- entrez.dict$Entrez_GeneID[which(entrez.dict$Symbol == short)]
	
	rna.up.entrez <- c(rna.up.entrez, new.entrez)
}

# 
rna.do.entrez <- c()
for (name in rna.do) {
	short <- dict$short.name[which(dict$mRNA_ID == name)]
	new.entrez <- entrez.dict$Entrez_GeneID[which(entrez.dict$Symbol == short)]
	
	rna.do.entrez <- c(rna.do.entrez, new.entrez)
}


cond1.p  <- {(A3.p - err > A2.p) & (A2.p + err > A1.p)}
cond2.p  <- {(A2.p - err > A1.p) & (A3.p + err > A2.p)}
cond2B.p <- {(A3.p - err > A1.p) & (A2.p + err > A1.p) & (A2.p - err < A1.p)}

cond3.p  <- {(A3.p + err < A2.p) & (A2.p < A1.p + err)}
cond4.p  <- {(A2.p + err < A1.p) & (A3.p < A2.p + err)}
cond4B.p <- {(A3.p + err < A1.p) & (A2.p + err < A1.p) & (A2.p - err > A1.p)}

pro.up <- pro.fit[cond1.p | cond2.p | cond2B.p , 1]
pro.do <- pro.fit[cond3.p | cond4.p | cond4B.p , 1]

sum(pro.up %in% pro.do)

pro.up.entrez <- c()
for (name in pro.up) {
	short <- dict$short.name[which(dict$Protein_id == name)]
	new.entrez <- entrez.dict$Entrez_GeneID[which(entrez.dict$Symbol == short)]
	
	pro.up.entrez <- c(pro.up.entrez, new.entrez)
}

# 
pro.do.entrez <- c()
for (name in pro.do) {
	short <- dict$short.name[which(dict$Protein_id == name)]
	new.entrez <- entrez.dict$Entrez_GeneID[which(entrez.dict$Symbol == short)]
	
	pro.do.entrez <- c(pro.do.entrez, new.entrez)
}

##########################
#### Query DAVID API ####
##########################

library(RDAVIDWebService)

# Add link to sign up
david1 <- DAVIDWebService(email="spencer.woody@utexas.edu", 
url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")

david2 <- DAVIDWebService(email="spencer.woody@utexas.edu", 
url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")

david3 <- DAVIDWebService(email="spencer.woody@utexas.edu", 
url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")

david4 <- DAVIDWebService(email="spencer.woody@utexas.edu", 
url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")

# Add lists to each DAVID
result <- addList(david1, rna.up.entrez, idType = "ENTREZ_GENE_ID", 
listName = "rna.up", listType = "Gene")

result <- addList(david2, rna.do.entrez, idType = "ENTREZ_GENE_ID", 
listName = "rna.do", listType = "Gene")

result <- addList(david3, pro.up.entrez, idType = "ENTREZ_GENE_ID", 
listName = "rna.up", listType = "Gene")

result <- addList(david4, pro.do.entrez, idType = "ENTREZ_GENE_ID", 
listName = "rna.up", listType = "Gene")

setAnnotationCategories(david1, "GOTERM_BP_FAT")
setAnnotationCategories(david2, "GOTERM_BP_FAT")
setAnnotationCategories(david3, "GOTERM_BP_FAT")
setAnnotationCategories(david4, "GOTERM_BP_FAT")

# Get functional clustering reports saved to R object
FuncAnnotChart.rna.up <- getFunctionalAnnotationChart(david1)
FuncAnnotChart.rna.do <- getFunctionalAnnotationChart(david2)
FuncAnnotChart.pro.up <- getFunctionalAnnotationChart(david3)
FuncAnnotChart.pro.do <- getFunctionalAnnotationChart(david4)

FuncAnnotChart.rna.up$Term
FuncAnnotChart.rna.do$Term
FuncAnnotChart.pro.up$Term
FuncAnnotChart.pro.do$Term

# Save functional clustering reports to .tsv files
getFunctionalAnnotationChartFile(david1, "Results/FuncAnnotChart_rna_up.tsv")
getFunctionalAnnotationChartFile(david2, "Results/FuncAnnotChart_rna_do.tsv")
getFunctionalAnnotationChartFile(david3, "Results/FuncAnnotChart_pro_up.tsv")
getFunctionalAnnotationChartFile(david4, "Results/FuncAnnotChart_pro_do.tsv")
