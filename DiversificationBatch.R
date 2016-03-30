#You can use code you wrote for the correlation exercise here.
source("DiversificationFunctions.R")
#setwd("~/Desktop/2016Spring/Phylometh/Diversification")
tree <- read.tree("Euc.final.dated.noNit.tre")
#First, let's look at a sister group comparison. Imagine you have one clade you think is especially noteworthy. 

library(caper)
ntax.focal.clade <- length(clade.members(30,tree)) # subgenus eucalyptus
ntax.sister.clade <- length(clade.members(42,tree)) # subgenus symphyomyrtus
depth.both <- branching.times(tree)[1] #time of the MRCA
actual.ratio <- min(c(ntax.focal.clade, ntax.sister.clade)) / max(c(ntax.focal.clade, ntax.sister.clade)) 

estimated.div.rate <- log(ntax.focal.clade + ntax.sister.clade)/depth.both # N(t) = N0 * exp(r*t)

# Simulate diversification for left and right clades with our estimated birth rate and the age of the tree:
nsim <- 10000
sim.ratios <- rep(NA, nsim)
for (i in sequence(nsim)) {
	left.clade <- sim.bd(b=estimated.div.rate, times=depth.both)[2,2] #get the number of taxa. We're assuming a pure birth model. This is dumb: if there's one thing we know about life, it's that extinction happens. But it's convenient for this case. This is known as a Yule model.
	right.clade <- sim.bd(b=estimated.div.rate, times=depth.both)[2,2] 
	sim.ratios[i] <- min(c(left.clade, right.clade)) / max(c(left.clade, right.clade))
	if(i%%500==0) {
		print(paste("Now", 100*i/nsim, "percent done"))	
	}
}

hist(sim.ratios, breaks=100, col="black", main=paste("Fraction of simulations with more disparity is", ecdf(sim.ratios)(actual.ratio)))
abline(v=actual.ratio, col="red")

#So, what does this mean about your observed result? What's the p-value?
# This means that over 90% of our simulations yielded clades that had more dissimilar numbers of taxa than our actual phylogeny. p-value is ecdf(sim.ratios)(actual.ratio)

#Now, try fitting different models for diversification.
div.results <- TryMultipleDivModels(tree)

best.model <- names(which(div.results[[5]][c(1,2,3)]==0))

# What are the parameters of the best model? What do you think they mean?
# the best model is the birth-death model:

paste("The maximum likelihood is", round(div.results[[2]]$LH, 2), "and the AIC score is", round(div.results[[2]]$aic, 2)) 
paste("The diversification rate (b-d) is", div.results[[2]]$r1, "per arbitrary unit of time")
paste("The extinction fraction (d/b) is", div.results[[2]]$a)


# Now try running BAMM. Use the tutorial at http://bamm-project.org/quickstart.html to do diversification analyses.


# check to see if the tree is ultrametric, fully bifurcating, and branch lengths greater than 0:
is.ultrametric(tree)
is.binary.tree(tree)
min(tree$edge.length)

# generate prior block:
library(BAMMtools)
priors <- setBAMMpriors(tree) # copy and paste into control file

# account for incomplete taxon sampling
source("isItMonophyletic.R") # the isItMonophyletic function determines whether your phylogeny is monophyletic compared to another bigger phylogeny. I will test my phylogeny against Zanne et al.'s:
zanneTree <- read.tree("bigtree.tre")
isItMonophyletic(tree, zanneTree) # there are 110 missing species according to Zanne et al.'s tree
samplingFraction <- 28/(28+110) # set globalSamplingFraction equal to 0.2028986 in            control file 

# download bamm file in current working directory and type the following in terminal:
tar -xzf bamm-2.5.0-MacOSX.tar.gz

# to run the control file:
./bamm -c myControlFile.txt




#### BAMM OUTPUT ####


	# import event_data.txt in to R
	edata <- getEventData(tree, eventdata = "event_data.txt", burnin=0.1)
	
	# look for convergence:
	mcmcout <- read.csv("mcmc_out.txt", header=T)
	plot(mcmcout$logLik ~ mcmcout$generation)
	
	# take out burn-in
	burnstart <- floor(0.1 * nrow(mcmcout))
	postburn <- mcmcout[burnstart:nrow(mcmcout), ]

	# and look at effective sample size:
	library(coda)
	effectiveSize(postburn$N_shifts)
	effectiveSize(postburn$logLik)
	# (In general, we want these to be at least 200 (and 200 is on the low side, but might be reasonable for very large datasets)

	# Posterior probs of models (The probability of model ‘0’ is the posterior probability of a model with just a single evolutionary rate dynamic, or no rate shifts):
	shift_probs <- summary(edata)

	# Bayes factors (many researchers consider values greater than 12 to be consistent with at least some effect)
	mcmc <- read.csv("mcmc_out.txt", stringsAsFactors = FALSE)
	expectedNumberOfShifts=1
	burnin=0.1
	
	# Calculate prior
	mcmc2 <- mcmc[floor(burnin * nrow(mcmc)):nrow(mcmc), ]; obsK <- seq(from = 0, to = max(mcmc2[, "N_shifts"]), by = 1)
	prob.k <- function(k, poissonRatePrior=1)	{
		Denom <- (poissonRatePrior + 1)^(k+1)
		Prob <- poissonRatePrior / Denom
		return(Prob)
	}; prior <- sapply(obsK, prob.k, poissonRatePrior = 1/expectedNumberOfShifts); prior <- data.frame(N_shifts = obsK, prob = prior)
    
    # Calculate posterior:
    posterior <- sapply(obsK, function(x) length(which(mcmc2[,"N_shifts"] == x)))/nrow(mcmc2); names(posterior) <- obsK; posterior <- data.frame(N_shifts = names(posterior), prob = posterior); posterior$prob[which(posterior$prob==0)] <- 1/1000

	# Bayes Factors relative to the null model (0 shifts)
	bayesFactors <- data.frame(N_shifts = obsK, ratios=(posterior$prob/prior$prob)/posterior$prob[1]/prior$prob[1])
	
	
	# Mean rate plot (model-averaged diversification rates)
	plot.bammdata(edata, lwd=2, legend=T)

	# identify the configurations with highest probabilities:
	css <- credibleShiftSet(edata, expectedNumberOfShifts=1, threshold=5, set.limit = 0.95)
	summary(css) # one distinct configuration
	plot.credibleshiftset(css)

	# extract the shift configuration that maximizes the marginal probability of rate shifts along individual branches
	msc.set <- maximumShiftCredibility(edata, maximize='product')
	msc.config <- subsetEventData(edata, index = msc.set$sampleindex)
	plot.bammdata(msc.config, lwd=2)
	addBAMMshifts(msc.config, cex = 2)

	# get overall mean rates 
	allrates <- getCladeRates(edata)
	mean(allrates$lambda) # overall speciation ra
	quantile(allrates$lambda, c(0.05, 0.95))

	# get mean rates for subgenus
	eucRates <- getCladeRates(edata, node=30)
	mean(eucRates $lambda)
	symphRates <- getCladeRates(edata, node=42)
	mean(symphRates $lambda)
	
	# rate through time plot
	plotRateThroughTime(edata, ratetype="speciation")

	# Macroevolutionary cohort analysis
	cmat <- getCohortMatrix(edata)
	cohorts(cmat, edata)


