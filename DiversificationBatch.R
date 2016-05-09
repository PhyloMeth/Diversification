#You can use code you wrote for the correlation exercise here.
source("DiversificationFunctions.R")
# tree <- read.tree("____PATH_TO_TREE_OR_SOME_OTHER_WAY_OF_GETTING_A_TREE____")
data(whales)
tree <- whales
#First, let's look at a sister group comparison. Imagine you have one clade you think is especially noteworthy. 

ntax.focal.clade <- 35
ntax.sister.clade <- 8
depth.both <- branching.times(tree)[as.character(getMRCA(tree, c("Orcinus_orca", "Phocoena_sinus")))] #time of the MRCA
actual.ratio <- min(c(ntax.focal.clade, ntax.sister.clade)) / max(c(ntax.focal.clade, ntax.sister.clade))

estimated.div.rate <- log(ntax.focal.clade + ntax.sister.clade)/depth.both #N(t) = N0 * exp(r*t)

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

#Now, try fitting different models for diversification.
div.results <- TryMultipleDivModels(tree)

best.model <- (result.list$AIC)
names(best.model) <- c("yule", "bd", "ddl")

best.model
# What are the parameters of the best model? What do you think they mean?


_____________________________________

# Now try running BAMM. Use the tutorial at http://bamm-project.org/quickstart.html to do diversification analyses.

