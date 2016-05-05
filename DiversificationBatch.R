#You can use code you wrote for the correlation exercise here.
source("/Users/oschwery/Documents/UTK-Orlando/Coursework/Phylometh/Diversification/DiversificationFunctions.R")
#Desktop
tree <- read.tree("/Users/oschwery/Documents/UTK-Orlando/MiscSync/dateddung2.tre")
#Laptop
tree <- read.tree("/Users/orlandoschwery/Documents/UTK/MiscSync/dateddung2.tre")

#First, let's look at a sister group comparison. Imagine you have one clade you think is especially noteworthy.
names1 <- tree$tip.label[c((getDescendants(tree, getMRCA(tree, tip=c("Macroderes_sp._dgi1", "Coprophanaeus_sp._dgi1")))))]
names2 <- tree$tip.label[c((getDescendants(tree, getMRCA(tree, tip=c("Demarziella_imitatrix", "Lepanus_glaber")))))]
ntax.focal.clade <- length(names1[!is.na(names1)])
ntax.sister.clade <- length(names2[!is.na(names2)])
depth.both <- branching.times(tree)[as.character(getMRCA(tree, tip=c("Macroderes_sp._dgi1","Demarziella_imitatrix")))] #time of the MRCA
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

# Fraction with more disparity was 0.3675. Actual ratio is 0.2307692. This means that under the bd model and assuming the same diversification rate, we would expect them to have a lower difference in a bit less than 2/3 of the cases and a higher difference in a bit more than 1/3 of the cases.
# We want to know whether the difference is significantly deviating from 1 (which means there's no difference between the clades). Since this is a probability density graph, we can see it as a one-sided test, and the p value is the ratio of cases that are further away from the null hypothesis than our actual data. This in turn is the area under the curve to the left of the red line, so 0.37 is our p-value and we have to assume that the difference between our clades is not significantly different from assuming 'both are equally big'. Given BD is the correct model of course.

#Now, try fitting different models for diversification.
div.results <- TryMultipleDivModels(tree)

best.model <- div.results[which(div.results$deltaAIC.vector %in% min(div.results$deltaAIC.vector))]

# What are the parameters of the best model? What do you think they mean?

# according to the AIC and deltaAIC, BD is the best model of the ones we tested. It gives us a net diversification rate r of 6.894714e-08 and an extinction fracton a of 1. r should be lambda-mu, and a shold be mu/lambda. (this doesn't really make sense to me, becuase if a is 1, then mu=lambda and then r should be 0... It is pretty low though, so I presume it's a rounding thing.)

# Now try running BAMM. Use the tutorial at http://bamm-project.org/quickstart.html to do diversification analyses.
