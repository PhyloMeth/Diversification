#You can use code you wrote for the correlation exercise here.
setwd("~/GitHub/Diversification")
source("DiversificationFunctions.R")

set.seed(123)
tree<-pbtree(n=20,scale=1,tip.label=(c("T.ochroleucum","T.palmeri","T.calocephalum","T.nanum","T.pallescens","T.repens","T.amabile","T.decorum","T.wigginsii","T.virginicum","T.vesiculosum","T.vernum","T.uniflorum","T.tomentosa","T.thalii","T.sylvatica","T.alpestre","T.argentinense","T.pratense","T.hybridum")))
plotTree(tree)
plot(nodelabels)
nodelabels()
#nodelabels look like a mess.  Not sure how to clean them up

#tree <- read.tree("____PATH_TO_TREE_OR_SOME_OTHER_WAY_OF_GETTING_A_TREE____")

#First, let's look at a sister group comparison. Imagine you have one clade you think is especially noteworthy. 

#assign these by node
ntax.focal.clade <- 24
ntax.sister.clade <- 30

depth.both <- 23 #time of the MRCA
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


#I had to seek out help on this function
#I still don't think this is really working correctly
best.model <- names(which(div.results[[5]][c(1,2,3)]==0))
#I can go into the results here and see that it is the 'bd' model
#it just doesn't display correctly



# What are the parameters of the best model? What do you think they mean?
paste("Max likelihood =", round(div.results[[2]]$LH, 2), "AIC", round(div.results[[2]]$aic, 2)) 
#LH: log-likelihood at maximum
paste("r =", div.results[[2]]$r1)
#r: net diversification rate giving the maximum log-likelihood
paste("a =", div.results[[2]]$a)
#a: extinction fraction giving the maximum log-likelihood



_____________________________________

# Now try running BAMM. Use the tutorial at http://bamm-project.org/quickstart.html to do diversification analyses.