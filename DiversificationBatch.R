setwd("~/Desktop/PhyloMeth/Diversification")
library(ape)
library(geiger) 
library(laser)
library(phytools)
library(rotl)

#You can use code you wrote for the correlation exercise here.
source("DiversificationFunctions.R")
#primate_id <- tnrs_match_names("Primates")
#primate_tree <- tol_subtree(ott_id = primate_id$ott_id[1])

primate.studies <- studies_find_studies(property="ot:focalCladeOTTTaxonName", value="Primates")
primate.ids<- unlist(primate.studies$study_ids)

primate.meta <- get_study_meta("pg_2656")
get_publication(primate.meta)
get_tree_ids(primate.meta)
candidate_for_synth(primate.meta)
#primate.dat<-get_study("pg_2656")
#summary(primate.dat)
#names(primate.dat)

primate.tree<- get_study_tree("pg_2656", "tree6185")
#time-calibrated,autocorrelated rates, soft-bounded constraints
print(primate.tree)
plot(primate.tree)



#First, let's look at a sister group comparison. Imagine you have one clade you think is especially noteworthy. 

ntax.focal.clade <- 25
ntax.sister.clade <- 113
depth.both <- 25.1 #time of the MRCA
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
div.results <- TryMultipleDivModels(primate.tree) #birth death model

best.model <- div.results[div.results[[5]]==0][[1]]

# What are the parameters of the best model? What do you think they mean?

LH: the log-likelihood at the maximum
r: the net diversification rate giving the maximum log-likelihood
a: the extinction fraction giving the maximum log-likelihood


# Now try running BAMM. Use the tutorial at http://bamm-project.org/quickstart.html to do diversification analyses.

