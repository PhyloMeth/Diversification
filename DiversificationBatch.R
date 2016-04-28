#You can use code you wrote for the correlation exercise here.
source("~/phyloMeth.packages/Diversification/DiversificationFunctions.R")
tree <- read.tree("~/R/Vascular_Plants_rooted.dated.tre")
ploidy.data <- read.csv(file="~/Documents/Data/ploidy.data.txt", stringsAsFactors=FALSE, row.names=1) #death to factors.

cleaned.discrete <- treedata(tree, ploidy.data, sort=TRUE)

VisualizeData(tree, cleaned.discrete)
dat <- cleaned.discrete$data[,1]
#First, let's use parsimony to look at ancestral states
cleaned.discrete.phyDat <- phyDat(data=dat, type="USER", levels=c("2", "3", "4", "6", "8", "12")) #phyDat is a data format used by phangorn

anc.p <- ancestral.pars(cleaned.discrete$phy, cleaned.discrete.phyDat, type="ACCTRAN")
plotAnc(tree=cleaned.discrete$phy, data=anc.p, i=1, cex.pie=0.05, cex=0.08, edge.width=0.4, pos = "topleft")

#Do you see any uncertainty? What does that meean for parsimony?

#now plot the likelihood reconstruction
anc.ml <- ancestral.pml(pml(cleaned.discrete$phy, cleaned.discrete.phyDat), type="ml")
plotAnc(cleaned.discrete$phy, anc.ml, 1)

#How does this differ from parsimony? 
#Why does it differ from parsimony?
#What does uncertainty mean?

#How many changes are there in your trait under parsimony? 
parsimony.score <- parsimony(cleaned.discrete$phy, cleaned.discrete.phyDat)
print(parsimony.score)

#Can you estimate the number of changes under a likelihood-based model? 

#Well, we could look at branches where the reconstructed state changed from one end to the other. But that's not really a great approach: at best, it will underestimate the number of changes (we could have a change on a branch, then a change back, for example). A better approach is to use stochastic character mapping.

estimated.histories <- make.simmap(cleaned.discrete$phy, cleaned.discrete$data[,1], model="ARD", nsim=5)

#always look to see if it seems reasonable
plotSimmap(estimated.histories)

counts <- countSimmap(estimated.histories)
print(counts)

#Depending on your biological question, investigate additional approaches:
#  As in the correlation week, where hypotheses were examined by constraining rate matrices, one can constrain rates to examine hypotheses. corHMM, ape, and other packages have ways to address this.
#  Rates change over time, and this could be relevant to a biological question: have rates sped up post KT, for example. Look at the models in geiger for ways to do this.
#  You might observe rates for one trait but it could be affected by some other trait: you only evolve wings once on land, for example. corHMM can help investigate this.

