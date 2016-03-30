library(ape)
library(geiger) 
library(laser) # Likelihood analysis of speciation/extinction rates from phylogenies; contains getBtimes, pureBirth, bd functions
library(phytools)
library(MuMIn)

TryMultipleDivModels <- function(tree) {
	# obtain branching times ordered from oldest to youngest: 
	tree.branching <- getBtimes(string=write.tree(tree)) 
	
	# fit pure birth (Yule) model to the branching times:
	yule.result <- pureBirth(tree.branching) # yields LH (max log-likelihood), aic, and r1 (speciation rate giving max likelihood)
	
	# fit birth-death model to the branching times:
	bd.result <- bd(tree.branching) # yields LH (max log-likelihood), aic, r (diversification rate giving max likelihood), and a (extinction fraction giving max likelihood)
	
	# fit logistic density-dependent speciation rate model
	ddl.result <- DDL(tree.branching) # yields LH (max log-likelihood), aic, and r1 (initial speciation rate), and kparam (carrying capacity)
	
	AIC.vector <- c(yule.result$aic, bd.result$aic, ddl.result$aic)
	names(AIC.vector) <- c("yule", "bd", "ddl")
	deltaAIC.vector <- min(AIC.vector)-AIC.vector
	names(deltaAIC.vector) <- c("yule", "bd", "ddl")
	AkaikeWeight.vector <- Weights(AIC.vector)
	names(AkaikeWeight.vector) <- c("yule", "bd", "ddl")
	result.list <- list(yule.result, bd.result, ddl.result, AIC.vector, deltaAIC.vector, AkaikeWeight.vector)
	return(result.list)
}