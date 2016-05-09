library(ape)
library(geiger) 
library(laser)
library(phytools)
library(BAMMtools)

TryMultipleDivModels <- function(tree) {
	tree.branching <- getBtimes(string=write.tree(tree))
	yule.result <- pureBirth(tree.branching)
	bd.result <- bd(tree.branching)
	ddl.result <- DDL(tree.branching)
	AIC.vector <- c(yule.result$aic, bd.result$aic, ddl.result$aic)
	deltaAIC.vector <- AIC.vector-(min(AIC.vector)) 
	rel.likelihood.vector <- exp(-.5 * deltaAIC.vector)
	AkaikeWeight.vector <- rel.likelihood.vector/sum(rel.likelihood.vector)
	result.list <- list(yule=yule.result, bd=bd.result, ddl=ddl.result, AIC=AIC.vector, deltAIC=deltaAIC.vector, Akaike=AkaikeWeight.vector)
	return(result.list)
}
