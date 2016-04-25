library(ape)
library(geiger) 
library(laser)
library(phytools)

TryMultipleDivModels <- function(tree) {
	tree.branching <- getBtimes(string=write.tree(tree))
	yule.result <- pureBirth(tree.branching)
	bd.result <- bd(tree.branching)
	ddl.result <- DDL(tree.branching)
	AIC.vector <- c(yule.result$aic, bd.result$aic, ddl.result$aic)
	deltaAIC.vector <- ________________
	rel.likelihood.vector <- 
	AkaikeWeight.vector <- ______________
	result.list <- list(bd=bd.result, yule=yule.result, ddl=ddl.result, AIC=AIC.vector, deltAIC=deltaAIC.vector, Akaike=AkaikeWeight.vector)
	return(result.list)
}