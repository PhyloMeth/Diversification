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
	deltaAIC.vector <- AIC.vector-min(AIC.vector)
	#AkaikeWeight.vector <- ______________
	result.list <- list(bd.result, yule.result, ddl.result, AIC.vector, deltaAIC.vector)#, AkaikeWeight.vector)
	return(result.list)
}