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
	deltaAIC.vector <- c((AIC.vector[1]-(min(AIC.vector))), (AIC.vector[2]-(min(AIC.vector))), (AIC.vector[3]-(min(AIC.vector))))
	rel.likelihood.vector <- c(exp(-0.5*deltaAIC.vector[1]), exp(-0.5*deltaAIC.vector[2]), exp(-0.5*deltaAIC.vector[3]))
	AkaikeWeight.vector <- c(rel.likelihood.vector[1]/sum(rel.likelihood.vector), rel.likelihood.vector[2]/sum(rel.likelihood.vector), rel.likelihood.vector[3]/sum(rel.likelihood.vector))
	result.list <- list(yule.result=yule.result, bd.result=bd.result, ddl.result=ddl.result, AIC.vector=AIC.vector, deltaAIC.vector=deltaAIC.vector, AkaikeWeight.vector=AkaikeWeight.vector)
	return(result.list)
}
