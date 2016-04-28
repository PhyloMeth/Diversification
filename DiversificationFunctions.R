library(ape)
library(geiger) 
library(corHMM)
library(phytools)
library(phangorn)

#You can use code you wrote for the correlation exercise here.


MergeData <- function(data1,data2) {
	merged.data<-merge(data1,data2,all=FALSE)
	merged.data.withrows<-as.data.frame(merged.data, row.names=merged.data[,1])
	merged.data.cleaned<-merged.data.withrows[,-1]
	print(merged.data.cleaned)
}

VisualizeData <- function(phy, data) {
	plot.phylo(phy)
	print(data)
	#Important here is to LOOK at your data before running it. Any weird values? Does it all make sense? What about your tree? Polytomies?
}


CleanData <- function(phy, data, sort) {
	phy.data<-treedata(phy,data,sort=sort)
	print(phy.data)
	#treedata() in Geiger is probably my favorite function in R.
}


