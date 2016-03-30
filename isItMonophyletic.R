
isItMonophyletic <- function(myTree, zanneTree){
	tipsInZanneTree <- which(zanneTree$tip.label %in% myTree$tip.label)
	prunedZanneTree <- drop.tip(zanneTree, seq(1,length(zanneTree$tip.label))[-tipsInZanneTree])
	tipsInYourTree <- which(myTree$tip.label %in% zanneTree$tip.label)
	print(paste(length(prunedZanneTree$tip.label), "of your species are in Zanne's phylogeny. Following is the list of your species that are not in Zanne's phylogeny:"))
	print(myTree$tip.label[-tipsInYourTree])
	
	distMat <- cophenetic(prunedZanneTree); rownames(distMat) <- prunedZanneTree$tip.label; colnames(distMat) <- prunedZanneTree$tip.label

	col=NA
	for(i in 1:dim(distMat)[1]){
		index <- which(distMat==max(distMat))[1] - (18*i)
		if( index < dim(distMat)[1] & index > 0){
		col=i
		}
	}
	row=NA
	for(i in 1:dim(distMat)[1]){
		index <- which(distMat==max(distMat))[1]
		if( index %in% c(i,i+seq(18,324,18))){
		row=i
		}
	}	
	tip1 <- rownames(distMat)[row]
	tip2 <- colnames(distMat)[col]
	
		mrcA <- fastMRCA(zanneTree, tip1, tip2)
	
	getDescendants<-function(tree,node,curr=NULL){
		if(is.null(curr)) curr<-vector()
		daughters<-tree$edge[which(tree$edge[,1]==node),2]
		curr<-c(curr,daughters)
		w<-which(daughters>=length(tree$tip))
		if(length(w)>0) for(i in 1:length(w)) 
		curr<-getDescendants(tree,daughters[w[i]],curr)
		return(curr)
		}
	descendantS <- na.omit(zanneTree$tip.label[getDescendants(zanneTree, mrcA)])
	if(length(descendantS) > length(prunedZanneTree$tip.label)){
		print(paste("Your tree is not monophyletic. The Zanne tree has found", length(descendantS)-length(prunedZanneTree$tip.label), "more species that are descendants of the common ancestor of your tree."))
	}
	if(length(descendantS) == length(prunedZanneTree$tip.label)){
		print(paste("Based on Zanne's tree, your tree is monophyletic. Congratulations!"))
	}
}
