#BAMM script

mcmc <- read.csv("C:\\Users\\Orlando\\Documents\\UZH\\Master\\Publication\\Analyses\\BAMM\\Ericac Results/Ericac_450_gen_tribe_200mio_1000/posterior_mcmc_out.txt", header=T)
dim(mcmc)

head(mcmc)

###################################################
######    Diagnosing convergence:            ######
###################################################

###################################################
#plot log likelihood trace

windows()
plot(mcmc$logLik ~ mcmc$generation)

windows()
plot(mcmc$N_shifts ~ mcmc$generation)

##################################################
#Discard some as burnin

burnstart <- floor(0.1 * nrow(mcmc))
mcmcpost <- mcmc[burnstart:nrow(mcmc), ]


##
library(coda)

windows()
plot(mcmc$logLik, type='l')

effectiveSize(mcmcpost$logLik)  ## calculates autocorrelation function
# effective size on N_shifts
effectiveSize(mcmcpost$N_shift)


## do this on pre-burnin MCMC

effectiveSize(mcmc$logLik)

# effective size on N_shifts
effectiveSize(mcmc$N_shift)

# compute posterior probabilities
post <- table(mcmc$N_shifts) / nrow(mcmcpost)
prior <- read.csv("C:\\Users\\Orlando\\Documents\\UZH\\Master\\Radiation Meeting Zurich 2014\\BAMM_Workshop\\Saxifragales\\prior_probs.txt")


###########
mcmc1chain <- read.csv("C:\\Users\\Orlando\\Documents\\UZH\\Master\\Radiation Meeting Zurich 2014\\BAMM_Workshop\\bamm-2.0.0\\metropolis coupled\\mcmc_out.txt")
mcmc1chainpost <- mcmc1chain[[100:1000,]

####################
library(BAMMtools)
tree <- read.tree("C:\\Users\\Orlando\\Documents\\UZH\\Master\\Publication\\Analyses\\BAMM\\Ericac_MCC_2-2_mean_fixed_450.tre")
ed <- getEventData(tree, "C:\\Users\\Orlando\\Documents\\UZH\\Master\\Publication\\Analyses\\BAMM\\Ericac Results/Ericac_450_gen_tribe_200mio_1000/event_data.txt", burnin=0.1, nsamples=5000)
#subsample event data file (good to get idea on how long it will take to get it in)


class(ed)
names(ed)

windows()
plot.bammdata(ed, legend=T)

windows()
plot.bammdata(ed, legend=T, spex='ex')

windows()
plot.bammdata(ed, legend=T, spex='se')

############
# macroevolutionary cohort analysis

cmat <- getCohortMatrix(ed)
dim(cmat)
 cmat[1:5, 1:5]

windows()
cohorts(cmat, ed)


#### this doesn't work, different output files########################3
prior <-  getBranchShiftPriors(tree, "C:\\Users\\Orlando\\Documents\\UZH\\Master\\Publication\\Analyses\\BAMM\\Ericac Results/Ericac_450_gen_tribe_200mio_1000\\prior_probs.txt")
windows()
plot(tree$edge.lenght,prior$edge.length)
# compute credible shift set (95%)
css <- credibleShiftSet(ed, prior)
summary(css)
windows()
plot.credibleshiftset(css)
windows()
plot(css)
####################################################################3




#############################################################################################
#looking at marginal shift probabilities on indiv. branches (marg.prob. that branch contains shift)
marg_probs <- marginalShiftProbsTree(ed) #copy of our tree, where branch lenghts will be transformed into marg.shift prob.
#quartz()
windows()
plot.phylo(marg_probs)
#cumulative shift probability for each branch (that shift occured somewhere between focal branch and root)

#OR in color; red = cumulative shift probability 0.95
cst <- cumulativeShiftProbsTree(ed)
edgecols <- rep('black', length(tree$edge.length))
is_highprobshift <- cst$edge.length >= 0.95
edgecols[ is_highprobshift ] <- "red"

#quartz()
#plot.phylo(mytree, edge.color = edgecols, cex=0.3)

quartz.options(height=40, width=10, dpi=72)
quartz()
windows()
plot.phylo(tree, edge.color = edgecols, edge.width=3, cex=0.6)

##OR maximum credibility shift configuration
#msc <- maximumShiftCredibility(edata) # BUNCH OF INFO TO BE ADDED
##select single representative and plot
#samp <- msc$sampleindex
#shiftnodes <- getShiftNodesFromIndex(edata, samp)
#quartz()
#plot(mytree)
#nodelabels(mytree, node = shiftnodes, cex=0.3, bg="red", pch=21)

### RATHER USE jetztree-script:

# Compute the marginal shift probabilities on branches
marg_tree <- marginalShiftProbsTree(ed);

# Compute the cumulative shift probabilities for branches:
cst <- cumulativeShiftProbsTree(ed);

# Compute the maximum shift credibility configuration:
# This is the joint distribution of shift events that maximizes the
#	marginal probability of the data.
msc <- maximumShiftCredibility(ed);

# get the relevant rate shift nodes from the
#	maximum shift credibility configuration
nodes <- getShiftNodesFromIndex(edata, index = msc$sampleindex);

######### get scaling and colors for shift nodes based on
#			their importance

library(gplots);
pal <- rich.colors(10);


cexmin <- 2;
cexmax <- 6;



# Get marginal probabilities associated with each node in the the
#		maximum shift credibility configuration
probvec <- numeric(length(nodes));

for (i in 1:length(nodes)){
	probvec[i] <- marg_probs$edge.length[marg_probs$edge[,2] == nodes[i] ];
}
size <- cexmin + probvec*(cexmax - cexmin);


## Get colors for nodes based on marginal shift probs
colvec <- rep(pal[length(pal)], length(nodes));

cuts <- seq(0, 1, length.out=length(pal)+1);

for (i in 1:(length(cuts) - 1)){
	for (k in 1:length(probvec)){
		if (probvec[k] > cuts[i] & probvec[k] <= cuts[i+1]){
			colvec[k] <- pal[i];
		}
	}
}


lmat <- matrix(c(1,1,1,1,1,1,1, 2), nrow=1);


edgewid <- rep(0.5, length(tree$edge.length));
edgecol <- rep('gray60', length(tree$edge.length));

edgewid[cst$edge.length > 0.9] <- 0.7;
edgecol[cst$edge.length > 0.9] <- pal[10];


## First, the marginal shift probability tree:
quartz.options(height=5, width=17);
#quartz()
windows()
plot.new();
layout(lmat);


plot.phylo(as.phylo.bammdata(ed), show.tip.label=F, edge.width=0.8, edge.color='gray40', no.margin=T, type = 'p', direction='upwards');

nodelabels(node=nodes, cex=size, pch=21, bg=colvec);

plot.new();
par(mar=c(3,0,3,6));
plot.window(xlim=c(0,4), ylim=c(-0.3, 1.3));
for (i in 1:(length(cuts)-1)){
	xco <- c(0,2,2, 0);
	yco <- c(cuts[i], cuts[i], cuts[i+1], cuts[i+1]);
	polygon(x=xco, y=yco, col=pal[i]);
}
axis(side=4,at=seq(0,1, by=0.2), cex.axis=2, font=2, las=1, pos=2.2)
mtext(text='Shift probability', side=4, line=1, cex=1.8)


##############################################
###### THe cumulative shift probability tree:
#
#  Highlight all branches (red) with a probability > 0.95 of
#	having rate dynamics different from the root of the tree.

edgecol <- rep('gray40', length(mytree$edge.length));
edgecol[cst$edge.length > 0.95] <- 'red';

quartz.options(height=5, width=17);
quartz()
plot.new();
layout(lmat);

plot.phylo(as.phylo.bammdata(edata), show.tip.label=F, edge.width=0.8, edge.color=edgecol, no.margin=T, type = 'p', direction='upwards');
nodelabels(node=nodes, cex=size, pch=21, bg=colvec);

plot.new();
par(mar=c(3,0,3,6));
plot.window(xlim=c(0,4), ylim=c(-0.3, 1.3));
for (i in 1:(length(cuts)-1)){
	xco <- c(0,2,2, 0);
	yco <- c(cuts[i], cuts[i], cuts[i+1], cuts[i+1]);
	polygon(x=xco, y=yco, col=pal[i]);
}
axis(side=4,at=seq(0,1, by=0.2), cex.axis=2, font=2, las=1, pos=2.2)
mtext(text='Shift probability', side=4, line=1, cex=1.8)


##########################33
######### clade specific rates

allrates <- getCladeRates(ed)
names(allrates)

mean(allrates$lambda)
quantile(allrates$lambda, c(0.05, 0.95))

mean(allrates$mu)
quantile(allrates$mu, c(0.05, 0.95))


allrates$netdiv<-allrates$lambda-allrates$mu

mean(allrates$netdiv)
quantile(allrates$netdiv, c(0.05, 0.95))


windows()
hist(wrates$lambda, breaks=50) #calculates time integrated average rates (KEEP IN MIND YOU PULL THIS OUT OF IT)

#dolphin's only
plot(tree, show.tip.label=T)
nodelabels()

Rhodorates <- getCladeRates(ed, node = 772)
Vaccrates <- getCladeRates(ed, node = 643)
Ericrates <- getCladeRates(ed, node = 842)
Gaulthrates <- getCladeRates(ed, node = 710)
Richrates <- getCladeRates(ed, node = 588)
Rhodonestrates <- getCladeRates(ed, node = 779)

Rhodorates$netdiv<-Rhodorates$lambda-Rhodorates$mu
Vaccrates$netdiv<-Vaccrates$lambda-Vaccrates$mu
Ericrates$netdiv<-Ericrates$lambda-Ericrates$mu
Gaulthrates$netdiv<-Gaulthrates$lambda-Gaulthrates$mu
Richrates$netdiv<-Richrates$lambda-Richrates$mu
Rhodonestrates$netdiv<-Rhodonestrates$lambda-Rhodonestrates$mu


mean(Rhodorates$lambda)
quantile(Rhodorates$lambda, c(0.05, 0.95))
mean(Vaccrates$lambda)
quantile(Vaccrates$lambda, c(0.05, 0.95))
mean(Ericrates$lambda)
quantile(Ericrates$lambda, c(0.05, 0.95))
mean(Gaulthrates$lambda)
quantile(Gaulthrates$lambda, c(0.05, 0.95))
mean(Richrates$lambda)
quantile(Richrates$lambda, c(0.05, 0.95))
mean(Rhodonestrates$lambda)
quantile(Rhodonestrates$lambda, c(0.05, 0.95))

mean(Rhodorates$mu)
quantile(Rhodorates$mu, c(0.05, 0.95))
mean(Vaccrates$mu)
quantile(Vaccrates$mu, c(0.05, 0.95))
mean(Ericrates$mu)
quantile(Ericrates$mu, c(0.05, 0.95))
mean(Gaulthrates$mu)
quantile(Gaulthrates$mu, c(0.05, 0.95))
mean(Richrates$mu)
quantile(Richrates$mu, c(0.05, 0.95))
mean(Rhodonestrates$mu)
quantile(Rhodonestrates$mu, c(0.05, 0.95))

mean(Rhodorates$netdiv)
quantile(Rhodorates$netdiv, c(0.05, 0.95))
mean(Vaccrates$netdiv)
quantile(Vaccrates$netdiv, c(0.05, 0.95))
mean(Ericrates$netdiv)
quantile(Ericrates$netdiv, c(0.05, 0.95))
mean(Gaulthrates$netdiv)
quantile(Gaulthrates$netdiv, c(0.05, 0.95))
mean(Richrates$netdiv)
quantile(Richrates$netdiv, c(0.05, 0.95))
mean(Rhodonestrates$netdiv)
quantile(Rhodonestrates$netdiv, c(0.05, 0.95))


nonRadrate <- getCladeRates(ed, node = c(772, 643, 842, 710, 588), nodetype = "exclude")
mean(nondolphinrate$lambda)
quantile(nondolphinrate$lambda, c(0.05, 0.95))

nonRadrate1 <- getCladeRates(ed, node = c(772), nodetype = "exclude")
nonRadrate2 <- getCladeRates(nonRadrate1, node = c(643), nodetype = "exclude")
nonRadrate3 <- getCladeRates(ed, node = c(842), nodetype = "exclude")
nonRadrate4 <- getCladeRates(ed, node = c(710), nodetype = "exclude")
nonRadrate5 <- getCladeRates(ed, node = c(588), nodetype = "exclude")


windows()
plot.new()
par(mfrow=c(2,1))
par(mar=c(1,1,1,1))
hist(wrates$lambda, breaks=50, xlim=c(0,0.5))
hist(drates$lambda, breaks=50, xlim=c(0,0.5))

non_drates <- getCladeRates(ed, node = 141, nodetype = 'exclude')
plot.new()
par(mfrow=c(2,1))
par(mar=c(1,1,1,1))
hist(wrates$lambda, breaks=50, xlim=c(0,0.5))
hist(non_drates$lambda, breaks=50, xlim=c(0,0.5))

#####
par(mfrow=c(1,2))
plot(whales, show.tip.label=F)

trates <- getTipRates(ed)
names(trates)
trates$lambda.avg
plot(trates$lambda.avg)   # goes in order of tree tips! (plot toghether to see!

plotRateThroughTime(ed)   #and with excluding and shizzle
