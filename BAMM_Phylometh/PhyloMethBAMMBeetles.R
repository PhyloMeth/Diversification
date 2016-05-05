# BAMM for beetles

library(BAMMtools)

#Convergence
mcmc <- read.csv("/Users/orlandoschwery/Documents/UT/Courses&Workshops/Spring16/Phylometh/Diversification/BAMM_Phylometh/mcmc_out.txt", header=T)  # load file
dim(mcmc)  # check dimensions
head(mcmc)  # check head

quartz()
plot(mcmc$logLik ~ mcmc$generation)
quartz()
plot(mcmc$N_shifts ~ mcmc$generation)

burnstart <- floor(0.1 * nrow(mcmc))
mcmcpost <- mcmc[burnstart:nrow(mcmc), ]

library(coda)
effectiveSize(mcmcpost$logLik)  # calculates autocorrelation function
effectiveSize(mcmcpost$N_shift)  # effective size on N_shifts

## do this on pre-burnin MCMC
effectiveSize(mcmc$logLik)
effectiveSize(mcmc$N_shift)

#Rate Shifts
library(BAMMtools)
tree <- read.tree("/Users/orlandoschwery/Documents/UT/Courses&Workshops/Spring16/Phylometh/Diversification/BAMM_Phylometh/dateddung2.tre")
ed <- getEventData(tree, "/Users/orlandoschwery/Documents/UT/Courses&Workshops/Spring16/Phylometh/Diversification/BAMM_Phylometh/event_data.txt", burnin=0.1)


#Phylorates Plots
quartz()
plot.bammdata(ed, legend=T, spex='netdiv')
quartz()
plot.bammdata(ed, legend=T, spex='s')
quartz()
plot.bammdata(ed, legend=T, spex='e')

summary(ed)

#BayesFactors for Shift Models
postfile <- mcmc  # since we already loaded this
#priorfile <- read.csv("/Users/orlandoschwery/Documents/UT/Courses&Workshops/Spring16/Phylometh/Diversification/BAMM_Phylometh/prior_probs.txt", header=T)
bfmat <- computeBayesFactors(postfile, expectedNumberOfShifts=3, burnin=0.1)

#Credible Shift Sets
#pset <- getBranchShiftPriors(tree, priorfile)
cset <- credibleShiftSet(ed, expectedNumberOfShifts=2, threshold=5, set.limit=0.95)
quartz()
plot.credibleshiftset(cset, lwd=2.5)
cset$number.distinct
summary(cset)
#Best Shift Set
best <- getBestShiftConfiguration(ed, expectedNumberOfShifts=2)
quartz()
plot.bammdata(best, lwd = 2)
addBAMMshifts(best, cex=2.5)

# Maximum Shift Credibility
msc.set <- maximumShiftCredibility(ed, maximize='product')
msc.config <- subsetEventData(ed, index = msc.set$sampleindex)
quartz()
plot.bammdata(msc.config, lwd=2)
addBAMMshifts(msc.config, cex = 2)

#Macroevo Cohort analysis
cmat <- getCohortMatrix(ed)
quartz()
cohorts(cmat, ed)
