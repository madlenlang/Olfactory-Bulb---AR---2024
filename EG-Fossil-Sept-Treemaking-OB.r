library(phytools) 
library(paleotree)
library(dispRity)
library(phangorn)

setwd("")

## READING THE MESQUITE TREE ##
	mytree <- read.nexus("mytree.nex")	

# Sorting out node ages for cal3 node.mins #
	nodetree<-makeNodeLabel(mytree, method = "number")	# assign nodelabels to the phy
	plot(nodetree, show.node.label = TRUE, cex = 0.5)	# plot the phy object
	nodeages<-read.csv("Node-Legend.csv",header=TRUE)
	nodeDates<-c(nodeages$Age1)
	
#### Make a calibrated tree 
#### bin_cal3TimePaleoPhy method 
#taxa range for age
	taxon.times<-read.csv("taxon_times_OB.csv",row.names=1)
	int.times<-read.csv("int_times_OB.csv",row.names=1)
	fossil.range<-list(int.times,taxon.times)

#bin_cal3TimePaleoPhy method
	likFun<-make_durationFreqDisc(fossil.range)
	spRes<-optim(parInit(likFun),likFun,lower=parLower(likFun),upper=parUpper(likFun),
		method="L-BFGS-B",control=list(maxit=1000000))

#sampling PROBABILITY per bin
	sProb <- spRes[[1]][2]

#calculate meanInt. We need to use an average int.length (intervals not the same duration)
	intLength<--apply(fossil.range[[1]],1,diff)
		hist(intLength)
	meanInt <--mean(apply(fossil.range[[1]], 1, diff)) #close to 1.8 Million years
	sRate <- sProb2sRate(sProb,int.length = meanInt)
	
# we also need extinction rate and branching rate (see above)
# need to divide by int.length...
	divRate <- spRes[[1]][1]/meanInt

#calibrated tree
	mytree1 <- bin_cal3TimePaleoPhy(mytree, 
		fossil.range, 
		FAD.only = TRUE,				# must be FALSE using "randObs" or "minMax" dateTreatment
		dateTreatment = "firstLast",	#  dates are first and last appearance, but that the time of observation is unknown
		root.max = 80.92,
		brRate = divRate, 
		extRate = divRate,
		sampRate = sRate, 
		ntrees = 100, 
		node.mins = nodeDates,	
		plot = FALSE
	)

#MCC tree 
	mcctree <- maxCladeCred(mytree1, tree = TRUE, part = NULL, rooted = TRUE)	#gives Maximum clade credibility tree (selects tree with highest clade support from list)
	
## Assigning node labels to tree (trying to make node.mins for cal3 function) ##
	mcctreenode<-makeNodeLabel(mcctree, method = "number") 
	plot(mcctreenode, show.node.label = TRUE, cex = 0.5)	# some extant clades will look like they have collapsed BUT it's just because branches are short
	
	treeage <- tree.age(mcctreenode, age=80.92, order="past", fossil=TRUE, digits=3)	# ages of the nodes and the root look correct 
	#write.nexus(mcctreenode, file = "timecaltree-OB.nex", translate = TRUE)
	#treetest <- read.nexus("timecaltree.nex")	










