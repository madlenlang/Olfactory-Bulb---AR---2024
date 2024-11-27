library(Arothron)
library(phytools)
library(phylolm)
library(geiger)
library(nlme)		# gls()
library(AICcmodavg)	# AIC testing
library(nortest)	# lillie test
library(RRPP)
library(ggpmisc)
library(ggplot2)

setwd("")

# Load Files #
	mytree <-read.nexus("mytree.nex")		# also used timecaltree-OB.nex for ancestral state reconstructions
	data <- read.csv("OB-Classifier-log.csv") 
	rownames(data) <- data$Sp

# Setting up data #
	tree_br <- compute.brlen(mytree, method = "Grafen", power = 1)								# compute branche lengths (tree has missing brlen)
	tree_node <- makeNodeLabel(tree_br, method = "number", prefix = "Node", nodeList = list())	# assign node labels
	tree_pr <- treedata(phy = tree_node, data = data, sort=TRUE)$phy							# prune tree - removes species not in dataset 
	data_pr <- data.frame(treedata(phy = tree_pr, data = data, sort=TRUE)$data)					# prune data - removes species not in tree; sorts data
	name.check(tree_pr, data_pr)																# check the names in the tree against the data

### VOLUME ANALYSES ###
# Setting up data - removing Plesiadapis cookei from volume analyses #
	datav <- data[!(row.names(data) %in% c("Plesiadapis_cookei_")),]			# drop P.c. from volume analyses
	tree_v <- treedata(phy = tree_node, data = datav, sort=TRUE)$phy			# prune P.c. from tree; tree used for PGLS volume analyses (below)
	name.check(tree_v, datav)
	tree_v$tip.label==datav$Sp
	
## NON-PHYLO REGRESSIONS ##
# ORDER - OLS model: OBV~ECV #
	modelv <- lm.rrpp(OBV~ECV, data = datav, print.progress = FALSE, RRPP = TRUE, iter = 9999)	
	summary(modelv)
	coef(modelv)
	residv <- residuals(modelv)

## Volume ANOVAS by groups ###
# Removing species with small sample size for ANOVAs #
	resid_vdf <- data.frame(residv, row.names = datav$Sp, vSp = datav$Sp, vOrd = datav$Ord, vC1 = datav$Clade1, vC2 = datav$Clade2)		# creating a df of the residuals 
	resid_vdf <- resid_vdf[!(row.names(resid_vdf) %in% 
		c("Cynocephalus_volans_","Galeopterus_variegatus_", "Carcinella_sigei_", "Anagale_gobiensis_","Labidolemur_kayi_")),]	# species to be removed due to small clade sample sizes

# Group = Clade1 (Hap/Strep seperated); non-phylo ANOVA  #	
	resid_vdf2 <- resid_vdf[!(row.names(resid_vdf) %in% 
		c("Ignacius_graybullianus_", "Microsyops_annectens_", "Niptomomys_cf__N__doreenae_", "Plesiadapis_tricuspidens_")),]	# additional species to be removed due to small clade sample sizes
	anova_vc1 <- lm.rrpp(residv~vC1, data = resid_vdf2, iter = 9999)
	anova(anova_vc1) 
	pw_vc1 <- pairwise(anova_vc1, groups = resid_vdf2$vC1)	
	summary(pw_vc1, stat.table = TRUE) 

# Group = Clade2 (Strep/Anthro/Tars seperated); non-phylo ANOVA # --> SD
	resid_vdf2 <- resid_vdf[!(row.names(resid_vdf) %in% 
		c("Ignacius_graybullianus_", "Microsyops_annectens_", "Niptomomys_cf__N__doreenae_", "Plesiadapis_tricuspidens_")),]	
	anova_vc2 <- lm.rrpp(residv~vC2, data = resid_vdf2, iter = 9999)
	anova(anova_vc2) 
	pw_vc1 <- pairwise(anova_vc2, groups = resid_vdf2$vC2)	
	summary(pw_vc1, stat.table = TRUE) 

## Specific groups regressions ##
# Haplorhini #
	hap.obv <- data$ECV[data$Clade1 == "Haplorhini"]
	hap.ecv <- data$OBV[data$Clade1 == "Haplorhini"]
	
	model.hap.v <- lm.rrpp(hap.obv~hap.ecv, print.progress = FALSE, RRPP = TRUE, iter = 9999)
	summary(model.hap.v)
	coef(model.hap.v)

# Strepsirrhini #
	str.obv <- data$ECV[data$Clade1 == "Strepsirrhini"]
	str.ecv <- data$OBV[data$Clade1 == "Strepsirrhini"]
	
	model.str.v <- lm.rrpp(str.obv~str.ecv, print.progress = FALSE, RRPP = TRUE, iter = 9999)
	summary(model.str.v)
	coef(model.str.v)
	
# Rodentia #
	rod.obv <- data$ECV[data$Ord == "Rodentia"]
	rod.ecv <- data$OBV[data$Ord == "Rodentia"]
	
	model.rod.v <- lm.rrpp(rod.obv~rod.ecv, print.progress = FALSE, RRPP = TRUE, iter = 9999)
	summary(model.rod.v)
	coef(model.rod.v)
	
# Lagomorpha #
	lag.obv <- data$ECV[data$Ord == "Lagomorpha"]
	lag.ecv <- data$OBV[data$Ord == "Lagomorpha"]

	model.lag.v <- lm.rrpp(lag.obv~lag.ecv, print.progress = FALSE, RRPP = TRUE, iter = 9999)
	summary(model.lag.v)
	coef(model.lag.v)

# Scandentia # 
	sca.obv <- data$ECV[data$Ord == "Scandentia"]
	sca.ecv <- data$OBV[data$Ord == "Scandentia"]
	
	model.sca.v <- lm.rrpp(sca.obv~sca.ecv, print.progress = FALSE, RRPP = TRUE, iter = 9999)
	summary(model.sca.v)
	coef(model.sca.v)

### MASS ANALYSES ###
# Setting up data - removing Plesiadapis tricuspidens from mass analyses #
	datam <- data[!(row.names(data) %in% c("Plesiadapis_tricuspidens_")),]		# drop P.t. from mass analyses
	tree_m <- treedata(phy = tree_node, data = datam, sort=TRUE)$phy			# prune P.t. from tree; tree used for PGLS mass analyses (below)
	name.check(tree_m, datam)	
	
## NON-PHYLO REGRESSIONS ##
# ORDER - OLS model: OBM~BM #
	modelm <- lm.rrpp(OBM~BM, data = datam, print.progress = FALSE, RRPP = TRUE, iter = 9999)	
	summary(modelm)
	coef(modelm)
	residm <- residuals(modelm)

## Mass ANOVAS by groups ###	
# Removing species with small sample size for ANOVAs #
	resid_mdf <- data.frame(residm, row.names = datam$Sp, mSp = datam$Sp, mOrd = datam$Ord, mC1 = datam$Clade1, mC2 = datam$Clade2)		# creating a df of the residuals 
	resid_mdf <- resid_mdf[!(row.names(resid_mdf) %in% 
		c("Cynocephalus_volans_","Galeopterus_variegatus_", "Carcinella_sigei_", "Anagale_gobiensis_","Labidolemur_kayi_")),]	# species to be removed due to small clade sample sizes

# Group = Clade1 (Hap/Step seperated); non-phylo ANOVA  #	
	resid_mdf2 <- resid_mdf[!(row.names(resid_mdf) %in% 
		c("Ignacius_graybullianus_", "Microsyops_annectens_", "Niptomomys_cf__N__doreenae_", "Plesiadapis_cookei_")),]	# additional species to be removed due to small clade sample sizes
	anova_mc1 <- lm.rrpp(residm~mC1, data = resid_mdf2, iter = 9999)
	anova(anova_mc1) 
	pw_mc1 <- pairwise(anova_mc1, groups = resid_mdf2$mC1)	
	summary(pw_mc1, stat.table = TRUE) 

# Group = Clade2 (Strep/Anthro/Tars seperated); non-phylo ANOVA #	--> SD
	resid_mdf2 <- resid_mdf[!(row.names(resid_mdf) %in% 
		c("Ignacius_graybullianus_", "Microsyops_annectens_", "Niptomomys_cf__N__doreenae_", "Plesiadapis_cookei_")),]	
	anova_mc2 <- lm.rrpp(residm~mC2, data = resid_mdf2, iter = 9999)
	anova(anova_mc2) 
	pw_mc1 <- pairwise(anova_mc2, groups = resid_mdf2$mC2)	
	summary(pw_mc1, stat.table = TRUE) 

## Specific groups regressions ##
# Haplorhini #
	hap.obm <- data$OBM[data$Clade1 == "Haplorhini"]
	hap.bm <- data$BM[data$Clade1 == "Haplorhini"]
	
	model.hap.m <- lm.rrpp(hap.obm~hap.bm, print.progress = FALSE, RRPP = TRUE, iter = 9999)	
	summary(model.hap.m)
	coef(model.hap.m)

# Strepsirrhini #
	str.obm <- data$OBM[data$Clade1 == "Strepsirrhini"]
	str.bm <- data$BM[data$Clade1 == "Strepsirrhini"]
	
	model.str.m <- lm.rrpp(str.obm~str.bm, print.progress = FALSE, RRPP = TRUE, iter = 9999)	
	summary(model.str.m)
	coef(model.str.m)
	
# Rodentia #
	rod.obm <- data$OBM[data$Ord == "Rodentia"]
	rod.bm <- data$BM[data$Ord == "Rodentia"]
	
	model.rod.m <- lm.rrpp(rod.obm~rod.bm, print.progress = FALSE, RRPP = TRUE, iter = 9999)	
	summary(model.rod.m)
	coef(model.rod.m)
	
# Lagomorpha #
	lag.obm <- data$OBM[data$Ord == "Lagomorpha"]
	lag.bm <- data$BM[data$Ord == "Lagomorpha"]
	
	model.lag.m <- lm.rrpp(lag.obm~lag.bm, print.progress = FALSE, RRPP = TRUE, iter = 9999)	
	summary(model.lag.m)
	coef(model.lag.m)

# Scandentia # 
	sca.obm <- data$OBM[data$Ord == "Scandentia"]
	sca.bm <- data$BM[data$Ord == "Scandentia"]
	
	model.sca.m <- lm.rrpp(sca.obm~sca.bm, print.progress = FALSE, RRPP = TRUE, iter = 9999)	
	summary(model.sca.m)
	coef(model.sca.m)
		
## Phylogenetic signals ##
	tree_v$tip.label==datav$Sp		# Making sure data and tip.labels are in the same order
	tree_m$tip.label==datam$Sp		
	
# Prune data can make values read as characters; this converts information back into numeric format #
	vECV <-as.numeric(as.character(datav$ECV))		# from df with P.c. removed 
		names(vECV) <- rownames(datav)
	vOBV <-as.numeric(as.character(datav$OBV))		# from df with P.c. removed
		names(vOBV) <- rownames(datav)
	vperv <-as.numeric(as.character(datav$perv))	# from df with P.c. removed
		names(vperv) <- rownames(datav)
	mBM <-as.numeric(as.character(datam$BM))		# from df with P.t. removed
		names(mBM) <- rownames(datam)
	mOBM <-as.numeric(as.character(datam$OBM))		# from df with P.t. removed
		names(mOBM) <- rownames(datam)
	mperm <-as.numeric(as.character(datam$perm))	# from df with P.t. removed
		names(mperm) <- rownames(datam)

# Pagel's lambda for different varaibles #
	phylosig(tree_v, vOBV, method = "lambda", test = TRUE, nsim = 9999)	 
	phylosig(tree_v, vECV, method = "lambda", test = TRUE, nsim = 9999)	
	phylosig(tree_m, mBM, method = "lambda", test = TRUE, nsim = 9999)	
	phylosig(tree_v, vperv, method = "lambda", test = TRUE, nsim = 9999)	
	phylosig(tree_m, mperm, method = "lambda", test = TRUE, nsim = 9999)	
	
	resid_v <- data.frame(residv, row.names = datav$Sp)
	phylosig(tree_v, resid_v$residv, method = "lambda", test = TRUE, nsim = 9999)
	resid_m <- data.frame(residm, row.names = datam$Sp)
	phylosig(tree_m, resid_m$residm, method = "lambda", test = TRUE, nsim = 9999)
	
# Clade Specific Phylogenetic signals #
	Species.der <- data_pr$Sp[data_pr$Clade1 == "Dermoptera"]		# assign species label used to match tip names
	Species.lag <- data_pr$Sp[data_pr$Clade1 == "Lagomorpha"]		
	Species.hap <- data_pr$Sp[data_pr$Clade1 == "Haplorhini"]
	Species.rod <- data_pr$Sp[data_pr$Clade1 == "Rodentia"]
	Species.sca <- data_pr$Sp[data_pr$Clade1 == "Scandentia"]
	Species.str <- data_pr$Sp[data_pr$Clade1 == "Strepsirrhini"]
	Species.pri <- data_pr$Sp[data_pr$Clade1 == "Primates"]
	Species.eua <- data_pr$Sp[data_pr$Clade1 == "Euarchontoglires"]
		
	lag.tree <- drop.tip(tree_pr, c(Species.hap, Species.rod, Species.sca, Species.str, Species.der, Species.pri, Species.eua))	# drop tips from tree for groups not of interest
	hap.tree <- drop.tip(tree_pr, c(Species.lag, Species.rod, Species.sca, Species.der, Species.str, Species.pri, Species.eua))
	rod.tree <- drop.tip(tree_pr, c(Species.lag, Species.hap, Species.sca, Species.str, Species.der, Species.pri, Species.eua))
	sca.tree <- drop.tip(tree_pr, c(Species.lag, Species.rod, Species.hap, Species.str, Species.der, Species.pri, Species.eua))
	str.tree <- drop.tip(tree_pr, c(Species.lag, Species.rod, Species.sca, Species.hap, Species.der, Species.pri, Species.eua))
	
# Haplorhini #
	hap.df <- datav[!(row.names(datav) %in% 
		c(Species.lag, Species.rod, Species.sca, Species.str, Species.der, Species.pri, Species.eua)),]	
	hap.tree$tip.label==hap.df$Sp	# Making sure data and tip.labels are in the same order
	phylosig(hap.tree, hap.df$OBV, method = "lambda", test = TRUE, nsim = 9999)	
	phylosig(hap.tree, hap.df$ECV, method = "lambda", test = TRUE, nsim = 9999)	
	phylosig(hap.tree, hap.df$BM, method = "lambda", test = TRUE, nsim = 9999)	
	phylosig(hap.tree, hap.df$perv, method = "lambda", test = TRUE, nsim = 9999)	
	phylosig(hap.tree, hap.df$perm, method = "lambda", test = TRUE, nsim = 9999)
	
# Strepsirrhini #
	str.df <- datav[!(row.names(datav) %in% 
		c(Species.lag, Species.rod, Species.sca, Species.hap, Species.der, Species.pri, Species.eua)),]	
	str.tree$tip.label==str.df$Sp
	phylosig(str.tree, str.df$OBV, method = "lambda", test = TRUE, nsim = 9999)		
	phylosig(str.tree, str.df$ECV, method = "lambda", test = TRUE, nsim = 9999)	
	phylosig(str.tree, str.df$BM, method = "lambda", test = TRUE, nsim = 9999)	
	phylosig(str.tree, str.df$perv, method = "lambda", test = TRUE, nsim = 9999)	
	phylosig(str.tree, str.df$perm, method = "lambda", test = TRUE, nsim = 9999)
	
# Rodentia #
	rod.df <- datav[!(row.names(datav) %in% 
		c(Species.lag, Species.hap, Species.sca, Species.str, Species.der, Species.pri, Species.eua)),]	
	rod.tree$tip.label==rod.df$Sp
	phylosig(rod.tree, rod.df$OBV, method = "lambda", test = TRUE, nsim = 9999)		
	phylosig(rod.tree, rod.df$ECV, method = "lambda", test = TRUE, nsim = 9999)	
	phylosig(rod.tree, rod.df$BM, method = "lambda", test = TRUE, nsim = 9999)	
	phylosig(rod.tree, rod.df$perv, method = "lambda", test = TRUE, nsim = 9999)	
	phylosig(rod.tree, rod.df$perm, method = "lambda", test = TRUE, nsim = 9999)

# Lagomorpha #
	lag.df <- datav[!(row.names(datav) %in% 
		c(Species.hap, Species.rod, Species.sca, Species.str, Species.der, Species.pri, Species.eua)),]	
	lag.tree$tip.label==lag.df$Sp
	phylosig(lag.tree, lag.df$OBV, method = "lambda", test = TRUE, nsim = 9999)		
	phylosig(lag.tree, lag.df$ECV, method = "lambda", test = TRUE, nsim = 9999)	
	phylosig(lag.tree, lag.df$BM, method = "lambda", test = TRUE, nsim = 9999)	
	phylosig(lag.tree, lag.df$perv, method = "lambda", test = TRUE, nsim = 9999)	
	phylosig(lag.tree, lag.df$perm, method = "lambda", test = TRUE, nsim = 9999)
	
# Scandentia #
	sca.df <- datav[!(row.names(datav) %in% 
		c(Species.lag, Species.rod, Species.hap, Species.str, Species.der, Species.pri, Species.eua)),]	
	sca.tree$tip.label==sca.df$Sp
	phylosig(sca.tree, sca.df$OBV, method = "lambda", test = TRUE, nsim = 9999)		
	phylosig(sca.tree, sca.df$ECV, method = "lambda", test = TRUE, nsim = 9999)	
	phylosig(sca.tree, sca.df$BM, method = "lambda", test = TRUE, nsim = 9999)	
	phylosig(sca.tree, sca.df$perv, method = "lambda", test = TRUE, nsim = 9999)	
	phylosig(sca.tree, sca.df$perm, method = "lambda", test = TRUE, nsim = 9999)

## PGLS REGRESSIONS - VOLUME ##
# PGLS model: OBV~ECV using Pagel's correltation structure; P.c. removed #
	OBV.pgls <- phylolm(vOBV~vECV, phy = tree_v, model = "lambda")
		summary(OBV.pgls)
	OBVres <- OBV.pgls$residuals

# Group = Order; phylo ANOVA  # - small groups sizes removed #
	ptree_v <- drop.tip(tree_v, c("Cynocephalus_volans_","Galeopterus_variegatus_", "Carcinella_sigei_", "Anagale_gobiensis_","Labidolemur_kayi_")) # dropping tips with low clade sample sizes
	pvres <- OBVres[-c(184, 183, 182, 16, 15)]		# removing species names from the species list (based on their number seq in the dataset)		
		names(pvres) <-rownames(resid_vdf)			# assign the species names to the residuals from the nonphylo residuals df
	name.check(ptree_v, pvres)

	phanova_v <- lm.rrpp(pvres~resid_vdf$vOrd, Cov = vcv.phylo(ptree_v), iter = 9999)
		anova(phanova_v) 
	phypw_v <- pairwise(phanova_v, groups = resid_vdf$vOrd)	
		summary(phypw_v, stat.table = TRUE) 
	
# Group = Clade1; phylo ANOVA  #	
	ptree_v2 <- drop.tip(ptree_v, c("Ignacius_graybullianus_", "Microsyops_annectens_", "Niptomomys_cf__N__doreenae_", "Plesiadapis_tricuspidens_")) # dropping tips with low clade sample sizes
	pvres2 <- OBVres[-c(205, 204, 203, 202, 184, 183, 182, 16, 15)]
		names(pvres2) <-rownames(resid_vdf2)	
	
	phanova_vc1 <- lm.rrpp(pvres2~resid_vdf2$vC1, Cov = vcv.phylo(ptree_v2), iter = 9999)
		anova(phanova_vc1) 
	phypw_vc1 <- pairwise(phanova_vc1, groups = resid_vdf2$vC1)	
		summary(phypw_vc1, stat.table = TRUE) 
	
# Group = Clade2; phylo ANOVA  #	
	phanova_vc2 <- lm.rrpp(pvres2~resid_vdf2$vC2, Cov = vcv.phylo(ptree_v2), iter = 9999)
		anova(phanova_vc2) 
	phypw_vc2 <- pairwise(phanova_vc2, groups = resid_vdf2$vC2)	
		summary(phypw_vc2, stat.table = TRUE) 

## PGLS REGRESSIONS - MASS ##
# PGLS model: OBM~BM using Pagel's correltation structure; P.t. removed #
	OBM.pgls <- phylolm(mOBM~mBM, phy = tree_m, model = "lambda")
		summary(OBM.pgls)
	OBMres <- OBM.pgls$residuals

# Group = Order; phylo ANOVA  # - small groups sizes removed #
	ptree_m <- drop.tip(tree_m, c("Cynocephalus_volans_","Galeopterus_variegatus_", "Carcinella_sigei_", "Anagale_gobiensis_","Labidolemur_kayi_")) 
	pmres <- OBMres[-c(184, 183, 182, 16, 15)]			
		names(pmres) <-rownames(resid_mdf)			
	name.check(ptree_m, pmres)

	phanova_m <- lm.rrpp(pmres~resid_mdf$mOrd, Cov = vcv.phylo(ptree_m), iter = 9999)
		anova(phanova_m) 
	phypw_m <- pairwise(phanova_m, groups = resid_mdf$mOrd)	
		summary(phypw_m, stat.table = TRUE) 

# Group = Clade1; phylo ANOVA  #	
	ptree_m2 <- drop.tip(ptree_v, c("Ignacius_graybullianus_", "Microsyops_annectens_", "Niptomomys_cf__N__doreenae_", "Plesiadapis_tricuspidens_")) # dropping tips with low clade sample sizes
	pmres2 <- OBMres[-c(205, 204, 203, 202, 184, 183, 182, 16, 15)]
		names(pmres2) <-rownames(resid_mdf2)	
	
	phanova_mc1 <- lm.rrpp(pmres2~resid_mdf2$mC1, Cov = vcv.phylo(ptree_m2), iter = 9999)
		anova(phanova_mc1) 
	phypw_mc1 <- pairwise(phanova_mc1, groups = resid_mdf2$mC1)	
		summary(phypw_mc1, stat.table = TRUE) 
	
# Group = Clade2; phylo ANOVA  #	
	phanova_mc2 <- lm.rrpp(pmres2~resid_mdf2$mC2, Cov = vcv.phylo(ptree_m2), iter = 9999)
		anova(phanova_mc2) 
	phypw_mc2 <- pairwise(phanova_mc2, groups = resid_mdf2$mC2)	
		summary(phypw_mc2, stat.table = TRUE) 

## PLOTS ##
## Ancestral state reconstruction - volume ##
	fit.v <- fastAnc(tree_v,vperv,vars=TRUE,CI=TRUE) #compute variances & 95% confidence intervals for each node
	fit.v$CI[1,]
	
	plot(tree_v,no.margin=TRUE,edge.width=2,cex=0.7)
	nodelabels(frame="none",adj=c(1.1,-0.4), cex=0.7)	# produce a tree with node labels to compare to ASR in ace_v

# projection of the reconstruction onto the edges of the tree
# mapping is accomplished by estimating states at internal nodes using ML with fastAnc, 
# and then interpolating the states along each edge using equation [2] of Felsenstein (1985).
	obj.v <- contMap(tree_v,vperv,plot=FALSE) #maps the observed and reconstructed phenotypic trait values onto the tree using a color gradient
	obj.v <- setMap(obj.v, colors=c("#660099", "#3399FF", "#66FF66", "#FF9900", "#FF3300"))
	plot(obj.v, 
		outline = FALSE, 
		res = 100,
		type = "fan", 
		legend = 1*max(nodeHeights(tree_v)),	
		fsize = c(0.5),	#font size
		ftype = "i",	#font type-"reg", "i" (italics), "b" (bold), or "bi" (bold-italics)	
		lwd = 3.9		#line width for branches
	)
	
# Ancestral state reconstruction - mass #
	fit.m <- fastAnc(tree_m,mperm,vars=TRUE,CI=TRUE) 
	fit.m$CI[1,]
	
	plot(tree_m,no.margin=TRUE,edge.width=2,cex=0.7)
	nodelabels(frame="none",adj=c(1.1,-0.4), cex=0.7)	

# projection of the reconstruction onto the edges of the tree
	obj.m <- contMap(tree_m,mperm,plot=FALSE) 
	obj.m <- setMap(obj.m, colors=c("#660099", "#3399FF", "#66FF66", "#FF9900", "#FF3300"))
	plot(obj.m, 
		outline = FALSE, 
		res = 100,
		type = "fan", 
		legend = 1*max(nodeHeights(tree_m)),	
		fsize = c(0.5),
		ftype = "i",	
		lwd = 3.9			
	)

# data frame with Ancestral State Reconstrictions ##
	ace <-data.frame(ASRvol = fit.v$ace, VARvol = fit.v$var, CI95vol = fit.v$CI95, ASRmass = fit.m$ace, VARmass = fit.m$var, CI95mass = fit.m$CI95) 
	#write.csv(ace, 'ASR.csv') 

## SCATTER PLOTS ##	
# Plot colours #

	cols <- c("#3399FF", "#FFCC00", "#339966", "#FF6600", "#330099", "#66CC00", "#CC0000", "#6600FF", "#FF33CC")

# VOLUME - Clade Specific LM regressions line #
	vdf <- data.frame(vECV, vOBV, vperv, datav$Clade1, datav$Fossil) 
	
	ggplot(vdf, aes(vECV, vOBV, color = datav$Clade2)) +
		geom_point(aes(color = datav$Clade2, shape=datav$Fossil), size = 4, show.legend = TRUE, alpha = 0.85) +
		geom_smooth(method="lm", se = FALSE, fullrange = FALSE, level = 0.95, size = 1, linetype = "solid") +	# add a regression line based on Clade2
		stat_smooth(method = "lm", col = "black", linetype = "longdash", se = FALSE, fullrange = FALSE, size = 1) +	
		scale_color_manual(values=cols) +
		scale_shape_manual(values=c(16, 13)) +
		theme(legend.position="top",
			axis.text.x = element_text(size = 15), 	# font size 13 on x and y axes
			axis.text.y = element_text(size = 15),
			legend.key = element_rect(fill = "white"),
			panel.background = element_rect(fill = "white"),
			panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid', colour = "gray90"), 
			panel.grid.minor = element_line(linewidth = 0.5, linetype = 'solid', colour = "gray90"),
			panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5)				
		) + 
		labs(x = "Log Endocranial Volume",
			y = "Log Olfactory Bulb Volume",
		) 

# MASS - Clade Specific LM regressions line #
	mdf <- data.frame(mBM, mOBM, mperm, datam$Clade2, datam$Fossil) 
	
	ggplot(mdf, aes(mBM, mOBM, color = datam$Clade2)) +
		geom_point(aes(color = datam$Clade2, shape=datam$Fossil), size = 4, show.legend = TRUE, alpha = 0.85) +
		geom_smooth(method="lm", se = FALSE, fullrange = FALSE, level = 0.95, size = 1, linetype = "solid") +	
		stat_smooth(method = "lm", col = "black", linetype = "longdash", se = FALSE, fullrange = FALSE, size = 1) +	
		scale_color_manual(values=cols) +
		scale_shape_manual(values=c(16, 13)) +
		theme(legend.position="top",
			axis.text.x = element_text(size = 15), 	
			axis.text.y = element_text(size = 15 ),
			legend.key = element_rect(fill = "white"),
			panel.background = element_rect(fill = "white"),
			panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid', colour = "gray90"), 
			panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid', colour = "gray90"),
			panel.border = element_rect(colour = "black", fill=NA, linewidth=1)	
		) + 
		labs(x = "Log Body Mass",
			y = "Log Olfactory Bulb Mass",
		) 

## BOX PLOTS ##
# Volume - Boxplot of residuals from non-phylo OBM~BM linear model ##
	resid_vbp <- data.frame(residv, row.names = datav$Sp, vC2 = datav$Clade2)		
	
	ggplot(resid_vbp, aes(vC2, residv, fill = vC2)) +
		geom_boxplot(lwd = 0.9, 
			alpha = 0.8,
			colour = "black") +
		scale_color_manual(values=cols) +
		scale_fill_manual(values=cols) +
		labs(x = "Order",
				y = "Residuals of Log Olfactory Bulb Volume") +
		theme(legend.position="top",
			axis.text.x = element_text(size = 12), 	
			axis.text.y = element_text(size = 12 ),
			legend.key = element_rect(fill = "white"),
			panel.background = element_rect(fill = "white"),
			panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid', colour = "gray90"), 
			panel.grid.minor = element_line(linewidth = 0.5, linetype = 'solid', colour = "gray90"),
			panel.border = element_rect(colour = "black", fill=NA, linewidth=1)				
		) 

# Mass - Boxplot of residuals from non-phylo OBM~BM linear model ##
	resid_mbp <- data.frame(residm, row.names = datam$Sp, mC2 = datam$Clade2)		

	ggplot(resid_mbp, aes(mC2, residm, fill = mC2)) +
		geom_boxplot(lwd = 0.9, 
			alpha = 0.8,
			colour = "black") +
		scale_color_manual(values=cols) +
		scale_fill_manual(values=cols) +
		labs(x = "Order",
				y = "Residuals of Log Olfactory Bulb Mass") +
		theme(legend.position="top",
			axis.text.x = element_text(size = 12), 
			axis.text.y = element_text(size = 12 ),
			legend.key = element_rect(fill = "white"),
			panel.background = element_rect(fill = "white"),
			panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid', colour = "gray90"), 
			panel.grid.minor = element_line(linewidth = 0.5, linetype = 'solid', colour = "gray90"),
			panel.border = element_rect(colour = "black", fill=NA, linewidth=1)				
		) 
		
# Volume - lightning plot - no crown #
	lp <- read.csv("Node-Fig-data-nocrown.csv") 
	cols1 <- lp$cols1

	ggplot(lp, aes(vASR, node_n, color = Node_Name)) +
		geom_point(aes(fill=Node_Name), size = 7, pch = 21, colour = "black", show.legend = TRUE, alpha = 1) +
		scale_fill_manual(values=cols1) +
		theme(legend.position="top",
			axis.text.x = element_text(size = 15), 
			axis.text.y = element_text(size = 15 ),
			legend.key = element_rect(fill = "white"),
			panel.background = element_rect(fill = "white"),
			panel.border = element_rect(colour = "black", fill=NA, linewidth=1)			
		) + 
		labs(x = "Percent Olfactory Bulb Volume",
			y = "Distance from Root",
		) 

# Mass - lightning plot - no crown #
	ggplot(lp, aes(mASR, node_n, color = Node_Name)) +
		geom_point(aes(fill = Node_Name), size = 7, pch = 21, colour = "black", show.legend = TRUE, alpha = 1) +
		scale_fill_manual(values=cols1) +
		theme(legend.position="top",
			axis.text.x = element_text(size = 15),
			axis.text.y = element_text(size = 15 ),
			legend.key = element_rect(fill = "white"),
			panel.background = element_rect(fill = "white"),
			panel.border = element_rect(colour = "black", fill=NA, linewidth=1)		
		) + 
		labs(x = "Percent Olfactory Bulb Mass",
			y = "Distance from Root",
		) 

# Volume - lightning plot - crown #
	lp <- read.csv("Node-Fig-data-crown.csv") 
	cols1 <- lp$cols1

	ggplot(lp, aes(vASR, node_n, color = Node_Name)) +
		geom_point(aes(fill = Node_Name), size = 8, pch = 21, colour = "black", show.legend = TRUE, alpha = 1) +
		scale_fill_manual(values=cols1) +
		theme(legend.position="top",
			axis.text.x = element_text(size = 15),
			axis.text.y = element_text(size = 15 ),
			legend.key = element_rect(fill = "white"),
			panel.background = element_rect(fill = "white"),
			panel.border = element_rect(colour = "black", fill=NA, linewidth=1)			
		) + 
		labs(x = "Percent Olfactory Bulb Volume",
			y = "Distance from Root",
		) 

# Mass - lightning plot - crown #
	ggplot(lp, aes(mASR, node_n, color = Node_Name)) +
		geom_point(aes(fill = Node_Name), size = 8, pch = 21, colour = "black", show.legend = TRUE, alpha = 1) +
		scale_fill_manual(values=cols1) +
		theme(legend.position="top",
			axis.text.x = element_text(size = 15),
			axis.text.y = element_text(size = 15 ),
			legend.key = element_rect(fill = "white"),
			panel.background = element_rect(fill = "white"),
			panel.border = element_rect(colour = "black", fill=NA, linewidth=1)		
		) + 
		labs(x = "Percent Olfactory Bulb Mass",
			y = "Distance from Root",
		)
		
## SUPP MAT - NORMALITY TESTING for PGLS models##	
# Prune data can make values read as characters (same issue in phylosig); this converts information back into numeric format #
	ECV <-as.numeric(as.character(data_pr$ECV))
		names(ECV) <- rownames(data_pr)
	OBV <-as.numeric(as.character(data_pr$OBV))
		names(OBV) <- rownames(data_pr)
	BM <-as.numeric(as.character(data_pr$BM))
		names(BM) <- rownames(data_pr)
	OBM <-as.numeric(as.character(data_pr$OBM))
		names(OBM) <- rownames(data_pr)
	perv <-as.numeric(as.character(data_pr$perv))
		names(perv) <- rownames(data_pr)
	perm <-as.numeric(as.character(data_pr$perm))
		names(perm) <- rownames(data_pr)
		
#Define correlation structures #  
	bm.tree<-corBrownian(1, phy=tree_pr)				# Brownian Correlation 
	pa.tree<-corPagel(1, phy=tree_pr, fixed = FALSE)	# Pagel's correcltion (derived from Brownian Motion model)
	ou.tree<-corMartins(1, phy=tree_pr, fixed = TRUE)	# Ornstein-Uhlenbeck Correlation (based on Martin's)
	gr.tree<-corGrafen(1, phy=tree_pr, fixed = FALSE)	# Grafen's Correlation 

## VOLUME ##
#PGLS regressions with different correlation structures #
	model.bm <- gls(OBV ~ ECV, correlation = bm.tree, method = "ML")
		summary(model.bm)
		coef(model.bm)
		resid.bm <- residuals(model.bm)	
		
	model.pa <- gls(OBV ~ ECV, correlation = pa.tree, method = "ML")
		summary(model.pa)
		coef(model.pa)
		resid.pa <- residuals(model.pa)	

	model.ou <- gls(OBV ~ ECV, correlation = ou.tree, method = "ML")
		summary(model.ou)
		coef(model.ou)
		resid.ou <- residuals(model.ou)	
		
# Residual error normality #
	lillie.test(chol(solve(vcv(tree_pr)))%*%resid.bm)	
	lillie.test(chol(solve(vcv(tree_pr)))%*%resid.ou)
	lillie.test(chol(solve(vcv(tree_pr)))%*%resid.pa)
		
	# Interpret lillie test using table - if D is less than critical value
	# you can accept the null hypothesis that there is a normal distribution 
	# See for more info - http://www.real-statistics.com/tests-normality-and-symmetry/statistical-tests-normality-symmetry/lilliefors-test-normality/
	# n = 222, critical number (0.5) is 0.895
		
# Model assessment using AICc #	
# AIC(c) Table #	
	Cand.models <- list()
		Cand.models[[1]] <- model.bm
		Cand.models[[2]] <- model.pa
		Cand.models[[3]] <- model.ou
	Modnames <- paste("mod", 1:length(Cand.models), sep = " ")
	print(aictab(cand.set = Cand.models, modnames =Modnames, sort = TRUE), digits = 4 , LL = TRUE)		# print table with AIC information, 'aictab' creates a ranked model selection table 
	
## MASS ##
#PGLS regressions with different correlation structures #
	model.bm <- gls(OBM ~ BM, correlation = bm.tree, method = "ML")
		summary(model.bm)
		coef(model.bm)
		resid.bm <- residuals(model.bm)	
		
	model.pa <- gls(OBM ~ BM, correlation = pa.tree, method = "ML")
		summary(model.pa)
		coef(model.pa)
		resid.pa <- residuals(model.pa)	

	model.ou <- gls(OBM ~ BM, correlation = ou.tree, method = "ML")
		summary(model.ou)
		coef(model.ou)
		resid.ou <- residuals(model.ou)	
		
# Residual error normality #
	lillie.test(chol(solve(vcv(tree_pr)))%*%resid.bm)	
	lillie.test(chol(solve(vcv(tree_pr)))%*%resid.ou)
	lillie.test(chol(solve(vcv(tree_pr)))%*%resid.pa)
		
# Model assessment using AICc #	
# AIC(c) Table #	
	Cand.models <- list()
		Cand.models[[1]] <- model.bm
		Cand.models[[2]] <- model.pa
		Cand.models[[3]] <- model.ou
	Modnames <- paste("mod", 1:length(Cand.models), sep = " ")
	print(aictab(cand.set = Cand.models, modnames =Modnames, sort = TRUE), digits = 4 , LL = TRUE)		# print table with AIC information, 'aictab' creates a ranked model selection table 
	
