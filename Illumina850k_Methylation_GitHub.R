## Author: Marcos Elizalde Horcada - marcos.elizaldeh@gmail.com
## Date: January 2021

## Based in A cross-package Bioconductor workflow for analysing methylation 
## array (Jovana Maksimovic, Belinda Phipson and Alicia Oshlack, 2020)

#--------------------------------------------------------------------#
#------------- Setting up the workspace and data loading ------------# 
#--------------------------------------------------------------------#

library(knitr)
library(limma)
library(minfi)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(IlluminaHumanMethylationEPICmanifest)
library(RColorBrewer)
library(missMethyl)
library(minfiData)
library(Gviz)
library(DMRcate)
library(stringr)
library(wateRmelon)
library(ChAMP)
library(readxl)
library(ggplot2)
library(reshape2)

# Get the Illumina 850k annotation data
ann850k <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
# Set up data directory
dataDirectory <- "IDAT"
# Read in the sample sheet for the experiment
targets <- read.metharray.sheet(dataDirectory, 
                                pattern="Samplesheet epilepsia.csv")

# Read in the raw data from the IDAT files
rgSet <- read.metharray.exp(targets=targets)
# Give the samples descriptive names
targets$ID <- paste(targets$Sample_Group,targets$Sample_Name,sep=".")
sampleNames(rgSet) <- targets$ID
rgSet

# Calculate the detection p-values
detP <- detectionP(rgSet)

# Examine detection p-values for all samples
pal <- brewer.pal(12,"Paired")
#par(mfrow=c(1,2)) for double plotting, LEGENDS ARE PROVISIONAL
barplot(colMeans(detP), col=pal[factor(targets$Sample_Group)], las=2, 
        cex.names=0.8, ylab="Mean detection p-values")
abline(h=0.05,col="red")
legend("topleft", legend=levels(factor(targets$Sample_Group)), fill=pal,
       bg="white", cex = 0.5)

barplot(colMeans(detP), col=pal[factor(targets$Sample_Group)], las=2, 
        cex.names=0.8, ylim=c(0,0.0015), ylab="Mean detection p-values")
abline(h=0.05,col="red")
legend("topleft", legend=levels(factor(targets$Sample_Group)), fill=pal, 
       bg="white", cex = 0.5)

# QC Report from the minfy package 
qcReport(rgSet, sampNames=targets$ID, sampGroups=targets$Sample_Group, 
         pdf="qcReport.pdf")

#--------------------------------------------------------------------#
#------------------------ Noob Normalization ------------------------# 
#--------------------------------------------------------------------#

# Remove poor quality samples
keep <- colMeans(detP) < 0.01
rgSet <- rgSet[,keep]
rgSet #1051943
# Remove poor quality samples from experimental design
targets <- targets[keep,]
targets
# Remove poor quality samples from detection p-value table
detP <- detP[,keep]
dim(detP)

## Normalization is performed with preprocessNoob() from minfy package,
## then normalized betas are extracted with getBeta() from minfy package.
## Later, these values are normalized again with BMIQ via wateRmelon through
## the champ.norm() function. Finally, different types of probes are removed.
## Data suggest that Noob + BMIQ improves reproducibility

# Normalization with Noob
mSetSq <- preprocessNoob(rgSet)
# Extract raw data for plotting
mSetRaw <- preprocessRaw(rgSet)

# Normalized data vs. Raw data
par(mfrow=c(1,2))
densityPlot(rgSet, sampGroups=targets$Sample_Group,main="Raw", legend=FALSE)
legend("top", legend = levels(factor(targets$Sample_Group)), 
       text.col=brewer.pal(12,"Paired"))
densityPlot(getBeta(mSetSq), sampGroups=targets$Sample_Group,
            main="Normalized", legend=FALSE)
legend("top", legend = levels(factor(targets$Sample_Group)), 
       text.col=brewer.pal(12,"Paired"))

# MDS plot 
plotMDS(getM(mSetSq), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Group)], xlim=(c(-4,6)),
        ylim=c(-2,2))
#legend("top", legend=levels(factor(targets$Sample_Group)), text.col=pal,
#bg="white", cex=0.7)

plotMDS(getM(mSetSq), top=1000, gene.selection="common", 
        col=pal[factor(targets$cellType)], dim=c(1,3))
#legend("top", legend=levels(factor(targets$Sample_Group)), text.col=pal, 
#cex=0.7, bg="white")

plotMDS(getM(mSetSq), top=1000, gene.selection="common", 
        col=pal[factor(targets$cellType)], dim=c(2,3),
        ylim=c(-1,1), xlim=c(-1,1))
#legend("topleft", legend=levels(factor(targets$Sample_Group)), text.col=pal,
#cex=0.7, bg="white")

plotMDS(getM(mSetSq), top=1000, gene.selection="common", 
        col=pal[factor(targets$cellType)], dim=c(3,4),
        xlim=c(-1,1))
#legend("topright", legend=levels(factor(targets$Sample_Group)), text.col=pal,
#cex=0.7, bg="white")

#--------------------------------------------------------------------#
#------------------------ BMIQ Normalization ------------------------# 
#--------------------------------------------------------------------#

# Extract beta values
bVals <- getBeta(mSetSq)
write.table(bVals, "betasNoob.txt")

#Normalize betasNoob with BMIQ via wateRmelon: con champ.norm()
betasNoobBMIQ <- champ.norm(beta=bVals,
                 rgSet=rgSet,
                 mset=mSetSq,
                 method="BMIQ",
                 plotBMIQ=FALSE,
                 arraytype="EPIC",
                 cores=3)

write.table(betasNoobBMIQ, "betasNoobBMIQ.txt")

#--------------------------------------------------------------------#
#----------------------------- Filtering ----------------------------# 
#--------------------------------------------------------------------#

# Filtering probes from betasNoobBMIQ:
betasNoobBMIQ <- as.data.frame(betasNoobBMIQ)
dim(betasNoobBMIQ) #866238

# Ensure probes are in the same order in the mSetSq and detP objects
detP <- detP[match(rownames(betasNoobBMIQ),rownames(detP)),]

# Remove any probes that have failed in one or more samples
keep <- rowSums(detP < 0.01) == ncol(mSetSq)
table(keep)
betasNoobBMIQ <- betasNoobBMIQ[keep,]
dim(betasNoobBMIQ) #858172

## If your data includes males and females, remove probes on the sex chromosomes
keep <- !(rownames(betasNoobBMIQ) %in% ann850k$Name[ann850k$chr %in%
                                                        c("chrX","chrY")])
table(keep)
betasNoobBMIQ <- betasNoobBMIQ[keep,]
dim(betasNoobBMIQ) #839254

## Exclude cross reactive probes (Chen et al 2013)
xReactiveProbes1 <- read_excel("48639-non-specific-probes-Illumina450k.xlsx",
                               sheet = 1)

keep <- !(rownames(betasNoobBMIQ) %in% xReactiveProbes1$TargetID)
table(keep)
betasNoobBMIQ <- betasNoobBMIQ[keep,] 
dim(betasNoobBMIQ) #812946

xReactiveProbes2 <- read_excel("48639-non-specific-probes-Illumina450k.xlsx",
                               sheet = 2)

keep <- !(rownames(betasNoobBMIQ) %in% xReactiveProbes2$TargetID)
table(keep)
betasNoobBMIQ <- betasNoobBMIQ[keep,] #811359
dim(betasNoobBMIQ)

# Exclude cross reactive probes (Pidsley et al.) Table S1
xReactiveProbes3 <- read.csv("13059_2016_1066_MOESM1_ESM.csv",
                             stringsAsFactors = FALSE)
keep <- !(rownames(betasNoobBMIQ) %in% xReactiveProbes3$TargetID)
table(keep)
betasNoobBMIQ <- betasNoobBMIQ[keep,] 
dim(betasNoobBMIQ) #796229

# Exclude SNPs.147 for CpG_maf >= 0.01
dataSNPs <- read.csv("SNPs.147CommonSingle.csv", stringsAsFactors = FALSE)
dataSNPs[is.na(dataSNPs)] = 0
SNPs.remove <- subset(dataSNPs, CpG_maf >= 0.01)
keep <- !(rownames(betasNoobBMIQ) %in% SNPs.remove$TargetID)
table(keep)
betasNoobBMIQ <- betasNoobBMIQ[keep,]
dim(betasNoobBMIQ) #772002

# Convert beta values to M values
mVals <- log2(betasNoobBMIQ/(1-betasNoobBMIQ))

# Plot M-values
par(mfrow=c(1,2))
plotMDS(mVals, top=1000, gene.selection="common",
        col=pal[factor(targets$Sample_Group)], cex=0.8)
legend("right", legend=levels(factor(targets$Sample_Group)), text.col=pal,
       cex=0.65, bg="white")

plotMDS(mVals, top=1000, gene.selection="common",
        col=pal[factor(targets$Tissue)])
legend("right", legend=levels(factor(targets$Sample_Source)), text.col=pal,
       cex=0.7, bg="white")

# Clean workspace and save R workspace image
rm(dataSNPs)
rm(detP)
rm(mSetRaw)
rm(mSetSq)
rm(probeInfoALL.lv)
rm(SNPs.remove)
rm(xReactiveProbes1)
rm(xReactiveProbes2)
rm(xReactiveProbes3)
rm(keep)

# WORKSPACE IMAGE CREATED HERE (can run HOUSEMAN part from this workspace)

#--------------------------------------------------------------------#
#------------ Probe-wise differential methylation analysis ----------# 
#--------------------------------------------------------------------#

betasNoobBMIQ <- read.csv("betasNoobBMIQ_filtered.csv")

# Choose the factor of interest
#cellType <- factor(targets$Sample_Group)
#cellType <- factor(targets$cellType)
cellType <- factor(targets$cellType)

# Use the above to create a design matrix
design <- model.matrix(~0+cellType, data=targets)
colnames(design) <- levels(cellType)
 
# Fit the linear model 
fit <- lmFit(mVals, design)
# Create a contrast matrix for specific comparisons
contMatrix <- makeContrasts(ap-ac,
                           cp-cc,
                           hp-hc,
                           sp-sc,
                           levels=design)
contMatrix

# fit the contrasts
fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)

# look at the numbers of DM CpGs at FDR < 0.05
summary(decideTests(fit2))

#       ap - ac cp - cc hp - hc sp - sc
#Down        23   24186      13   13115
#NotSig  771944  696120  771969  753018
#Up          35   51696      20    5869

# get the table of results for the contrast
ann850kSub <- ann850k[match(rownames(mVals),ann850k$Name),
                      c(1:4,12:19,24:ncol(ann850k))]

DMPs_a <- topTable(fit2, num=Inf, coef=1, genelist=ann850kSub)
head(DMPs_a)
write.csv(DMPs_a, "DMPs_ac_ap.csv")

DMPs_c <- topTable(fit2, num=Inf, coef=2, genelist=ann850kSub)
head(DMPs_c)
write.csv(DMPs_a, "DMPs_cc_cp.csv")

DMPs_h <- topTable(fit2, num=Inf, coef=3, genelist=ann850kSub)
head(DMPs_h)
write.csv(DMPs_a, "DMPs_hc_hp.csv")

DMPs_s <- topTable(fit2, num=Inf, coef=4, genelist=ann850kSub)
head(DMPs_s)
write.csv(DMPs_a, "DMPs_sc_sp.csv")

# plot the top 4 most significantly differentially methylated CpGs/tissue
par(mfrow=c(2,2))
sapply(rownames(DMPs_a)[1:4], function(cpg){
  plotCpg(betasNoobBMIQ, cpg=cpg, pheno=targets$Sample_Group, ylab = "Beta values")
})

par(mfrow=c(2,2))
sapply(rownames(DMPs_c)[1:4], function(cpg){
  plotCpg(betasNoobBMIQ, cpg=cpg, pheno=targets$Sample_Group, ylab = "Beta values")
})

par(mfrow=c(2,2))
sapply(rownames(DMPs_h)[1:4], function(cpg){
  plotCpg(betanorm6, cpg=cpg, pheno=targets$Sample_Group, ylab = "Beta values")
})

par(mfrow=c(2,2))
sapply(rownames(DMPs_s)[1:4], function(cpg){
  plotCpg(betanorm6, cpg=cpg, pheno=targets$Sample_Group, ylab = "Beta values")
})

## Compute the betamean and beta difference for each tissue

##### Cortex
cellType <- factor(targets$cellType)
cp <- targets[targets$cellType == "cp",] #Dataframe containing just CP samples
cc <- targets[targets$cellType == "cc",]

betacc <- betasNoobBMIQ[, colnames(betasNoobBMIQ) %in% cc$ID]
betacp <- betasNoobBMIQ[, colnames(betasNoobBMIQ) %in% cp$ID]

betacpmean <- rowMeans(betacp) # Compute mean of every probe accross samples
betaccmean <- rowMeans(betacc)

# Calculate the difference between patients and controls
difference <- (betacpmean-betaccmean)
betaincrease <- cbind(betacpmean, betaccmean, difference)
betaincrease <- as.data.frame(betaincrease)

write.csv(betaincrease, "betaincrease_CORTEX_2.csv")

DMPs_ctotal <- cbind(DMPs_c,betaincrease)
write.csv(DMPs_ctotal, "DMPs_ctotal.csv")

##### Hippocampus
cellType <- factor(targets$cellType)
hp <- targets[targets$cellType == "hp",] #Dataframe containing just HP samples
hc <- targets[targets$cellType == "hc",]

betahp <- betasNoobBMIQ[, colnames(betasNoobBMIQ)%in% hp$ID]
betahc <- betasNoobBMIQ[, colnames(betasNoobBMIQ)%in% hc$ID]

betahpmean <- rowMeans(betahp) # Compute mean of every probe accross samples
betahcmean <- rowMeans(betahc)

# Calculate the difference between patients and controls
difference <- (betahpmean-betahcmean)
betaincrease <- cbind(betahpmean, betahcmean, difference)

write.csv(betaincrease, "betaincrease_HIPPOCAMPUS.csv")

DMPs_htotal <- cbind(DMPs_h, betaincrease)
write.csv(DMPs_htotal, "DMPs_htotal.csv")

#### Amygdala
cellType <- factor(targets$cellType)
ap <- targets[targets$cellType == "ap",]
ac <- targets[targets$cellType == "ac",]

betaap <- betasNoobBMIQ[, colnames(betasNoobBMIQ)%in% ap$ID]
betaac <- betasNoobBMIQ[, colnames(betasNoobBMIQ)%in% ac$ID]

betaapmean <- rowMeans(betaap)
betaacmean <- rowMeans(betaac)

# Calculate the difference between patients and controls
difference <- (betaapmean - betaacmean)
betaincrease <- cbind(betaapmean, betaacmean, difference)

write.csv(betaincrease, "betaincrease_AMYGDALA.csv")

DMPs_atotal <- cbind(DMPs_a, betaincrease)
write.csv(DMPs_atotal, "DMPs_atotal.csv")

##### Blood
cellType <- factor(targets$cellType)
sp <- targets[targets$cellType == "sp",]
sc <- targets[targets$cellType == "sc",]

betasp <- betasNoobBMIQ[, colnames(betasNoobBMIQ)%in% sp$ID]
betasc <- betasNoobBMIQ[, colnames(betasNoobBMIQ)%in% sc$ID]

betaspmean <- rowMeans(betasp)
betascmean <- rowMeans(betasc)

# Calculate the difference between patients and controls
difference <- (betaspmean - betascmean)
betaincrease <- cbind(betaspmean, betascmean, difference)

write.csv(betaincrease, "betaincrease_BLOOD.csv")

DMPs_stotal <- cbind(DMPs_s, betaincrease)
write.csv(DMPs_stotal, "DMPs_stotal.csv")

#--------------------------------------------------------------------#
#------------ Differential methylation analysis of regions ----------# 
#--------------------------------------------------------------------#

## Extract Differential Methylated Regions (DMRs) for each tissue
mVals <- as.matrix(mVals)
##### Amygdala
myAnnotation <- cpg.annotate(object = mVals, datatype = "array", what = "M", 
                             analysis.type = "differential", design = design, 
                             contrasts = TRUE, cont.matrix = contMatrix, 
                             coef = "ap - ac", arraytype = "EPIC")

str(myAnnotation)

DMRs_a <- dmrcate(myAnnotation, lambda=1000, C=2, pcutoff = "fdr")
results.ranges <- extractRanges(DMRs_a)
results.ranges

write.csv(as.data.frame(results.ranges), "DMRs_a.csv")

##### Cortex
myAnnotation <- cpg.annotate(object = mVals, datatype = "array", what = "M", 
                             analysis.type = "differential", design = design, 
                             contrasts = TRUE, cont.matrix = contMatrix, 
                             coef = "cp - cc", arraytype = "EPIC")

str(myAnnotation)

DMRs_c <- dmrcate(myAnnotation, lambda=1000, C=2, pcutoff = "fdr")
results.ranges <- extractRanges(DMRs_c)
results.ranges

write.csv(as.data.frame(results.ranges), "DMRs_c.csv")

##### Hypoccampus
myAnnotation <- cpg.annotate(object = mVals, datatype = "array", what = "M", 
                             analysis.type = "differential", design = design, 
                             contrasts = TRUE, cont.matrix = contMatrix, 
                             coef = "hp - hc", arraytype = "EPIC")

str(myAnnotation)

DMRs_h <- dmrcate(myAnnotation, lambda=1000, C=2, pcutoff = "fdr")
results.ranges <- extractRanges(DMRs_h)
results.ranges

write.csv(as.data.frame(results.ranges), "DMRs_h.csv")

##### Blood
myAnnotation <- cpg.annotate(object = mVals, datatype = "array", what = "M", 
                             analysis.type = "differential", design = design, 
                             contrasts = TRUE, cont.matrix = contMatrix, 
                             coef = "sp - sc", arraytype = "EPIC")

str(myAnnotation)

DMRs_s <- dmrcate(myAnnotation, lambda=1000, C=2, pcutoff = "fdr")
results.ranges <- extractRanges(DMRs_s)
results.ranges

write.csv(as.data.frame(results.ranges), "DMRs_s.csv")

## Finally we check for blood cell type proportions
## FlowSorted.Blood.EPIC for estimate cell type composition on blood samples
library("FlowSorted.Blood.EPIC")
library(ExperimentHub)  
hub <- ExperimentHub()  
query(hub, "FlowSorted.Blood.EPIC")  
FlowSorted.Blood.EPIC <- hub[["EH1136"]]  
head(IDOLOptimizedCpGs) #OptimizedCpGs
FlowSorted.Blood.EPIC 

# Take the blood samples
RGsetTargets <- rgSet[, rgSet$Tissue == "blood"] 
RGsetTarget

# Estimate cell counts
if (memory.limit()>8000){  
  countsEPIC<-estimateCellCounts2(RGsetTargets, compositeCellType = "Blood",   
                                  processMethod = "preprocessNoob",  
                                  probeSelect = "IDOL",  
                                  cellTypes = c("CD8T", "CD4T", "NK", "Bcell",  
                                                "Mono", "Neu"),  
                                  referencePlatform =   
                                    "IlluminaHumanMethylationEPIC",  
                                  referenceset = NULL,  
                                  IDOLOptimizedCpGs =IDOLOptimizedCpGs,   
                                  returnAll = FALSE)  
  
  head(countsEPIC$counts)  
}

# Save results as a data frame
write.csv(as.data.frame(countsEPIC$counts), "cell_distro.csv")

# Plot 
controles <- read.csv("Results/Blood/cells_controles.csv")
rownames(controles) <- controles$X
controles <- controles[,-1]
treatment <- read.csv("Results/Blood/cells_treatment.csv")
rownames(treatment) <- treatment$X
treatment <- treatment[,-1]

controles$group <- rep('control', nrow(controles))
controles <- melt(controles, id.vars = 'group')
treatment$group <- rep('treatment', nrow(treatment))
treatment <- melt(treatment, id.vars = 'group')

allData <- rbind(controles, treatment)

ggplot(allData, aes(x=variable, y=value, color = group)) +
  geom_boxplot() +
  labs(title = "Differences in blood cell types", x="Cell Type", y="Proportion")

# Clean workspace and save Workspace image
rm(ac)
rm(allData)
rm(ap)
rm(betaac)
rm(betaap)
rm(betacc)
rm(betacp)
rm(betahc)
rm(betahp)
rm(betasc)
rm(betasp)
rm(bloodbetas)
rm(cc)
rm(CellTypeMeans450K)
rm(contMatrix)
rm(controles)
rm(countsEPIC)
rm(CorrectedBeta)
rm(cp)
rm(design)
rm(DMPs_a)
rm(DMPs_atotal)
rm(DMPs_c)
rm(DMPs_ctotal)
rm(DMPs_h)
rm(DMPs_htotal)
rm(DMPs_s)
rm(DMPs_stotal)
rm(DMRs_a)
rm(DMRs_c)
rm(DMRs_h)
rm(DMRs_s)
rm(fit)
rm(fit2)
rm(FlowSorted.Blood.EPIC)
rm(hc)
rm(hp)
rm(hub)
rm(sc)
rm(sp)
rm(betaacmean)
rm(betaapmean)
rm(betaccmean)
rm(betacpmean)
rm(betahcmean)
rm(betahpmean)
rm(betascmean)
rm(betaspmean)
rm(betaincrease)
rm(myAnnotation)
rm(results.ranges)
rm(RGsetTargets)
rm(treatment)
rm(difference)
rm(cellType)

# Can run Houseman correction method from saved workspace image

#--------------------------------------------------------------------#
#----------------- Houseman method for blood samples ----------------# 
#--------------------------------------------------------------------#

# Get beta values from blood samples
blood <- targets[targets$Tissue == "blood",] 
betablood <- betasNoobBMIQ[, colnames(betasNoobBMIQ) %in% blood$ID]
write.csv(betablood, "betablood.csv")

# Perform Houseman correction using myRef function from ChAMP library
myRefbase <- champ.refbase(beta=betablood, arraytype="EPIC")
CorrectedBeta <- myRefbase$CorrectedBeta
write.csv(CorrectedBeta, "CorrectedBeta_BLOOD.csv")

# Convert beta values to M values in blood samples
mValsb <- log2(CorrectedBeta / (1 - CorrectedBeta))
mValsb <- replace(mValsb, is.na(mValsb), 0) # Replace NAs for 0s
write.csv(mValsb, "Mvalues_BLOOD_Houseman.csv")

# Get dataframe of betasc + betasp
bloodbetas <- cbind(betasc, betasp)
write.csv(bloodbetas, "bloodbetas.csv")

# Limma
cellType <- factor(blood$cellType)

# Use the above to create a design matrix
design2 <- model.matrix(~0+cellType, data=blood)
colnames(design2) <- levels(cellType)

# Fit the linear model
fit <- lmFit(mValsb, design2)
# Create a contrast matrix for specific comparisons
contMatrix2 <- makeContrasts(sp-sc,
                             levels=design2)
contMatrix2

# fit the contrasts
fit2 <- contrasts.fit(fit, contMatrix2)
fit2 <- eBayes(fit2)

# look at the numbers of DM CpGs at FDR < 0.05
summary(decideTests(fit2))

#sp - sc
#Down         0
#NotSig  772002
#Up           0

# get the table of results for the contrast
ann850kSub <- ann850k[match(rownames(mValsb),ann850k$Name),
                      c(1:4,12:19,24:ncol(ann850k))]

DMPs_shouse <- topTable(fit2, num=Inf, coef=1, sort.by = "B", genelist=ann850kSub)
head(DMPs_shouse)
write.csv(DMPs_shouse, "DMPs_shouse.csv")
