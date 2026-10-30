## Author: Marcos Elizalde Horcada - marcos.elizaldeh@gmail.com
## Date: January 2021
## Revised: October 2026

## Based in A cross-package Bioconductor workflow for analysing methylation
## array (Jovana Maksimovic, Belinda Phipson and Alicia Oshlack, 2020)

## Revision notes (October 2026)
## Bug fixes only: normalization (Noob + BMIQ), probe filters, limma
## contrasts, DMRcate parameters, cell counts (estimateCellCounts2 with IDOL)
## and the Houseman correction (champ.refbase) are unchanged.
##  - DMPs_*total.csv: the beta differences are now joined to the topTable()
##    results by probe ID. cbind() paired them by position, but topTable()
##    sorts probes by significance, so each row got the delta-beta of
##    another probe.
##  - DMPs_cc_cp.csv, DMPs_hc_hp.csv and DMPs_sc_sp.csv saved the amygdala
##    results (DMPs_a); each file now holds its own contrast.
##  - Undefined objects: betanorm6 -> betasNoobBMIQ, RGsetTarget ->
##    RGsetTargets.
##  - Removed the rm() blocks (they deleted betasc/betasp before Houseman)
##    and the bloodbetas object (same data as betablood.csv).
##  - The filtered betas are written to betasNoobBMIQ_filtered.csv and saved
##    in checkpoint_filtered.RData, instead of reading files the script never
##    created. The blood cell-type boxplot is built from countsEPIC$counts
##    instead of Results/Blood/cells_*.csv.
##  - MDS legends now match the colouring variable.
##  - Removed memory.limit() (defunct since R 4.2) and unused packages
##    (knitr, minfiData, Gviz, missMethyl, stringr).
##  - DMPs, delta-betas and DMRs are computed in loops over the tissues;
##    output file names are unchanged.

#--------------------------------------------------------------------#
#------------- Setting up the workspace and data loading ------------#
#--------------------------------------------------------------------#

library(limma)
library(minfi)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(IlluminaHumanMethylationEPICmanifest)
library(RColorBrewer)
library(DMRcate)
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
#legend("top", legend=levels(factor(targets$cellType)), text.col=pal,
#cex=0.7, bg="white")

plotMDS(getM(mSetSq), top=1000, gene.selection="common",
        col=pal[factor(targets$cellType)], dim=c(2,3),
        ylim=c(-1,1), xlim=c(-1,1))
#legend("topleft", legend=levels(factor(targets$cellType)), text.col=pal,
#cex=0.7, bg="white")

plotMDS(getM(mSetSq), top=1000, gene.selection="common",
        col=pal[factor(targets$cellType)], dim=c(3,4),
        xlim=c(-1,1))
#legend("topright", legend=levels(factor(targets$cellType)), text.col=pal,
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

# Back to a matrix (same values) for limma, plotCpg() and DMRcate
betasNoobBMIQ <- as.matrix(betasNoobBMIQ)
write.csv(betasNoobBMIQ, "betasNoobBMIQ_filtered.csv")

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
legend("right", legend=levels(factor(targets$Tissue)), text.col=pal,
       cex=0.7, bg="white")

# Checkpoint: the rest of the script (DMPs, DMRs, cell counts and Houseman)
# can be run from here in a new R session
save(targets, rgSet, betasNoobBMIQ, mVals, file="checkpoint_filtered.RData")

#--------------------------------------------------------------------#
#------------ Probe-wise differential methylation analysis ----------#
#--------------------------------------------------------------------#

# To resume in a new R session: load the libraries and ann850k (first
# section), then run from here
if (!exists("mVals")) load("checkpoint_filtered.RData")

# Tissue codes used in targets$cellType: <code>p = patients, <code>c = controls
tissues <- c(a = "amygdala", c = "cortex", h = "hippocampus", s = "blood")

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

# get the table of results for each contrast
ann850kSub <- ann850k[match(rownames(mVals),ann850k$Name),
                      c(1:4,12:19,24:ncol(ann850k))]

DMPs <- list()
for (tissue in names(tissues)) {
  contrast <- paste0(tissue, "p - ", tissue, "c") # e.g. "ap - ac"
  DMPs[[tissue]] <- topTable(fit2, num=Inf, coef=contrast, genelist=ann850kSub)
  print(head(DMPs[[tissue]]))
  write.csv(DMPs[[tissue]], paste0("DMPs_", tissue, "c_", tissue, "p.csv"))

  # plot the top 4 most significantly differentially methylated CpGs
  par(mfrow=c(2,2))
  plotCpg(betasNoobBMIQ, cpg=rownames(DMPs[[tissue]])[1:4],
          pheno=targets$Sample_Group, ylab = "Beta values")
}

## Compute the betamean and beta difference for each tissue
betaincrease_files <- c(a = "betaincrease_AMYGDALA.csv",
                        c = "betaincrease_CORTEX_2.csv",
                        h = "betaincrease_HIPPOCAMPUS.csv",
                        s = "betaincrease_BLOOD.csv")

betaincrease <- list()
DMPs_total <- list()
for (tissue in names(tissues)) {
  patients <- targets[targets$cellType == paste0(tissue, "p"),]
  controls <- targets[targets$cellType == paste0(tissue, "c"),]

  betap <- betasNoobBMIQ[, colnames(betasNoobBMIQ) %in% patients$ID, drop=FALSE]
  betac <- betasNoobBMIQ[, colnames(betasNoobBMIQ) %in% controls$ID, drop=FALSE]

  betapmean <- rowMeans(betap) # Compute mean of every probe accross samples
  betacmean <- rowMeans(betac)

  # Calculate the difference between patients and controls
  difference <- (betapmean - betacmean)
  betaincrease[[tissue]] <- data.frame(betapmean, betacmean, difference)
  # Same column names as before, e.g. betacpmean, betaccmean, difference
  colnames(betaincrease[[tissue]])[1:2] <- paste0("beta", tissue,
                                                  c("p", "c"), "mean")
  write.csv(betaincrease[[tissue]], betaincrease_files[[tissue]])

  # Join by probe ID: topTable() sorts the probes by significance, so the
  # rows are not in the same order as in betaincrease
  stopifnot(setequal(rownames(DMPs[[tissue]]), rownames(betaincrease[[tissue]])))
  DMPs_total[[tissue]] <- cbind(DMPs[[tissue]],
                                betaincrease[[tissue]][rownames(DMPs[[tissue]]),])
  write.csv(DMPs_total[[tissue]], paste0("DMPs_", tissue, "total.csv"))
}

#--------------------------------------------------------------------#
#------------ Differential methylation analysis of regions ----------#
#--------------------------------------------------------------------#

## Extract Differential Methylated Regions (DMRs) for each tissue
## Note: recent DMRcate versions (with EPIC v2 support) may require
## arraytype = "EPICv1" instead of "EPIC" for these 850k (EPIC v1) arrays
DMRs <- list()
for (tissue in names(tissues)) {
  myAnnotation <- cpg.annotate(object = mVals, datatype = "array", what = "M",
                               analysis.type = "differential", design = design,
                               contrasts = TRUE, cont.matrix = contMatrix,
                               coef = paste0(tissue, "p - ", tissue, "c"),
                               arraytype = "EPIC")

  str(myAnnotation)

  DMRs[[tissue]] <- dmrcate(myAnnotation, lambda=1000, C=2, pcutoff = "fdr")
  results.ranges <- extractRanges(DMRs[[tissue]])
  print(results.ranges)

  write.csv(as.data.frame(results.ranges), paste0("DMRs_", tissue, ".csv"))
}

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
RGsetTargets

# Estimate cell counts (needs about 8 GB of RAM)
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

# Save results as a data frame
write.csv(as.data.frame(countsEPIC$counts), "cell_distro.csv")

# Plot cell type proportions in patients (sp) vs controls (sc)
cellCounts <- as.data.frame(countsEPIC$counts)
cellCounts$group <- factor(targets$cellType[match(rownames(cellCounts),
                                                  targets$ID)],
                           levels = c("sc", "sp"),
                           labels = c("controls", "patients"))
stopifnot(!anyNA(cellCounts$group))
cellCounts <- melt(cellCounts, id.vars = 'group')

print(ggplot(cellCounts, aes(x=variable, y=value, color = group)) +
        geom_boxplot() +
        labs(title = "Differences in blood cell types", x="Cell Type",
             y="Proportion"))

#--------------------------------------------------------------------#
#----------------- Houseman method for blood samples ----------------#
#--------------------------------------------------------------------#

# Can also be run on its own from the checkpoint (see above)
if (!exists("betasNoobBMIQ")) load("checkpoint_filtered.RData")

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
