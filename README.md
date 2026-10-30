# Illumina EPIC (850k) DNA methylation analysis

R analysis script for Illumina Infinium MethylationEPIC (850k) arrays, developed for an epigenome-wide study of drug-resistant temporal lobe epilepsy at Hospital Universitario de La Princesa (Madrid). The study compares patients and controls in four tissues: hippocampus, cortex, amygdala and blood.

The workflow follows the cross-package Bioconductor workflow for methylation arrays by Maksimovic, Phipson and Oshlack, extended with BMIQ normalization and cell-composition correction for blood.

## Workflow

1. **Quality control:** detection p-values, minfi QC report, and removal of samples with mean detection p-value ≥ 0.01.
2. **Normalization:** Noob (minfi) followed by BMIQ (wateRmelon, via ChAMP).
3. **Probe filtering:** removal of probes that failed detection in any sample, probes on sex chromosomes, cross-reactive probes (Chen et al. 2013; Pidsley et al. 2016) and probes with a SNP at the CpG site (dbSNP 147, MAF ≥ 0.01).
4. **Differentially methylated positions (DMPs):** limma on M-values, patients vs. controls in each tissue (FDR < 0.05), plus mean beta values and Δβ per tissue.
5. **Differentially methylated regions (DMRs):** DMRcate (lambda = 1000, C = 2).
6. **Blood cell composition:** cell-type proportions (CD8T, CD4T, NK, B cells, monocytes, neutrophils) estimated with `estimateCellCounts2` and the FlowSorted.Blood.EPIC IDOL reference.
7. **Cell-composition correction in blood:** reference-based (Houseman) correction with ChAMP, followed by a new limma DMP analysis.
8. **Visualization:** detection p-value barplots, density plots, MDS plots, top CpG plots and cell-proportion boxplots.

## Requirements

R with Bioconductor and CRAN packages, mainly: minfi, limma, DMRcate, ChAMP, wateRmelon, FlowSorted.Blood.EPIC, ExperimentHub, IlluminaHumanMethylationEPICanno.ilm10b2.hg19, IlluminaHumanMethylationEPICmanifest, ggplot2, reshape2, readxl and RColorBrewer.

## Input

- IDAT files and sample sheet (`IDAT/` folder)
- Cross-reactive probe lists from Chen et al. (2013) and Pidsley et al. (2016)
- SNP annotation table (`SNPs.147CommonSingle.csv`)

Patient data are not included. The script is organized in sections meant to be run step by step; it is not a packaged pipeline.

## Publication

Sánchez-Jiménez P, Elizalde-Horcada M, et al. DNA Methylation Description of Hippocampus, Cortex, Amygdala, and Blood of Drug-Resistant Temporal Lobe Epilepsy. *Mol Neurobiol*. 2023. https://doi.org/10.1007/s12035-022-03180-z
