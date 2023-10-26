# TODO: edit to match the style guidelines

# ==============================================================================
# Author(s) : Heini M Natri, hnatri@asu.edu
# Date : July 2019
# Description: Analyzing and visualizing the ICGC HCC Japanese data
# ==============================================================================

# ======================================
# Load libraries
# ======================================

library(RColorBrewer)
library(matrixStats)
library(ggplot2)
library(edgeR)
library(DESeq2)
library(limma)
library(doParallel)
library(variancePartition)
library(org.Hs.eg.db)
library(clusterProfiler)
library(GOSemSim)
library(biomaRt)
library(UpSetR)
library(VennDiagram)
library(ggrepel)
library(dplyr)
library(stringr)
library(forcats)

# ======================================
# Environment parameters
# ======================================

# Working directory
setwd("~/3.0 Hasting Research/Liver Cancer Analysis")

# Defining colors
viralPalette <- brewer.pal(8, "Set1")
hbvColor <- viralPalette[1]
hcvColor <- viralPalette[2]
bothColor <- viralPalette[3]
neitherColor <- viralPalette[4]

sexTissuePalette <- brewer.pal(12, "Paired")
maleTumorColor <- sexTissuePalette[4]
maleAdjacentColor <- sexTissuePalette[3]
femaleTumorColor <- sexTissuePalette[6]
femaleAdjacentColor <- sexTissuePalette[5]

# ======================================
# Read in data
# ======================================

metadata <- read.table("metadata_for_de.csv", row.names=1,header=TRUE, sep=",")
tumorAdjacentExp <- read.table("japan_all_samples_salmon_expression_counts.txt", row.names = 1, header=TRUE)
colnames(tumorAdjacentExp) <- gsub("\\.", "-", colnames(tumorAdjacentExp))

# Importing gene annotations
#genes <- read.table("gencode.v25.chr_patch_hapl_scaff.annotation.bed", header=FALSE, sep="\t")
genes <- read.table("~/3.0 Hasting Research/Novel Model/Cohen Melanoma/Penalized_Regression_Data/gencodeTranscripts.txt", header=TRUE, sep="\t")
genes <- data.frame(genes)
tumorAdjacentExp <- tumorAdjacentExp[rownames(tumorAdjacentExp) %in% genes$GENEID ,]
genes <- genes[match(rownames(tumorAdjacentExp), genes$GENEID),]
# Calculating gene length, this is needed for calculating the FPKM values
genes$length <- with(genes, end - start)

# Removing RK023 due to low quality
metadata <- metadata[!(metadata$ID == "RK023") , ]

# Subsetting and ordering metadata to match the count matrix
tumorAdjacentExpSubset <- tumorAdjacentExp[,colnames(tumorAdjacentExp) %in% metadata$sampleid]
metadataSubset <- metadata[metadata$sampleid %in% colnames(tumorAdjacentExpSubset),]
metadataSubset <- metadataSubset[match(colnames(tumorAdjacentExpSubset), metadataSubset$sampleid),]
rownames(metadataSubset) <- metadataSubset$sampleid

# Adding tissue type, converting categorical variables to factors
metadataSubset$tumor <- as.numeric(grepl('tumor', metadataSubset$sampleid, ignore.case=T))
metadataSubset$gender_tissue <- paste(metadataSubset$Gender, metadataSubset$tumor, sep="_")
metadataSubset$gender_tissue_viral <- paste(metadataSubset$gender_tissue, metadataSubset$Virus_infection, sep="_")
metadataSubset$library_type <- metadataSubset$strandedness
metadataSubset$library_type <- factor(metadataSubset$library_type)
metadataSubset$tumor <- factor(metadataSubset$tumor)
metadataSubset$Ta <- factor(metadataSubset$Ta)
metadataSubset$Portal_vein_invasion <- factor(metadataSubset$Portal_vein_invasion)
metadataSubset$Hepatic_vein_invasion <- factor(metadataSubset$Hepatic_vein_invasion)
metadataSubset$Bile_duct_invasion <- factor(metadataSubset$Bile_duct_invasion)
metadataSubset$Liver_fibrosisc <- factor(metadataSubset$Liver_fibrosisc)
metadataSubset$Prognosis <- factor(metadataSubset$Prognosis)


# Creating the DGEList object
dge <- DGEList(counts=tumorAdjacentExpSubset, genes=genes)
colnames(dge) <- colnames(tumorAdjacentExpSubset)
dge$samples$sex <- metadataSubset$Gender
dge$samples$viral <- factor(metadataSubset$Virus_infection)
dge$samples$ID <- metadataSubset$ID
dge$samples$tumor <- metadataSubset$tumor
dge$samples$gender_tissue <- metadataSubset$gender_tissue
dge$samples$gender_tissue_viral <- metadataSubset$gender_tissue_viral
dge$samples$library_type <- metadataSubset$library_type
dge$samples$edmonson_grade <- metadataSubset$Edmondson_grade
dge$samples$Ta <- metadataSubset$Ta
dge$samples$survival <- metadataSubset$Overall_survival_month

# Inspecting the N of samples in each group
table(dge$samples$gender_tissue_viral)

# ======================================
# Filtering expression data
# ======================================

# Keeping genes that have a mean FPKM of at least 0.5 in at least one of the
# groups under investigation and at least 6 reads in at least 10 samples
fpkm <- rpkm(dge, gene.length=dge$genes$length)

M_1_HBV_mean_fpkm <- apply(as.data.frame(fpkm)[(dge$samples$gender_tissue_viral=="M_1_HBV")],1,mean,na.rm=TRUE)
M_0_HBV_mean_fpkm <- apply(as.data.frame(fpkm)[(dge$samples$gender_tissue_viral=="M_0_HBV")],1,mean,na.rm=TRUE)
M_1_HCV_mean_fpkm <- apply(as.data.frame(fpkm)[(dge$samples$gender_tissue_viral=="M_1_HCV")],1,mean,na.rm=TRUE)
M_0_HCV_mean_fpkm <- apply(as.data.frame(fpkm)[(dge$samples$gender_tissue_viral=="M_0_HCV")],1,mean,na.rm=TRUE)
M_1_HBVHCV_mean_fpkm <- apply(as.data.frame(fpkm)[(dge$samples$gender_tissue_viral=="M_1_both")],1,mean,na.rm=TRUE)
M_0_HBVHCV_mean_fpkm <- apply(as.data.frame(fpkm)[(dge$samples$gender_tissue_viral=="M_0_both")],1,mean,na.rm=TRUE)
M_1_NBNC_mean_fpkm <- apply(as.data.frame(fpkm)[(dge$samples$gender_tissue_viral=="M_1_NBNC")],1,mean,na.rm=TRUE)
M_0_NBNC_mean_fpkm <- apply(as.data.frame(fpkm)[(dge$samples$gender_tissue_viral=="M_0_NBNC")],1,mean,na.rm=TRUE)

F_1_HBV_mean_fpkm <- apply(as.data.frame(fpkm)[(dge$samples$gender_tissue_viral=="F_1_HBV")],1,mean,na.rm=TRUE)
F_0_HBV_mean_fpkm <- apply(as.data.frame(fpkm)[(dge$samples$gender_tissue_viral=="F_0_HBV")],1,mean,na.rm=TRUE)
F_1_HCV_mean_fpkm <- apply(as.data.frame(fpkm)[(dge$samples$gender_tissue_viral=="F_1_HCV")],1,mean,na.rm=TRUE)
F_0_HCV_mean_fpkm <- apply(as.data.frame(fpkm)[(dge$samples$gender_tissue_viral=="F_0_HCV")],1,mean,na.rm=TRUE)
F_1_NBNC_mean_fpkm <- apply(as.data.frame(fpkm)[(dge$samples$gender_tissue_viral=="F_1_NBNC")],1,mean,na.rm=TRUE)
F_0_NBNC_mean_fpkm <- apply(as.data.frame(fpkm)[(dge$samples$gender_tissue_viral=="F_0_NBNC")],1,mean,na.rm=TRUE)

keep <- (M_1_HBV_mean_fpkm > 0.5 | M_0_HBV_mean_fpkm > 0.5 | 
           M_1_HCV_mean_fpkm > 0.5 | M_0_HCV_mean_fpkm > 0.5 |
           M_1_HBVHCV_mean_fpkm > 0.5 | M_0_HBVHCV_mean_fpkm > 0.5 |
           M_1_NBNC_mean_fpkm > 0.5 | M_0_NBNC_mean_fpkm > 0.5 |
           F_1_HBV_mean_fpkm > 0.5 | F_0_HBV_mean_fpkm > 0.5 |
           F_1_HCV_mean_fpkm > 0.5 | F_0_HCV_mean_fpkm > 0.5 |
           F_1_NBNC_mean_fpkm > 0.5 | F_0_NBNC_mean_fpkm > 0.5)

dge <- dge[keep,,keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge, method="TMM")
keep <- rowSums(dge$counts > 6) >= 10
dge <- dge[keep,,keep.lib.size=FALSE]
dge <- calcNormFactors(dge, method="TMM")

# N of genes retained after filtering
dim(dge$genes)

# ===========================================================
# ===========================================================
# Analysis of all tumor vs. tumor-adjacent regardless of sex
# ===========================================================
# ===========================================================

# Creating a new design model matrix with the variable of interest and the
# library type
design <- model.matrix(~0+dge$samples$tumor+dge$samples$library_type+dge$samples$Ta)
colnames(design) <- gsub("dge\\$samples\\$tumor", "tumor", colnames(design))
colnames(design) <- gsub("dge\\$samples\\$library_typeunstranded", "library_type", colnames(design))
colnames(design) <- gsub("dge\\$samples\\$Ta2", "Ta2", colnames(design))
colnames(design) <- gsub("dge\\$samples\\$Ta3", "Ta3", colnames(design))
colnames(design) <- gsub("dge\\$samples\\$Ta4", "Ta4", colnames(design))
head(design)

# Running voom again with the new design matrix.
v <- voomWithQualityWeights(dge, design, plot=TRUE)

# ======================================
# Multi-dimensional scaling plot
# ======================================

# MDS plot. The dimensions are in the mds object.
pdf("~/3.0 Hasting Research/Liver Cancer Analysis/DEG FIgures/MDS_figure_comp1_2.pdf")

mds <- plotMDS(v, top = 50, dim.plot = c(1,2), plot=TRUE, cex=2,
               pch=ifelse(v$targets$viral %in% c("HBV"), 17,
                          ifelse(v$targets$viral %in% c("HCV"), 15,
                                 ifelse(v$targets$viral %in% c("both"), 16,  3))),
               col=ifelse(v$targets$gender_tissue=="M_1", maleTumorColor,
                          ifelse(v$targets$gender_tissue=="M_0", maleAdjacentColor,
                                 ifelse(v$targets$gender_tissue=="F_1", femaleTumorColor,
                                        ifelse(v$targets$gender_tissue=="F_0", femaleAdjacentColor, "azure3")))),
               gene.selection = "common") #, gene.selection = "pairwise", labels=tissue,
legend("topleft", pch=c(15), 
       col=c(maleTumorColor, maleAdjacentColor, femaleTumorColor, femaleAdjacentColor),
       legend=c("Male tumor", "Male adjacent", "Female tumor", "Female adjacent"))
legend("topright", pch=c(17, 15, 16, 3),
       col=c("black"),
       legend=c("HBV", "HCV", "Both", "Neither"))

dev.off()

pdf("~/3.0 Hasting Research/Liver Cancer Analysis/DEG FIgures/MDS_figure_comp2_3_colored_by_lib_type.pdf")
mds <- plotMDS(v, top = 50, dim.plot = c(1,2), plot=TRUE, cex=2,
               pch=ifelse(v$targets$viral %in% c("HBV"), 17,
                          ifelse(v$targets$viral %in% c("HCV"), 15,
                                 ifelse(v$targets$viral %in% c("both"), 16,  3))),
               col=ifelse(v$targets$library_type=="stranded", "red", 
                          ifelse(v$targets$library_type=="unstranded", "blue", "black")),
               gene.selection = "common")
legend("topleft", pch=c(15), 
       col=c(maleTumorColor, maleAdjacentColor, femaleTumorColor, femaleAdjacentColor),
       legend=c("Male tumor", "Male adjacent", "Female tumor", "Female adjacent"))
legend("topright", pch=c(17, 15, 16, 3),
       col=c("black"),
       legend=c("HBV", "HCV", "Both", "Neither"))

dev.off()


pdf("~/3.0 Hasting Research/Liver Cancer Analysis/DEG FIgures/MDS_figure_comp2_3.pdf")
mds <- plotMDS(v, top = 50, dim.plot = c(2,3), plot=TRUE, cex=2,
               pch=ifelse(v$targets$viral %in% c("HBV"), 17,
                          ifelse(v$targets$viral %in% c("HCV"), 15,
                                 ifelse(v$targets$viral %in% c("both"), 16,  3))),
               col=ifelse(v$targets$gender_tissue=="M_1", maleTumorColor,
                          ifelse(v$targets$gender_tissue=="M_0", maleAdjacentColor,
                                 ifelse(v$targets$gender_tissue=="F_1", femaleTumorColor,
                                        ifelse(v$targets$gender_tissue=="F_0", femaleAdjacentColor, "azure3")))),
               gene.selection = "common") #, gene.selection = "pairwise", labels=tissue,
legend("topleft", pch=c(15), 
       col=c(maleTumorColor, maleAdjacentColor, femaleTumorColor, femaleAdjacentColor),
       legend=c("Male tumor", "Male adjacent", "Female tumor", "Female adjacent"))
legend("topright", pch=c(17, 15, 16, 3),
       col=c("black"),
       legend=c("HBV", "HCV", "Both", "Neither"))

dev.off()



# Removing batch effects
vCorrectLibtype <- removeBatchEffect(v, batch=v$targets$library_type)


#pdf("~/3.0 Hasting Research/Liver Cancer Analysis/DEG FIgures/MDS_figure_comp1_2_afteradj.pdf")
mds <- plotMDS(vCorrectLibtype, top = 50, ndim = 10, dim.plot = c(1,2), plot=TRUE, cex=2,
               pch=ifelse(v$targets$viral %in% c("HBV"), 17,
                          ifelse(v$targets$viral %in% c("HCV"), 15,
                                 ifelse(v$targets$viral %in% c("both"), 16,  3))),
               col=ifelse(v$targets$gender_tissue=="M_1",
                          maleTumorColor,
                          ifelse(v$targets$gender_tissue=="M_0",
                                 maleAdjacentColor,
                                 ifelse(v$targets$gender_tissue=="F_1",
                                        femaleTumorColor, femaleAdjacentColor))),
               gene.selection = "common")
legend("topleft", pch=c(15), 
       col=c(maleTumorColor, maleAdjacentColor, femaleTumorColor, femaleAdjacentColor),
       legend=c("Male tumor", "Male adjacent", "Female tumor", "Female adjacent"))
legend("topright", pch=c(17, 15, 16, 3),
       col=c("black"),
       legend=c("HBV", "HCV", "Both", "Neither"))
#dev.off()


# ======================================
# PCA
# ======================================

# Select most variable genes based on coefficient of variance (mean scaled)
# Voom transformed counts
voomCounts <- v$E
voomCountsMatrix <- data.matrix(voomCounts, rownames.force = NA)
ntop = length(dge$genes$TXNAME)

means <- rowMeans(voomCountsMatrix)
Pvars <- rowVars(voomCountsMatrix)
cv2 <- Pvars/means^2
select <- order(cv2, decreasing = TRUE)[seq_len(min(ntop, length(cv2)))]
head(select)
highly_variable_exp <- ((voomCountsMatrix)[select, ])
dim(highly_variable_exp)

# Running PCA
pca_exp <- prcomp(t(highly_variable_exp), scale=T, center=T)
head(pca_exp$x)

# Dataframe with the first 10 PCs
dim1_10 <- data.frame(pca_exp$x[,1:10])
# Adding metadata
pcaWithMetadata <- merge(dim1_10, metadataSubset, by=0, all=TRUE)
pcaWithMetadata$Virus.infection <- factor(pcaWithMetadata$Virus_infection,
                                          levels=c("HBV", "HCV", "both", "NBNC", NA))
pcaWithMetadata$gender_tissue <- factor(pcaWithMetadata$gender_tissue,
                                        levels=c("M_1", "M_0", "F_1", "F_0"))

# Plotting
pdf("~/3.0 Hasting Research/Liver Cancer Analysis/DEG FIgures/PCA_with_metadata.pdf", width=10, height=10)
ggplot(data=pcaWithMetadata, aes(x = PC1, y = PC2, shape = Virus.infection,
                                 color = gender_tissue)) +
  geom_point(size = 8) +
  theme_bw() +
  scale_color_manual(values=c(maleTumorColor, maleAdjacentColor,
                              femaleTumorColor, femaleAdjacentColor,
                              "azure3")) +
  theme(plot.title = element_text(size = 12, face = "bold"),
        legend.title=element_text(size = 30),
        legend.text=element_text(size = 30),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        axis.title.x = element_text(size = 22),
        axis.title.y = element_text(size = 22)) +
  guides(color = guide_legend(order = 2),
         shape = guide_legend(order = 1)) +
  theme(legend.title = element_blank()) +
  xlab("PC1") +
  ylab("PC2")

dev.off()
# ======================================
# Variance explained
# ======================================

# Variance explained by variancePartition
metadataSubset$Alcohol_intake <- as.factor(metadataSubset$Alcohol_intake)
metadataSubset$Smoking <- as.factor(metadataSubset$Smoking)
# Specifying variables to consider
form <- ~ 
  Tumor_size_mm + 
  Overall_survival_month + 
  (1|ID) +
  (1|tumor) + 
  (1|Gender) + 
  (1|Virus_infection) +
  (1|Ta) + 
  (1|Edmondson_grade) + 
  (1|Portal_vein_invasion) + 
  (1|Hepatic_vein_invasion) + 
  (1|Bile_duct_invasion) + 
  (1|Liver_fibrosisc) + 
  (1|Alcohol_intake) + 
  (1|Smoking) + 
  (1|Prognosis)

# Running parallel
cl <- makeCluster(6)
registerDoParallel(cl)

# Fitting the linear model
library(variancePartition)

varPart <- fitExtractVarPartModel(voomCounts, form, metadataSubset)

# Sorting variables by median fraction of variance explained
vp <- sortCols(varPart)

# Violin plot of contribution of each variable to total variance
pdf("~/3.0 Hasting Research/Liver Cancer Analysis/DEG FIgures/VarPart.pdf", width=10, height=10)
plotVarPart(vp)
dev.off()

# Stop parallel
registerDoSEQ()
stopCluster(cl)

# ======================================
# Checking male vs. female grade and stage to see if these should be added as cofactors
# ======================================
metadataSubset_female <- metadataSubset[which(metadataSubset$Gender=="F"),]
metadataSubset_male <- metadataSubset[which(metadataSubset$Gender=="M"),]

#metadataSubset_male$Ta <- log(as.numeric(metadataSubset_male$Ta), 10)
#metadataSubset_female$Ta <- log(as.numeric(metadataSubset_female$Ta),10)
boxplot(metadataSubset_female$Ta, metadataSubset_male$Ta)
#t.test(as.numeric(metadataSubset_female$Ta), as.numeric(metadataSubset_male$Ta), alternative = "two.sided", paired=FALSE)
wilcox.test(as.numeric(metadataSubset_female$Ta), as.numeric(metadataSubset_male$Ta))
pdf("~/3.0 Hasting Research/Liver Cancer Analysis/Survival/Tumor_stage_boxplot.pdf")
boxplot(metadataSubset_female$Ta, metadataSubset_male$Ta, names=c("Female", "Male"))
dev.off()

metadataSubset_female$Edmondson_grade <- gsub("1~2", "1.5", metadataSubset_female$Edmondson_grade)
metadataSubset_female$Edmondson_grade <- gsub("2~3", "2.5", metadataSubset_female$Edmondson_grade)
metadataSubset_female$Edmondson_grade <- as.numeric(metadataSubset_female$Edmondson_grade)

metadataSubset_male$Edmondson_grade <- gsub("1~2", "1.5", metadataSubset_male$Edmondson_grade)
metadataSubset_male$Edmondson_grade <- gsub("2~3", "2.5", metadataSubset_male$Edmondson_grade)
metadataSubset_male$Edmondson_grade <- as.numeric(metadataSubset_male$Edmondson_grade)
boxplot(metadataSubset_female$Edmondson_grade, metadataSubset_male$Edmondson_grade)
t.test(metadataSubset_female$Edmondson_grade, metadataSubset_male$Edmondson_grade, alternative = "two.sided", paired=FALSE)
pdf("~/3.0 Hasting Research/Liver Cancer Analysis/Survival/Edmonson_grade_boxplot.pdf")
boxplot(metadataSubset_female$Edmondson_grade, metadataSubset_male$Edmondson_grade, names=c("Female", "Male"))
dev.off()


# ============================================================================================
# Differential expression analysis with limma - all male tumor adjacent vs. non-tumor-adjacent
# ===========================================================================================

# Block design for individual. This is used in tumor-normal comparisons with
# paired samples.
corfit <- duplicateCorrelation(v, design, block = v$targets$ID)
# This should give a positive correlation value. It represents the
# correlation between measurements made on the same person.
corfit$consensus

# Fitting the linear model with limma.
# If using paired samples, the within-patient correlation and a block design
# for patient is used to account for pairwise samples
fit <- lmFit(v, design, block = v$targets$ID, correlation = corfit$consensus)

# Contrast design for differential expression
# Defining pairwise comparisons
contrasts <- makeContrasts(Adjacent_vs_Tumor = tumor1 - tumor0,
                           levels=colnames(design))

head(contrasts)

# Assigning all comparisons to a vector for later
allComparisons <- colnames(contrasts)

# Running contrast analysis
vfit <- contrasts.fit(fit, contrasts = contrasts)
# Looking at N of DEGs with adj. p <0.01 and log2FC>2
summary(decideTests(vfit, adjust.method = "BH", p.value = 0.01, lfc = 2))

# Computing differential expression based on the empirical Bayes moderation of
# the standard errors towards a common value. Robust = should the estimation of
# the empirical Bayes prior parameters be robustified against outlier sample
# variances?
veBayesFit <- eBayes(vfit, robust=TRUE)
plotSA(veBayesFit, main = "Final model: Mean-variance trend")

vTopTable <- topTable(veBayesFit, n=Inf, p.value=1, lfc=0)

DEGs <- topTable(veBayesFit, n=Inf, p.value=0.01, lfc=2)

# ===========================================
#Volcano plot of tumor vs tumor-adjacent
# ===========================================

df <- data.frame(vTopTable$adj.P.Val, vTopTable$logFC, vTopTable$chr, vTopTable$GENEID, vTopTable$gene_name)
colnames(df) <- c("adj.P.Val", "logFC", "chr", "id", "name")
dfSig <- df[(abs(df$logFC) >= 2 & df$adj.P.Val <= 0.01),]$id


dfAnons <- subset(df, chr != "chrX" & chr != "chrY" & !(id %in% dfSig))
dfAnons <- cbind(dfAnons , rep(1, nrow(dfAnons)))
colnames(dfAnons)[6] <- "Color"

dfXnons <- subset(df, chr == "chrX" & !(id %in% dfSig))
dfXnons <- cbind(dfXnons, rep(2, nrow(dfXnons)))
colnames(dfXnons)[6] <- "Color"

dfYnons <- subset(df, chr == "chrY" & !(id %in% dfSig))
dfYnons <- cbind(dfYnons, rep(3, nrow(dfYnons)))
colnames(dfYnons)[6] <- "Color"

dfA <- subset(df, chr != "chrX" & chr != "chrY" & id %in% dfSig)
dfA <- cbind(dfA, rep(4, nrow(dfA)))
colnames(dfA)[6] <- "Color"

dfX <- subset(df, chr == "chrX" & id %in% dfSig)
dfX <- cbind(dfX, rep(5, nrow(dfX)))
colnames(dfX)[6] <- "Color"

dfY <- subset(df, chr == "chrY" & id %in% dfSig)
dfY <- cbind(dfY, rep(6, nrow(dfY)))
colnames(dfY)[6] <- "Color"

dfPlot <- rbind(dfAnons, dfXnons, dfYnons, dfA, dfX, dfY)
dfPlot$Color <- as.factor(dfPlot$Color)

# Constructing the plot object. The colors will not work the way they should if
# one of the groups doesn't have any genes
p <- ggplot(data = dfPlot, aes(x = logFC, y = -log10(adj.P.Val), color=Color )) +
  geom_point(alpha = 0.5, size = 8) +
  theme_bw() +
  theme(legend.position = "none") +
  xlim(c(-15, 15)) + ylim(c(0, 130)) +
  scale_color_manual(values = c("azure3", "pink", "seagreen2", "black", "mediumvioletred", "springgreen")) +
  labs(x=expression(log[2](FC)),
       y=expression(-log[10] ~ "(FDR-adjusted " ~ italic("p") ~ "-value)")) +
  theme(axis.title.x=element_text(size=12), 
        axis.text.x=element_text(size=10)) +
  theme(axis.title.y=element_text(size=12),
        axis.text.y=element_text(size=10))
p

forLabel <- subset(dfPlot, adj.P.Val<=0.01 & abs(logFC)>=2)

# Adding lines for significance thresholds
pdf("~/3.0 Hasting Research/Liver Cancer Analysis/DEG Figures/Volcano_plot_tumor_tumoradjacent_DEGs.pdf", width=12, height=12)
p + geom_hline(yintercept = 2, colour="#000000", linetype="dashed"
) + geom_vline(xintercept = 2, colour="#000000", linetype="dashed"
) + geom_vline(xintercept = -2, colour="#000000", linetype="dashed"
) + geom_text_repel(data=forLabel, max.iter=1000, box.padding = 0.25, force=1, aes(x = logFC, y = -log10(adj.P.Val), label=name, color=Color), size=8)
dev.off()

# =================
# Enriched pathways
# =================

degResult_genes <- DEGs
degResult_genes$hgnc_symbol <- degResult_genes$name
ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
biomart <- getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol', 'entrezgene_id', 'kegg_enzyme'), filters = 'hgnc_symbol', values = degResult_genes$gene_name, mart = ensembl)
colnames(degResult_genes) <- c("chr", "start", "end", "TXNAME", "GENEID", "hgnc_symbol", "length", 
                               "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")
degResult_genes_biomart <- merge(degResult_genes, biomart, by="hgnc_symbol")
head(degResult_genes_biomart)
degResult <- DEGs
degResult$hgnc_symbol <- degResult$gene_name
degResult_biomart <- merge(degResult, biomart, by="hgnc_symbol")
degResult_biomart <- degResult_biomart[complete.cases(degResult_biomart), ]
geneList <- degResult_biomart$entrezgene_id
geneList_quant <- degResult_biomart$logFC
ego <- enrichGO(geneList, OrgDb = "org.Hs.eg.db", ont="BP", readable=TRUE)
pdf("~/3.0 Hasting Research/Liver Cancer Analysis/GO_KEGG Figures/dotplot_all_tumor_vs_tumor_adjacent.pdf", width=12, height=12)
dotplot(ego, showCategory=30)
dev.off()


# ===========================================================
# ===========================================================
# Analysis of tumor vs. non-tumor differentiated by sex
# ===========================================================
# ===========================================================

# ======================================
# voom transformation
# ======================================

# Creating a design model matrix with the variable of interest
design <- model.matrix(~0+dge$samples$gender_tissue_viral)
colnames(design) <- gsub("dge\\$samples\\$gender_tissue_viral", "", colnames(design))
head(design)

# Running voom with quality weights. Normalizes expression intensities so that
# the log-ratios have similar distributions across a set of samples.
# To quantile normalize, add normalize.method="quantile"
# Running parallel

v <- voomWithQualityWeights(dge, design, plot=TRUE)


# ============================================================================================
# Differential expression analysis with limma - all male tumor adjacent vs. non-tumor-adjacent
# ===========================================================================================

# Creating a new design model matrix with the variable of interest and the
# first dimension of the MDS
design <- model.matrix(~0+v$targets$gender_tissue+v$targets$library_type+v$targets$Ta)
colnames(design) <- gsub("v\\$targets\\$gender_tissue", "", colnames(design))
colnames(design) <- gsub("v\\$targets\\$library_typeunstranded", "library_type", colnames(design))
colnames(design) <- gsub("v\\$targets\\$Ta2", "Ta2", colnames(design))
colnames(design) <- gsub("v\\$targets\\$Ta3", "Ta3", colnames(design))
colnames(design) <- gsub("v\\$targets\\$Ta4", "Ta4", colnames(design))
head(design)

# Running voom again with the new design matrix.
v <- voomWithQualityWeights(dge, design, plot=TRUE)

# Block design for individual. This is used in tumor-normal comparisons with
# paired samples.
corfit <- duplicateCorrelation(v, design, block = v$targets$ID)
# This should give a positive correlation value. It represents the
# correlation between measurements made on the same person.
corfit$consensus

# Fitting the linear model with limma.
# If using paired samples, the within-patient correlation and a block design
# for patient is used to account for pairwise samples
fit <- lmFit(v, design, block = v$targets$ID, correlation = corfit$consensus)

# Contrast design for differential expression
# Defining pairwise comparisons
contrasts <- makeContrasts(Adjacent_M_vs_Tumor_M = M_1 - M_0,
                           Adjacent_F_vs_Tumor_F = F_1 - F_0,
                           levels=colnames(design))

head(contrasts)

# Assigning all comparisons to a vector for later
allComparisons <- colnames(contrasts)

# Running contrast analysis
vfit <- contrasts.fit(fit, contrasts = contrasts)
# Looking at N of DEGs with adj. p <0.01 and log2FC>2
summary(decideTests(vfit, adjust.method = "BH", p.value = 0.01, lfc = 2))

# Computing differential expression based on the empirical Bayes moderation of
# the standard errors towards a common value. Robust = should the estimation of
# the empirical Bayes prior parameters be robustified against outlier sample
# variances?
veBayesFit <- eBayes(vfit, robust=TRUE)
plotSA(veBayesFit, main = "Final model: Mean-variance trend")

vTopTable_M <- topTable(veBayesFit, coef=1, n=Inf, p.value=1, lfc=0)
vTopTable_F <- topTable(veBayesFit, coef=2, n=Inf, p.value=1, lfc=0)

DEGs_M <- topTable(veBayesFit, coef=1, n=Inf, p.value=0.01, lfc=2)
DEGs_F <- topTable(veBayesFit, coef=2, n=Inf, p.value=0.01, lfc=2)
DEGs_F_relax_p <- topTable(veBayesFit, coef=2, n=Inf, p.value=0.1, lfc=2)

# ===========================================
#Volcano plot of male tumor vs tumor-adjacent
# ===========================================

df <- data.frame(vTopTable_M$adj.P.Val, vTopTable_M$logFC, vTopTable_M$chr, vTopTable_M$GENEID, vTopTable_M$gene_name)
colnames(df) <- c("adj.P.Val", "logFC", "chr", "id", "name")
dfSig <- df[(abs(df$logFC) >= 2 & df$adj.P.Val <= 0.01),]$id


dfAnons <- subset(df, chr != "chrX" & chr != "chrY" & !(id %in% dfSig))
dfAnons <- cbind(dfAnons , rep(1, nrow(dfAnons)))
colnames(dfAnons)[6] <- "Color"

dfXnons <- subset(df, chr == "chrX" & !(id %in% dfSig))
dfXnons <- cbind(dfXnons, rep(2, nrow(dfXnons)))
colnames(dfXnons)[6] <- "Color"

dfYnons <- subset(df, chr == "chrY" & !(id %in% dfSig))
dfYnons <- cbind(dfYnons, rep(3, nrow(dfYnons)))
colnames(dfYnons)[6] <- "Color"

dfA <- subset(df, chr != "chrX" & chr != "chrY" & id %in% dfSig)
dfA <- cbind(dfA, rep(4, nrow(dfA)))
colnames(dfA)[6] <- "Color"

dfX <- subset(df, chr == "chrX" & id %in% dfSig)
dfX <- cbind(dfX, rep(5, nrow(dfX)))
colnames(dfX)[6] <- "Color"

dfY <- subset(df, chr == "chrY" & id %in% dfSig)
dfY <- cbind(dfY, rep(6, nrow(dfY)))
colnames(dfY)[6] <- "Color"

dfPlot <- rbind(dfAnons, dfXnons, dfYnons, dfA, dfX, dfY)
dfPlot$Color <- as.factor(dfPlot$Color)

# Constructing the plot object. The colors will not work the way they should if
# one of the groups doesn't have any genes
p <- ggplot(data = dfPlot, aes(x = logFC, y = -log10(adj.P.Val), color=Color )) +
  geom_point(alpha = 0.5, size = 8) +
  theme_bw() +
  theme(legend.position = "none") +
  xlim(c(-15, 15)) + ylim(c(0, 120)) +
  scale_color_manual(values = c("azure3", "pink", "seagreen2", "black", "mediumvioletred", "springgreen")) +
  labs(x=expression(log[2](FC)),
       y=expression(-log[10] ~ "(FDR-adjusted " ~ italic("p") ~ "-value)")) +
  theme(axis.title.x=element_text(size=12), 
        axis.text.x=element_text(size=10)) +
  theme(axis.title.y=element_text(size=12),
        axis.text.y=element_text(size=10))
p

forLabel <- subset(dfPlot, adj.P.Val<=0.01 & abs(logFC)>=2)

# Adding lines for significance thresholds
pdf("~/3.0 Hasting Research/Liver Cancer Analysis/DEG Figures/Volcano_plot_male_DEGs.pdf", width=12, height=12)
p + geom_hline(yintercept = 2, colour="#000000", linetype="dashed"
) + geom_vline(xintercept = 2, colour="#000000", linetype="dashed"
) + geom_vline(xintercept = -2, colour="#000000", linetype="dashed"
) + geom_text_repel(data=forLabel, max.iter=1000, box.padding = 0.25, force=1, aes(x = logFC, y = -log10(adj.P.Val), label=name, color=Color), size=8)
dev.off()


# ===========================================
#Volcano plot of female tumor vs tumor-adjacent
# ===========================================

df <- data.frame(vTopTable_F$adj.P.Val, vTopTable_F$logFC, vTopTable_F$chr, vTopTable_F$GENEID, vTopTable_F$gene_name)
colnames(df) <- c("adj.P.Val", "logFC", "chr", "id", "name")
dfSig <- df[(abs(df$logFC) >= 2 & df$adj.P.Val <= 0.01),]$id


dfAnons <- subset(df, chr != "chrX" & chr != "chrY" & !(id %in% dfSig))
dfAnons <- cbind(dfAnons , rep(1, nrow(dfAnons)))
colnames(dfAnons)[6] <- "Color"

dfXnons <- subset(df, chr == "chrX" & !(id %in% dfSig))
dfXnons <- cbind(dfXnons, rep(2, nrow(dfXnons)))
colnames(dfXnons)[6] <- "Color"

dfYnons <- subset(df, chr == "chrY" & !(id %in% dfSig))
dfYnons <- cbind(dfYnons, rep(3, nrow(dfYnons)))
colnames(dfYnons)[6] <- "Color"

dfA <- subset(df, chr != "chrX" & chr != "chrY" & id %in% dfSig)
dfA <- cbind(dfA, rep(4, nrow(dfA)))
colnames(dfA)[6] <- "Color"

dfX <- subset(df, chr == "chrX" & id %in% dfSig)
dfX <- cbind(dfX, rep(5, nrow(dfX)))
colnames(dfX)[6] <- "Color"

dfY <- subset(df, chr == "chrY" & id %in% dfSig)
dfY <- cbind(dfY, rep(6, nrow(dfY)))
colnames(dfY)[6] <- "Color"

dfPlot <- rbind(dfAnons, dfXnons, dfYnons, dfA, dfX, dfY)
dfPlot$Color <- as.factor(dfPlot$Color)

# Constructing the plot object. The colors will not work the way they should if
# one of the groups doesn't have any genes
p <- ggplot(data = dfPlot, aes(x = logFC, y = -log10(adj.P.Val), color=Color )) +
  geom_point(alpha = 0.5, size = 8) +
  theme_bw() +
  theme(legend.position = "none") +
  xlim(c(-15, 15)) + ylim(c(0, 50)) +
  scale_color_manual(values = c("azure3", "pink", "seagreen2", "black", "mediumvioletred", "springgreen")) +
  labs(x=expression(log[2](FC)),
       y=expression(-log[10] ~ "(FDR-adjusted " ~ italic("p") ~ "-value)")) +
  theme(axis.title.x=element_text(size=12), 
        axis.text.x=element_text(size=10)) +
  theme(axis.title.y=element_text(size=12),
        axis.text.y=element_text(size=10))
p

forLabel <- subset(dfPlot, adj.P.Val<=0.01 & abs(logFC)>=2)

# Adding lines for significance thresholds
pdf("~/3.0 Hasting Research/Liver Cancer Analysis/DEG Figures/Volcano_plot_female_DEGs.pdf", width=12, height=12)
p + geom_hline(yintercept = 2, colour="#000000", linetype="dashed"
) + geom_vline(xintercept = 2, colour="#000000", linetype="dashed"
) + geom_vline(xintercept = -2, colour="#000000", linetype="dashed"
) + geom_text_repel(data=forLabel, max.iter=1000, box.padding = 0.25, force=1, aes(x = logFC, y = -log10(adj.P.Val), label=name, color=Color), size=8)
dev.off()

library(VennDiagram)
venn.diagram(List("Female"=DEGs_F$gene_name, "Male"=DEGs_M$gene_name),  filename = "~/3.0 Hasting Research/Liver Cancer Analysis/DEG Figures/Overlap_DEGs.png")
venn.diagram(List("Female"=DEGs_F_relax_p$gene_name, "Male"=DEGs_M$gene_name),  filename = "~/3.0 Hasting Research/Liver Cancer Analysis/DEG Figures/Overlap_DEGs_Frelaxp.png")

write.csv(DEGs_F, "~/3.0 Hasting Research/Liver Cancer Analysis/DEG Figures/Female_all_DEGs.csv")
write.csv(DEGs_M, "~/3.0 Hasting Research/Liver Cancer Analysis/DEG Figures/Male_all_DEGs.csv")

# ===================================
# Enriched pathways males vs. females
# ===================================

path <- "~/3.0 Hasting Research/Liver Cancer Analysis/DEGs/DEGs_newcounts_Salmon_batchMD1cov_fpkm05_HBVAdjacent_M_vs_HBVAdjacent_F_fdr1_lfc0.txt"
degResult_genes <- read.table(path)
degResult_genes$hgnc_symbol <- degResult_genes$name
ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
biomart <- getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol', 'entrezgene_id', 'kegg_enzyme'), filters = 'hgnc_symbol', values = degResult_genes$gene_name, mart = ensembl)
colnames(degResult_genes) <- c("chr", "start", "end", "TXNAME", "GENEID", "hgnc_symbol", "length", 
                               "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")
degResult_genes_biomart <- merge(degResult_genes, biomart, by="hgnc_symbol")
head(degResult_genes_biomart)
degResult <- DEGs_F
degResult$hgnc_symbol <- degResult$gene_name
degResult_biomart <- merge(degResult, biomart, by="hgnc_symbol")
degResult_biomart <- degResult_biomart[complete.cases(degResult_biomart), ]
geneList <- degResult_biomart$entrezgene_id
geneList_quant <- degResult_biomart$logFC
ego <- enrichGO(geneList, OrgDb = "org.Hs.eg.db", ont="BP", readable=TRUE)
pdf("~/3.0 Hasting Research/Liver Cancer Analysis/GO_KEGG Figures/Goplot_all_female.pdf", width=12, height=12)
goplot(ego)
dev.off()
pdf("~/3.0 Hasting Research/Liver Cancer Analysis/GO_KEGG Figures/dotplot_all_female.pdf", width=12, height=12)
dotplot(ego, showCategory=30)
dev.off()

degResult <- DEGs_M
degResult$hgnc_symbol <- degResult$gene_name
degResult_biomart <- merge(degResult, biomart, by="hgnc_symbol")
degResult_biomart <- degResult_biomart[complete.cases(degResult_biomart), ]
geneList <- degResult_biomart$entrezgene_id
geneList_quant <- degResult_biomart$logFC
ego <- enrichGO(geneList, OrgDb = "org.Hs.eg.db", ont="BP", readable=TRUE)
pdf("~/3.0 Hasting Research/Liver Cancer Analysis/GO_KEGG Figures/Goplot_all_male.pdf", width=12, height=12)
goplot(ego)
dev.off()
pdf("~/3.0 Hasting Research/Liver Cancer Analysis/GO_KEGG Figures/dotplot_all_male.pdf", width=12, height=12)
dotplot(ego, showCategory=30)
dev.off()

# =========================================================================================================
# Differential expression analysis with limma - tumor adjacent vs. non-tumor-adjacent seperated by etiology
# =========================================================================================================

# Creating a new design model matrix with the variable of interest and the
# first dimension of the MDS
design <- model.matrix(~0+v$targets$gender_tissue_viral+v$targets$library_type+v$targets$Ta)
colnames(design) <- gsub("v\\$targets\\$gender_tissue_viral", "", colnames(design))
colnames(design) <- gsub("v\\$targets\\$library_typeunstranded", "library_type", colnames(design))
colnames(design) <- gsub("v\\$targets\\$Ta2", "Ta2", colnames(design))
colnames(design) <- gsub("v\\$targets\\$Ta3", "Ta3", colnames(design))
colnames(design) <- gsub("v\\$targets\\$Ta4", "Ta4", colnames(design))
head(design)

# Running voom again with the new design matrix.
v <- voomWithQualityWeights(dge, design, plot=TRUE)

# Block design for individual. This is used in tumor-normal comparisons with
# paired samples.
corfit <- duplicateCorrelation(v, design, block = v$targets$ID)
# This should give a positive correlation value. It represents the
# correlation between measurements made on the same person.
corfit$consensus

# Fitting the linear model with limma.
# If using paired samples, the within-patient correlation and a block design
# for patient is used to account for pairwise samples
fit <- lmFit(v, design, block = v$targets$ID, correlation = corfit$consensus)

# Contrast design for differential expression
# Defining pairwise comparisons
contrasts <- makeContrasts(HBVAdjacent_M_vs_HBVTumor_M = M_1_HBV - M_0_HBV,
                           HCVAdjacent_M_vs_HCVTumor_M = M_1_HCV - M_0_HCV,
                           NeitherAdjacent_M_vs_NeitherTumor_M = M_1_NBNC - M_0_NBNC,
                           HBVAdjacent_F_vs_HBVTumor_F = F_1_HBV - F_0_HBV,
                           HCVAdjacent_F_vs_HCVTumor_F = F_1_HCV - F_0_HCV,
                           NeitherAdjacent_F_vs_NeitherTumor_F = F_1_NBNC - F_0_NBNC,
                           levels=colnames(design))

head(contrasts)

# Assigning all comparisons to a vector for later
allComparisons <- colnames(contrasts)

# Running contrast analysis
vfit <- contrasts.fit(fit, contrasts = contrasts)
# Looking at N of DEGs with adj. p <0.01 and log2FC>2
summary(decideTests(vfit, adjust.method = "BH", p.value = 0.01, lfc = 2))

# Computing differential expression based on the empirical Bayes moderation of
# the standard errors towards a common value. Robust = should the estimation of
# the empirical Bayes prior parameters be robustified against outlier sample
# variances?
veBayesFit <- eBayes(vfit, robust=TRUE)
plotSA(veBayesFit, main = "Final model: Mean-variance trend")

vTopTable_M_HBV <- topTable(veBayesFit, coef=1, n=Inf, p.value=1, lfc=0)
vTopTable_M_HCV <- topTable(veBayesFit, coef=2, n=Inf, p.value=1, lfc=0)
vTopTable_M_Neither <- topTable(veBayesFit, coef=3, n=Inf, p.value=1, lfc=0)
vTopTable_F_HBV <- topTable(veBayesFit, coef=4, n=Inf, p.value=1, lfc=0)
vTopTable_F_HCV <- topTable(veBayesFit, coef=5, n=Inf, p.value=1, lfc=0)
vTopTable_F_Neither <- topTable(veBayesFit, coef=6, n=Inf, p.value=1, lfc=0)

DEGs_M_HBV <- topTable(veBayesFit, coef=1, n=Inf, p.value=0.01, lfc=2)
DEGs_M_HCV <- topTable(veBayesFit, coef=2, n=Inf, p.value=0.01, lfc=2)
DEGs_M_Neither <- topTable(veBayesFit, coef=3, n=Inf, p.value=0.01, lfc=2)
DEGs_F_HBV <- topTable(veBayesFit, coef=4, n=Inf, p.value=0.01, lfc=2)
DEGs_F_HBV_relax_p <- topTable(veBayesFit, coef=4, n=Inf, p.value=0.1, lfc=2)
DEGs_F_HCV <- topTable(veBayesFit, coef=5, n=Inf, p.value=0.01, lfc=2)
DEGs_F_HCV_relax_p <- topTable(veBayesFit, coef=5, n=Inf, p.value=0.1, lfc=2)
DEGs_F_Neither <- topTable(veBayesFit, coef=6, n=Inf, p.value=0.01, lfc=2)
DEGs_F_Neither_relax_p <- topTable(veBayesFit, coef=6, n=Inf, p.value=0.1, lfc=2)

# ===============================================
#Volcano plot of male HBV tumor vs tumor-adjacent
# ===============================================

df <- data.frame(vTopTable_M_HBV$adj.P.Val, vTopTable_M_HBV$logFC, vTopTable_M_HBV$chr, vTopTable_M_HBV$GENEID, vTopTable_M_HBV$gene_name)
colnames(df) <- c("adj.P.Val", "logFC", "chr", "id", "name")
dfSig <- df[(abs(df$logFC) >= 2 & df$adj.P.Val <= 0.01),]$id


dfAnons <- subset(df, chr != "chrX" & chr != "chrY" & !(id %in% dfSig))
dfAnons <- cbind(dfAnons , rep(1, nrow(dfAnons)))
colnames(dfAnons)[6] <- "Color"

dfXnons <- subset(df, chr == "chrX" & !(id %in% dfSig))
dfXnons <- cbind(dfXnons, rep(2, nrow(dfXnons)))
colnames(dfXnons)[6] <- "Color"

dfYnons <- subset(df, chr == "chrY" & !(id %in% dfSig))
dfYnons <- cbind(dfYnons, rep(3, nrow(dfYnons)))
colnames(dfYnons)[6] <- "Color"

dfA <- subset(df, chr != "chrX" & chr != "chrY" & id %in% dfSig)
dfA <- cbind(dfA, rep(4, nrow(dfA)))
colnames(dfA)[6] <- "Color"

dfX <- subset(df, chr == "chrX" & id %in% dfSig)
dfX <- cbind(dfX, rep(5, nrow(dfX)))
colnames(dfX)[6] <- "Color"

dfY <- subset(df, chr == "chrY" & id %in% dfSig)
dfY <- cbind(dfY, rep(6, nrow(dfY)))
colnames(dfY)[6] <- "Color"

dfPlot <- rbind(dfAnons, dfXnons, dfYnons, dfA, dfX, dfY)
dfPlot$Color <- as.factor(dfPlot$Color)

# Constructing the plot object. The colors will not work the way they should if
# one of the groups doesn't have any genes
p <- ggplot(data = dfPlot, aes(x = logFC, y = -log10(adj.P.Val), color=Color )) +
  geom_point(alpha = 0.5, size = 8) +
  theme_bw() +
  theme(legend.position = "none") +
  xlim(c(-15, 15)) + ylim(c(0, 50)) +
  scale_color_manual(values = c("azure3", "pink", "seagreen2", "black", "mediumvioletred", "springgreen")) +
  labs(x=expression(log[2](FC)),
       y=expression(-log[10] ~ "(FDR-adjusted " ~ italic("p") ~ "-value)")) +
  theme(axis.title.x=element_text(size=12), 
        axis.text.x=element_text(size=10)) +
  theme(axis.title.y=element_text(size=12),
        axis.text.y=element_text(size=10))
p

forLabel <- subset(dfPlot, adj.P.Val<=0.01 & abs(logFC)>=2)

# Adding lines for significance thresholds
pdf("~/3.0 Hasting Research/Liver Cancer Analysis/DEG Figures/Volcano_plot_male_HBV_DEGs.pdf", width=12, height=12)
p + geom_hline(yintercept = 2, colour="#000000", linetype="dashed"
) + geom_vline(xintercept = 2, colour="#000000", linetype="dashed"
) + geom_vline(xintercept = -2, colour="#000000", linetype="dashed"
) + geom_text_repel(data=forLabel, max.iter=1000, box.padding = 0.25, force=1, aes(x = logFC, y = -log10(adj.P.Val), label=name, color=Color), size=8)
dev.off()


# ===========================================
#Volcano plot of female HBV tumor vs tumor-adjacent
# ===========================================

df <- data.frame(vTopTable_F_HBV$adj.P.Val, vTopTable_F_HBV$logFC, vTopTable_F_HBV$chr, vTopTable_F_HBV$GENEID, vTopTable_F_HBV$gene_name)
colnames(df) <- c("adj.P.Val", "logFC", "chr", "id", "name")
dfSig <- df[(abs(df$logFC) >= 2 & df$adj.P.Val <= 0.01),]$id


dfAnons <- subset(df, chr != "chrX" & chr != "chrY" & !(id %in% dfSig))
dfAnons <- cbind(dfAnons , rep(1, nrow(dfAnons)))
colnames(dfAnons)[6] <- "Color"

dfXnons <- subset(df, chr == "chrX" & !(id %in% dfSig))
dfXnons <- cbind(dfXnons, rep(2, nrow(dfXnons)))
colnames(dfXnons)[6] <- "Color"

dfYnons <- subset(df, chr == "chrY" & !(id %in% dfSig))
dfYnons <- cbind(dfYnons, rep(3, nrow(dfYnons)))
colnames(dfYnons)[6] <- "Color"

dfA <- subset(df, chr != "chrX" & chr != "chrY" & id %in% dfSig)
dfA <- cbind(dfA, rep(4, nrow(dfA)))
colnames(dfA)[6] <- "Color"

dfX <- subset(df, chr == "chrX" & id %in% dfSig)
dfX <- cbind(dfX, rep(5, nrow(dfX)))
colnames(dfX)[6] <- "Color"

dfY <- subset(df, chr == "chrY" & id %in% dfSig)
dfY <- cbind(dfY, rep(6, nrow(dfY)))
colnames(dfY)[6] <- "Color"

dfPlot <- rbind(dfAnons, dfXnons, dfYnons, dfA, dfX, dfY)
dfPlot$Color <- as.factor(dfPlot$Color)

# Constructing the plot object. The colors will not work the way they should if
# one of the groups doesn't have any genes
p <- ggplot(data = dfPlot, aes(x = logFC, y = -log10(adj.P.Val), color=Color )) +
  geom_point(alpha = 0.5, size = 8) +
  theme_bw() +
  theme(legend.position = "none") +
  xlim(c(-15, 15)) + ylim(c(0, 30)) +
  scale_color_manual(values = c("azure3", "pink", "seagreen2", "black", "mediumvioletred", "springgreen")) +
  labs(x=expression(log[2](FC)),
       y=expression(-log[10] ~ "(FDR-adjusted " ~ italic("p") ~ "-value)")) +
  theme(axis.title.x=element_text(size=12), 
        axis.text.x=element_text(size=10)) +
  theme(axis.title.y=element_text(size=12),
        axis.text.y=element_text(size=10))
p

forLabel <- subset(dfPlot, adj.P.Val<=0.01 & abs(logFC)>=2)

# Adding lines for significance thresholds
pdf("~/3.0 Hasting Research/Liver Cancer Analysis/DEG Figures/Volcano_plot_female_HBV_DEGs.pdf", width=12, height=12)
p + geom_hline(yintercept = 2, colour="#000000", linetype="dashed"
) + geom_vline(xintercept = 2, colour="#000000", linetype="dashed"
) + geom_vline(xintercept = -2, colour="#000000", linetype="dashed"
) + geom_text_repel(data=forLabel, max.iter=1000, box.padding = 0.25, force=1, aes(x = logFC, y = -log10(adj.P.Val), label=name, color=Color), size=8)
dev.off()

library(VennDiagram)
venn.diagram(List("Female"=DEGs_F_HBV$gene_name, "Male"=DEGs_M_HBV$gene_name),  filename = "~/3.0 Hasting Research/Liver Cancer Analysis/DEG Figures/Overlap_DEGs_HBV.png")
venn.diagram(List("Female"=DEGs_F_HBV_relax_p$gene_name, "Male"=DEGs_M_HBV$gene_name),  filename = "~/3.0 Hasting Research/Liver Cancer Analysis/DEG Figures/Overlap_DEGs_HBV_Frelaxp.png")


# ===============================================
#Volcano plot of male HCV tumor vs tumor-adjacent
# ===============================================

df <- data.frame(vTopTable_M_HCV$adj.P.Val, vTopTable_M_HCV$logFC, vTopTable_M_HCV$chr, vTopTable_M_HCV$GENEID, vTopTable_M_HCV$gene_name)
colnames(df) <- c("adj.P.Val", "logFC", "chr", "id", "name")
dfSig <- df[(abs(df$logFC) >= 2 & df$adj.P.Val <= 0.01),]$id


dfAnons <- subset(df, chr != "chrX" & chr != "chrY" & !(id %in% dfSig))
dfAnons <- cbind(dfAnons , rep(1, nrow(dfAnons)))
colnames(dfAnons)[6] <- "Color"

dfXnons <- subset(df, chr == "chrX" & !(id %in% dfSig))
dfXnons <- cbind(dfXnons, rep(2, nrow(dfXnons)))
colnames(dfXnons)[6] <- "Color"

dfYnons <- subset(df, chr == "chrY" & !(id %in% dfSig))
dfYnons <- cbind(dfYnons, rep(3, nrow(dfYnons)))
colnames(dfYnons)[6] <- "Color"

dfA <- subset(df, chr != "chrX" & chr != "chrY" & id %in% dfSig)
dfA <- cbind(dfA, rep(4, nrow(dfA)))
colnames(dfA)[6] <- "Color"

dfX <- subset(df, chr == "chrX" & id %in% dfSig)
dfX <- cbind(dfX, rep(5, nrow(dfX)))
colnames(dfX)[6] <- "Color"

dfY <- subset(df, chr == "chrY" & id %in% dfSig)
dfY <- cbind(dfY, rep(6, nrow(dfY)))
colnames(dfY)[6] <- "Color"

dfPlot <- rbind(dfAnons, dfXnons, dfYnons, dfA, dfX, dfY)
dfPlot$Color <- as.factor(dfPlot$Color)

# Constructing the plot object. The colors will not work the way they should if
# one of the groups doesn't have any genes
p <- ggplot(data = dfPlot, aes(x = logFC, y = -log10(adj.P.Val), color=Color )) +
  geom_point(alpha = 0.5, size = 8) +
  theme_bw() +
  theme(legend.position = "none") +
  xlim(c(-15, 15)) + ylim(c(0, 70)) +
  scale_color_manual(values = c("azure3", "pink", "seagreen2", "black", "mediumvioletred", "springgreen")) +
  labs(x=expression(log[2](FC)),
       y=expression(-log[10] ~ "(FDR-adjusted " ~ italic("p") ~ "-value)")) +
  theme(axis.title.x=element_text(size=12), 
        axis.text.x=element_text(size=10)) +
  theme(axis.title.y=element_text(size=12),
        axis.text.y=element_text(size=10))
p

forLabel <- subset(dfPlot, adj.P.Val<=0.01 & abs(logFC)>=2)

# Adding lines for significance thresholds
pdf("~/3.0 Hasting Research/Liver Cancer Analysis/DEG Figures/Volcano_plot_male_HCV_DEGs.pdf", width=12, height=12)
p + geom_hline(yintercept = 2, colour="#000000", linetype="dashed"
) + geom_vline(xintercept = 2, colour="#000000", linetype="dashed"
) + geom_vline(xintercept = -2, colour="#000000", linetype="dashed"
) + geom_text_repel(data=forLabel, max.iter=1000, box.padding = 0.25, force=1, aes(x = logFC, y = -log10(adj.P.Val), label=name, color=Color), size=8)
dev.off()


# ===========================================
#Volcano plot of female HCV tumor vs tumor-adjacent
# ===========================================

df <- data.frame(vTopTable_F_HCV$adj.P.Val, vTopTable_F_HCV$logFC, vTopTable_F_HCV$chr, vTopTable_F_HCV$GENEID, vTopTable_F_HCV$gene_name)
colnames(df) <- c("adj.P.Val", "logFC", "chr", "id", "name")
dfSig <- df[(abs(df$logFC) >= 2 & df$adj.P.Val <= 0.01),]$id


dfAnons <- subset(df, chr != "chrX" & chr != "chrY" & !(id %in% dfSig))
dfAnons <- cbind(dfAnons , rep(1, nrow(dfAnons)))
colnames(dfAnons)[6] <- "Color"

dfXnons <- subset(df, chr == "chrX" & !(id %in% dfSig))
dfXnons <- cbind(dfXnons, rep(2, nrow(dfXnons)))
colnames(dfXnons)[6] <- "Color"

dfYnons <- subset(df, chr == "chrY" & !(id %in% dfSig))
dfYnons <- cbind(dfYnons, rep(3, nrow(dfYnons)))
colnames(dfYnons)[6] <- "Color"

dfA <- subset(df, chr != "chrX" & chr != "chrY" & id %in% dfSig)
dfA <- cbind(dfA, rep(4, nrow(dfA)))
colnames(dfA)[6] <- "Color"

dfX <- subset(df, chr == "chrX" & id %in% dfSig)
dfX <- cbind(dfX, rep(5, nrow(dfX)))
colnames(dfX)[6] <- "Color"

dfY <- subset(df, chr == "chrY" & id %in% dfSig)
dfY <- cbind(dfY, rep(6, nrow(dfY)))
colnames(dfY)[6] <- "Color"

dfPlot <- rbind(dfAnons, dfXnons, dfYnons, dfA, dfX, dfY)
dfPlot$Color <- as.factor(dfPlot$Color)

# Constructing the plot object. The colors will not work the way they should if
# one of the groups doesn't have any genes
p <- ggplot(data = dfPlot, aes(x = logFC, y = -log10(adj.P.Val), color=Color )) +
  geom_point(alpha = 0.5, size = 8) +
  theme_bw() +
  theme(legend.position = "none") +
  xlim(c(-15, 15)) + ylim(c(0, 40)) +
  scale_color_manual(values = c("azure3", "pink", "seagreen2", "black", "mediumvioletred", "springgreen")) +
  labs(x=expression(log[2](FC)),
       y=expression(-log[10] ~ "(FDR-adjusted " ~ italic("p") ~ "-value)")) +
  theme(axis.title.x=element_text(size=12), 
        axis.text.x=element_text(size=10)) +
  theme(axis.title.y=element_text(size=12),
        axis.text.y=element_text(size=10))
p

forLabel <- subset(dfPlot, adj.P.Val<=0.01 & abs(logFC)>=2)

# Adding lines for significance thresholds
pdf("~/3.0 Hasting Research/Liver Cancer Analysis/DEG Figures/Volcano_plot_female_HCV_DEGs.pdf", width=12, height=12)
p + geom_hline(yintercept = 2, colour="#000000", linetype="dashed"
) + geom_vline(xintercept = 2, colour="#000000", linetype="dashed"
) + geom_vline(xintercept = -2, colour="#000000", linetype="dashed"
) + geom_text_repel(data=forLabel, max.iter=1000, box.padding = 0.25, force=1, aes(x = logFC, y = -log10(adj.P.Val), label=name, color=Color), size=8)
dev.off()

library(VennDiagram)
venn.diagram(List("Female"=DEGs_F_HCV$gene_name, "Male"=DEGs_M_HCV$gene_name),  filename = "~/3.0 Hasting Research/Liver Cancer Analysis/DEG Figures/Overlap_DEGs_HCV.png")
venn.diagram(List("Female"=DEGs_F_HCV_relax_p$gene_name, "Male"=DEGs_M_HCV$gene_name),  filename = "~/3.0 Hasting Research/Liver Cancer Analysis/DEG Figures/Overlap_DEGs_HCV_Frelaxp.png")

# ===============================================
#Volcano plot of male non-viral tumor vs tumor-adjacent
# ===============================================

df <- data.frame(vTopTable_M_Neither$adj.P.Val, vTopTable_M_Neither$logFC, vTopTable_M_Neither$chr, vTopTable_M_Neither$GENEID, vTopTable_M_Neither$gene_name)
colnames(df) <- c("adj.P.Val", "logFC", "chr", "id", "name")
dfSig <- df[(abs(df$logFC) >= 2 & df$adj.P.Val <= 0.01),]$id


dfAnons <- subset(df, chr != "chrX" & chr != "chrY" & !(id %in% dfSig))
dfAnons <- cbind(dfAnons , rep(1, nrow(dfAnons)))
colnames(dfAnons)[6] <- "Color"

dfXnons <- subset(df, chr == "chrX" & !(id %in% dfSig))
dfXnons <- cbind(dfXnons, rep(2, nrow(dfXnons)))
colnames(dfXnons)[6] <- "Color"

dfYnons <- subset(df, chr == "chrY" & !(id %in% dfSig))
dfYnons <- cbind(dfYnons, rep(3, nrow(dfYnons)))
colnames(dfYnons)[6] <- "Color"

dfA <- subset(df, chr != "chrX" & chr != "chrY" & id %in% dfSig)
dfA <- cbind(dfA, rep(4, nrow(dfA)))
colnames(dfA)[6] <- "Color"

dfX <- subset(df, chr == "chrX" & id %in% dfSig)
dfX <- cbind(dfX, rep(5, nrow(dfX)))
colnames(dfX)[6] <- "Color"

dfY <- subset(df, chr == "chrY" & id %in% dfSig)
dfY <- cbind(dfY, rep(6, nrow(dfY)))
colnames(dfY)[6] <- "Color"

dfPlot <- rbind(dfAnons, dfXnons, dfYnons, dfA, dfX, dfY)
dfPlot$Color <- as.factor(dfPlot$Color)

# Constructing the plot object. The colors will not work the way they should if
# one of the groups doesn't have any genes
p <- ggplot(data = dfPlot, aes(x = logFC, y = -log10(adj.P.Val), color=Color )) +
  geom_point(alpha = 0.5, size = 8) +
  theme_bw() +
  theme(legend.position = "none") +
  xlim(c(-15, 15)) + ylim(c(0, 40)) +
  scale_color_manual(values = c("azure3", "pink", "seagreen2", "black", "mediumvioletred", "springgreen")) +
  labs(x=expression(log[2](FC)),
       y=expression(-log[10] ~ "(FDR-adjusted " ~ italic("p") ~ "-value)")) +
  theme(axis.title.x=element_text(size=12), 
        axis.text.x=element_text(size=10)) +
  theme(axis.title.y=element_text(size=12),
        axis.text.y=element_text(size=10))
p

forLabel <- subset(dfPlot, adj.P.Val<=0.01 & abs(logFC)>=2)

# Adding lines for significance thresholds
pdf("~/3.0 Hasting Research/Liver Cancer Analysis/DEG Figures/Volcano_plot_male_Neither_DEGs.pdf", width=12, height=12)
p + geom_hline(yintercept = 2, colour="#000000", linetype="dashed"
) + geom_vline(xintercept = 2, colour="#000000", linetype="dashed"
) + geom_vline(xintercept = -2, colour="#000000", linetype="dashed"
) + geom_text_repel(data=forLabel, max.iter=1000, box.padding = 0.25, force=1, aes(x = logFC, y = -log10(adj.P.Val), label=name, color=Color), size=8)
dev.off()


# ===========================================
#Volcano plot of female non-viral tumor vs tumor-adjacent
# ===========================================

df <- data.frame(vTopTable_F_Neither$adj.P.Val, vTopTable_F_Neither$logFC, vTopTable_F_Neither$chr, vTopTable_F_Neither$GENEID, vTopTable_F_Neither$gene_name)
colnames(df) <- c("adj.P.Val", "logFC", "chr", "id", "name")
dfSig <- df[(abs(df$logFC) >= 2 & df$adj.P.Val <= 0.01),]$id


dfAnons <- subset(df, chr != "chrX" & chr != "chrY" & !(id %in% dfSig))
dfAnons <- cbind(dfAnons , rep(1, nrow(dfAnons)))
colnames(dfAnons)[6] <- "Color"

dfXnons <- subset(df, chr == "chrX" & !(id %in% dfSig))
dfXnons <- cbind(dfXnons, rep(2, nrow(dfXnons)))
colnames(dfXnons)[6] <- "Color"

dfYnons <- subset(df, chr == "chrY" & !(id %in% dfSig))
dfYnons <- cbind(dfYnons, rep(3, nrow(dfYnons)))
colnames(dfYnons)[6] <- "Color"

dfA <- subset(df, chr != "chrX" & chr != "chrY" & id %in% dfSig)
dfA <- cbind(dfA, rep(4, nrow(dfA)))
colnames(dfA)[6] <- "Color"

dfX <- subset(df, chr == "chrX" & id %in% dfSig)
dfX <- cbind(dfX, rep(5, nrow(dfX)))
colnames(dfX)[6] <- "Color"

dfY <- subset(df, chr == "chrY" & id %in% dfSig)
dfY <- cbind(dfY, rep(6, nrow(dfY)))
colnames(dfY)[6] <- "Color"

dfPlot <- rbind(dfAnons, dfXnons, dfYnons, dfA, dfX, dfY)
dfPlot$Color <- as.factor(dfPlot$Color)

# Constructing the plot object. The colors will not work the way they should if
# one of the groups doesn't have any genes
p <- ggplot(data = dfPlot, aes(x = logFC, y = -log10(adj.P.Val), color=Color )) +
  geom_point(alpha = 0.5, size = 8) +
  theme_bw() +
  theme(legend.position = "none") +
  xlim(c(-15, 15)) + ylim(c(0, 7.5)) +
  scale_color_manual(values = c("azure3", "pink", "seagreen2", "black", "mediumvioletred", "springgreen")) +
  labs(x=expression(log[2](FC)),
       y=expression(-log[10] ~ "(FDR-adjusted " ~ italic("p") ~ "-value)")) +
  theme(axis.title.x=element_text(size=12), 
        axis.text.x=element_text(size=10)) +
  theme(axis.title.y=element_text(size=12),
        axis.text.y=element_text(size=10))
p

forLabel <- subset(dfPlot, adj.P.Val<=0.01 & abs(logFC)>=2)

# Adding lines for significance thresholds
pdf("~/3.0 Hasting Research/Liver Cancer Analysis/DEG Figures/Volcano_plot_female_Neither_DEGs.pdf", width=12, height=12)
p + geom_hline(yintercept = 2, colour="#000000", linetype="dashed"
) + geom_vline(xintercept = 2, colour="#000000", linetype="dashed"
) + geom_vline(xintercept = -2, colour="#000000", linetype="dashed"
) + geom_text_repel(data=forLabel, max.iter=1000, box.padding = 0.25, force=1, aes(x = logFC, y = -log10(adj.P.Val), label=name, color=Color), size=8)
dev.off()

library(VennDiagram)
venn.diagram(List("Female"=DEGs_F_Neither$gene_name, "Male"=DEGs_M_Neither$gene_name),  filename = "~/3.0 Hasting Research/Liver Cancer Analysis/DEG Figures/Overlap_DEGs_Neither.png")
venn.diagram(List("Female"=DEGs_F_Neither_relax_p$gene_name, "Male"=DEGs_M_Neither$gene_name),  filename = "~/3.0 Hasting Research/Liver Cancer Analysis/DEG Figures/Overlap_DEGs_Neither_Frelaxp.png")
