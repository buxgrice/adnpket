## ---------------------------
##
## Script name: Pre-Processesing and Differential Expression Analysis for ADNP Ketamine Stuy 
##
## Purpose of script: Starting from a raw count matrix, we will sort the columns by trial timepoint, 
## remove those that failed QC, then filter and voom nomalize the resulting matrix. 
## We will then run Differential Expression Analysis to compare expression profiles across the trail 
## timecourse in response to ketamine treatment
##
## Author: Ariela Buxbaum Grice
##
## Date Created: 2023-06-14
##
## 
## Email: buxgrice@gmail.com
##
##
## ---------------------------

#||----------------------------------------------------------------------------||-----------
#||----------------------------------------------------------------------------||-----------
#||                                                                            ||
#||                  ||LOAD PACKAGES & INPUTS, SET DIRECTORIES||               ||
#||                                                                            ||
#||----------------------------------------------------------------------------||-----------
#||----------------------------------------------------------------------------||-----------

library(easypackages)

#processing 
libraries("edgeR","stringr","biomaRt")
# normalization and DEG analysis: 
libraries("doParallel", "variancePartition", "ape", "limma")

#||----------------------------------------------------------------------------||-----------
#||                       ||SET BASIC WORKING DIRECTORY||                      ||
#||----------------------------------------------------------------------------||-----------

setwd ("/Users/arielabuxbaumgrice/Desktop/projects_new/ADNP_ketamine")
getwd()

#||----------------------------------------------------------------------------||-----------
#||     SET UP PARAMETERS AND CALLS TO STREAMLINE GRAPHING, TABLES, ANALYSES   ||
#||----------------------------------------------------------------------------||-----------

input.dir = "/Users/arielabuxbaumgrice/Desktop/projects_new/ADNP_ketamine/input.files/"
output.dir = "/Users/arielabuxbaumgrice/Desktop/projects_new/ADNP_ketamine/outputs/"

options(stringsAsFactors = FALSE)

#||----------------------------------------------------------------------------||-----------
#||----------------------------------------------------------------------------||-----------
#||                                ||LOAD DATA||                               ||
#||----------------------------------------------------------------------------||-----------
#||----------------------------------------------------------------------------||-----------

#||-------------------------------------------------------||-----
#|| Subject Metadata                                      ||
#||-------------------------------------------------------||------

meta = read.csv("Final_MetadataCombined.csv", row.names = 1) # a complete file with Subject metadata, estimated cell fractions, and QC metrics

#||-------------------------------------------------------||-----
#|| Raw Counts                                            ||
#||-------------------------------------------------------||-----

raw = read.delim(paste0(input.dir,"ADNP_counts.txt"))

#||----------------------------------------------------------------------------||-----------
#||----------------------------------------------------------------------------||-----------
#||                          ||CLEAN AND FORMAT SAMPLES||                      ||
#||----------------------------------------------------------------------------||-----------
#||----------------------------------------------------------------------------||-----------

# CC1329.202.Post is a OUTLIER; 
# CC1668.201.Post and CC1628.201.Wk4 are swapped and identity could not be determined

counts = raw[,!(names(raw) %in% c("CC1329.202.Post", "CC1668.201.Post", "CC1628.201.Wk4"))]
dim(counts) #[1] 58929    52

#Re-order count matrix to the trial timecourse
genes = counts[1]
pre = counts[grepl("Pre", colnames(counts))]
post = counts[grepl("Post", colnames(counts))]
D1 = counts[grepl("D1", colnames(counts))]
Wk1 = counts[grepl("Wk1", colnames(counts))]
Wk2 = counts[grepl("Wk2", colnames(counts))]
Wk4 = counts[grepl("Wk4", colnames(counts))]
counts = data.frame(genes, pre, post, D1, Wk1, Wk2, Wk4)

table(counts[1] == raw[1]) # TRUE
rownames(counts) = counts$Gene
counts = counts[,-1]

#Save generated ordered raw count matrix of the final sample group
write.table(counts, "Final_CountMatrix_ADNP.txt", sep = "\t") 

#||----------------------------------------------------------------------------||-----------
#||----------------------------------------------------------------------------||-----------
#||                         ||FILTER AND VOOM NORMALIZE||                      ||
#||----------------------------------------------------------------------------||-----------
#||----------------------------------------------------------------------------||-----------

filter=ncol(counts)/3 #4
DGE = DGEList(counts=counts[,1:51], genes = rownames(counts))

dim(DGE) #[1] 58929    51
keep = rowSums(cpm(DGE)>1) >= filter
DGE = DGE[keep,] 
dim(DGE) #[1] 17218    51

DGE = calcNormFactors(DGE)
DGE$samples

#VOOM NORMALIZATION
VST.mat = voom(DGE, design = NULL, plot=TRUE)
voom.mat = cbind(VST.mat$genes, VST.mat$E)
rownames(voom.mat)=voom.mat$genes
voom.mat = voom.mat[,-1]

#Save a voom-normalized filtered expression matrix
write.table(voom.mat, "VSTmat_ADNP.txt", sep="\t") 

#||----------------------------------------------------------------------------||-----------
#||----------------------------------------------------------------------------||-----------
#||                                                                            ||
#||                                ||DEG ANALYSIS||                            ||
#||                                                                            ||
#||----------------------------------------------------------------------------||-----------
#||----------------------------------------------------------------------------||-----------

#||-------------------------------------------------------||-----
#|| Load Voom Matrix                                      ||
#||-------------------------------------------------------||-----

#voom.mat = read.delim("VSTmat_ADNP.txt", check.names=TRUE, 
#                      stringsAsFactors=FALSE, row.names=1, header=TRUE) 

#||-------------------------------------------------------||-----
#|| Isolate Metadata for Covariates only                  ||
#||-------------------------------------------------------||------

covariates = meta[,c("CCID", "Timepoint", "Sex", "MutClass", "RIN")]
table(rownames(covariates)==colnames(voom.mat))

#||-------------------------------------------------------||-----
#|| MDs and Clustering                                    ||
#||-------------------------------------------------------||-----

time = factor(covariates$Timepoint)
time.cols = c("darkgoldenrod", "plum4", "darkcyan", 
              "cornsilk4", "mediumvioletred", "midnightblue")

subjects = colnames(voom.mat)
subjects = str_replace_all(subjects,"\\..*","")
plot.mat = voom.mat # make a copy with the abbreviated sample IDs
colnames(plot.mat) = subjects

CorMat = 1 - cor(plot.mat)
distMat = as.dist(CorMat)
MDS = cmdscale( distMat, eig=TRUE )
clust = stats::hclust( distMat, "ward.D" )

par(mfrow=c(1,1))
plotMDS(plot.mat, col=time.cols[time], cex.axis=0.9, las=1, cex=0.8)
legend(x = "bottom",inset = c(0, .5),legend = c("Baseline", "Post",
                                                "Day1", "Week1", 
                                                "Week2", "Week4"),  
       col = c("darkgoldenrod", "plum4", "darkcyan", 
               "cornsilk4", "mediumvioletred", "midnightblue"), lwd = 2,bty = "n")

par(mfrow=c(1,1))
plot(as.phylo(clust), type="phylogram", tip.color = time.cols[time],
     cex=0.7, lwd=2, main="Clustering by Donor and Timepoint")

#||----------------------------------------------------------------------------||-----------
#||----------------------------------------------------------------------------||-----------
#||                          ||VARIANCE PARTITION||                            ||
#||----------------------------------------------------------------------------||-----------
#||----------------------------------------------------------------------------||-----------

clust = makeCluster(4)
registerDoParallel(clust)

Time = factor(covariates$Timepoint)
Subj = factor(covariates$CCID)
Sex = factor(covariates$Sex)
Class = factor(covariates$MutClass, levels = c("Class I", "Class II"), 
               labels = c("A", "B"))
RIN = covariates$RIN

form = ~ RIN + (1|Sex) + (1|Subj) + (1|Class) + (1|Time) 

varPart = fitExtractVarPartModel(voom.mat, form, covariates)
vp = sortCols( varPart )
plotVarPart(vp,label.angle=50)

#||-------------------------------------------------------||----
#|| Design Matrix                                         ||
#||-------------------------------------------------------||----

design = model.matrix(~0 + Time + Sex)

colnames(design)
covariates$CCID = factor(covariates$CCID)
IDdupcorr = duplicateCorrelation(voom.mat, design, block = covariates$CCID)
IDdupcorr$consensus.correlation # 0.3843199
fit1 = lmFit(voom.mat, design, block = covariates$CCID,
             correlation=IDdupcorr$consensus)

contrast.matrix = makeContrasts(
  T1 =TimeDosingDay0-TimeBaseline,
  T2 =TimeDay1_post-TimeBaseline,
  T3 =TimeWeek1_post-TimeBaseline,
  T4 =TimeWeek2_post-TimeBaseline,
  T5 =TimeWeek4_post-TimeBaseline,
  levels=design)

#apply linear model
fit2 = contrasts.fit(fit1, contrast.matrix)
fit2 = eBayes(fit2)

#||-------------------------------------------------------||-----
#|| Save DEGs                                             ||
#||-------------------------------------------------------||-----

# Pre vs Post
topTable(fit2, coef="T1", adjust="BH")
DEGs_T1 = topTable(fit2, coef="T1", n=nrow(voom.mat))

# If you want to save your full list of differentially expressed genes 
# > write.table(DEGs_T1, "DEGs_Time_PREvPOST.txt", sep="\t")
# If you want to make a quick volcano plot to visualize
# > plot(DEGs_T1$logFC, -log10(DEGs_T1$adj.P.Val))
# > abline(h=-log10(0.05))

# Pre vs Day 1
topTable(fit2, coef="T2", adjust="BH")
DEGs_T2 = topTable(fit2, coef="T2", n=nrow(voom.mat))

# Pre vs Week 1
topTable(fit2, coef="T3", adjust="BH")
DEGs_T3 = topTable(fit2, coef="T3", n=nrow(voom.mat))

# Pre vs Week 2
topTable(fit2, coef="T4", adjust="BH")
DEGs_T4 = topTable(fit2, coef="T4", n=nrow(voom.mat))

# Pre vs Week 4
topTable(fit2, coef="T5", adjust="BH")
DEGs_T5 = topTable(fit2, coef="T5", n=nrow(voom.mat))

# Isolate only those with pval < 0.05
for (i in c(1:5)) {
  deg_list = get((paste0("DEGs_T", i)))
  keep = deg_list$adj.P.Val < 0.05
  if (any(keep)   == TRUE) {
    sig = deg_list[keep,]
    list = data.frame(genes = rownames(sig), comparison = paste0("Time",i))
    assign((paste0("sigGenes", i)), list)}}

allDEGs = rbind(sigGenes1, sigGenes2, sigGenes4, sigGenes5) # No genes from T3 (Baseline vs. One Week) were signficant

# If you want to save your full list of differentially expressed genes by timepoint
# > write.table(allDEGs, "EnsIDs_allDEGS_TimeComparisons.txt", sep="\t")


