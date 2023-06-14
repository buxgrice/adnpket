## ---------------------------
##
## Script name: WCGNA Analysis for ADNP Ketamine Study 
##
## Purpose of script: Starting from a normalized count matrix, and a metadata file,
## we will format the selected covariates for WGCNA. We will select the soft threshold 
## to make the gene dendogram, and then generate modules using a few methods. 
## We use the dynamic cut modules and then associate trial timepoint with each module to test for signficance
## Significantly dynamically regulated modules associated with time are plotted and their genes 
## extreacted for downstream enrichment analyses.
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
libraries("WGCNA","cluster")

options(stringsAsFactors  =  FALSE)
allowWGCNAThreads(n=10)

#||----------------------------------------------------------------------------||-----------
#||                       ||SET BASIC WORKING DIRECTORY||                      ||
#||----------------------------------------------------------------------------||-----------

setwd ("/Users/arielabuxbaumgrice/Desktop/projects_new/ADNP_ketamine")
getwd()

#||----------------------------------------------------------------------------||-----------
#||----------------------------------------------------------------------------||-----------
#||                          ||CLEAN AND FORMAT METADATA||                     ||
#||                                 ||FOR WGCNA||                              ||
#||----------------------------------------------------------------------------||-----------
#||----------------------------------------------------------------------------||-----------

#||-------------------------------------------------------||-----
#|| Expression Data                                       ||
#||-------------------------------------------------------||-----

# The voom-normalized filtered expression matrix generated in "Processing_DifferentialExpression.R"
voom.mat = read.delim("VSTmat_ADNP.txt", check.names=TRUE, stringsAsFactors=FALSE, row.names=1, header=TRUE) 

expr = as.data.frame(t(voom.mat))
dim(expr) #[1]  51 17218

gsg=goodSamplesGenes(expr,verbose=3)
gsg$allOK

nGenes = ncol(expr)
nSamples = nrow(expr)

#||-------------------------------------------------------||-----
#|| Subject Metadata                                      ||
#||-------------------------------------------------------||------

# a complete file with Subject metadata, estimated cell fractions, and QC metrics
meta = read.csv("Final_MetadataCombined.csv", row.names = 1) 
meta_wgcna = data.frame(matrix(nrow=ncol(voom.mat),ncol=4))
rownames(meta_wgcna) = colnames(voom.mat)
colnames(meta_wgcna) = c("Donor", "Time", "Sex", "MutClass")

# > Format metadata for WGNCA 

meta_wgcna$Donor  = str_sub(rownames(meta),  1, 10)
meta_wgcna$Time = str_sub(rownames(meta), 12)
meta_wgcna$Time = factor((meta_wgcna$Time), levels = c("Pre", "Post", "D1", "Wk1", "Wk2", "Wk4"),
                   labels = c("0", "1", "2", "3", "4", "5"))

males = intersect(rownames(meta[(meta$Sex == "M"),]), (rownames(meta_wgcna)))
females = intersect(rownames(meta[(meta$Sex == "F"),]), (rownames(meta_wgcna)))
any(males %in% females) # should be false
meta_wgcna[males,]$Sex = "0"
meta_wgcna[females,]$Sex = "1"

class1 = intersect(rownames(meta[(meta$MutClass == "Class I"),]), (rownames(meta_wgcna)))
class2 = intersect(rownames(meta[(meta$MutClass == "Class II"),]), (rownames(meta_wgcna)))
any(class1 %in% class2) # should be false
meta_wgcna[class1,]$MutClass = "0"
meta_wgcna[class2,]$MutClass = "1"

table(rownames(meta_wgcna)==rownames(expr))

#Save a metadata file designed specifically for WGCNA analysis  
write.table(meta_wgcna, "ADNP_WGCNA_metadata.txt", sep = "\t")

#||----------------------------------------------------------------------------||-----------
#||----------------------------------------------------------------------------||-----------
#||                                                                            ||
#||                                 ||WGCNA||                                  ||
#||                                                                            ||
#||----------------------------------------------------------------------------||-----------
#||----------------------------------------------------------------------------||-----------

#||----------------------------------------------------------------------------||-----------
#||----------------------------------------------------------------------------||-----------
#||                             ||THRESHOLDING||                               ||
#||----------------------------------------------------------------------------||-----------
#||----------------------------------------------------------------------------||-----------

powers = c(1:30)
thresh = pickSoftThreshold(expr, powerVector=powers, networkType="signed")

par(mfrow=c(2,2))
plot(thresh$fitIndices[,1],-sign(thresh$fitIndices[,3])*thresh$fitIndices[,2], xlab="Soft Threshold (power)",ylab="SFT, signed R^2",type="n",main=paste("Scale independence - NGN2"))
text(thresh$fitIndices[,1],-sign(thresh$fitIndices[,3])*thresh$fitIndices[,2],labels=powers,col="steelblue4")
abline(h=0.80,col="red")  
plot(thresh$fitIndices[,1],thresh$fitIndices[,5],type="n",
     xlab="Soft Threshold (power)",ylab="Mean Connectivity",main=paste("Mean connectivity"))
text(thresh$fitIndices[,1],thresh$fitIndices[,5],labels=powers,col="steelblue4")

#||----------------------------------------------------------------------------||-----------
#||----------------------------------------------------------------------------||-----------
#||                               ||GENE TREE||                                ||
#||----------------------------------------------------------------------------||-----------
#||----------------------------------------------------------------------------||-----------

allowWGCNAThreads(n=10)

# > select your power based on the thresholding tree; we used 11

adjacencyPre = adjacency((expr),power=11,type="signed") 
diag(adjacencyPre)=0
dissTOMPre   = 1-TOMsimilarity(adjacencyPre, TOMType="signed")
geneTreePre  = hclust(as.dist(dissTOMPre), method="average")

#||-------------------------------------------------------||-----
#|| HYBRID CUT                                            ||
#||-------------------------------------------------------||-----

mColorh=NULL
for (ds in 0:3){
  tree = cutreeHybrid(dendro = geneTreePre, pamStage=FALSE,
                      minClusterSize = (50), cutHeight = 0.99999999999, 
                      deepSplit = ds, distM = dissTOMPre)
  mColorh=cbind(mColorh,labels2colors(tree$labels));
}

par(mfrow=c(1,1))
plotDendroAndColors(geneTreePre, mColorh, paste("dpSplt =",0:4), main = "Hybrid tree cut",dendroLabels=FALSE);

#||-------------------------------------------------------||-----
#|| DYNAMIC CUT                                           ||
#||-------------------------------------------------------||-----

DetectedColors0 = NULL;
DetectedColors0 = 
  cbind(DetectedColors0,labels2colors(cutreeDynamic(dendro = geneTreePre,
                                                    cutHeight = 0.99999999999, minClusterSize = 50,
                                                    method = "tree", deepSplit = 0)));
DetectedColors1 = NULL;
DetectedColors1 = 
  cbind(DetectedColors1,labels2colors(cutreeDynamic(dendro = geneTreePre,
                                                    cutHeight = 0.99999999999, minClusterSize = 50,
                                                    method = "tree", deepSplit = 1)));
DetectedColors =cbind(DetectedColors0,DetectedColors1);

plotDendroAndColors(geneTreePre, mColorh, paste("dpSplt =",0:4), main = "Hybrid tree cut",dendroLabels=FALSE);
plotDendroAndColors(geneTreePre, DetectedColors, paste("dpSplt =",0:1), main = "Dynamic tree cut",dendroLabels=FALSE);\

#||-------------------------------------------------------||-----
#|| Deep Split                                            ||
#||-------------------------------------------------------||-----

modulesPRE =  DetectedColors[,1]

#||-------------------------------------------------------||-----
#|| Look at Original Modules                              ||
#||-------------------------------------------------------||-----

table(modulesPRE)

MEs0  =  moduleEigengenes(expr,  modulesPRE)$eigengenes
MEs = orderMEs(MEs0)

# > Save original module eigenvalues, associated with gene expression
write.table(MEs, "Orig.DynamicCut_WGCNA_ExprMEs.txt", sep = "\t")

#||-------------------------------------------------------||-----
#|| Module Visualization                                  ||
#||-------------------------------------------------------||-----
#|
distMEs = 1-abs(cor(MEs0,use="p"))
distMEs = ifelse(is.na(distMEs), 0, distMEs)
Dist_Tree = hclust(as.dist(distMEs),method="average") 
Dist_ME   = cmdscale(as.dist(distMEs),2)
colors = names(table(modulesPRE))
names = row.names((expr))

pdf("Dynamic_Module_Visualization.pdf",height=5,width=9)
par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 3) + 0.1, cex=1)
plot(Dist_Tree, xlab="",ylab="",main="",sub="")
plot(Dist_ME, col= colors,  main="MDS plot", cex=2, pch=19)

for (which.module in names(table(modulesPRE)))
{
  par(mfrow=c(2,1), mar=c(4, 4.1, 4.1, 2))
  plotMat(t(scale(expr[,modulesPRE==which.module])),cex.axis=2,nrgcols=100,rlabels=F,tck=0, rcols=which.module,main=paste("Heatmap",which.module,"Module"))
  
  ME = MEs0[, paste("ME",which.module, sep="")] 
  n<- barplot(ME, col=which.module, cex.main=1, ylab="Eigengene Expression",xlab="")
  axis(1,at=n, labels=row.names(expr), las=2, cex.axis=0.5, font=2)
};
dev.off();

#||----------------------------------------------------------------------------||-----------
#||----------------------------------------------------------------------------||-----------
#||                           ||ADD TIME DATA||                                ||
#||----------------------------------------------------------------------------||-----------
#||----------------------------------------------------------------------------||-----------
#|
#||-------------------------------------------------------||-----
#|| Dunnett's Test                                        ||
#||-------------------------------------------------------||-----

MEs = orderMEs(MEs0)
dunnetts = data.frame(matrix(nrow=(ncol(MEs)), ncol=5))
colnames(dunnetts) = c("Post_Pre", "D1_Pre", "Wk1_Pre", "Wk2_Pre", "Wk4_Pre")
rownames(dunnetts) = colnames(MEs)

for (i in (1:(ncol(MEs)))) {
  split = data.frame(MEs[,i])
  rownames(split) = rownames(MEs)
  colnames(split) = "Value"
  split$Time = substring((rownames(split)), 12)
  split$Time = factor(split$Time, levels = c("Pre", "Post", "D1", "Wk1", "Wk2", "Wk4"))
  dun = DunnettTest(x=split$Value, g=split$Time)
  pvals = dun$Pre[,4]
  if(any(pvals < 0.1)) {
    cat(paste((colnames(MEs)[i])), ": ", 
        names(pvals[which(pvals < 0.1)]),
        pvals[which(pvals < 0.1)], fill = TRUE)}
  return = t(dun$Pre[,4])
  dunnetts[(paste0((colnames(MEs)[i]))),] = return
  assign((paste0(colnames(MEs)[i])), dun) }

# > Save Dunnetts test output
write.csv(dunnetts, "DunnettsTest_timepoint_moduleCorr.csv")

#||-------------------------------------------------------||-----
#|| Module Boxplots                                       ||
#|| of those significantly associated with time           ||
#||-------------------------------------------------------||-----

Time = meta_wgcna$Time

par(mfrow = c(3,2))
boxplot(MEs$MEtan~Time,col="tan",las=1, cex.axis=0.6, outline=F);abline(v=7.5);abline(h=0, lty=2)
boxplot(MEs$MEviolet~Time,col="violet",las=1, cex.axis=0.6, outline=F);abline(v=7.5);abline(h=0, lty=2)
boxplot(MEs$MEsteelblue~Time,col="steelblue",las=1, cex.axis=0.6, outline=F);abline(v=7.5);abline(h=0, lty=2)
boxplot(MEs$MEyellow~Time,col="yellow",las=1, cex.axis=0.6, outline=F);abline(v=7.5);abline(h=0, lty=2)
boxplot(MEs$MEdarkgreen~Time,col="darkgreen",las=1, cex.axis=0.6, outline=F);abline(v=7.5);abline(h=0, lty=2)
boxplot(MEs$MEpurple~Time,col="purple",las=1, cex.axis=0.6, outline=F);abline(v=7.5);abline(h=0, lty=2)

#||-------------------------------------------------------||-----
#|| Extract Genes in Time-Associated Modules              ||
#||-------------------------------------------------------||-----

modNames = substring(names(MEs), 3)
nGenes = ncol(expr) 
nSamples = nrow(expr) 
geneModuleMembership = as.data.frame(cor(expr,MEs,use = "p")) #PEARSON CORRELATION MEMBERSHIP OF EACH GENE TO A MODULE 
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples)) #P VALUE SIGNIFICANCE FOR EACH GENE IN EACH MODULE

MM_names = names(geneModuleMembership)
names(geneModuleMembership) = paste("MM.",MM_names,sep="")
MM_names = substring(MM_names,3,length(MM_names))
names(MMPvalue) = paste("Mp.",MM_names,sep="")

Colors=modulesPRE
tempout<-cbind(geneModuleMembership,MMPvalue,Colors)
sortedout<-tempout[,sort(names(tempout))]
geneinfo<-sortedout[order(sortedout[,1]),]

# > Save full network output of all modules
write.csv(geneinfo,file="WGCNA_DynamicMods_Network_output.csv")

module_genes = split(geneinfo, geneinfo$Colors)

for (i in (1:(ncol(MEs)))) {
  genes_column = rownames(module_genes[[i]])
  genes = data.frame(genes_column)
  genes$Module = names(module_genes[i])
  assign((paste0(names(module_genes[i]))), genes)}

mod_assign = rbind(tan,violet,steelblue,yellow,darkgreen,purple)

mod_assign = apply(mod_assign,2,as.character)

# > Save Time-Associated Modules Gene lists
write.table(mod_assign, "EnsIDGene_TimeModules_Genes.txt",  
            quote = FALSE, sep = "\t", row.names = FALSE)

mod_assign = rbind(darkred, green, pink,tan,violet,steelblue,saddlebrown, 
                   yellow,darkturquoise, midnightblue, turquoise, darkolivegreen, 
                   paleturquoise, grey60, lightyellow, skyblue, white, darkgrey, 
                   salmon, black, greenyellow, lightgreen, darkorange, blue, 
                   royalblue, lightcyan,  magenta, cyan, brown, darkgreen,purple,
                   orange, red, grey)
# > Save All Modules Gene lists
write.table(mod_assign, "EnsIDGene_AllModules_Genes.txt",  
            quote = FALSE, sep = "\t", row.names = FALSE)

