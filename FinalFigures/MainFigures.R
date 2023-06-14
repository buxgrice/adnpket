## ---------------------------
##
## Script name:Recreating Final Figures in this Manuscript  
##
## Purpose of script: Generating main figures 1-3 of this paper
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
#||                  ||FINAL FIGURES FOR KETAMINE ADNP PAPER||                 ||
#||                                                                            ||
#||----------------------------------------------------------------------------||-----------
#||----------------------------------------------------------------------------||-----------

library(easypackages)

# > Figure 1
libraries("trackViewer","stringr")
# > Figure 2
libraries("reshape","ggplot2", "biomaRt", "ggrepel", 
          "dplyr", "pheatmap", "corrplot")
# > Figure 3
libraries("ggplot2") 


#||----------------------------------------------------------------------------||-----------
#||                       ||SET BASIC WORKING DIRECTORY||                      ||
#||----------------------------------------------------------------------------||-----------

setwd ("/Users/arielabuxbaumgrice/Desktop/projects_new/ADNP_ketamine")
getwd()


#||----------------------------------------------------------------------------||-----------
#||----------------------------------------------------------------------------||-----------
#||                               ||FIGURE 1A||                                ||
#||----------------------------------------------------------------------------||-----------
#||----------------------------------------------------------------------------||-----------
#Clinical trial workflow and data collected across the timecourse. 

#||----------------------------------------------------------------------------||-----------
#||----------------------------------------------------------------------------||-----------
#||                               ||FIGURE 1B||                                ||
#||----------------------------------------------------------------------------||-----------
#||----------------------------------------------------------------------------||-----------

#Pathogenic variants their locations along the ADNP gene locus, with percent mutated allele shows

# ADNP : 50888918-50931437, length = 42519, 1102 Amino acids
#Zinc Finger: 
#  AA: 74-97, 107-129, 165-188, 221-244, 447-469, 
#489-510, 512-535, 622-647, 662-686
#DNA Binding: 
#  AA: 754-814
#NAP Region: (NAPVSIPQ)
#  AA: 354-361

pile = read.csv("MPileUp_Results.csv", row.names = 1)
pile$Time = str_sub(rownames(pile), start=19)
pile$Time = factor(pile$Time, levels = c("Pre","Post","D1", "Wk1", "Wk2", "Wk4"))
dim(pile) #[1] 51 15
pile$Donor = str_sub(rownames(pile), start=8, end = 13)

mutations = pile[(pile$Time == "Pre"),]
# missing Pre: Sample_CC1680, Sample_CC1360: using the 4 week timepoint for these
mutations[9,] = pile[c("Sample_CC1360.201-Wk4"),]
mutations[10,] = pile[c("Sample_CC1680-202-Wk4"),]

mutations$pM = 1-(as.numeric(mutations$pR))
mutations$start = NA

mutations["Sample_CC1368-202-Pre","start"] = as.numeric(50892557) 
mutations["Sample_CC1420-204-Pre","start"] = as.numeric(50892557)
mutations["Sample_CC1566-203-Pre","start"] = as.numeric(50892501)
mutations["Sample_CC1327-202-Pre","start"] = as.numeric(50892215)
mutations["Sample_CC1329-202-Pre","start"] = as.numeric(50894352)
mutations["Sample_CC1409-203-Pre","start"] = as.numeric(50893606)
mutations["Sample_CC1628-201-Pre","start"] = as.numeric(50894374)
mutations["Sample_CC1668-201-Pre","start"] =  as.numeric(50893895)
mutations["Sample_CC1360.201-Wk4","start"] =  as.numeric(50894375)
mutations["Sample_CC1680-202-Wk4","start"] =  as.numeric(50894063)

features = GRanges("chr20", IRanges(c(1, 74, 107, 165, 221, 354, 447, 489, 512, 622, 662, 754), 
                                    # ^ amino acid start locations of ADNP, then Zn_Fn 1-4,  NAP, Zn_Fn 5-9, & then homeobox domain 
                                     width=c(1102, 23, 22, 23, 23, 7, 22, 21, 23, 25, 24, 60), 
                                    # ^ lengths of ADNP, then Zn_Fn 1-4,  NAP, Zn_Fn 5-9, & then homeobox domain
                                     height = c(0.075), names = c("ADNP", "Zn_Fn1", "Zn_Fn2", "Zn_Fn3", "Zn_Fn4", 
                                                                  "NAP Region", "Zn_Fn5", "Zn_Fn6", "Zn_Fn7", 
                                                                  "Zn_Fn8", "Zn_Fn9", "DNA-binding homeobox domain"), 
                                     fill = c("grey54", "thistle4", "thistle4",  "thistle4", "thistle4", 
                                              "gold","thistle4","thistle4", "thistle4", 
                                              "thistle4", "thistle4", "paleturquoise3")))

SNP = c(113, 114, 121, 218, 274, 369, 719, 719, 738, 832) # start positions of mutations
sample.gr = GRanges("chr20", IRanges(SNP, width=1, names= c("p.Thr113Serfs*47", "p.Phe114Serfs*47",
                                                             "p.Leu121Glyfs*5",  "p.Glu218*", 
                                                             "p.Lys274Asnfs*31", "p.Leu369Serfs*30", 
                                                             "p.Tyr719*", "p.Tyr719*",
                                                             "p.Ser738*", "p.Asn832Lysfs*81")),
                     color = c("darkgoldenrod", "darkgoldenrod", "darkgoldenrod","darkgoldenrod",
                               "darkgoldenrod", "darkgoldenrod", "coral3", "coral3", 
                               "darkgoldenrod", "darkgoldenrod"))
#lolliplot(sample.gr, features)

# ADD PIE CHARTS FOR ALLELIC FREQUENCIES 
sample.gr$score = NULL 
sample.gr$label = NULL
sample.gr$node.label.col = NULL
sample.gr$pR = (mutations$pR)*100
sample.gr$pM = 100 - sample.gr$pR
sample.gr$color = rep(list(c("goldenrod3", 'coral3')), length(SNP))
sample.gr$border = "gray30"
sample.gr$label.parameter.rot = 45

# > Returns the final figure < 
lolliplot(sample.gr, features, type="pie", yaxis=FALSE)


#||----------------------------------------------------------------------------||-----------
#||----------------------------------------------------------------------------||-----------
#||                               ||FIGURE 2A||                                ||
#||----------------------------------------------------------------------------||-----------
#||----------------------------------------------------------------------------||-----------

counts = read.csv("SignificantDEGCounts_Figure2A.csv")
counts$comp = c("PrePost", "PreD1", "PreWk1", "PreWk2", "PreWk4")
counts$comp = factor(counts$comp, 
                     levels = c("PrePost", "PreD1", "PreWk1", "PreWk2", "PreWk4"),
                     labels = c(" vs. Post", " vs. Day 1", "vs. Week 1", 
                                "vs. Week 2", "vs. Week 4"))
plot = melt(counts)

# > Returns the final figure < 
ggplot(plot, aes(x=comp, y=value, fill = variable)) + 
  geom_bar(stat="identity", position=position_dodge()) + 
  geom_point(position=position_jitterdodge()) +
  scale_fill_manual(values=c("indianred3", "steelblue4")) + 
  theme_classic()


#||----------------------------------------------------------------------------||-----------
#||----------------------------------------------------------------------------||-----------
#||                               ||FIGURE 2B||                                ||
#||----------------------------------------------------------------------------||-----------
#||----------------------------------------------------------------------------||-----------

DEGs_T1 = read.delim("DEGs_Time_PREvPOST.txt", row.names = 1)
results.df = data.frame(DEGs_T1)

padj_cutoff = 0.05

results.df = results.df[order(results.df$adj.P.Val),]
results.df = results.df[complete.cases(results.df[,6]),]
up = results.df[(results.df$logFC > log(2)),]
down = results.df[(results.df$logFC < -log(2)),]
top20.up_cutoff = sort(up$adj.P.Val)[20] # the 20th smallest value of res$padj # there are actually 22 here, as some of the pvals are the same
top20.down_cutoff = sort(down$adj.P.Val)[20] #the 20th smallest value of res$padj

results.df$labels = ifelse(
  (((results.df$logFC > log(2)) & (results.df$adj.P.Val <= top20.up_cutoff)) | 
     ((results.df$logFC < -log(2)) & (results.df$adj.P.Val <= top20.down_cutoff))), rownames(results.df), NA) 

labs = rownames(results.df[!is.na(results.df$labels),])

# Convert EnsIDs to gene symbols to label the top genes

ensembl = c('') 

for (i in (1:length(labs))) {
  gene = labs[i]
  ss = unlist(strsplit((results.df[gene,]$label),"[.]"))
  ensembl[i] = ss[1] }

labels = data.frame(ShortID = ensembl, FullID = rownames(results.df[!is.na(results.df$labels),]))
mart = useDataset("hsapiens_gene_ensembl", useMart("ensembl")) 
genes = c(labels$ShortID)
gene_list = getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"), values = genes, mart = mart)
labels$Symbol = gene_list[order(match(gene_list[,1],labels[,1])),]$hgnc_symbol
results.df$labels_symb = NA
rownames(labels) = labels$FullID

rows = c(which(rownames(results.df) %in% rownames(labels)))
names = c(which(rownames(labels) %in% rownames(results.df)))
results.df[c(rows),]$labels_symb = labels[c(names),]$Symbol
results.df["ENSG00000101126.17",]$labels_symb = "ADNP" #Label for ADNP, as well 

t = results.df %>% 
  mutate(
    Significance = case_when(
      logFC > 0 & adj.P.Val <= 0.05  ~ "+ LFC & p.adj 0.05", 
      logFC < 0 & adj.P.Val <= 0.05  ~ "- LFC & p.adj 0.05", 
      TRUE ~ "NS"))

# > Returns the final figure < 
ggplot(t) +
  geom_point(aes(x = logFC, y = -log10(adj.P.Val), color = Significance), 
             alpha = 0.8, size = 1.6) +
  geom_label_repel(max.overlaps = 50, box.padding = 0.6,
                   segment.alpha = 0.3,
                   mapping = aes(logFC, -log10(adj.P.Val), label = labels_symb),
                   size = 2)  + scale_y_continuous(limits=c(0,3)) + 
  scale_color_manual(values=c("steelblue2", "tomato1", "grey70")) + 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = -log(2), linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = log(2), linetype = "dashed", color = "grey40") + 
  theme_classic()

#||----------------------------------------------------------------------------||-----------
#||----------------------------------------------------------------------------||-----------
#||                               ||FIGURE 2C||                                ||
#||----------------------------------------------------------------------------||-----------
#||----------------------------------------------------------------------------||-----------

plot_data = read.csv("DEGintersection_CellCovariate.csv", row.names =1)
plot_data = plot_data[-c(1,4,8,10),]

plot_data2 = cbind(
  100-plot_data$pct_DEGT1orig_DEGT1cell,
  100-plot_data$pct_DEGT2orig_DEGT2cell, 
  100- plot_data$pct_DEGT4orig_DEGT4cell,
  100-plot_data$pct_DEGT5orig_DEGT5cell)
rownames(plot_data2) = rownames(plot_data)
colnames(plot_data2) = c("Pre vs. Post", "Pre vs. Day 1", "Pre vs. Week 2", "Pre vs. Week 4")

cols = colorRampPalette(c("white", "firebrick3"))(30)

# > Returns the final figure < 
pheatmap(plot_data2, color = cols, cluster_rows = FALSE)

#||----------------------------------------------------------------------------||-----------
#||----------------------------------------------------------------------------||-----------
#||                               ||FIGURE 2D||                                ||
#||----------------------------------------------------------------------------||-----------
#||----------------------------------------------------------------------------||-----------

data = read.delim("logFC_all_forValidation.txt", row.names=1)
plot = data[,c(1:5, 11:13)]

colnames(plot) = c("Pre vs. Post", "Pre vs. Day1", "Pre vs. Wk1", 
                   "Pre vs. Wk2", "Pre vs. Wk4", "Ho E2 + Ket",
                   "Ho Ket", "Cathomas Pre vs. Post")

# > Returns the final figure < 
corrplot(cor(plot), method = "circle", is.corr = FALSE, col = COL1('OrRd'))

#||----------------------------------------------------------------------------||-----------
#||----------------------------------------------------------------------------||-----------
#||                               ||FIGURE 2E||                                ||
#||----------------------------------------------------------------------------||-----------
#||----------------------------------------------------------------------------||-----------

cathomas = read.delim("logFC_Cathomas&OurData_forValidation.txt", row.names= 1)
ho = read.delim("logFC_Ho2019&OurData_forValidation.txt", row.names = 1)

# > Returns the final figure < 
ggplot(cathomas, aes(x = PrePost, y =Cathomas)) + 
  stat_density_2d(aes(fill = ..level..), geom = "polygon", color = "white") + 
  scale_fill_distiller(palette = "Reds", direction = 1) +
  geom_smooth(method='lm',color = "grey70", linewidth = .4) + 
  theme_classic() + stat_poly_line() + stat_poly_eq() # R2 = 0.12  

# > Returns the final figure < 
ggplot(ho, aes(x = PrePost, y =Ho.E2.Ket)) + 
  stat_density_2d(aes(fill = ..level..), geom = "polygon", color = "white") + 
  scale_fill_distiller(palette = "Reds", direction = 1) +
  geom_smooth(method='lm',color = "grey70", linewidth = .4) + 
  theme_classic() + 
  stat_poly_line() + stat_poly_eq() # R2 = 0.25

# > Returns the final figure < 
ggplot(ho, aes(x = PrePost, y =Ho.Ket)) + 
  stat_density_2d(aes(fill = ..level..), geom = "polygon", color = "white") + 
  scale_fill_distiller(palette = "Reds", direction = 1) +
  geom_smooth(method='lm',color = "grey70", linewidth = .4) + 
  theme_classic() + 
  stat_poly_line() + stat_poly_eq() # R2 = 0.19

#||----------------------------------------------------------------------------||-----------
#||----------------------------------------------------------------------------||-----------
#||                               ||FIGURE 3A||                                ||
#||----------------------------------------------------------------------------||-----------
#||----------------------------------------------------------------------------||-----------
# See WGCNA_Code

#||----------------------------------------------------------------------------||-----------
#||----------------------------------------------------------------------------||-----------
#||                               ||FIGURE 3B||                                ||
#||----------------------------------------------------------------------------||-----------
#||----------------------------------------------------------------------------||-----------

MEs = read.delim("Orig.DynamicCut_WGCNA_ExprMEs.txt")
meta_wgcna = read.delim("ADNP_WGCNA_metadata.txt")

Time = meta_wgcna$Time

# > Returns the final figure < 
par(mfrow = c(3,2))
boxplot(MEs$MEtan~Time,col="tan",las=1, cex.axis=0.6, outline=F);abline(v=7.5);abline(h=0, lty=2)
boxplot(MEs$MEviolet~Time,col="violet",las=1, cex.axis=0.6, outline=F);abline(v=7.5);abline(h=0, lty=2)
boxplot(MEs$MEsteelblue~Time,col="steelblue",las=1, cex.axis=0.6, outline=F);abline(v=7.5);abline(h=0, lty=2)
boxplot(MEs$MEyellow~Time,col="yellow",las=1, cex.axis=0.6, outline=F);abline(v=7.5);abline(h=0, lty=2)
boxplot(MEs$MEdarkgreen~Time,col="darkgreen",las=1, cex.axis=0.6, outline=F);abline(v=7.5);abline(h=0, lty=2)
boxplot(MEs$MEpurple~Time,col="purple",las=1, cex.axis=0.6, outline=F);abline(v=7.5);abline(h=0, lty=2)

#||----------------------------------------------------------------------------||-----------
#||----------------------------------------------------------------------------||-----------
#||                               ||FIGURE 3C||                                ||
#||----------------------------------------------------------------------------||-----------
#||----------------------------------------------------------------------------||-----------
\