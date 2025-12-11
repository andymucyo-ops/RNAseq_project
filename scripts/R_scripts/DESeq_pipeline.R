#load libraries 
library("DESeq2")
library("clusterProfiler")
library("pheatmap")
library("ggrepel")
library("ggplot2")
library("EnhancedVolcano")


# ------------------------------------------------------------------------------
# 1) Data formatting before analysis
# ------------------------------------------------------------------------------
# set file path
sample_counts_path <- "~/R-projects/RNAseq_project/data/count.txt"


# save the feature count file as a table
countdata <- read.table(sample_counts_path, header = TRUE)
head(countdata)

 
# changing the count into gene id
rownames(countdata) <- countdata$Geneid 
countdata <- countdata[, -c(1)] #remove the gene id column
head(countdata)


# Create a coldata data frame containing the data information according to
# README in reads_Blood, acting as the metadata of the counts
coldata <- data.frame( Sample = c("SRR7821949", "SRR7821950",
"SRR7821951", "SRR7821952", "SRR7821953", "SRR7821968", "SRR7821969",
"SRR7821970", "SRR7821954", "SRR7821955", "SRR7821956", "SRR7821957",
"SRR7821971", "SRR7821972", "SRR7821973"), Genotype = c(rep("WT", 8),
rep("DKO", 7)), Condition = c(rep("Case", 5), rep("Control", 3),
rep("Case", 4), rep("Control", 3)), row.names = NULL )

head(coldata, n=10L)

# converting sample into factor to match DESeq2 input format, and order by
# the sample names to match column names of the count data table
coldata$Sample <- factor(coldata$Sample, levels =
sort(unique(coldata$Sample)))

# saving the sorted meta data as a new data frame
sorted_coldata <- coldata[order(coldata$Sample),]

# setting sample names as rownames
rownames(sorted_coldata) <- sorted_coldata$Sample

# remonving sample names column
sorted_coldata <-sorted_coldata[,-c(1)]

head(sorted_coldata, n=10L)

# final check before launching DESeq
all(rownames(sorted_coldata) == colnames(countdata))

#-------------------------------------------------------------------------------
# 2) start Differential gene expression analysis
#-------------------------------------------------------------------------------
# create DESeq data set from the feature counts data and the counts meta data
dds <- DESeqDataSetFromMatrix( countdata, sorted_coldata, design = ~
Genotype + Condition + Genotype:Condition )

#transform the meta data values into factor, to match DESeq expected data types, and set references for downstream data analysis 
dds$Genotype <- relevel(dds$Genotype, ref = "WT")

dds$Condition <- relevel(dds$Condition, ref = "Control")

# run differential expression analysis
dds <- DESeq(dds)

#------------------------------------------------------------------------------- 
# 3) Results and data visualization
#------------------------------------------------------------------------------- 
# extraction of transformed values, and principal component plot of the samples
vsd <- vst(dds, blind=FALSE) # extracts the corrected counts for each entry of the differential expression matrix 

png(filename = "./results/R_plots/PCAplot.png",
    width = 3000,
    height = 2400,
    res = 300)

plotPCA(vsd,
  intgroup = c("Genotype","Condition"),
  ntop = 500,
  ) # visualization plot of the clustering of the data

dev.off()

# get results coefficient to use for highlight in the following result step
resultsNames(dds)

# get results comparison with 3 different methods 
res_Genotype <- DESeq2::results(dds, name = "Genotype_DKO_vs_WT", alpha = 0.05) 
# highlights the effect of the Genotype (DKO vs WT) in control condition, enables observation of the genes affected by DKO, regardless of the condition
# expected in heatmap: clear difference between 

res_Condition <- DESeq2::results(dds, name = "Condition_Case_vs_Control", alpha = 0.05)
# highlights the response of the gene to the condition in WT (control vs case), and enables for comparison of the response in DKO 

res_Interaction <- DESeq2::results(dds, name = "GenotypeDKO.ConditionCase", alpha = 0.05)
# highlight the response of the genes to the condition in WT and DKO
# enables for verification of the genes for which the response to the condition depends on the Genotype


# get summary of the result data output for each comparison method
summary_Genotype <- summary(res_Genotype)
summary_Condition <- summary(res_Condition)
summary_Interaction <- summary(res_Interaction)


# get the significant genes for each result method, then reduce it to the top 50 differentially expressed genes
sig_genes_Genotype <- rownames(res_Genotype)[which(!is.na(res_Genotype$padj) & 
                                        res_Genotype$padj < 0.05 & 
                                        abs(res_Genotype$log2FoldChange) > 1)]

top50_genes_Genotype <- sig_genes_Genotype[1:50]


sig_genes_Condition <- rownames(res_Condition)[which(!is.na(res_Condition$padj) & 
                                        res_Condition$padj < 0.05 & 
                                        abs(res_Condition$log2FoldChange) > 1)]

top50_genes_Condition <- sig_genes_Condition[1:50]


sig_genes_Interaction <- rownames(res_Interaction)[which(!is.na(res_Interaction$padj) & 
                                        res_Interaction$padj < 0.05 & 
                                        abs(res_Interaction$log2FoldChange) > 1)]

top50_genes_Interaction <- sig_genes_Interaction[1:50]

 
# extreact the transformed count matching the to 50 genes for each condition
mat_Genotype <- assay(vsd)[top50_genes_Genotype, ]
mat_Condition <- assay(vsd)[top50_genes_Condition, ]
mat_Interaction <- assay(vsd)[top50_genes_Interaction, ] 


# heatmap for different highlighted result data frame
pheatmap(mat_Genotype,
         scale = "row",
         main = "Genotype Effect (DKO vs WT in Control)",
         annotation_col = as.data.frame(colData(dds)[, c("Genotype", "Condition")]), 
         show_rownames = FALSE,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         fontsize = 10, 
         fontsize_row = 8,
         filename = "./results/R_plots/heatmap_Genotype.png"
         )

pheatmap(mat_Condition,
         scale = "row",
         main = "Condition Effect (Case vs Control in WT)",
         annotation_col = as.data.frame(colData(dds)[, c("Genotype", "Condition")]), 
         show_rownames = FALSE, 
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         fontsize = 10, 
         fontsize_row = 8,
         filename = "./results/R_plots/heatmap_Condition.png"
         )

pheatmap(mat_Interaction,
         scale = "row",
         main = "Interaction Effect",
         annotation_col = as.data.frame(colData(dds)[,c("Genotype", "Condition")]), 
         show_rownames = FALSE,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         fontsize = 10, 
         fontsize_row = 8,
         filename = "./results/R_plots/heatmap_Interaction.png"
         )

# volcano plot for each result set and save them as png files 

# Genotype DKO vs WT volcano plot
png(filename = "./results/R_plots/volcano_plot_Genotype.png",
    width = 3000,
    height = 2400,
    res = 300)

EnhancedVolcano(
  res_Genotype,
  lab = rownames(res_Genotype),
  x = 'log2FoldChange',
  y = 'padj',
  title = 'Genotype Effect: DKO vs WT',
  subtitle = 'In Control condition'
  )

dev.off()

# Condition case vs control volcano plot
png(filename = "./results/R_plots/volcano_plot_Condition.png",
    width = 3000,
    height = 2400,
    res = 300)

EnhancedVolcano(
  res_Condition,
  lab = rownames(res_Condition),
  x = 'log2FoldChange',
  y = 'padj',
  title = 'Condition Effect: Case vs Control',
  subtitle = 'In WT genotype'
  )

dev.off()

# Interaction volcano plot
png(filename = "./results/R_plots/volcano_plot_Interaction.png",
    width = 3000,
    height = 2400,
    res = 300)

EnhancedVolcano(
  res_Interaction,
  lab = rownames(res_Interaction),
  x = 'log2FoldChange',
  y = 'padj',
  title = 'Interaction Effect',
  subtitle = 'Differential response to Case in DKO vs WT'
  )

dev.off()

# ------------------------------------------------------------------------------


