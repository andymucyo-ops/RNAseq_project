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
res_Genotype <- results(dds, name = "Genotype_DKO_vs_WT", alpha = 0.05) 
# highlights the effect of the Genotype (DKO vs WT) in control condition, enables observation of the genes affected by DKO, regardless of the condition
# expected in heatmap: clear difference between 

res_Condition <- results(dds, name = "Condition_Case_vs_Control", alpha = 0.05 )
# highlights the response of the gene to the condition in WT (control vs case), and enables for comparison of the response in DKO 

res_Interaction <- results(dds, name = "GenotypeDKO.ConditionCase", alpha = 0.05)
# highlight the response of the genes to the condition in WT and DKO
# enables for verification of the genes for which the response to the condition depends on the Genotype

# get summary of the result data output for each comparison method
summary_Genotype <- summary(res_Genotype)
summary_Condition <- summary(res_Condition)
summary_Interaction <- summary(res_Interaction)


#transform results data into data frame
res_Genotype_df <- as.data.frame(res_Genotype)
res_Condition_df <- as.data.frame(res_Condition)
res_Interaction_df <- as.data.frame(res_Interaction)

#remove NA values from data frame 
res_Genotype_df <- res_Genotype_df[!is.na(res_Genotype_df$padj),]
res_Condition_df <- res_Condition_df[!is.na(res_Condition_df$padj),]
res_Interaction_df <- res_Interaction_df[!is.na(res_Interaction_df$padj),]

# selection of the top 50 genes for each result data frame 
# first we create a vector containing the to 50 differentialy expressed genes in each result table based on padj values 
# then create a matrix containg each column of the transformed count matrix but only the rows of the most differentialy expressed genes
top_genes_Genotype <- head(order(res_Genotype_df$padj), 30)
mat_Genotype <- assay(vsd)[top_genes_Genotype, ]

top_genes_Condition <- head(order(res_Condition_df$padj), 30)
mat_Condition <- assay(vsd)[top_genes_Condition, ]

top_genes_Interaction <- head(order(res_Interaction_df$padj), 30)
mat_Interaction <- assay(vsd)[top_genes_Interaction, ] 


clean_matrix_for_heatmap <- function(mat) {
  # 1. Supprimer les lignes avec des NA
  mat <- mat[complete.cases(mat), ]
  
  # 2. Supprimer les lignes avec variance nulle
  row_vars <- apply(mat, 1, var)
  mat <- mat[row_vars > 0, ]
  
  # 3. Supprimer les valeurs infinies
  mat <- mat[is.finite(rowSums(mat)), ]
  
  return(mat)
}

mat_Genotype_clean <-  clean_matrix_for_heatmap(mat_Genotype)
  
mat_Condition_clean <- clean_matrix_for_heatmap(mat_Condition)
  
mat_Interaction_clean <- clean_matrix_for_heatmap(mat_Interaction)
# heatmap for different highlighted result data frame
pheatmap(mat_Genotype_clean,
         #scale = "row",
         main = "Genotype Effect (DKO vs WT in Control)",
         annotation_col = as.data.frame(colData(dds)[, c("Genotype", "Condition")]), 
         show_rownames = TRUE,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         fontsize = 10, 
         fontsize_row = 8,
         filename = "./results/R_plots/heatmap_Genotype.png"
         )

pheatmap(mat_Condition_clean,
         #scale = "row",
         main = "Condition Effect (Case vs Control in WT)",
         annotation_col = as.data.frame(colData(dds)[, c("Genotype", "Condition")]), 
         show_rownames = TRUE, 
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         fontsize = 10, 
         fontsize_row = 8,
         filename = "./results/R_plots/heatmap_Condition.png"
         )

pheatmap(mat_Interaction_clean,
         #scale = "row",
         main = "Interaction Effect",
         annotation_col = as.data.frame(colData(dds)[,c("Genotype", "Condition")]), 
         show_rownames = TRUE,
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
  res_Genotype_df,
  lab = rownames(res_Genotype_df),
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
  res_Condition_df,
  lab = rownames(res_Condition_df),
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
  res_Interaction_df,
  lab = rownames(res_Interaction_df),
  x = 'log2FoldChange',
  y = 'padj',
  title = 'Interaction Effect',
  subtitle = 'Differential response to Case in DKO vs WT'
  )

dev.off()

# ------------------------------------------------------------------------------


