#load libraries

library("DESeq2")
library("clusterProfiler")
library("pheatmap")
library("ggplot2")
library("ggrepel")
library("EnhancedVolcano")
library("org.Mm.eg.db")
library("enrichplot")

#set wroking directory
setwd("~/Projects/R/RNAseq_project/")

# ---------------------------------------------------------------------
# 1) Data formatting before analysis
# ---------------------------------------------------------------------

# set file path
sample_counts_path <- "~/Projects/R/RNAseq_project/data/count.txt"

# save the feature count file as a table
counts <- read.table(sample_counts_path, header = TRUE)

# changing rownames to gene id
rownames(counts) <- counts$Geneid

#remove Geneid column
counts$Geneid <- NULL
head(counts)


# Create a sample_info data frame containing the Data information according to README in reads_Blood, acting as the metadata of the counts

sample_info_path <- "~/Projects/R/RNAseq_project/data/sample_info.txt"

sample_info <- read.table(sample_info_path, header = TRUE)

#order sample_info acording to sample name
sample_info <- sample_info[order(sample_info$Sample), ]

#change rownames to sample names
rownames(sample_info) <- sample_info$Sample

#remove sample column
sample_info$Sample <- NULL
head(sample_info)

# final check before launching DESeq (Must return TRUE to perform DESeq analysis!)
all(rownames(sample_info) == colnames(counts))


#-----------------------------------------------------------------------
# 2) start Differential gene expression analysis
#-----------------------------------------------------------------------

# create DESeq data set from the feature counts data and the counts meta data(sample_info)
dds <- DESeqDataSetFromMatrix(
  counts,
  sample_info,
  design = ~Group
)

# change Group data type to match expected type in DESeq() function
dds$Group <- as.factor(dds$Group)


# run differential expression analysis
dds <- DESeq(dds)

#-----------------------------------------------------------------------
# 3) Results and data visualization
#-----------------------------------------------------------------------

# extraction of transformed values (corrected counts) of counts into new matrix to nromalize counts and compare them
vsd <- vst(dds, blind = FALSE)


#plot PCA (principal component analysis) to check for clustering of the data, and save as imageb file
png(
  filename = "./results/R_plots/PCAplot.png",
  width = 3000,
  height = 2400,
  res = 300
)

plotPCA(
  vsd,
  intgroup = c("Group"),
  ntop = 500,
)

dev.off()

#-----------------------------------------------------------------------
# get results comparison with 4 different methods
#-----------------------------------------------------------------------
#Highlight differences of expression between case and control condition in WT genotype
res_WT <- results(
  dds,
  contrast = c(
    "Group",
    "Blood_WT_Case",
    "Blood_WT_Control"
  ),
  alpha = 0.05
)

# highlight of the difference of expression between case and control condition in DKO genotype
res_DKO <- results(
  dds,
  contrast = c(
    "Group",
    "Blood_DKO_Case",
    "Blood_DKO_Control"
  ),
  alpha = 0.05
)

# highlight the difference of expression in control condition between WT and DKO in Case condition
res_DKOvsWT_case <- results(
  dds,
  contrast = c(
    "Group",
    "Blood_DKO_Case",
    "Blood_WT_Case"
  ),
  alpha = 0.05
)

# highlight the difference of expression in control condition between WT and DKO in control condition
res_DKOvsWT_control <- results(
  dds,
  contrast = c(
    "Group",
    "Blood_DKO_Control",
    "Blood_WT_Control"
  ),
  alpha = 0.05
)

# get summary of the result data output for each comparison method
summary_WT <- summary(res_WT)
summary_DKO <- summary(res_DKO)
summary_DKOvsWT_case <- summary(res_DKOvsWT_case)
summary_DKOvsWT_control <- summary(res_DKOvsWT_control)


#Create function to retrieve top 50 genes for each result sets
get_top_genes <- function(res) {
  #set filtering criteria for the results with padj defined and absolute value of log2FoldChange > 1
  filtering_criteria <- which(
    !is.na(res$padj) &
      abs(res$log2FoldChange) > 1
  )

  # get significant genes order by padj value
  sig_genes <- rownames(res)[filtering_criteria][order(res$padj[
    filtering_criteria
  ])]

  top50_genes <- sig_genes[1:50]

  return(top50_genes)
}


#-----------------------------------------------------------------------
#get to 50 genes for each comparison
#-----------------------------------------------------------------------
top50_WT <- get_top_genes(res_WT)

top50_DKO <- get_top_genes(res_DKO)

top50_DKOvsWT_case <- get_top_genes(res_DKOvsWT_case)

top50_DKOvsWT_control <- get_top_genes(res_DKOvsWT_control)

# extract the corrected count matrix matching the top 50 genes for each results highlight
mat_WT <- assay(vsd)[top50_WT, ]
mat_DKO <- assay(vsd)[top50_DKO, ]
mat_DKOvsWT_case <- assay(vsd)[top50_DKOvsWT_case, ]
mat_DKOvsWT_control <- assay(vsd)[top50_DKOvsWT_control, ]

#---------------------------------------------------------------------
# heatmap for different highlighted result data frame
#---------------------------------------------------------------------
pheatmap(
  mat_WT,
  scale = "row",
  main = "Case vs Control condition comparison in WT genotype",
  annotation_col = data.frame(
    Group = colData(dds)$Group,
    row.names = colnames(mat_WT)
  ),
  show_rownames = FALSE,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize = 10,
  fontsize_row = 8,
  filename = "./results/R_plots/heatmap_WT.png"
)

pheatmap(
  mat_DKO,
  scale = "row",
  main = "Case vs Control condition comparison in DKO genotype",
  annotation_col = data.frame(
    Group = colData(dds)$Group,
    row.names = colnames(mat_DKO)
  ),
  show_rownames = FALSE,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize = 10,
  fontsize_row = 8,
  filename = "./results/R_plots/heatmap_DKO.png"
)

pheatmap(
  mat_DKOvsWT_case,
  scale = "row",
  main = "DKO vs WT in case condition",
  annotation_col = data.frame(
    Group = colData(dds)$Group,
    row.names = colnames(mat_DKOvsWT_case)
  ),
  show_rownames = FALSE,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize = 10,
  fontsize_row = 8,
  filename = "./results/R_plots/heatmap_DKOvsWT_case.png"
)

pheatmap(
  mat_DKOvsWT_control,
  scale = "row",
  main = "DKO vs WT in control condition",
  annotation_col = data.frame(
    Group = colData(dds)$Group,
    row.names = colnames(mat_DKOvsWT_control)
  ),
  show_rownames = FALSE,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize = 10,
  fontsize_row = 8,
  filename = "./results/R_plots/heatmap_DKOvsWT_control.png"
)

#-----------------------------------------------------------------------
# volcano plot for each result set and save them as png files
#-----------------------------------------------------------------------

# Case vs Control condition comparison in WT genotype
png(
  filename = "./results/R_plots/volcano_plot_WT.png",
  width = 3000,
  height = 2400,
  res = 300
)

EnhancedVolcano(
  res_WT,
  lab = rownames(res_WT),
  x = 'log2FoldChange',
  y = 'padj',
  title = 'Genotype effect: Case vs Control',
  subtitle = 'Differential response in WT'
)

dev.off()

# Case vs Control condition comparicon in DKO Genotype
png(
  filename = "./results/R_plots/volcano_plot_DKO.png",
  width = 3000,
  height = 2400,
  res = 300
)

EnhancedVolcano(
  res_DKO,
  lab = rownames(res_DKO),
  x = 'log2FoldChange',
  y = 'padj',
  title = 'Genotype effect: Case vs Control',
  subtitle = 'Differential response in DKO'
)

dev.off()

# DKO vs WT genotype comparison in Case condition
png(
  filename = "./results/R_plots/volcano_plot_DKOvsWT_case.png",
  width = 3000,
  height = 2400,
  res = 300
)

EnhancedVolcano(
  res_DKOvsWT_case,
  lab = rownames(res_DKOvsWT_case),
  x = 'log2FoldChange',
  y = 'padj',
  title = 'Condition effect: DKO vs WT',
  subtitle = 'Differential response in Case condition'
)

dev.off()


# Interaction volcano plot
png(
  filename = "./results/R_plots/volcano_plot_DKOvsWT_control.png",
  width = 3000,
  height = 2400,
  res = 300
)

EnhancedVolcano(
  res_DKOvsWT_case,
  lab = rownames(res_DKOvsWT_control),
  x = 'log2FoldChange',
  y = 'padj',
  title = 'Condition effect: DKO vs WT',
  subtitle = 'Differential response in Control condition'
)

dev.off()

# -----------------------------------------------------------------------
# 4) enrichGO analysis
# -----------------------------------------------------------------------

# save gene list of differentialy expressed genes, for each results data
gene_list_WT <- top50_WT

gene_list_DKO <- top50_DKO

gene_list_DKOvsWT_case <- top50_DKOvsWT_case

gene_list_DKOvsWT_control <- top50_DKOvsWT_control
# set universe list (ensemble of all genes), cf. summary(res)

universe_WT <- rownames(res_WT)[which(!is.na(res_WT$padj))]

universe_DKO <- rownames(res_DKO)[which(!is.na(res_DKO$padj))]

universe_DKOvsWT_case <- rownames(res_DKOvsWT_case)[which(
  !is.na(res_DKOvsWT_case$padj)
)]

universe_DKOvsWT_control <- rownames(res_DKOvsWT_control)[which(
  !is.na(res_DKOvsWT_control$padj)
)]


# create enrivhGO (ego) object for each result

ego_WT <- enrichGO(
  gene = gene_list_WT,
  universe = universe_WT,
  OrgDb = org.Mm.eg.db,
  ont = "BP",
  keyType = "ENSEMBL",
  readable = TRUE
)

head(ego_WT, 3)

ego_DKO <- enrichGO(
  gene = gene_list_DKO,
  universe = universe_DKO,
  OrgDb = org.Mm.eg.db,
  ont = "BP",
  keyType = "ENSEMBL",
  readable = TRUE
)

head(ego_DKO, 3)

ego_DKOvsWT_case <- enrichGO(
  gene = gene_list_DKOvsWT_case,
  universe = universe_DKOvsWT_case,
  OrgDb = org.Mm.eg.db,
  ont = "BP",
  keyType = "ENSEMBL",
  readable = TRUE
)

head(ego_DKOvsWT_case, 3)

ego_DKOvsWT_control <- enrichGO(
  gene = gene_list_DKOvsWT_control,
  universe = universe_DKOvsWT_control,
  OrgDb = org.Mm.eg.db,
  ont = "BP",
  keyType = "ENSEMBL",
  readable = TRUE
)

head(ego_DKOvsWT_control, 3)

#-----------------------------------------------------------------------
# ego's bar plots
#-----------------------------------------------------------------------
png(
  filename = "./results/R_plots/barplot_enrichGO_WT.png",
  width = 3000,
  height = 2400,
  res = 300
)
barplot(ego_WT, x = "Count") +
  ggtitle(
    "Enriched terms bar plot",
    "In WT: case vs control comparison"
  )
dev.off()


png(
  filename = "./results/R_plots/barplot_enrichGO_DKO.png",
  width = 3000,
  height = 2400,
  res = 300
)
barplot(ego_DKO, x = "Count") +
  ggtitle("Enriched terms bar plot", "In DKO: case vs control comparison")
dev.off()


png(
  filename = "./results/R_plots/barplot_enrichGO_DKOvsWT_case.png",
  width = 3000,
  height = 2400,
  res = 300
)
barplot(ego_DKOvsWT_case, x = "Count") +
  ggtitle(
    "Enriched terms bar plot",
    "In Case condition: DKO vs WT comparison"
  )
dev.off()


png(
  filename = "./results/R_plots/barplot_enrichGO_DKOvsWT_control.png",
  width = 3000,
  height = 2400,
  res = 300
)
barplot(ego_DKOvsWT_control, x = "Count") +
  ggtitle(
    "Enriched terms bar plot",
    "In Control condition: DKO vs WT comparison"
  )
dev.off()
