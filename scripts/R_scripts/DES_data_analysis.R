#set wroking directory
setwd("~/R-projects/RNAseq_project/")

#load libraries

library("DESeq2")
library("clusterProfiler")
library("pheatmap")
library("ggplot2")
library("ggrepel")
library("EnhancedVolcano")
library("org.Mm.eg.db")
library("enrichplot")

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

coldata <- data.frame(
	Sample = c(
		"SRR7821949",
		"SRR7821950",
		"SRR7821951",
		"SRR7821952",
		"SRR7821953",
		"SRR7821968",
		"SRR7821969",
		"SRR7821970",
		"SRR7821954",
		"SRR7821955",
		"SRR7821956",
		"SRR7821957",
		"SRR7821971",
		"SRR7821972",
		"SRR7821973"
	),
	Genotype = c(rep("WT", 8), rep("DKO", 7)),
	Condition = c(
		rep("Case", 5),
		rep("Control", 3),
		rep("Case", 4),
		rep("Control", 3)
	),
	row.names = NULL
)

head(coldata, n = 10L)

# converting sample into factor to match DESeq2 input format, and order by
# the sample names to match column names of the count data table
coldata$Sample <- factor(coldata$Sample, levels = sort(unique(coldata$Sample)))

# saving the sorted meta data as a new data frame
sorted_coldata <- coldata[order(coldata$Sample), ]

# setting sample names as rownames
rownames(sorted_coldata) <- sorted_coldata$Sample

# remonving sample names column
sorted_coldata <- sorted_coldata[, -c(1)]

head(sorted_coldata, n = 10L)

# final check before launching DESeq
all(rownames(sorted_coldata) == colnames(countdata))

#-------------------------------------------------------------------------------
# 2) start Differential gene expression analysis
#-------------------------------------------------------------------------------
# create DESeq data set from the feature counts data and the counts meta data

dds <- DESeqDataSetFromMatrix(
	countdata,
	sorted_coldata,
	design = ~ Genotype + Condition + Genotype:Condition
)

#transform the meta data values into factor, to match DESeq expected data types, and set references for downstream data analysis

dds$Genotype <- relevel(dds$Genotype, ref = "WT")
dds$Condition <- relevel(dds$Condition, ref = "Control")

# run differential expression analysis
dds <- DESeq(dds)

#-------------------------------------------------------------------------------
# 3) Results and data visualization
#-------------------------------------------------------------------------------
# extraction of transformed values, and principal component plot of the samples
vsd <- vst(dds, blind = FALSE) # extracts the corrected counts for each entry of the differential expression matrix

png(
	filename = "./results/R_plots/PCAplot.png",
	width = 3000,
	height = 2400,
	res = 300
)

plotPCA(vsd, intgroup = c("Genotype", "Condition"), ntop = 500, ) # visualization plot of the clustering of the data

dev.off()

# get results coefficient to use for highlight in the following result step
resultsNames(dds)

# get results comparison with 3 different methods
res_WT <- results(dds, name = "Condition_Case_vs_Control", alpha = 0.05)
#Highlight differences of expression between control and case in WT

res_DKO <- results(
	dds,
	contrast = list(c(
		"GenotypeDKO.ConditionCase",
		"Condition_Case_vs_Control"
	)),
	alpha = 0.05
)
# highlight of the difference of expression between control and case in DKO

res_DKOvsWT_case <- results(
	dds,
	name = "GenotypeDKO.ConditionCase",
	alpha = 0.05
)

# highlight the response of the genes in control condition between WT and DKO
# enables for verification of the genes for which the response to the condition depends on the Genotype

# get summary of the result data output for each comparison method

summary_WT <- summary(res_WT)
summary_DKO <- summary(res_DKO)
summary_DKOvsWT_case <- summary(res_DKOvsWT_case)


# get the significant genes for each result method, then reduce it to the top 50 differentialy expressed genes

#set filtering criteria for each result set
filtering_criteria <- which(
  !is.na(res_WT$padj) &
    res_WT$padj < 0.05 &
    abs(res_WT$log2FoldChange) > 1
)

# WT
sig_genes_WT <- rownames(res_WT)[filtering_criteria][order(res_WT$padj[filtering_criteria])]
filtered_counts_WT <- length(sig_genes_WT) # nmbr of genes after filtering: 6768
top50_genes_WT <- sig_genes_WT[1:50]

# DKO

sig_genes_DKO <- rownames(res_DKO)[filtering_criteria][order(res_DKO$padj[filtering_criteria])]
filtered_counts_DKO <- length(sig_genes_DKO) # nmbr of genes after filtering: 5450
top50_genes_DKO <- sig_genes_DKO[1:50]

# DKO vs WT
sig_genes_DKOvsWT_case <- rownames(res_DKOvsWT_case)[filtering_criteria][order(res_DKOvsWT_case$padj[filtering_criteria])]
filtered_counts_DKOvsWT <- length(sig_genes_DKOvsWT_case) # nmbr of genes after filtering: 3745
top50_genes_DKOvsWT_case <- sig_genes_DKOvsWT_case[1:50]


# extract the transformed count matrix matching the top 50 genes for each results highlight
mat_WT <- assay(vsd)[top50_genes_WT, ]
mat_DKO <- assay(vsd)[top50_genes_DKO, ]
mat_DKOvsWT_case <- assay(vsd)[top50_genes_DKOvsWT_case, ]


# heatmap for different highlighted result data frame
pheatmap(
	mat_WT,
	scale = "row",
	main = "Case vs Control condition comparison in WT genotype",
	annotation_col = as.data.frame(colData(dds)[, c(
		"Genotype",
		"Condition"
	)]),
	# clustering_distance_cols = "binary",
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
	annotation_col = as.data.frame(colData(dds)[, c(
		"Genotype",
		"Condition"
	)]),
	# clustering_distance_cols = "binary",
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
	annotation_col = as.data.frame(colData(dds)[, c(
		"Genotype",
		"Condition"
	)]),
	# clustering_distance_cols = "binary",
	show_rownames = FALSE,
	cluster_rows = TRUE,
	cluster_cols = TRUE,
	fontsize = 10,
	fontsize_row = 8,
	filename = "./results/R_plots/heatmap_DKOvsWT_case.png"
)

# volcano plot for each result set and save them as png files

# Genotype DKO vs WT volcano plot
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

# Condition case vs control volcano plot
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

# Interaction volcano plot
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
	subtitle = 'Differential response in Control control condition'
)

dev.off()

# ------------------------------------------------------------------------------

#start enrichGO analysis

# save gene list of differentialy expressed genes, for each results data
gene_list_WT <- sig_genes_WT

gene_list_DKO <- sig_genes_DKO

gene_list_DKOvsWT_case <- sig_genes_DKOvsWT_case


# set universe list (ensemble of all genes), cf. summary(res)

universe_WT <- rownames(res_WT)[which(!is.na(res_WT$padj))]

universe_DKO <- rownames(res_DKO)[which(!is.na(res_DKO$padj))]

universe_DKOvsWT <- rownames(res_DKOvsWT_case)[which(
	!is.na(res_DKOvsWT_case$padj)
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

ego_DKOvsWT <- enrichGO(
	gene = gene_list_DKOvsWT_case,
	universe = universe_DKOvsWT,
	OrgDb = org.Mm.eg.db,
	ont = "BP",
	keyType = "ENSEMBL",
	readable = TRUE
)

head(ego_DKOvsWT, 3)

# ego's bar plots

png(
	filename = "./results/R_plots/barplot_enrichGO_WT.png",
	width = 3000,
	height = 2400,
	res = 300
)
barplot(ego_WT, showCategory = 10) +
	ggtitle(
		"Enriched terms bar plot",
		"In WT: in case vs control comparison"
	)
dev.off()


png(
	filename = "./results/R_plots/barplot_enrichGO_DKO.png",
	width = 3000,
	height = 2400,
	res = 300
)
barplot(ego_DKO, showCategory = 8) +
	ggtitle("Enriched terms bar plot", "In DKO: case vs control comparison")
dev.off()


png(
	filename = "./results/R_plots/barplot_enrichGO_DKOvsWT.png",
	width = 3000,
	height = 2400,
	res = 300
)
barplot(ego_DKOvsWT, showCategory = 10) +
	ggtitle(
		"Enriched terms bar plot",
		"In Control condition: DKO vs WT comparison"
	)
dev.off()
