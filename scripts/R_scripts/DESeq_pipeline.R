#' ---
#' editor_options: 
#'   markdown: 
#'     wrap: 72
#' ---
#' 
#' install all required packages for analysis
#' 
#' note: use github version of enhance volcano if not able to install using
#' BiocManager otherwise can be installed using:
#' `BiocManager::install("EnhancedVolcano")`
#' 
## ---------------------------------------------------------------------------------------------
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
devtools::install_github("kevinblighe/EnhancedVolcano")

#load libraries 
library("DESeq2")
library("pheatmap")
library("ggrepel")
library("ggplot2")
library("EnhancedVolcano")

#' 
#' # 1) Data formatting before analysis
#' 
#' set file path
#' 
## ---------------------------------------------------------------------------------------------
sample_counts_path <- "~/R-projects/RNAseq_project/data/count.txt"

#' 
#' save the feature count file as a table
#' 
## ---------------------------------------------------------------------------------------------
countdata <- read.table(sample_counts_path, header = TRUE)
head(countdata)

#' 
#' changing the count into gene id
#' 
## ---------------------------------------------------------------------------------------------
rownames(countdata) <- countdata$Geneid 
countdata <- countdata[, -c(1)] #remove the gene id column
head(countdata)

#' 
#' Create a coldata data frame containing the data information according to
#' README in reads_Blood, acting as the metadata of the counts
#' 
## ---------------------------------------------------------------------------------------------
coldata <- data.frame( Sample = c("SRR7821949", "SRR7821950",
"SRR7821951", "SRR7821952", "SRR7821953", "SRR7821968", "SRR7821969",
"SRR7821970", "SRR7821954", "SRR7821955", "SRR7821956", "SRR7821957",
"SRR7821971", "SRR7821972", "SRR7821973"), Genotype = c(rep("WT", 8),
rep("DKO", 7)), Condition = c(rep("Case", 5), rep("Control", 3),
rep("Case", 4), rep("Control", 3)), row.names = NULL )

head(coldata, n=10L)

#' 
#' converting sample into factor to match DESeq2 input format, and order by
#' the sample names to match column names of the count data table
#' 
## ---------------------------------------------------------------------------------------------
coldata$Sample <- factor(coldata$Sample, levels =
sort(unique(coldata$Sample)))

# saving the sorted meta data as a new data frame
sorted_coldata <- coldata[order(coldata$Sample),]

#setting sample names as rownames
rownames(sorted_coldata) <- sorted_coldata$Sample

#remonving sample names column
sorted_coldata <-sorted_coldata[,-c(1)]

head(sorted_coldata, n=10L)

#final check before launching DESeq
all(rownames(sorted_coldata) == colnames(countdata))

#' 
#' # 2) start Differential gene expression analysis
#' 
#' create DESeq data set from the feature counts data and the counts meta
#' data
#' 
## ---------------------------------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix( countdata, sorted_coldata, design = ~
Genotype + Condition + Genotype:Condition )

#transform the meta data values into factor, to match DESeq expected data types, and set references for downstream data analysis 
dds$Genotype <- relevel(dds$Genotype, ref = "WT")

dds$Condition <- relevel(dds$Condition, ref = "Control")

#' 
#' then run differential expression analysis
#' 
## ---------------------------------------------------------------------------------------------
dds <- DESeq(dds)

#' 
#' # 3) Results and data visualization
#' 
#' extraction of transformed values, and principal component plot of the
#' samples
#' 
## ---------------------------------------------------------------------------------------------
vsd <- vst(dds, blind=FALSE)

plotPCA(vsd, intgroup = c("Genotype","Condition"))

#' 
#' get results coefficient to highlight in the following result step
#' 
## ---------------------------------------------------------------------------------------------
resultsNames(dds)

#' 
#' get results comparation of case vs control form both genotype (WT and
#' DKO)
#' 
## ---------------------------------------------------------------------------------------------
res_WT <- results(dds, name = "Condition_Case_vs_Control", alpha = 0.05)

res_DKO <- results(dds, list(
c("Condition_Case_vs_Control","GenotypeDKO.ConditionCase") ), alpha = 0.05 )

summary_WT <- summary(res_WT)
  
summary_DKO <- summary(res_DKO)

#' 
#' transform into data frame, and remove NA values
#' 
## ---------------------------------------------------------------------------------------------
#transform into dataframe
res_WT_df <- as.data.frame(res_WT)

res_DKO_df <- as.data.frame(res_DKO)


#remove NA values from data frame 
res_WT_df <- res_WT_df[!is.na(res_WT_df$padj),]

res_DKO_df <- res_DKO_df[!is.na(res_DKO_df$padj),]

#' 
#' heatmap of the results
#' 
## ---------------------------------------------------------------------------------------------
select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:50] 
df <- as.data.frame(colData(dds)[,c("Genotype","Condition")])
pheatmap(assay(vsd)[select,],cluster_rows = FALSE, cluster_cols = FALSE, annotation_col=df)

# heatmap clustering by similarity
## ---------------------------------------------------------------------------------------------
library("RColorBrewer")
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- rownames(sorted_coldata)
colnames(sampleDistMatrix) <- rownames(sorted_coldata)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

#' 
#' volcano plot for each result set
#' 
## ---------------------------------------------------------------------------------------------
EnhancedVolcano(
  res_WT_df,
  lab = rownames(res_WT_df),
  x = 'log2FoldChange',
  y = 'pvalue',
  title = 'WT case vs control'
)

#' 
## ---------------------------------------------------------------------------------------------
EnhancedVolcano(
  res_DKO_df,
  lab = rownames(res_DKO_df),
  x = 'log2FoldChange',
  y = 'pvalue',
  title = 'DKO case vs control'
)

