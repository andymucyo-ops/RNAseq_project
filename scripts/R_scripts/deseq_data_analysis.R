library(DESeq2)

#-------------------------------------------------------------------------------
# data formating before analysis
#-------------------------------------------------------------------------------

# seting file path
sample_counts_path <- "/Users/mucyo/R-projects/RNAseq_project/data/count.txt"

# save the feature count file as a table
countData <- read.table(sample_counts_path, header = TRUE)

# changing the count into gene id
rownames(countData) <- countData$Geneid
#remove the gene id column
countData <- countData[, -c(1)]

# create a colData data frame containing the data information according to README in reads_Blood
colData <- data.frame(
  Sample = c("SRR7821949", "SRR7821950", "SRR7821951", "SRR7821952", "SRR7821953",
             "SRR7821968", "SRR7821969", "SRR7821970",
             "SRR7821954", "SRR7821955", "SRR7821956", "SRR7821957",
             "SRR7821971", "SRR7821972", "SRR7821973"),
  Genotype = c(rep("WT", 8), rep("DKO", 7)),
  Condition = c(rep("Case", 5), rep("Control", 3),
                rep("Case", 4), rep("Control", 3)),
  row.names = NULL
)

# converting sample into factor to sort them
colData$Sample <- factor(colData$Sample, levels = sort(unique(colData$Sample)))

# sorting data per sample
sorted_colData <- colData[order(colData$Sample),]

#-------------------------------------------------------------------------------
# start Differential gene expression analysis
#-------------------------------------------------------------------------------

# create DESeq2 dataset from the feature counts data and the sample information 
dds <- DESeqDataSetFromMatrix(
  countData,
  sorted_colData,
  design = ~ Genotype + Condition + Genotype:Condition
)

# set reference for Conditions and Genotype 
dds$Genotype <- relevel(dds$Genotype, ref = "WT")

dds$Condition <- relevel(dds$Condition, ref = "Control")

# differential expression analysis on the imported data 
dds <- DESeq(dds)


#-------------------------------------------------------------------------------
# results and data visualization
#-------------------------------------------------------------------------------

# extraction of transformed values, and principal component plot of the samples 
vsd <- vst(dds, blind=FALSE)

plotPCA(vsd, intgroup = c("Genotype","Condition"))

# result output 

dds$Group <- factor(paste0(dds$Genotype, dds$Condition))

design(dds) <- ~ Group

resultsNames(dds)

# get results case vs control for WT and DKO
res_WT <- results(dds, name = "Condition_Case_vs_Control")

res_DKO <- results(dds,
                   list( c("Condition_Case_vs_Control","GenotypeDKO.ConditionCase") )
                  )

# transform into data frame 
res_WT_df <- as.data.frame(res_WT)

res_DKO_df <- as.data.frame(res_DKO)

#remove NA values from data frame 
res_WT_df <- res_WT_df[!is.na(res_WT_df$padj),]

res_DKO_df <- res_DKO_df[!is.na(res_DKO_df$padj),]
