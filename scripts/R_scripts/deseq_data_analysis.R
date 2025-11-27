library(DESeq2)

# seting file path
sample_counts_path <- "/Users/mucyo/R-projects/RNAseq_project/data/count.txt"


# save the feature count file as a table
countData <- read.table(sample_counts_path, header = TRUE)

# format the data input to match expected DESeq format
rownames(countData) <- countData$Geneid

countData <- countData[, -c(1)]

# create a colData dataframe containing the data informations according to README in reads_Blood

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
sorted_colData <- colData[order(colData$Sample), ]


# create DESeq2 dataset from the feature counts data and the sample information 
dds <- DESeqDataSetFromMatrix(
  countData,
  sorted_colData,
  design = ~ Genotype + Condition
)

# differential expression analysis on the imported data 
dds <- DESeq(dds)


# extraction of transformed values, and principal component plot of the samples 
vsd <- vst(dds, blind=FALSE)

plotPCA(vsd, intgroup = c("Genotype","Condition"))

