##03_Pseudobulk

library(Seurat)
library(DESeq2)

TEC <- readRDS("SeuratObjectFile")

##Add extra metadata to set conditions later
TEC$barcode2 <- paste0(TEC$Experiment, "_", TEC$barcode)
TEC$barcode2 <- paste0(TEC$labels, "_", TEC$barcode)

##Generate Counts Matrix - sample level---------------------------------------
##aggregate counts to sample 
Counts <- AggregateExpression(TEC,
                              group.by = c("labels", "barcode2"),
                              assays = "RNA",
                              return.seurat = FALSE)
Counts <- AggregateExpression(TEC,
                              group.by = c("meta_clusters", "barcode2"),
                              assays = "RNA",
                              return.seurat = FALSE)


Counts <- Counts$RNA
##transpose
Counts.t <- t(Counts)

##Convert to dataframe
Counts.t <- as.data.frame(Counts.t)

##Split rows to get cell names
splitRows <- gsub('_.*', '', rownames(Counts.t))

##Split 
Counts.split <- split.data.frame(Counts.t, f = factor(splitRows))

##Fix colnames and re-transpose
Counts.split.mod <- lapply(Counts.split, function(x){
  rownames(x) <- gsub('.*_(.*)', '\\1', rownames(x))
  t(x)
})

##Get Counts matrix from specific cluster---------------------------------------
Counts_cTEC <- Counts.split.mod$cTEC


X <- Counts_cTEC
##Generate sample level metadata------------------------------------------------
colData <- data.frame(samples = colnames(X))
#condition - age
colData <- colData %>%
  mutate(condition = ifelse(grepl('Adult', samples), 'Adult', 'Weaning')) %>%
  column_to_rownames(var = 'samples')
#condition - microbial status
colData <- colData %>%
  mutate(condition = ifelse(grepl('SPF', samples), 'SPF',
                            ifelse(grepl('GF-FMT', samples), 'GF_FMT', 'GF'))) %>%
  column_to_rownames(var = 'samples')

#Subset Analysis
colData <- (column_to_rownames(colData, var = 'samples'))
colData$age <- ifelse(grepl('Adult', rownames(colData)), 'Adult', 'Weaning')
colData$label <- gsub("([^-]*-[^-]*).*", "\\1", rownames(colData))
colData$sex <- ifelse(grepl('Male', rownames(colData)), 'Male', 'Female')
colData$samples <- gsub("^([^-]*-[^-]*)-", "", rownames(colData)) 


##Create comparison design------------------------------------------------------
colData$condition <- factor(colData$condition, levels = c("Weaning", "Adult"))
colData$condition <- factor(colData$condition, levels = c("SPF", "GF", "GF_FMT"))
colData$contrast1 <- factor(colData$condition, levels = c("Other", "cTEC-1"))


##Filter Matrix for gene counts >= 10 in at least 5 samples---------------------
keep <- rowSums(X >= 10) >= 5
X_filt <- X[keep,]

##Create DESeq2 object----------------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = X_filt,
                              colData = colData,
                              design = ~ condition)


##Run DESeq2--------------------------------------------------------------------
dds <- DESeq(dds)

##Extract results---------------------------------------------------------------
res <- results(dds, contrast = c("condition", "SPF", "GF")) #SPFvsGF
res <- results(dds, contrast = c("condition", "SPF", "GF_FMT")) #SPFvsGF_FMT
res <- results(dds, contrast = c("condition", "GF", "GF_FMT")) #GFvsGF_FMT
res <- results(dds, name = "condition_Weaning_vs_Adult") ##WeaningvsAdult
res <- results(dds, name = "contrast1_cTEC.1_vs_Other") ##cTEC subsets 