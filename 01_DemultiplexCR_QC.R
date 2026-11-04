##01_Demultiplexing Cell Ranger Filtered Matrix Output

library(Seurat)
library(deMULTIplex2)

##Read in counts data
counts <- Read10X("path_to_10Xdata")

##Extract HTO and GEX 
hto_counts <- counts$`Antibody Capture`
gex_counts <- counts$`Gene Expression`

##Create Seurat Object and CLR Normalize 
SO <- CreateSeuratObject(counts = gex_counts)
SO[["HTO"]] <- CreateAssayObject(counts = hto_counts)

SO <- NormalizeData(SO, 
                    assay = "HTO", 
                    normalization.method = "CLR",
                    margin = 2)

RidgePlot(SO, assay = "HTO", 
          features = rownames(SO[["HTO"]]))

hto_df <- as.data.frame(t(hto_counts))

##Demultiplex
res <- demultiplexTags(as.matrix(hto_df),
                       init.cos.cut = c(), 
                       prob.cut = 0.5,
                       plot.diagnostics = TRUE,
                       plot.name = "demux")

##Dataset Parameters
#SPF----------------------------------------------------------------------------
#WeaningFemale > init.cos.cut = 0.5
#WeaningMale > init.cos.cut = c(0.5,0.8,0.5)
#AdultMale/AdultFemale > init.cos.cut = c(0.5, 0.85, 0.5)
#GermFree-----------------------------------------------------------------------
#Adult Male > init.cos.cut = c(0.5, 0.65, 0.5), prob.cut = 0.1
#Adult Female > init.cos.cut = c(0.5, 0.6, 0.5), prob.cut = 0.1
#Weaning Male, Weaning Female > init.cos.cut = 0.5
##GermFree + FMT----------------------------------------------------------------
#Adult Female > init.cos.cut = c(0.1, 0.99, 0.1)
#Adult Male > init.cos.cut = c(0.5, 0.75, 0.5), prob.cut = 0.1
#Weaning Female > init.cos.cut = c(0.1, 0.95, 0.5)
#Weaning Male > init.cos.cut = c(0.5, 0.9)

##Inspect assignments 
table(res$final_assign)

##Inspect UMAP
umap_df <- as.data.frame(res$umap)
umap_df$hashtag1_prob <- res$prob_mtx[,"Ms.Hashtag_1"]
umap_df$hashtag2_prob <- res$prob_mtx[,"Ms.Hashtag_2"]
umap_df$hashtag3_prob <- res$prob_mtx[,"Ms.Hashtag_3"]

ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = hashtag1_prob)) +
  geom_point(size = 0.5, alpha = 0.5) +
  theme_classic() 

##Add assignments to seurat object
SO$demux_assign <- "negative"
matching_barcodes <- intersect(colnames(SO), names(res$final_assign))
SO$demux_assign[matching_barcodes] <- res$final_assign[matching_barcodes]

table(SO$demux_assign)

##Drop multiplets/missing cells
SO <- subset(SO,
             subset = demux_assign %in% c("Ms.Hashtag_1",
                                          "Ms.Hashtag_2",
                                          "Ms.Hashtag_3"))

table(SO$demux_assign)

###Filter by QC thresholds
##QC SO
SO[["percent.mt"]] <- PercentageFeatureSet(SO, pattern = "^mt-")

VlnPlot(SO, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        group.by = "demux_assign")
table(SO$demux_assign)

SO <- subset(SO, 
             subset = nFeature_RNA > 500 & 
               nFeature_RNA < 10000 & 
               nCount_RNA > 1000 &
               percent.mt < 10)

saveRDS(SO, "Filename")

