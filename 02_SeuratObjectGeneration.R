##02_Seurat Object Generation

library(Seurat)

##Read in Replicates------------------------------------------------------------
#SPF
AdultMale_SPF <- readRDS("filepath")
AdultFemale_SPF <- readRDS("filepath")
WeaningMale_SPF <- readRDS("filepath")
WeaningFemale_SPF <- readRDS("filepath")
#GF
AdultMale_GF <- readRDS("filepath")
AdultFemale_GF <- readRDS("filepath")
WeaningMale_GF <- readRDS("filepath")
WeaningFemale_GF <- readRDS("filepath")
#GF + FMT
AdultMale_GF_F <- readRDS("filepath")
AdultFemale_GF_F <- readRDS("filepath")
WeaningMale_GF_F <- readRDS("filepath")
WeaningFemale_GF_F <- readRDS("filepath")

##Set up metadata---------------------------------------------------------------
#Experiment
AdultMale_SPF$Experiment <- "SPF"
AdultFemale_SPF$Experiment <- "SPF"
WeaningMale_SPF$Experiment <- "SPF"
WeaningFemale_SPF$Experiment <- "SPF"

AdultMale_GF$Experiment <- "GF"
AdultFemale_GF$Experiment <- "GF"
WeaningMale_GF$Experiment <- "GF"
WeaningFemale_GF$Experiment <- "GF"

AdultMale_GF_F$Experiment <- "GF_FMT"
AdultFemale_GF_F$Experiment <- "GF_FMT"
WeaningMale_GF_F$Experiment <- "GF_FMT"
WeaningFemale_GF_F$Experiment <- "GF_FMT"

AdultMale_SPF$barcode <- paste0("AdultMale", "_", AdultMale_SPF$demux_assign)
AdultMale_SPF$barcode <- gsub("Ms.Hashtag_1", "Rep1", AdultMale_SPF$barcode)
AdultMale_SPF$barcode <- gsub("Ms.Hashtag_2", "Rep2", AdultMale_SPF$barcode)
AdultMale_SPF$barcode <- gsub("Ms.Hashtag_3", "Rep3", AdultMale_SPF$barcode)

AdultFemale_SPF$barcode <- paste0("AdultFemale", "_", AdultFemale_SPF$demux_assign)
AdultFemale_SPF$barcode <- gsub("Ms.Hashtag_1", "Rep1", AdultFemale_SPF$barcode)
AdultFemale_SPF$barcode <- gsub("Ms.Hashtag_2", "Rep2", AdultFemale_SPF$barcode)
AdultFemale_SPF$barcode <- gsub("Ms.Hashtag_3", "Rep3", AdultFemale_SPF$barcode)

WeaningMale_SPF$barcode <- paste0("WeaningMale", "_", WeaningMale_SPF$demux_assign)
WeaningMale_SPF$barcode <- gsub("Ms.Hashtag_1", "Rep1", WeaningMale_SPF$barcode)
WeaningMale_SPF$barcode <- gsub("Ms.Hashtag_2", "Rep2", WeaningMale_SPF$barcode)
WeaningMale_SPF$barcode <- gsub("Ms.Hashtag_3", "Rep3", WeaningMale_SPF$barcode)

WeaningFemale_SPF$barcode <- paste0("WeaningFemale", "_", WeaningFemale_SPF$demux_assign)
WeaningFemale_SPF$barcode <- gsub("Ms.Hashtag_1", "Rep1", WeaningFemale_SPF$barcode)
WeaningFemale_SPF$barcode <- gsub("Ms.Hashtag_2", "Rep2", WeaningFemale_SPF$barcode)
WeaningFemale_SPF$barcode <- gsub("Ms.Hashtag_3", "Rep3", WeaningFemale_SPF$barcode)

AdultMale_GF$barcode <- paste0("AdultMale", "_", AdultMale_GF$demux_assign)
AdultMale_GF$barcode <- gsub("Ms.Hashtag_1", "Rep1", AdultMale_GF$barcode)
AdultMale_GF$barcode <- gsub("Ms.Hashtag_2", "Rep2", AdultMale_GF$barcode)
AdultMale_GF$barcode <- gsub("Ms.Hashtag_3", "Rep3", AdultMale_GF$barcode)

AdultFemale_GF$barcode <- paste0("AdultFemale", "_", AdultFemale_GF$demux_assign)
AdultFemale_GF$barcode <- gsub("Ms.Hashtag_1", "Rep1", AdultFemale_GF$barcode)
AdultFemale_GF$barcode <- gsub("Ms.Hashtag_2", "Rep2", AdultFemale_GF$barcode)
AdultFemale_GF$barcode <- gsub("Ms.Hashtag_3", "Rep3", AdultFemale_GF$barcode)

WeaningMale_GF$barcode <- paste0("WeaningMale", "_", WeaningMale_GF$demux_assign)
WeaningMale_GF$barcode <- gsub("Ms.Hashtag_1", "Rep1", WeaningMale_GF$barcode)
WeaningMale_GF$barcode <- gsub("Ms.Hashtag_2", "Rep2", WeaningMale_GF$barcode)
WeaningMale_GF$barcode <- gsub("Ms.Hashtag_3", "Rep3", WeaningMale_GF$barcode)

WeaningFemale_GF$barcode <- paste0("WeaningFemale", "_", WeaningFemale_GF$demux_assign)
WeaningFemale_GF$barcode <- gsub("Ms.Hashtag_1", "Rep1", WeaningFemale_GF$barcode)
WeaningFemale_GF$barcode <- gsub("Ms.Hashtag_2", "Rep2", WeaningFemale_GF$barcode)
WeaningFemale_GF$barcode <- gsub("Ms.Hashtag_3", "Rep3", WeaningFemale_GF$barcode)

AdultMale_GF_F$barcode <- paste0("AdultMale", "_", AdultMale_GF_F$demux_assign)
AdultMale_GF_F$barcode <- gsub("Ms.Hashtag_1", "Rep1", AdultMale_GF_F$barcode)
AdultMale_GF_F$barcode <- gsub("Ms.Hashtag_2", "Rep2", AdultMale_GF_F$barcode)
AdultMale_GF_F$barcode <- gsub("Ms.Hashtag_3", "Rep3", AdultMale_GF_F$barcode)

AdultFemale_GF_F$barcode <- paste0("AdultFemale", "_", AdultFemale_GF_F$demux_assign)
AdultFemale_GF_F$barcode <- gsub("Ms.Hashtag_1", "Rep1", AdultFemale_GF_F$barcode)
AdultFemale_GF_F$barcode <- gsub("Ms.Hashtag_2", "Rep2", AdultFemale_GF_F$barcode)
AdultFemale_GF_F$barcode <- gsub("Ms.Hashtag_3", "Rep3", AdultFemale_GF_F$barcode)

WeaningMale_GF_F$barcode <- paste0("WeaningMale", "_", WeaningMale_GF_F$demux_assign)
WeaningMale_GF_F$barcode <- gsub("Ms.Hashtag_1", "Rep1", WeaningMale_GF_F$barcode)
WeaningMale_GF_F$barcode <- gsub("Ms.Hashtag_2", "Rep2", WeaningMale_GF_F$barcode)

WeaningFemale_GF_F$barcode <- paste0("WeaningFemale", "_", WeaningFemale_GF_F$demux_assign)
WeaningFemale_GF_F$barcode <- gsub("Ms.Hashtag_1", "Rep1", WeaningFemale_GF_F$barcode)
WeaningFemale_GF_F$barcode <- gsub("Ms.Hashtag_2", "Rep2", WeaningFemale_GF_F$barcode)
WeaningFemale_GF_F$barcode <- gsub("Ms.Hashtag_3", "Rep3", WeaningFemale_GF_F$barcode)

##Merge Replicates/Datasets-----------------------------------------------------
TEC <- merge(AdultMale_SPF, y = list(AdultFemale_SPF, WeaningMale_SPF, WeaningFemale_SPF,
                                     AdultMale_GF, AdultFemale_GF, WeaningMale_GF, WeaningFemale_GF,
                                     AdultMale_GF_F, AdultFemale_GF_F, WeaningMale_GF_F, WeaningFemale_GF_F), 
             add.cell.ids = c("AdultMale_SPF", "AdultFemale_SPF", "WeaningMale_SPF", "WeaningFemale_SPF",
                              "AdultMale_GF", "AdultFemale_GF", "WeaningMale_GF", "WeaningFemale_GF",
                              "AdultMale_GF_F", "AdultFemale_GF_F", "WeaningMale_GF_F", "WeaningFemale_GF_F"), project = "TEC")

##Add additional metadata-------------------------------------------------------
TEC$age <- ifelse(grepl("Adult", TEC$barcode), "Adult", "Weaning")
TEC$sex <- ifelse(grepl("Male", TEC$barcode), "Male", "Female")
TEC$dataset <- "TEC"

#Analyze without integration----------------------------------------------------
TEC <- NormalizeData(TEC, normalization.method = "LogNormalize", scale.factor = 10000)
TEC <- FindVariableFeatures(TEC, selection.method = "vst", nfeatures = 10000)
TEC <- ScaleData(TEC, features = VariableFeatures(TEC))
TEC <- RunPCA(TEC)
ElbowPlot(TEC)
TEC <- FindNeighbors(TEC, dims = 20)
TEC <- FindClusters(TEC, resolution = 0.5, cluster.name = "unintegrated_clusters")
TEC <- RunUMAP(TEC, dims = 1:20, reduction = "pca", reduction.name = "umap.unintegrated")

DimPlot(TEC, reduction = "umap.unintegrated", group.by = "Experiment")

#Perform Integration (method = Harmony)-----------------------------------------
TEC <- IntegrateLayers(object = TEC, method = HarmonyIntegration,
                       orig.reduction = "pca", new.reduction = "harmony",
                       verbose = FALSE)

##Visualize different integrations and compare to unintegrated umaps
TEC <- FindNeighbors(TEC, reduction = "harmony", dims = 1:20)
TEC <- FindClusters(TEC, resolution = 0.5)
TEC <- RunUMAP(TEC, dims = 1:20, reduction = "harmony")
p1 <- DimPlot(TEC, reduction = "umap", group.by = "Experiment", raster = FALSE)
p2 <- DimPlot(TEC, reduction = "umap", group.by = "age", raster = FALSE)
plot <- wrap_plots(c(p1,p2), ncol = 2, byrow = F)
plot + plot_annotation(title = "Harmony Integration")


##Join RNA layers---------------------------------------------------------------
TEC[["RNA"]] <- JoinLayers(TEC[["RNA"]])

##Remove unwanted cells (CD45+)--------------------------------------------------
DimPlot(TEC, reduction = "umap", label = TRUE, raster = FALSE)
FeaturePlot(TEC, features = c("Epcam", "Ptprc", "Cd3e", "Psmb11", "Prss16", "H2-Aa"), reduction = "umap", raster = FALSE)
##Endothelial
FeaturePlot(TEC, features = c("Pecam1", "Kdr", "Cldn5"), reduction = "umap", raster = FALSE)
##Fibroblasts
FeaturePlot(TEC, features = c("Pdgfra", "Upk3b", "Lrrn4"), reduction = "umap", raster = FALSE)



Idents(TEC) <- TEC$seurat_clusters
VlnPlot(TEC, features = c("percent.mt", "nCount_RNA", "nFeature_RNA"), raster = FALSE)

#Filtered cells: 12-TCells, 20-DC/MACS, 19-Fibroblasts/Endothelial, 13 is mix of nurse cTEC and T cells, 10- low quality cluster)

TEC <- subset(TEC, idents = setdiff(levels(TEC), c("12", "20", "19", "10")))

TEC <- FindNeighbors(TEC, reduction = "harmony", dims = 1:20)
TEC <- FindClusters(TEC, resolution = 0.5)
TEC <- RunUMAP(TEC, dims = 1:20, reduction = "harmony")
DimPlot(TEC, reduction = "umap", label = TRUE, raster = FALSE)


##Add cluster annotations and analyze-------------------------------------------
##Identify Clusters and label
#cTEC
features <- c("Prss16", "Psmb11", "Enpep", "Ly75", "Cxcl12", "Pax1", "Cd3e", "Rag1", "Ccl21a")
VlnPlot(TEC, features = features, stack = T)
##15 = Nurse_cTEC/Cycling_cTEC
##0,5,13 = cTEC
##12 = Intertypical TEC

##mTEC
features <- c("Cxcl12","Ccl21a", "Cd80", "Cd52", "Aire", "Top2a", "Mki67", "Spink5", "Pigr", "Plet1", "Avil", "Pou2f3")
VlnPlot(TEC, features = features, stack = T)
##1,2,3,6= mTECI
##8,9,10,4 = mTECII
##7 = mTECIII
##11 = mTECIV

##Specialized
features <- c("Gp2", "Spib", "Ccl20", "Snap25", "Stxbp5l", "Dnah12", "Myog", "Myl1", "Col1a1", "Dcn", "Fbn1")
VlnPlot(TEC, features = features, stack = T)
##18 = Ciliated
##17 = Neuroendocrine
##16 = Muscle-like / Myoid
##14 = Microfold

##Set labels
labels <- c("0" = "cTEC-2", "1" = "mTECI-1", "2" = "mTECI-2", "3" = "mTECI-4", "4" = "mTECII-4", "5" = "cTEC-1", "6" = "mTECI-3", "7" = "mTECIII", "8" = "mTECII-1", "9" = "mTECII-2", "10" = "mTECII-3",
            "11" = "mTECIV", "12" = "Intertypical_TEC", "13" = "cTEC-3", "14" = "Microfold_TEC", "15" = "Nurse_cTEC", "16" = "Muscle-like_TEC", "17" = "Neuroendocrine_TEC", "18" = "Ciliated_TEC")
TEC <- RenameIdents(TEC, labels)

##Annotate meta clusters
Idents(TEC) <- TEC$seurat_clusters
meta_clusters <- c("0" = "cTEC", "1" = "mTECI", "2" = "mTECI", "3" = "mTECI", "4" = "mTECII", "5" = "cTEC", "6" = "mTECI", "7" = "mTECIII", "8" = "mTECII", "9" = "mTECII", "10" = "mTECII",
                   "11" = "mTECIV", "12" = "Intertypical_TEC", "13" = "cTEC", "14" = "Microfold_TEC", "15" = "Nurse_cTEC", "16" = "Muscle-like_TEC", "17" = "Neuroendocrine_TEC", "18" = "Ciliated_TEC")
TEC <- RenameIdents(TEC, meta_clusters)
TEC$meta_clusters <- Idents(TEC)


#Cell Cycle Phase
cc.genes <- Seurat::cc.genes
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

TEC <- CellCycleScoring(
  TEC,
  s.features = s.genes,
  g2m.features = g2m.genes,
  set.ident = TRUE
)
DimPlot(TEC, group.by = "Phase", reduction = "umap")
FeaturePlot(TEC, features = c("S.Score", "G2M.Score"), reduction = "umap")

cluster_cols <- c(
  "#003F5C",   #ncTEC
  "#08306B",   #cTEC
  "#5DADE2",   #cTEC
  "#9ECAE1",   #cTEC
  "#5B6770",   #Intertypical 
  "#B22222",   #mTECI 
  "#D95F02",   #mTECI  
  "#F28E2B",   #mTECI  
  "#FDBF6F",   #mTECI  
  "#238B45",   #mTECII 
  "#66C2A4",   #mTECII
  "#A1D99B",   #mTECII 
  "#E5F5E0",   #mTECII 
  "#7570B3",   #mTECIII
  "#E7298A",   #mTECIV
  "#5C4033",   #Other
  "#8C6D31",   #Other
  "#A89F91",   #Other
  "#C2B280"    #Other
)

cluster_cols <- c(
  "#003F5C",
  "#2171B5",
  "#5B6770",
  "#E8750A",
  "#41A85F",
  "#7570B3",
  "#E7298A",
  "#5C4033",   
  "#8C6D31",   
  "#A89F91",  
  "#C2B280")

