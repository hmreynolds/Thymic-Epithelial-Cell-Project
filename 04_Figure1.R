##04_Figure1

library(Seurat)
##1B-----------------------------------------------------------------------------
TEC <- readRDS("SeuratObjectfile.rds")

cluster_cols <- c(
  "#003F5C",   #n cTEC
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

##UMAP--------------------------------------------------------------------------
Idents(TEC) <- TEC$labels_num

p1 <- DimPlot(TEC, reduction = "umap", label = T, raster = F, 
              cols = cluster_cols, label.box = T)
p1 <- p1 +  theme(
  axis.title = element_text(size = 12),
  axis.text = element_text(size = 12),
  legend.title = element_blank(),
  legend.text = element_blank(),
  legend.key.size = unit(0.3, "cm"),
  legend.spacing.y = unit(0.05, "cm"),
  plot.margin = margin(5.5, 80, 5.5, 5.5)
)
p1


Idents(TEC) <- TEC$labelsandnum
my_cluster_order <- c("1:Nurse_cTEC", "2:cTEC-1", "3:cTEC-2", "4:cTEC-3", "5:Intertypical_TEC", "6:mTECI-1", "7:mTECI-2", "8:mTECI-3", "9:mTECI-4", 
                      "10:mTECII-1", "11:mTECII-2", "12:mTECII-3", "13:mTECII-4", "14:mTECIII", "15:mTECIV", "16:Muscle-like_TEC", "17:Microfold_TEC", "18:Neuroendocrine_TEC", "19:Ciliated_TEC")
Idents(TEC) <- factor(Idents(TEC), levels = my_cluster_order)
p1_1 <- DimPlot(TEC, reduction = "umap", label = T, raster = F, 
                cols = cluster_cols, label.box = T)
p1_1 <- p1_1 + theme(legend.text = element_text(size = 10),
                     axis.title = element_blank(),
                     legend.key.size = unit(0.3, "cm"),
                     legend.spacing.y = unit(0.05, "cm"))

##FeaturePlot-------------------------------------------------------------------
Idents(TEC) <- TEC$labels
my_cluster_order <- c("Nurse_cTEC", "cTEC-1", "cTEC-2", "cTEC-3", "Intertypical_TEC", "mTECI-1", "mTECI-2", "mTECI-3", "mTECI-4", 
                      "mTECII-1", "mTECII-2", "mTECII-3", "mTECII-4", "mTECIII", "mTECIV", "Muscle-like_TEC", "Microfold_TEC", "Neuroendocrine_TEC", "Ciliated_TEC")
Idents(TEC) <- factor(Idents(TEC), levels = my_cluster_order)
features <- c("Foxn1", "Epcam", "Krt8", "Krt14", "Krt5", "Cd3e", "Psmb11", "Enpep", "Cxcl12", "Ccl21a", 
              "Cd40", "Cd80", "Aire", "Cd52", "Spink5", "Avil", "Dclk1", "Myog", "Myl1", "Spib", "Gp2", "Snap25", "Stxbp5l", "Pigr", "Dnah12")
VlnPlot(TEC, features = features, stack = T) 
p2 <- DotPlot(TEC, features = features) + 
  coord_flip() +
  scale_color_gradient2(low = "darkblue", mid = "white", high = "darkred", midpoint = 0) +
  geom_point(
    aes(size = pct.exp),
    shape = 21, stroke = 0.5, color = "black") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "right",        # or "bottom", "top", "left", or c(x, y) e.g. c(1.1, 0.5)
    legend.title = element_text(size = 10, hjust = 0.5, ),
    legend.box = "vertical"
  )
p2

##ClusterFrequency---------------------------------------------------------------
TEC <- readRDS("/Volumes/Bettini_Lab/DATA/Hayley/Hayley/TECscRNAseq/10XTEC_Analysis/01_CRtoSO_v10/07_SPFAnalysis/TEC_SPFsubset.rds")
##All Clusters
my_cluster_order <- c("Nurse_cTEC", "cTEC-1", "cTEC-2", "cTEC-3", "Intertypical_TEC", "mTECI-1", "mTECI-2", "mTECI-3", "mTECI-4", 
                      "mTECII-1", "mTECII-2", "mTECII-3", "mTECII-4", "mTECIII", "mTECIV", "Muscle-like_TEC", "Microfold_TEC", "Neuroendocrine_TEC", "Ciliated_TEC")
ClusterFreq <- table(TEC$labels, TEC$age)
ClusterFreq <- as.data.frame(ClusterFreq)
ClusterFreq$Var1 <- as.character(ClusterFreq$Var1) 
ClusterFreq$Var1 <- factor(ClusterFreq$Var1, levels = my_cluster_order)
ClusterFreq$Var2 <- factor(ClusterFreq$Var2, levels = c("Weaning", "Adult"))

p3 <- ggplot(ClusterFreq, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 1.0, color = "black") +
  scale_fill_manual(values = cluster_cols) +
  facet_grid(cols = vars(Var2), scales = "free",
             labeller = labeller(
               Var3 = c(
                 Weaning = "Weaning",
                 Adult = "Adult"
               )
             )) +
  theme(axis.title.x = element_blank(), 
        panel.grid = element_blank(), 
        legend.title = element_blank(), 
        panel.border = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 12, color = "black"),
        strip.text.x = element_text(size =15),
        legend.text = element_text(size = 12),
        legend.key.size = unit(0.5, "cm"),
        legend.spacing.y = unit(0.05, "cm")
  ) + labs(y = "Proportion of TECs") 
p3

##Metaclusters
my_cluster_order <- c("Nurse_cTEC", "cTEC", "Intertypical_TEC", "mTECI", "mTECII", "mTECIII", "mTECIV", "Muscle-like_TEC", "Microfold_TEC", "Neuroendocrine_TEC", "Ciliated_TEC")
ClusterFreq <- table(TEC$meta_clusters, TEC$age, TEC$Experiment)
ClusterFreq <- as.data.frame(ClusterFreq)
ClusterFreq$Var1 <- as.character(ClusterFreq$Var1) 
ClusterFreq$Var1 <- factor(ClusterFreq$Var1, levels = my_cluster_order)
ClusterFreq$Var2 <- factor(ClusterFreq$Var2, levels = c("Weaning", "Adult"))
cluster_cols <- c(
  "#003F5C",
  "#2171B5",
  "#5B6770",
  "#E8750A",
  "#41A85F",
  "#7570B3",
  "#E7298A",
  "#5C4033",   #Other
  "#8C6D31",   #Other
  "#A89F91",   #Other
  "#C2B280")
p4 <- ggplot(ClusterFreq, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 1.0, color = "black") +
  scale_fill_manual(values = cluster_cols) +
  facet_grid(cols = vars(Var2), scales = "free",
             labeller = labeller(
               Var3 = c(
                 Weaning = "Weaning",
                 Adult = "Adult"
               )
             )) +
  theme(axis.title.x = element_blank(), 
        panel.grid = element_blank(), 
        legend.title = element_blank(), 
        panel.border = element_blank(), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.ticks.x = element_blank(),
        strip.text.x = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.key.size = unit(0.5, "cm"),
        legend.spacing.y = unit(0.05, "cm")
  ) + labs(y = "Proportion of TECs") 
p4

##Per-Cell density ratio score---------------------------------------------------

##Calculate a local density ratio for each cell — how enriched is its neighborhood for Dataset1 vs Dataset2 — and store it as metadata for plotting.
umap_df <- as.data.frame(Embeddings(TEC, "umap"))
umapSPF <- umap_df[grepl("SPF", rownames(umap_df)),]
umapSPF_weaning <- umapSPF[grepl("Weaning", rownames(umapSPF)),]
umapSPF_adult <- umapSPF[grepl("Adult", rownames(umapSPF)),]
dataset_label <- TEC$age

library(RANN)
k <- 50
nn <- nn2(umapSPF, k = k +1)$nn.idx
prop_d1 <- apply(nn[, -1], 1, function(neighbors) {
  mean(dataset_label[neighbors] == "Weaning")
})
TEC$density_ratio <- log2((prop_d1 + 0.01) / (1 - prop_d1 + 0.01))
p5 <- FeaturePlot(TEC, features = "density_ratio", reduction = "umap") +
  scale_color_gradient2(low = "#00008B", mid = "lightgrey", high = "#8B0000",
                        midpoint = 0,
                        name = "log2(Weaning/Adult)\ndensity ratio",
                        limits = c(-3,3),
                        oob = scales::squish) +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12)
  )
p5 

##Differential Cluster Proportions----------------------------------------------
library(speckle)
library(limma)
##Plot cluster proportions with sex as variable
Idents(TEC) <- TEC$meta_clusters
TEC$barcode2 <- paste0(TEC$Experiment, "_", TEC$sex, "_", TEC$barcode)
cluster_props <- TEC@meta.data %>%
  count(barcode2, cluster = Idents(TEC)) %>%
  group_by(barcode2) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()

cluster_props$age <- ifelse(grepl("Adult", cluster_props$barcode2), 'Adult', 'Weaning')
cluster_props$dataset <- ifelse(grepl("SPF", cluster_props$barcode2), 'SPF',
                                ifelse(grepl('GF_FMT', cluster_props$barcode2), 'GF_FMT', 'GF'))
cluster_props$sex <- ifelse(grepl("Male", cluster_props$barcode2), 'Male', 'Female')
cluster_props$dataset <- factor(cluster_props$dataset, levels = c("SPF", "GF", "GF_FMT"))
cluster_props$cluster <- factor(cluster_props$cluster, levels = my_cluster_order)
cluster_props$age <- factor(cluster_props$age, levels = c("Weaning", "Adult"))

cluster_props <- cluster_props %>% filter(cluster %in% c("cTEC", "mTECI", "Intertypical_TEC", "mTECII", "mTECIII", "mTECIV"))
cluster_props <- cluster_props %>% filter(dataset %in% c("SPF"))

#Feed Seurat Object into propeller
props <- getTransformedProps(
  clusters = Idents(TEC),
  sample = TEC@meta.data$barcode2,
  transform = "asin"
)

sample_meta <- TEC@meta.data %>%
  distinct(barcode2, dataset, age) %>%
  arrange(match(barcode2, colnames(props$Proportions)))  # match column order

# Build design matrix for the factorial model
design <- model.matrix(~ age, data = sample_meta)

# Fit limma model across clusters
fit <- lmFit(props$TransformedProps, design)
fit <- eBayes(fit)

# Pull out results — e.g., interaction term (dataset:age)
topTable(fit, coef = "ageWeaning", sort.by = "P", number = Inf)

stat_results <- topTable(fit, coef = "ageWeaning", number = Inf) %>%
  tibble::rownames_to_column("cluster") %>%
  mutate(
    p.signif = case_when(
      adj.P.Val < 0.001 ~ "***",
      adj.P.Val < 0.01  ~ "**",
      adj.P.Val < 0.05  ~ "*",
      TRUE ~ "ns"
    ),
    # position label near top of each facet's y-range — adjust per your data
    y.position = NA  # fill in based on max proportion per cluster, see below
  )

# get a sensible y position per cluster (e.g., just above the max proportion)
y_positions <- cluster_props %>%
  group_by(cluster) %>%
  summarise(y.position = max(proportion) * 1.2)

stat_results <- stat_results %>%
  filter(cluster %in% c("cTEC", "mTECI", "Intertypical_TEC", "mTECII", "mTECIII", "mTECIV"))

stat_results <- left_join(stat_results, y_positions, by = "cluster")

p6 <- ggplot(cluster_props, aes(x = dataset, y = proportion, fill = age)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, position = position_dodge(width = 0.75)) +
  geom_jitter(
    aes(color = age, shape = sex),
    size = 2,
    alpha = 0.8,
    position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75)
  ) +
  geom_text(
    data = stat_results,
    aes(x = 1.0, y = (y.position), label = p.signif),
    inherit.aes = FALSE,
    size = 5,
    fontface = "bold"
  ) +
  scale_fill_manual(values = c("#2171B5", "#C75C00")) +
  scale_color_manual(values = c("#2171B5", "#C75C00")) +
  scale_shape_manual(values = c(16,17)) +
  facet_wrap(~cluster, scales = "free_y") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 10, color = "black"),
        strip.text.x = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.key.size = unit(0.3, "cm"),
        legend.spacing.y = unit(0.05, "cm")) +
  ylab("Cluster proportion") +
  xlab("")
p6
##Variance Plot- SPF------------------------------------------------------------
library(variancePartition)
library(reformulas)
library(limma)

##load DDS objects, SPF WvsA design---------------------------------------------
#cTEC DDS
load("ddsobject.rda")
dds_cTEC <- dds
#Intertypical
load("ddsobject.rda")
dds_Intertypical <- dds
#mTECI
load("ddsobject.rda")
dds_mTECI <- dds
#mTECII
load("ddsobject.rda")
dds_mTECII <- dds
#mTECIII
load("ddsobject.rda")
dds_mTECIII <- dds
#mTECIV 
load("ddsobject.rda")
dds_mTECIV <- dds

dds_list <- list(
  "cTEC" = dds_cTEC,
  "Intertypical" = dds_Intertypical,
  "mTECI" = dds_mTECI,
  "mTECII" = dds_mTECII,
  "mTECIII" = dds_mTECIII,
  "mTECIV" = dds_mTECIV
)

##Look at variance 
extract_condition_variance <- function(dds_list) {
  library(DESeq2)
  library(variancePartition)
  
  results <- list()
  
  for (Name in names(dds_list)) {
    dds <- dds_list[[Name]]
    
    # VST normalization
    vsd <- vst(dds)
    
    form <- ~ condition + sex
    
    
    # Fit variance partition models
    varPart <- fitExtractVarPartModel(
      assay(vsd),
      form,
      colData(vsd)
    )
    
    vp_mean <- colMeans(varPart, na.rm = TRUE)
    
    df <- data.frame(
      factor   = names(vp_mean),
      variance = as.numeric(vp_mean),
      subset   = Name
    )
    
    results[[Name]] <- df
  }
  
  
  return(results)
}
results <- extract_condition_variance(dds_list)

# Combine into one data frame
combined <- do.call(rbind, results)
rownames(combined) <- NULL
combined

combined$factor <- factor(combined$factor, levels = c("Residuals", "sex", "condition"))

p <- ggplot(combined, aes(x = subset, y = variance, fill = factor)) +
  geom_col(width = 0.8, colour = "#000000") +
  scale_fill_manual(values = c("lightgrey", "#FFFFFF", "#000000"), labels = c("Residuals", "Sex", "Age")) +
  coord_flip() +
  theme_classic() +
  labs(x = NULL, y = "Mean proportion of variance", fill = "Factor") +
  ggtitle(label = paste0("SPF"))

p <- p + theme(plot.title = element_blank(),
               legend.text = element_text(size = 12),
               axis.text = element_text(size = 10),
               axis.title.x = element_text(size = 12),
               legend.title = element_blank())




