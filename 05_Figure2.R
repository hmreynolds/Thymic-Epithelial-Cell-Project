##05_Figure2
library(Seurat)
library(ComplexHeatmap)

##cTEC UMAP and cluster frequencies---------------------------------------------
TEC <- readRDS("SeuratObject.rds")
#Subset cTECs
Idents(TEC) <- TEC$meta_clusters
cTEC <- subset(TEC, idents = "cTEC")
#Run cTEC UMAP
Idents(cTEC) <- cTEC$labels
cTEC <- RunUMAP(cTEC, dims = 1:20, reduction = "harmony", reduction.name = "umap.cTEC")
order <- c("cTEC-1", "cTEC-2", "cTEC-3")
Idents(cTEC) <- factor(Idents(cTEC), levels = order)

#cTEC
my_cluster_order <- c("cTEC-1", "cTEC-2", "cTEC-3")
cluster_cols <- c("#08306B", "#5DADE2", "#9ECAE1")

Name <- "cTEC"
barcode_order <-c("WeaningMale_Rep1", "WeaningMale_Rep2", "WeaningMale_Rep3", 
                  "WeaningFemale_Rep1", "WeaningFemale_Rep2", "WeaningFemale_Rep3",
                  "AdultMale_Rep1", "AdultMale_Rep2", "AdultMale_Rep3", 
                  "AdultFemale_Rep1", "AdultFemale_Rep2", "AdultFemale_Rep3")

ClusterFreq <- table(cTEC$labels, cTEC$barcode)
ClusterFreq <- as.data.frame(ClusterFreq)
ClusterFreq$Var1 <- as.character(ClusterFreq$Var1) 
ClusterFreq$Var1 <- factor(ClusterFreq$Var1, levels = my_cluster_order)
ClusterFreq$Var2 <- factor(ClusterFreq$Var2, levels = barcode_order)

ClusterFreq$Var2 <- factor(ClusterFreq$Var2,
                           levels = c("WeaningMale_Rep1", "WeaningMale_Rep2", "WeaningMale_Rep3",
                                      "WeaningFemale_Rep1", "WeaningFemale_Rep2", "WeaningFemale_Rep3",
                                      "AdultMale_Rep1", "AdultMale_Rep2", "AdultMale_Rep3",
                                      "AdultFemale_Rep1", "AdultFemale_Rep2", "AdultFemale_Rep3"),
                           labels = c("WM-1", "WM-2", "WM-3",
                                      "WF-1", "WF-2", "WF-3",
                                      "AM-1", "AM-2", "AM-3",
                                      "AF-1", "AF-2", "AF-3"))

p2 <- ggplot(ClusterFreq, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 1.0, color = "black") +
  scale_fill_manual(values = cluster_cols) +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(), 
        legend.position = "none", 
        panel.border = element_blank(), 
        axis.line = element_line(size = 0.8, colour = "black"),
        axis.ticks = element_line(size = 0.8),
        axis.ticks.length = unit(0.25, "cm"),
        
  ) +
  labs(y = "Proportion of cTECs",
       x = c("M-1", "M-2", "M-3", "F-1", "F-2", "F-3",
             "M-1", "M-2", "M-3", "F-1", "F-2", "F-3")) 
p2


p3 <- DimPlot(cTEC, reduction = "umap.cTEC", cols = cluster_cols)
p3 <- p3 +  theme(
  axis.title = element_text(size = 10),
  axis.text = element_text(size = 10),
  legend.key.size = unit(0.3, "cm"),
  legend.spacing.y = unit(0.05, "cm")
)
p3


##FeaturePlots------------------------------------------------------------------
p4 <- FeaturePlot(cTEC, features = c("Cxcl12", "Ly6a", "Cd83", "Syngr1"), reduction = "umap.cTEC", ncol = 2) &
  theme(
    plot.title = element_text(size = 10),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.3, "cm"),
    legend.spacing.y = unit(0.05, "cm")
  )
p4

##Pathways----------------------------------------------------------------------
load("SPFcTECWvsA.rda")
res <- results(dds, name = "condition_Weaning_vs_Adult")

##Create GSEA objects
res_df <- data.frame(res@listData)
rownames(res_df) <- rownames(res)

res_GEA <- data.frame(res@listData)
rownames(res_GEA) <- rownames(res)
res_GEA <- res_GEA[!is.na(res_GEA$stat),]
ranks <- res_GEA$stat
names <- AnnotationDbi::select(org.Mm.eg.db,
                               keys = rownames(res_GEA),
                               columns = "ENTREZID",
                               keytype = "SYMBOL")
names(ranks) <- names$ENTREZID
ranks <- sort(ranks, decreasing = TRUE)
ranks <- ranks[!is.na(names(ranks))]

gsea_react <- gsePathway(
  geneList = ranks,
  organism = "mouse",
  pvalueCutoff = 0.05
)

gsea_go <- gseGO(
  geneList = ranks,
  pvalueCutoff = 0.05,
  OrgDb = org.Mm.eg.db,
  ont = "BP",
  keyType = "ENTREZID"
)


library(stringr)
library(ggplot2)

selected_GO <- c("GO:0002495", "GO:0019886", "GO:0032635", "GO:0043567", "GO:0030510", "GO:0071559",
                 "GO:0008610", "GO:0006694", "GO:0042060")
selected_Reactome <- c("R-MMU-9907900", "R-MMU-983168", "R-MMU-9912633", "R-MMU-1236974","R-MMU-1236975", 
                       "R-MMU-1236978", "R-MMU-2428928", "R-MMU-1474244")

GO_df <- as.data.frame(gsea_go) %>%
  filter(ID %in% selected_GO) %>%
  mutate(Database = "GO")

Reactome_df <- as.data.frame(gsea_react) %>%
  filter(ID %in% selected_Reactome) %>%
  mutate(Database = "Reactome")

combined_df <- bind_rows(GO_df, Reactome_df)

combined_df <- combined_df %>%
  mutate(Direction = case_when(
    NES > 1 ~ "Up in Weaning",
    NES < -1 ~ "Down in Weaning")
  )

combined_df$Direction <- factor(combined_df$Direction, levels = c("Up in Weaning", "Down in Weaning"))

library(tidytext)

Title = "GO+Reactome GSEA Pathways"
Name = "SPF cTEC Weaning vs Adult"
p5 <- ggplot(combined_df, aes(x = reorder_within(Description, NES, Direction))) +
  geom_bar(aes(y = NES), stat = "identity", fill = "grey85", width = 0.8, color = "darkgrey") +
  geom_point(aes(y = NES, size = -log10(p.adjust), fill = NES), shape = 21, color = "black") +
  facet_wrap(~Direction, ncol = 1, scales = "free_y", 
             labeller = as_labeller(c("Up in Weaning" = "Antigen Processing/Presentation",
                                      "Down in Weaning" = "Homeostatic Remodeling"))) +
  scale_fill_gradient2(low = "darkblue", mid = "white", high = "darkred", name = "NES") +
  scale_size_continuous(name = "-log10(p.adjust)") +
  scale_x_reordered(labels = function(x) str_wrap(str_remove(x, "___.+$"), width = 50)) +
  coord_flip() +
  labs(
    #title = paste0(Name, "\n" ,Title),
    x = "",
    y = "NES"
  ) +
  theme_bw() +
  theme(
    plot.subtitle = element_text(hjust = 0.5, face = "bold", size = 10),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    strip.text.x = element_text(size = 10, color = "black"),
    legend.position = "right", 
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10)
  )
p5



##Gene-level Heatmaps-----------------------------------------------------------
MHCI <- c("Tap1", "Tap2", "B2m", "Calr", "Canx", "Pdia3", 
          "Psmb9", "Psmb10", "Psmb11", "Psmb8", "Psmb5", "Psmb6", "Psmb7", "Nlrc5",
          "H2-K2",	  "H2-K1",	  "H2-D1",	 
          "H2-Q1",	  "H2-Q2",	  "H2-Q4","H2-Q5",	  "H2-Q6",	  "H2-Q7",	  "H2-Q10",	  "H2-T24",	  "H2-T23",	 
          "H2-T22",	  "H2-T10",	  "H2-T3",	  "H2-M10.2",	  "H2-M10.1",	  "H2-M10.3",	 
          "H2-M10.4",	  "H2-M11",	  "H2-M9",	  "H2-M1",	  "H2-M10.5",	  "H2-M10.6",	
          "H2-M5",	  "H2-M3",	  "H2-M2")
MHCII <- c("Ctsl", "Ifi30", "Cd74", "Prss16", "Ciita",   "H2-Oa",	  "H2-DMa",	 
           "H2-DMb2",	  "H2-DMb1",	  "H2-Ob",	  "H2-Ab1",	  "H2-Aa",	  "H2-Eb1",	
           "H2-Eb2",	  "H2-Ea",	  "H2-Bl", "Ifi30")
genes <- c(MHCI, MHCII)

load("SPFcTECWvsA.rda")
res <- results(dds, name = "condition_Weaning_vs_Adult")

res_df <- data.frame(res@listData)
rownames(res_df) <- rownames(res)
res_df <- res_df %>% 
  filter(padj < 0.05)

genes <- genes[genes %in% rownames(res_df)]

sample_order <- c("WeaningMale-Rep1",   "WeaningMale-Rep2",   "WeaningMale-Rep3",
                  "WeaningFemale-Rep1", "WeaningFemale-Rep2", "WeaningFemale-Rep3",
                  "AdultMale-Rep1", "AdultMale-Rep2",     "AdultMale-Rep3",
                  "AdultFemale-Rep1",   "AdultFemale-Rep2",   "AdultFemale-Rep3")

vsd <- vst(dds, blind = FALSE)
mat <- assay(vsd)
table(rownames(mat) %in% genes)
genes <- intersect(genes, rownames(mat))
mat <- mat[genes, ]
mat <- mat[,sample_order]
mat_scaled <- t(scale(t(mat)))

library(ComplexHeatmap)
annotation_row <- data.frame(
  log2FC_FMTvsGF = res_df[genes, "log2FoldChange"]
)
rownames(annotation_row) <- genes
annotation_row$log2FC_cat <- cut(annotation_row$log2FC_FMTvsGF,
                                 breaks = c(-Inf, -0.5, 0.5, Inf),
                                 labels = c("Down", "|<0.5|", "Up"))

annotation_row <- annotation_row[, c("log2FC_cat"), drop = FALSE]
colnames(annotation_row) <- "log2FC"

annotation_row$group <- case_when(
  rownames(annotation_row) %in% MHCI ~ "MHCI",
  rownames(annotation_row) %in% MHCII ~ "MHCII"
)

annotation_row$group <- factor(annotation_row$group, levels = c("MHCI", "MHCII"))
annotation_row$log2FC <- factor(annotation_row$log2FC, levels = c("Up", "|<0.5|", "Down"))

sample_order %in% colnames(mat_scaled)
mat_scaled <- mat_scaled[, sample_order]
p6 <-  pheatmap(
  mat_scaled,
  annotation_row = annotation_row[, "log2FC", drop = FALSE],
  name = "Expression",
  annotation_colors = list(
    log2FC = c(Down = "grey", "|<0.5|" = "white", Up = "#5E3C86")
  ),
  labels_col = c("M-1", "M-2", "M-3", "F-1", "F-2", "F-3",
                 "M-1", "M-2", "M-3", "F-1", "F-2", "F-3"),
  show_rownames = TRUE,
  show_colnames = TRUE,
  cluster_cols = FALSE,
  angle_col = "45",
  fontsize_row = 10,
  row_gap = unit(8, "mm"),
  row_title = NULL,
  row_split = annotation_row$group,
  cluster_row_slices = FALSE,
  annotation_names_row = FALSE,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100)
)
p6

