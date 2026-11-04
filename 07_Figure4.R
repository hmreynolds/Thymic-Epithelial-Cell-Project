##07_Figure4
library(decoupleR)
library(dorothea)
library(Seurat)
library(dplyr)
library(tibble)
library(dplyr)
library(tidyr)
library(rstatix)

##RestoredPathways---------------------------------------------------------------

##Read in GSEA files
SPFvsGF <- readRDS("cTEC_Weaning_SPFvsGF_GSEA_GO.rds")
SPFvsGF_FMT <- readRDS("cTEC_Weaning_SPFvsGF_FMT_GSEA_GO.rds")
GFvsGF_FMT <- readRDS("cTEC_Weaning_GFvsGF_FMT_GSEA_GO.rds")

##Create List of GSEA results (ranked)
gseaG_results_list <- list(
  "SPFvsGF"            = SPFvsGF,
  "SPFvsGF+FMT"        = SPFvsGF_FMT,
  "GFvsGF+FMT"         = GFvsGF_FMT
)

#Extract significant pathways per cluster into a single dataframe
pathway_df <- lapply(names(gseaG_results_list), function(cl) {
  res <- as.data.frame(gseaG_results_list[[cl]])
  res %>%
    filter(p.adjust < 0.05) %>%
    mutate(
      cluster   = cl,
      direction = ifelse(NES > 0, "UP", "DOWN")
    ) %>%
    dplyr::select(ID, Description, NES, p.adjust, setSize, cluster, direction)
}) %>%
  bind_rows() %>%
  mutate(cluster = paste(cluster))


#Compute semantic similarity matrix with rrvgo
all_ids <- unique(pathway_df$ID)


sim_matrix_G <- calculateSimMatrix(
  all_ids,
  orgdb     = "org.Mm.eg.db",
  ont       = "BP",              # rrvgo uses GO-style
  method    = "Rel"              # relevance-based similarity (recommended)
)

#Reduce redundant pathways
scores <- pathway_df %>%
  group_by(ID) %>%
  summarise(score = mean(-log10(p.adjust))) %>%
  deframe()

# Align scores to sim_matrix rows
scores_aligned <- scores[rownames(sim_matrix_G)]
scores_aligned[is.na(scores_aligned)] <- min(scores_aligned, na.rm = TRUE)

reduced_G <- reduceSimMatrix(
  sim_matrix_G,
  scores    = scores_aligned,
  threshold = 0.8,    # higher = more aggressive reduction; 0.7 is a good default
  orgdb     = "org.Mm.eg.db"
)

##Filter pathway_df on reduced terms 
# get representative terms only
parent_terms <- reduced_G %>%
  filter(term == parentTerm) %>%
  dplyr::select(go, term) %>%
  distinct()

#Add parent terms into pathway_df
pathway_df_parent <- pathway_df %>%
  left_join(
    reduced_G %>%
      dplyr::select(go, parentTerm),
    by = c("ID" = "go")
  )

# filter pathway_df to representatives only
pathway_df_reduced <- pathway_df_parent %>%
  filter(parentTerm %in% parent_terms$term)

pathway_df_reduced$cluster <- factor(pathway_df_reduced$cluster, levels = c("SPFvsGF", "SPFvsGF+FMT", "GFvsGF+FMT"))

##look at shared pathways- possible rescue (SPFvsGF & GFvsGF+FMT)
FMT <- pathway_df %>% filter(cluster == "GFvsGF+FMT" & NES < 0)
SPF <- pathway_df %>% filter(cluster == "SPFvsGF" & NES > 0)

shared <- FMT %>% filter(Description %in% SPF$Description)
shared <- pathway_df_reduced %>% filter(Description %in% shared$Description)

shared_id <- unique(shared$ID)

plot_df <- shared %>%
  group_by(cluster, parentTerm, direction) %>%
  summarise(
    NES = NES[which.min(p.adjust)],        # use most significant child term
    p.adjust = min(p.adjust),
    setSize = max(setSize),
    n_terms = n(),
    .groups = "drop"
  )

p1 <- plot_df %>%
  filter(p.adjust < 0.05) %>%
  ggplot(aes(x = cluster, y = parentTerm, 
             size = -log10(p.adjust), fill = NES)) +
  geom_point(shape = 21, stroke = 0.5, color = "black") +
  scale_fill_gradient2(low = "darkblue", mid = "white", high = "darkred", midpoint = 0) +
  #facet_grid(rows = vars(direction), scales = "free_y", space = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title.x = element_blank())  
p1 <- p1 +  scale_y_discrete(
  labels = function(x) str_wrap(x, width = 50))
p1

##Response to Interferon-beta---------------------------------------------------

pathway_df_reduced <- pathway_df %>%
  filter(grepl("interferon", Description, ignore.case = TRUE))

plot_df <- pathway_df_reduced %>%
  group_by(cluster, Description, direction) %>%
  summarise(
    NES = NES,      
    p.adjust = p.adjust,
    setSize = max(setSize),
    n_terms = n(),
    .groups = "drop"
  )
plot_df$cluster <- factor(plot_df$cluster, levels = c("SPFvsGF","GFvsGF+FMT"))

p2 <- plot_df %>%
  filter(p.adjust < 0.05) %>%
  ggplot(aes(x = cluster, y = Description, 
             size = -log10(p.adjust), fill = NES)) +
  geom_point(shape = 21, stroke = 0.5, color = "black") +
  scale_fill_gradient2(low = "darkblue", mid = "white", high = "darkred", midpoint = 0) +
  #facet_grid(rows = vars(direction), scales = "free_y", space = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title.x = element_blank())  
p2 <- p2 +  scale_y_discrete(
  labels = function(x) str_wrap(x, width = 50))
p2


##Interferon Regulons-----------------------------------------------------------
TEC <- readRDS("SeuratObject.rds")
TEC <- subset(TEC, Experiment == "SPF" & meta_cluster == "cTEC")

##Load regulons
data(dorothea_mm, package = "dorothea")

regulons <- dorothea_mm %>%
  filter(confidence %in% c("A", "B", "C"))

## Get matrix from Seurat Object
mat <- GetAssayData(
  TEC,
  assay = "RNA",
  layer = "data"
)

## Run decoupleR ULM
tf_acts <- run_ulm(
  mat = mat,
  network = regulons,
  .source = "tf",
  .target = "target",
  .mor = "mor"
)

#Store as Seurat assay object and scale
tf_mat <- tf_acts %>%
  select(score, source, condition) %>%
  tidyr::pivot_wider(
    names_from = condition,
    values_from = score
  ) %>%
  tibble::column_to_rownames("source") %>%
  as.matrix()

##Look at gene set size per regulon
setSize <- regulons %>%
  group_by(tf) %>%
  summarise(setSize = n(), .groups = "drop")

meta <- TEC@meta.data
meta$group <- paste(meta$meta_clusters, meta$Experiment, meta$barcode, sep = "/")

library(Matrix)

##Pull TF avg scores for each meta_cluster and Experiment (SPF/GF/GF+FMT) by replicate
meta_use <- meta[rownames(meta) %in% colnames(tf_mat), , drop = FALSE]

tf_avg <- sapply(unique(meta_use$group), function(g) {
  cells <- rownames(meta_use)[meta_use$group == g]
  Matrix::rowMeans(tf_mat[, cells, drop = FALSE])
})


tf_avg <- as.data.frame(tf_avg)

df <- tf_avg %>%
  tibble::rownames_to_column("TF") %>%
  pivot_longer(-TF, names_to = "group", values_to = "score")

df <- df %>%
  separate(group, into = c("cluster", "Experiment", "Replicate"), sep = "/")

df_summary <- df %>%
  group_by(TF, cluster, Experiment, Replicate) %>%
  summarise(mean_score = mean(score), .groups = "drop")

df_summary$Sample <- paste(df_summary$Experiment, df_summary$Replicate, sep = "_")

##Stats
library(rstatix)

run_regulon_stats <- function(df, cluster_name) {
  #subset cluster
  cluster_df <- df_summary %>%
    filter(cluster == cluster_name)
  
  #Kruskal-Wallis across conditions
  kw_results <- cluster_df %>%
    group_by(TF) %>%
    kruskal_test(mean_score ~ Experiment) %>%
    adjust_pvalue(method = "BH") %>%
    arrange(p.adj)
  
  #Significant TFs
  sig_TFs <- kw_results %>%
    filter(p.adj < 0.05) %>%
    pull(TF)
  
  #Post-hoc pairwise tests
  posthoc_results <- cluster_df %>%
    filter(TF %in% sig_TFs) %>%
    group_by(TF) %>%
    pairwise_wilcox_test(
      mean_score ~ Experiment,
      p.adjust.method = "BH"
    ) %>%
    add_significance()
  
  return(
    list(
      KW = kw_results, 
      significant_TFs = sig_TFs,
      pairwise = posthoc_results
    )
  )
}

stats <- run_regulon_stats(df, "cTEC")
statsPW <- stats$pairwise

ifn_tfs <- c("Irf3", "Irf9", "Irf1", "Irf2", "Irf7", "Stat1", "Stat2")

statsIFN <- statsPW %>% filter(TF %in% ifn_tfs)

##Create Heatmap
wide <- df_summary %>%
  pivot_wider(names_from = Experiment, values_from = mean_score)


cTECrep <- df_summary %>% 
  filter(cluster == "cTEC") %>%
  filter(TF %in% ifn_tfs)


heatmap_df <- cTECrep %>%
  group_by(TF, Experiment) %>%
  summarise(
    mean_score = mean(mean_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Experiment,
    values_from = mean_score
  )

mat <- heatmap_df %>%
  tibble::column_to_rownames("TF") %>%
  as.matrix()

condition_order <- c("SPF", "GF", "GF_FMT")

mat <- mat[, condition_order, drop = FALSE]

statsIFN_SPF <- statsIFN %>% filter(group1 == "GF" & group2 == "SPF") %>% select(TF, p.adj.signif)
colnames(statsIFN_SPF) <- c("TF", "SPFvsGF")
statsIFN_FMT <- statsIFN %>% filter(group1 == "GF" & group2 == "GF_FMT") %>% select(TF, p.adj.signif)
colnames(statsIFN_FMT) <- c("TF", "GF_FMTvsGF")
statsIFN_SvsF <- statsIFN %>% filter(group1 == "GF_FMT" & group2 == "SPF") %>% select(TF, p.adj.signif)
colnames(statsIFN_SvsF) <- c("TF", "SPFvsGF_FMT")

annotation_row <- statsIFN_SPF %>% 
  full_join(statsIFN_FMT, by = "TF") %>%
  full_join(statsIFN_SvsF, by = "TF") %>%
  column_to_rownames("TF")

library(ComplexHeatmap)
p3 <- pheatmap(
  mat,
  name = "Row-scaled Z-score",
  scale = "row",
  cluster_cols = FALSE,
  labels_col = c("SPF", "GF", "GF+FMT"),
  annotation_row = annotation_row,
  fontsize_row = 10,
  fontsize_col = 10,
  annotation_names_row = FALSE,
  annotation_colors = list(
    `SPFvsGF` = c("ns" = "white", "*" = "#BDAFCB", "**" = "#5E3C86"),
    `GF_FMTvsGF` = c("ns" = "white", "*" = "#A3D6A3", "**" = "#228B22"),
    `SPFvsGF_FMT` = c("ns" = "white", "*" = "#F0D99A", "**" = "#D99A00")
  ),
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  angle_col = "0"
)
p3


##IRSE/GAS Gene Heatmap---------------------------------------------------------

##load dds
load("cTEC_Weaning_MicrobialComparisons.rda")

##extract results for each comparison
res <- results(dds, contrast = c("condition", "SPF", "GF")) #SPFvsGF
res_df <- data.frame(res@listData)
rownames(res_df) <- rownames(res)
res_df <- res_df %>% filter(padj < 0.05)
res_df_SPF <- res_df


res <- results(dds, contrast = c("condition", "GF_FMT", "GF")) #FMTvsGF
res_df <- data.frame(res@listData)
rownames(res_df) <- rownames(res)
res_df <- res_df %>% filter(padj < 0.05)
res_df_FMT <- res_df 

res <- results(dds, contrast = c("condition", "SPF", "GF_FMT")) #SPFvsFMT
res_df <- data.frame(res@listData)
rownames(res_df) <- rownames(res)
res_df <- res_df %>% filter(padj < 0.05)
res_df_SPF_FMT <- res_df

##ISRE vs GAS
ISRE <- c("Mx1", "Mx2", "Isg15", "Oas1a", "Oas1b", "Oas2", "Ifi44", "Ifit1", "Ifit2", "Ifit3", 
          "Irf7", "Irf9", "Stat2", "Usp18")

GAS <- c("Cxcl9", "Cxcl10", "Gbp2", "Gbp5", "Stat1", "Irf1", "Ciita",
         "Psmb8", "Psmb9", "Psmb10", "Tap1", "Tap2", "Socs1", "Iigp1", "Igtp")

genes <- c(ISRE, GAS)

##Filter genes signficant against the GF group
genes <- genes[genes %in% rownames(res_df_SPF) | genes %in% rownames(res_df_FMT)]

##Plot
Name <- "cTEC Weaning"
Descrip <- "IRSE vs GAS Associated Genes"
vsd <- vst(dds, blind = FALSE)
mat <- assay(vsd)
table(rownames(mat) %in% genes)
genes <- intersect(genes, rownames(mat))
mat <- mat [genes, ]
mat_scaled <- t(scale(t(mat)))

library(pheatmap)
library(ComplexHeatmap)

annotation_row <- data.frame(
  log2FC_SPFvsGF = res_df_SPF[genes, "log2FoldChange"],
  log2FC_GF_FMTvsGF = res_df_FMT[genes, "log2FoldChange"],
  log2FC_SPFvsGF_FMT = res_df_SPF_FMT[genes, "log2FoldChange"]
)
rownames(annotation_row) <- genes
annotation_row$log2FC_SPF_cat <- cut(annotation_row$log2FC_SPFvsGF,
                                     breaks = c(-Inf, -0.05, 0.5, Inf),
                                     labels = c("Down", "|<0.5|", "Up"))

annotation_row$log2FC_FMT_cat <- cut(annotation_row$log2FC_GF_FMTvsGF,
                                     breaks = c(-Inf, -0.05, 0.5, Inf),
                                     labels = c("Down", "|<0.5|", "Up"))

annotation_row$log2FC_SPF_FMT_cat <- cut(annotation_row$log2FC_SPFvsGF_FMT,
                                         breaks = c(-Inf, -0.05, 0.5, Inf),
                                         labels = c("Down", "|<0.5|", "Up"))

annotation_row$group <- case_when(
  rownames(annotation_row) %in% ISRE ~ "ISRE",
  rownames(annotation_row) %in% GAS ~ "GAS"
)

annotation_row$group <- factor(annotation_row$group, levels = c("ISRE", "GAS"))

colnames(annotation_row) <- c("log2FC_SPF_all", "log2FC_GF_FMT_all", "log2FC_SPFvsGF_FMT_all", "log2FC_SPFvsGF", "log2FC_GF_FMTvsGF", "log2FC_SPFvsGF_FMT", "group")


rownames(annotation_row) <- genes
sample_order <- c("SPF-WeaningMale-Rep1", "SPF-WeaningMale-Rep2", "SPF-WeaningMale-Rep3", 
                  "SPF-WeaningFemale-Rep1", "SPF-WeaningFemale-Rep2", "SPF-WeaningFemale-Rep3",
                  "GF-WeaningMale-Rep1", "GF-WeaningMale-Rep2", "GF-WeaningMale-Rep3", 
                  "GF-WeaningFemale-Rep1", "GF-WeaningFemale-Rep2", "GF-WeaningFemale-Rep3",
                  "GF-FMT-WeaningMale-Rep1", "GF-FMT-WeaningMale-Rep2", 
                  "GF-FMT-WeaningFemale-Rep1", "GF-FMT-WeaningFemale-Rep2", "GF-FMT-WeaningFemale-Rep3")
sample_order %in% colnames(mat_scaled)
mat_scaled <- mat_scaled[, sample_order]
library(ComplexHeatmap)
p4 <- ComplexHeatmap::pheatmap(
  mat_scaled,
  name = "Expression",
  labels_col = c("SPF-WM1", "SPF-WM2", "SPF-WM3", "SPF-WF1", "SPF-WF2", "SPF-WF3",
                 "GF-WM1", "GF-WM2", "GF-WM3", "GF-WF1", "GF-WF2", "GF-WF3",
                 "FMT-WM1", "FMT-WM2", "FMT-WF1", "FMT-WF2", "FMT-WF3"),
  annotation_row = annotation_row[, c("log2FC_SPFvsGF", "log2FC_GF_FMTvsGF", "log2FC_SPFvsGF_FMT"), drop = FALSE],
  annotation_colors = list(
    log2FC_SPFvsGF = c(Down = "grey", "|<0.5|" = "white", Up = "#5E3C86"),
    log2FC_GF_FMTvsGF = c(Down = "grey", "|<0.5|" = "white", Up = "#228B22"),
    log2FC_SPFvsGF_FMT = c(Down = "grey", "|<0.5|" = "white", Up = "#D99A00")
  ),
  annotation_names_row = FALSE,
  show_rownames = TRUE,
  row_title = NULL,
  cluster_cols = FALSE,
  angle_col = "45",
  row_split = annotation_row$group,
  row_gap = unit(8, "mm"),
  fontsize_row = 10,
  cluster_row_slices = FALSE,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100)
)

p4


##cTEC Function DEG-------------------------------------------------------------
load("cTEC_Weaning_MicrobialComparisons.rda")

##Set gene targets
chemokines <- c("Cxcl1",	"Cxcl10",	"Cxcl9")
MHCI <- c("Pdia3",	"Nlrc5",	"Psmb10",	"Psmb6",	"Psmb5",	"Psmb11",	"H2-K1",
          "Tap2",	"H2-D1",	"H2-Q2",	"H2-Q4",	"H2-Q6",	"H2-Q7",	"H2-T23",	"H2-M3")
MHCII <- c("Ifi30",	"Prss16",	"Ciita",	"H2-DMa",	"H2-DMb1",	"H2-Eb1")
proteasome <- c("Psmd1",	"Psmd5",	"Psmd14",	"Psmd4",	"Psmb2",	"Psmc2",	"Psmc4",
                "Psmd7",	"Psma4",	"Psme4",	"Psmc5",
                "Psmd12",	"Psma6",	"Psma3",	"Psmg2")
genes <- c(chemokines, MHCI, MHCII, proteasome)

##Extract dds results and filter for signficance 
resultsNames(dds)
Name <- "cTEC Weaning"
res <- results(dds, contrast = c("condition", "GF_FMT", "GF")) #GF+FMTvsGF
sample_order <- c("GF-WeaningMale-Rep1", "GF-WeaningMale-Rep2", "GF-WeaningMale-Rep3", 
                  "GF-WeaningFemale-Rep1", "GF-WeaningFemale-Rep2", "GF-WeaningFemale-Rep3",
                  "GF-FMT-WeaningMale-Rep1", "GF-FMT-WeaningMale-Rep2", 
                  "GF-FMT-WeaningFemale-Rep1", "GF-FMT-WeaningFemale-Rep2", "GF-FMT-WeaningFemale-Rep3")

res_df <- data.frame(res@listData)
rownames(res_df) <- rownames(res)
res_df <- res_df %>% filter(padj < 0.05)

genes <- genes[genes %in% rownames(res_df)]

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

annotation_row$group <- case_when(
  rownames(annotation_row) %in% chemokines ~ "Chemokine",
  rownames(annotation_row) %in% proteasome ~ "Proteasome",
  rownames(annotation_row) %in% MHCI ~ "MHCI",
  rownames(annotation_row) %in% MHCII ~ "MHCII")
)

colnames(annotation_row) <- c("log2FC", "group")
annotation_row$group <- factor(annotation_row$group, levels = c("MHCI", "MHCII", "Proteasome", "Chemokine"))
library(ComplexHeatmap)

sample_order %in% colnames(mat_scaled)
mat_scaled <- mat_scaled[, sample_order]

p5 <- ComplexHeatmap::pheatmap(
  t(mat_scaled),
  annotation_col = annotation_row[, "log2FC", drop = FALSE],
  name = "Expression",
  annotation_colors = list(
    log2FC = c(Down = "grey", "|<0.5|" = "white", Up = "#5E3C86")
  ),
  labels_row = c("GF-WM1", "GF-WM2", "GF-WM3", "GF-WF1", "GF-WF2", "GF-WF3",
                 "FMT-WM1", "FMT-WM2", "FMT-WF1", "FMT-WF2", "FMT-WF3"),
  show_rownames = TRUE,
  show_colnames = TRUE,
  cluster_rows = FALSE,
  fontsize_row = 10,
  angle_col = "45",
  column_gap = unit(4, "mm"),
  column_title = NULL,
  annotation_names_col = FALSE,
  column_split = annotation_row$group,
  cluster_column_slices = FALSE,
  annotation_names_row = FALSE,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
)
p5

