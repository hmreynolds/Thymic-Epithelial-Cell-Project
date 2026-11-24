##06_Figure3

library(Seurat)
library(dplyr)
library(decoupleR)
library(dorothea)
library(tibble)
library(dplyr)
library(tidyr)
library(rstatix)
library(igraph)
library(ggraph)
library(ggplot2)
library(tibble)
library(scatterpie)
library(ggrepel)
library(rrvgo)
library(DESeq2)

##Differential age-associated interactions by cluster---------------------------
TEC <- readRDS("SeuratObject.rds")

##Exclude GF+FMT data
TEC <- subset(TEC, Experiment != "GF_FMT")

my_cluster_order <- c("Nurse_cTEC", "cTEC", "Intertypical_TEC", "mTECI", "mTECII", "mTECIII", "mTECIV", "Muscle-like_TEC", "Microfold_TEC", "Neuroendocrine_TEC", "Ciliated_TEC")

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
cluster_props <- cluster_props %>% filter(dataset %in% c("SPF", "GF"))

#Feed Seurat Object into propeller
library(speckle)
library(limma)
props <- getTransformedProps(
  clusters = Idents(TEC),
  sample = TEC@meta.data$barcode2,
  transform = "asin"
)

sample_meta <- TEC@meta.data %>%
  distinct(barcode2, Experiment, age) %>%
  arrange(match(barcode2, colnames(props$Proportions)))  # match column order

sample_meta$dataset <- sample_meta$Experiment

# Build design matrix for the factorial model
design <- model.matrix(~ dataset * age, data = sample_meta)

# Fit limma model across clusters
fit <- lmFit(props$TransformedProps, design)
fit <- eBayes(fit)

# Pull out results — e.g., interaction term (dataset:age)
table <- topTable(fit, coef = "datasetSPF:ageWeaning", sort.by = "P", number = Inf)

stat_results <- topTable(fit, coef = "datasetSPF:ageWeaning", number = Inf) %>%
  tibble::rownames_to_column("cluster") %>%
  mutate(
    p.signif = case_when(
      adj.P.Val < 0.001 ~ "***",
      adj.P.Val < 0.01  ~ "**",
      adj.P.Val < 0.05  ~ "*",
      TRUE ~ paste0("p=", signif(adj.P.Val, 2))
    ),
    # position label near top of each facet's y-range — adjust per your data
    y.position = NA  # fill in based on max proportion per cluster, see below
  )

# get a sensible y position per cluster (e.g., just above the max proportion)
y_positions <- cluster_props %>%
  group_by(cluster) %>%
  summarise(y.position = max(proportion) * 1.1)

stat_results <- stat_results %>%
  filter(cluster %in% c("cTEC", "mTECI", "Intertypical_TEC", "mTECII", "mTECIII", "mTECIV"))

stat_results <- left_join(stat_results, y_positions, by = "cluster")

cluster_props$age <- factor(cluster_props$age, levels = c("Weaning", "Adult"))

p1 <- ggplot(cluster_props, aes(x = dataset, y = proportion, fill = age)) +
  geom_boxplot(
    aes(group = interaction(dataset, age)),
    outlier.shape = NA, 
    alpha = 0.7,
    position = position_dodge(width = 0.75)
  ) +
  geom_point(
    aes(color = age, 
        shape = sex,
        group = interaction(dataset, age)),
    position = position_jitterdodge(
      jitter.width = 0.5,
      dodge.width = 0.75),
    size = 2,
    alpha = 0.8
  ) +
  geom_text(
    data = stat_results,
    aes(
      x = 1.5, 
      y = y.position.y, 
      label = paste0("interaction: ", p.signif),
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  scale_fill_manual(values = c("#0072B2", "#D55E00")) +
  scale_color_manual(values = c("#0072B2", "#D55E00")) +
  scale_shape_manual(values = c(16,17)) +
  facet_wrap(~cluster, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title.y = element_text(size = 12),
        strip.text.x = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.key.size = unit(0.3, "cm"),
        legend.spacing.y = unit(0.05, "cm")) +
  ylab("Cluster proportion") +
  xlab("")
p1

##Network Plot------------------------------------------------------------------

Name <- "cTEC"
##Read in GSEA files
SPFvsGF_W <- readRDS("SPFvsGF_Weaning_GSEA_GO.rds")
SPFvsGF_A <- readRDS("SPFvsGF_Adult_GSEA_GO.rds")
##Create List of GSEA results (ranked)
gseaG_results_list <- list(
  "SPFvsGF_W"            = SPFvsGF_W,
  "SPFvsGF_A"            = SPFvsGF_A
)

cluster_colors <- c(
  "SPFvsGF_W_UP"       = "#78C679",
  "SPFvsGF_W_DOWN"     = "#7570B3",
  "SPFvsGF_A_UP"       = "#00416a",
  "SPFvsGF_A_DOWN"     = "#ff7256"
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
  mutate(cluster = paste(cluster, ifelse(NES > 0, "UP", "DOWN"), sep = "_"))


cluster_names <- unique(pathway_df$cluster)
cat("Total significant pathway-cluster rows:", nrow(pathway_df), "\n")
cat("Unique GO terms:", length(unique(pathway_df$ID)), "\n")

#Compute semantic similarity matrix with rrvgo

all_ids <- unique(pathway_df$ID)


sim_matrix_G <- calculateSimMatrix(
  all_ids,
  orgdb     = "org.Mm.eg.db",
  ont       = "BP",              # rrvgo uses GO-style
  method    = "Rel"              # relevance-based similarity (recommended)
)

cat("Similarity matrix dimensions:", dim(sim_matrix_G), "\n")


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


#Extract coordinates from scatterPlot via ggplot_build


sp <- scatterPlot(sim_matrix_G, reduced_G, algorithm = "umap", addLabel = FALSE)
sp_data <- ggplot_build(sp)$data[[1]]   # first layer = points

# sp_data has x, y, and PANEL/group but not GO IDs
# Row order matches rownames(sim_matrix) — bind directly
coords_df <- data.frame(
  go = rownames(sim_matrix_G),
  x  = sp_data$x,
  y  = sp_data$y
)

#Build node attributes — one row per parent term

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


# Pre-filter: only keep parent/representative terms as nodes
parent_lookup <- reduced_G %>%
  filter(term == parentTerm) %>%          # a term is its own parent = representative
  dplyr::select(go, term) %>%
  distinct() %>%
  mutate(
    go        = as.character(go),
    short_name = ifelse(nchar(term) > 45, paste0(substr(term, 1, 42), "..."), term)
  )

parent_ids <- parent_lookup$go

# Mean setSize and best p.adjust per parent group
node_stats <- reduced_G %>%
  mutate(go = as.character(go)) %>%
  left_join(
    pathway_df %>% dplyr::select(ID, setSize, p.adjust, NES, cluster),
    by = c("go" = "ID") 
  ) %>%
  group_by(parentTerm) %>%
  summarise(
    node_size  = mean(setSize,   na.rm = TRUE),
    best_padj  = min(p.adjust,   na.rm = TRUE),
    mean_NES   = mean(NES,       na.rm = TRUE),
    .groups    = "drop"
  ) %>%
  mutate(parentTerm = as.character(parentTerm))

# Cluster proportions per parent node
cluster_props <- reduced_G %>%
  mutate(go = as.character(go)) %>%
  left_join(
    pathway_df %>%
      dplyr::select(ID, pathway_cluster = cluster),
    by = c("go" = "ID")
  ) %>%
  filter(!is.na(pathway_cluster)) %>%
  group_by(parentTerm, pathway_cluster) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(parentTerm) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  pivot_wider(
    names_from = pathway_cluster,
    values_from = prop,
    values_fill = 0
  )

cluster_props_collapsed <- cluster_props %>%
  group_by(parentTerm) %>%
  summarise(
    n = sum(n),
    across(-n, sum),
    .groups = "drop"
  )
# UMAP coordinates for parent nodes only
parent_coords <- coords_df %>%
  filter(coords_df$go %in% parent_ids)

# Assemble nodes — start from parent_ids as base to prevent row loss
library(scales)
nodes <- data.frame(id = parent_ids) %>%
  left_join(parent_lookup           %>% dplyr::rename(id = go), by = "id") %>%
  left_join(cluster_props_collapsed  %>% dplyr::rename(term = parentTerm), by = "term") %>%
  left_join(node_stats               %>% dplyr::rename(term = parentTerm), by = "term") %>%
  left_join(parent_coords            %>% dplyr::rename(id = go), by = "id") %>%
  mutate(
    node_size  = rescale(node_size, to = c(2, 10)),
    node_size  = ifelse(is.na(node_size), 2, node_size)
  )

# Fill any missing cluster columns with 0
for (cl in cluster_names) {
  if (!cl %in% colnames(nodes)) {
    nodes[[cl]] <- 0
  } else {
    nodes[[cl]][is.na(nodes[[cl]])] <- 0
  }
}


##Repel overlapping nodes

library(ggforce)   # install.packages("ggforce")

# After building nodes, repel overlapping x/y positions
repel_nodes <- function(nodes, min_dist = 0.2, iterations = 100) {
  xy <- as.matrix(nodes[, c("x", "y")])
  r  <- -log10(nodes$best_padj) / 8   # same radius as pie plot
  
  for (i in seq_len(iterations)) {
    for (a in seq_len(nrow(xy))) {
      for (b in seq_len(nrow(xy))) {
        if (a == b) next
        dx   <- xy[a, 1] - xy[b, 1]
        dy   <- xy[a, 2] - xy[b, 2]
        dist <- sqrt(dx^2 + dy^2)
        min_d <- r[a] + r[b] + min_dist
        if (dist < min_d && dist > 0) {
          push  <- (min_d - dist) / 2
          xy[a, 1] <- xy[a, 1] + push * dx / dist
          xy[a, 2] <- xy[a, 2] + push * dy / dist
          xy[b, 1] <- xy[b, 1] - push * dx / dist
          xy[b, 2] <- xy[b, 2] - push * dy / dist
        }
      }
    }
  }
  nodes$x <- xy[, 1]
  nodes$y <- xy[, 2]
  nodes
}

nodes <- repel_nodes(nodes, min_dist = 1.0)


#Build edge list from similarity matrix

# Subset sim_matrix to parent nodes only
sim_sub <- sim_matrix_G[parent_ids, parent_ids]

edge_cutoff <- 0.3

edge_df <- as.data.frame(as.table(sim_sub)) %>%
  setNames(c("from", "to", "similarity")) %>%
  mutate(
    from       = as.character(from),
    to         = as.character(to),
    similarity = as.numeric(similarity)
  ) %>%
  filter(
    from < to,                     # no duplicate/self edges
    similarity > edge_cutoff,
    from %in% nodes$id,
    to   %in% nodes$id
  )

cat("Edges after filtering:", nrow(edge_df), "\n")

#Join coordinates to edges

coord_lookup <- nodes %>% dplyr::select(id, x, y)
coord_lookup <- nodes %>%
  distinct(id, .keep_all = TRUE) %>%
  dplyr::select(id, x, y)

plot_edges <- edge_df %>%
  left_join(coord_lookup, by = c("from" = "id")) %>%
  dplyr::rename(x_from = x, y_from = y) %>%
  left_join(coord_lookup, by = c("to" = "id")) %>%
  dplyr::rename(x_to = x, y_to = y)

cat("NA in edge coords:", sum(is.na(plot_edges$x_from) | is.na(plot_edges$x_to)), "\n")

##Winsorize coordinates or scale factor
#Winsorize
xlim <- quantile(nodes$x, c(0.1, 0.99))
ylim <- quantile(nodes$y, c(0.1, 0.99))

nodes <- nodes %>%
  mutate(
    x = pmin(pmax(x, xlim[1]), xlim[2]),
    y = pmin(pmax(y, ylim[1]), ylim[2])
  )

plot_edges <- plot_edges %>%
  mutate(
    x_to = pmin(pmax(x_to, xlim[1]), xlim[2]),
    x_from = pmin(pmax(x_from, xlim[1]), xlim[2]),
    y_to = pmin(pmax(y_to, ylim[1]), ylim[2]),
    y_from = pmin(pmax(y_from, ylim[1]), ylim[2]),
  )
##Scale Factor
scale_factor <- -1
nodes <- nodes %>%
  mutate(x = x * scale_factor, y = y * scale_factor)
plot_edges <- plot_edges %>%
  mutate(x_from = x_from * scale_factor, y_from = y_from * scale_factor,
         x_to = x_to * scale_factor, y_to = y_to * scale_factor)

#Plot
order <- c("SPFvsGF_W_UP", "SPFvsGF_W_DOWN", "SPFvsGF_A_UP", "SPFvsGF_A_DOWN")

##Plot with node size = best padj 
p <- ggplot() +
  # Edges
  geom_segment(
    data = plot_edges,
    aes(x = x_from, y = y_from, xend = x_to, yend = y_to,
        alpha = similarity),
    color     = "grey60",
    linewidth = 0.5,
    show.legend = FALSE
  ) +
  scale_alpha_continuous(range = c(0.05, 0.5)) +
  # Pie chart nodes
  geom_scatterpie(
    data      = nodes,
    aes(x = x, y = y, r = (-log10(best_padj))/8),
    cols      = cluster_names,
    color     = "white",
    linewidth = 0.25
  ) +
  scale_fill_manual(values = cluster_colors, name = "Comparisons", breaks = order,
                    labels = c(
                      "SPFvsGF_W_UP" = "SPF Weaning",
                      "SPFvsGF_W_DOWN" = "GF Weaning",
                      "SPFvsGF_A_UP" = "SPF Adult",
                      "SPFvsGF_A_DOWN" = "GF Adult"
                    )) +
  # Pathway labels — parent terms only
  geom_text_repel(
    data          = nodes,
    aes(x = x, y = y, label = short_name),
    size          = 3.0,
    fontface      = "bold",
    color         = "grey10",
    bg.color      = "white",
    bg.r          = 0.1,
    box.padding   = 0.5,
    point.padding = 0.3,
    max.overlaps  = 25,
    segment.color = "black",
    segment.size  = 0.25,
    min.segment.length = 0.2
  ) +
  labs(
    title    = "GSEA GO network cTEC Weaning and Adult SPF vs GF",
    subtitle = paste0(
      "Nodes = representative GO terms  |  ",
      "Pie = cluster proportion  |  ",
      "Edge = semantic similarity > ", edge_cutoff
    ),
    fill = "Cluster"
  )  +
  theme_void(base_size = 10) +
  theme(
    plot.title      = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle   = element_text(size = 9, color = "grey40", hjust = 0.5),
    legend.position = "right",
    legend.title    = element_text(face = "bold", size = 10),
    plot.margin     = margin(5, 5, 5, 5)
  ) +
  coord_equal() +
  geom_scatterpie_legend(
    -log10(nodes$best_padj)/8,          # must match the r= argument used in geom_scatterpie
    x     = min(nodes$x) - 6,   # position — adjust to avoid overlap
    y     = min(nodes$y) - 1,
    n     = 3,                     # number of circles to show in legend
    labeller = function(r) {       # convert radius back to original setSize
      round(rescale(r * 20,
                    from = c(2, 10),
                    to   = range(-log10(node_stats$best_padj), na.rm = TRUE)))
    }
  )
p

##Interferon pathways-----------------------------------------------------------

##extract interferon related pathways from reduced pathways
interferon_cTEC <- pathway_df_reduced %>% filter(grepl("interferon", cTECweaning$Description))
cols <- c("SPFvsGF_W_UP" = "#78C679", "SPFvsGF_W_DOWN"  = "#7570B3", "SPFvsGF_A_UP" = "#00416a", "SPFvsGF_A_DOWN" = "#ff7256")
interferon_cTEC$comparison <- interferon_cTEC$cluster
order <- c("SPFvsGF_W_UP", "SPFvsGF_W_DOWN", "SPFvsGF_A_UP")
p2 <- interferon_cTEC %>%
  ggplot(aes(x = NES, 
             y = reorder(Description,NES), 
             size = -log10(p.adjust), 
             color = comparison)) +
  geom_point() +
  scale_color_manual(values = cols, breaks = order,
                     labels = c(
                       "SPFvsGF_W_UP" = "SPF Weaning",
                       "SPFvsGF_W_DOWN" = "GF Weaning",
                       "SPFvsGF_A_UP" = "SPF Adult"
                     ))+
  theme_bw() +
  theme(axis.text.x = element_text(size = 8, color = "black"), 
        plot.title = element_text(hjust = 0, size = 12),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.title.y = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10))    +
  labs(
    title    = str_wrap("Response to Interferon-beta", width = 30) 
  )
p2 <- p2 + scale_y_discrete(labels = function(x) str_wrap(x, width = 30))
p2


##DeltaRegulon-------------------------------------------------------------------
TEC <- readRDS("SeuratObject.rds")
##Extract weaning
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

##Compare Groups
df_summary <- df %>%
  group_by(TF, cluster, Experiment) %>%
  summarise(
    mean_score = mean(score, na.rm = TRUE),
    sd_score = sd(score, na.rm = TRUE),
    n = sum(!is.na(mean_score)),
    .groups = "drop"
  )

wide <- df_summary %>%
  tidyr::pivot_wider(names_from = Experiment, values_from = c(mean_score, sd_score, n))


SPFvsGF <- wide %>%
  mutate(delta = mean_score_SPF - mean_score_GF,
         delta_sd = sqrt(sd_score_SPF^2 + sd_score_GF^2),
         delta_se = sqrt(sd_score_SPF^2 / n_SPF +
                           sd_score_GF^2 / n_GF))

data <- SPFvsGF

Subset <- data %>% filter(cluster == "cTEC")
topSubset <- Subset %>% 
  arrange(desc(delta)) %>% 
  slice_head(n=20)
bottomSubset <- Subset %>% 
  arrange(delta) %>% 
  slice_head(n=20)

cool <- rbind(topSubset, bottomSubset)
cool <- Subset %>% filter(TF %in% cool$TF)
cool$direction <- ifelse(cool$delta > 0, "SPF", "GF")
cool$direction <- factor(cool$direction, levels = c("SPF", "GF"))

cool$setSize <- cool %>%
  left_join(setSize, by = c("TF" = "tf"))
cool$setSize <- cool$setSize$setSize

p <- cool %>% ggplot(aes(
  x = delta, 
  y = reorder(TF, delta), 
  color = delta, 
  size = setSize)) +
  geom_point() +
  geom_errorbar(
    aes(
      xmin = delta - delta_se,
      xmax = delta + delta_se),
    width = 0.05,
    linewidth = 0.3
  ) +
  scale_color_gradient2(low = "darkblue", mid = "white", 
                        high = "darkred", midpoint = 0) +
  facet_grid(rows = vars(direction), scales = "free_y", space = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8, color = "black"),
        axis.text.y = element_text(size = 8, color = "black"))  +
  labs(
    x = paste0("Delta Regulon Score"),
    y = "Transcrition Factor")
p

