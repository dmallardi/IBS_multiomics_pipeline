#############################################################################
#                       IBS MULTI-OMICS PIPELINE
#############################################################################

#=============================================================================
#STEP 1 - Import and merge patient metadata
#=============================================================================

#------------------------------------------------------------------------------
# Clean environment
#------------------------------------------------------------------------------

rm(list = ls())

#------------------------------------------------------------------------------
# Load required libraries
#------------------------------------------------------------------------------

library(dplyr)
library(readr)
library(tidyr)
library(phyloseq)
library(MintTea)
library(ggplot2)
library(igraph)
library(ggraph)
library(tidygraph)
library(pheatmap)
library(pROC)
library(glmnet)
library(caret)
library(car)
library(reticulate)

#------------------------------------------------------------------------------
# Helper function
#------------------------------------------------------------------------------

# Standardize patient identifiers
clean_patient_id <- function(x) {
  trimws(toupper(x))
}

#------------------------------------------------------------------------------
# Import datasets
#------------------------------------------------------------------------------

patients <- read_csv(
  "C:/Users/utente/Desktop/RICERCA/IBS Mintea/datasets/PRJNA812699/Patients.csv",
  show_col_types = FALSE
)

metadata <- read_delim(
  "C:/Users/utente/Desktop/RICERCA/IBS Mintea/datasets/PRJNA812699/metadata1.csv",
  show_col_types = FALSE
)

#------------------------------------------------------------------------------
# Harmonize patient identifiers
#------------------------------------------------------------------------------

patients$Patient <- clean_patient_id(patients$Patient)
metadata$Patient <- clean_patient_id(metadata$Patient)

#------------------------------------------------------------------------------
# Check for duplicated patient IDs
#------------------------------------------------------------------------------

stopifnot(!anyDuplicated(patients$Patient))
stopifnot(!anyDuplicated(metadata$Patient))

#------------------------------------------------------------------------------
# Merge patient information with metadata
#------------------------------------------------------------------------------

final_data <- patients %>%
  left_join(metadata, by = "Patient")

#------------------------------------------------------------------------------
# Identify missing metadata after merging
#------------------------------------------------------------------------------

missing_patients <- final_data %>%
  filter(is.na(Group))

#------------------------------------------------------------------------------
# Recover missing patient(s) from supplementary metadata
#------------------------------------------------------------------------------

metadata2 <- read_delim(
  "C:/Users/utente/Desktop/RICERCA/IBS Mintea/datasets/PRJNA812699/metadata2.csv",
  show_col_types = FALSE
)

metadata2$Patient <- clean_patient_id(metadata2$Patient)

# Remove the first column (not required)
metadata2 <- metadata2[, -1]

# Retrieve missing patient(s)
missing_row <- metadata2 %>%
  filter(Patient %in% missing_patients$Patient)

#------------------------------------------------------------------------------
# Rebuild complete metadata table
#------------------------------------------------------------------------------

metadata_complete <- bind_rows(metadata, missing_row)

# Verify absence of duplicated IDs
stopifnot(!anyDuplicated(metadata_complete$Patient))

#------------------------------------------------------------------------------
# Final merge
#------------------------------------------------------------------------------

final_data <- patients %>%
  left_join(metadata_complete, by = "Patient")

#------------------------------------------------------------------------------
# Quality control
#------------------------------------------------------------------------------

stopifnot(sum(is.na(final_data$Group)) == 0)

cat("Metadata successfully merged.\n")
cat("Number of patients:", nrow(final_data), "\n")

#------------------------------------------------------------------------------
# Optional: save merged metadata
#------------------------------------------------------------------------------

# write.csv(
#   final_data,
#   "C:/Users/utente/Desktop/RICERCA/IBS Mintea/datasets/PRJNA812699/dati_finali.csv",
#   row.names = FALSE
# )
# =============================================================================
# STEP 2 - Taxonomic data preparation
# =============================================================================

#------------------------------------------------------------------------------
# Import ASV abundance table and taxonomy
#------------------------------------------------------------------------------

asv <- read.csv(
  "C:/Users/utente/Desktop/RICERCA/IBS Mintea/datasets/PRJNA812699/data_request/16s_data.csv",
  row.names = 1
)

tax <- read.csv(
  "C:/Users/utente/Desktop/RICERCA/IBS Mintea/datasets/PRJNA812699/data_request/16s_taxonomy_table_updated.csv",
  row.names = 1
)

#------------------------------------------------------------------------------
# Quality control
#------------------------------------------------------------------------------

stopifnot(all(rownames(asv) %in% rownames(tax)))

#------------------------------------------------------------------------------
# Create phyloseq object
#------------------------------------------------------------------------------

OTU <- otu_table(as.matrix(asv), taxa_are_rows = TRUE)
TAX <- tax_table(as.matrix(tax))

physeq <- phyloseq(OTU, TAX)

#------------------------------------------------------------------------------
# Replace missing genus annotations
#------------------------------------------------------------------------------

tax_table(physeq)[, "Genus"][
  is.na(tax_table(physeq)[, "Genus"])
] <- "Unknown"

#------------------------------------------------------------------------------
# Aggregate ASVs at genus level
#------------------------------------------------------------------------------

ps_genus <- tax_glom(
  physeq,
  taxrank = "Genus"
)

otu <- as.data.frame(otu_table(ps_genus))
tax <- as.data.frame(tax_table(ps_genus))

otu$Genus <- tax$Genus

genus_agg <- otu %>%
  group_by(Genus) %>%
  summarise(across(where(is.numeric), sum))

genus_mat <- as.data.frame(genus_agg)

rownames(genus_mat) <- genus_mat$Genus
genus_mat <- genus_mat[, -1]

# Samples as rows, genera as columns
genus_mat <- t(genus_mat)

#------------------------------------------------------------------------------
# Retain only patients included in the clinical dataset
#------------------------------------------------------------------------------

ids <- intersect(
  patients$Patient,
  rownames(genus_mat)
)

genus_filtered <- genus_mat[ids, ]

stopifnot(
  nrow(genus_filtered) == length(ids)
)

#------------------------------------------------------------------------------
# Prefix taxonomic features
#------------------------------------------------------------------------------

colnames(genus_filtered) <- ifelse(
  grepl("^T__", colnames(genus_filtered)),
  colnames(genus_filtered),
  paste0("T__", colnames(genus_filtered))
)

#------------------------------------------------------------------------------
# Remove genera with zero variance
#------------------------------------------------------------------------------

genus_filtered <- genus_filtered[
  ,
  apply(genus_filtered, 2, var) != 0
]

#------------------------------------------------------------------------------
# Remove low-prevalence genera
# Present in >15% of samples
#------------------------------------------------------------------------------

genus_filtered <- genus_filtered[
  ,
  colSums(genus_filtered > 0) >
    (0.15 *nrow(genus_filtered))
]

#------------------------------------------------------------------------------
# Create TAXA matrix
#------------------------------------------------------------------------------

TAXA <- genus_filtered

# Remove unclassified genus
TAXA <- TAXA[
  ,
  colnames(TAXA) != "T__Unknown"
]

#------------------------------------------------------------------------------
# Total Sum Scaling (TSS) normalization
#------------------------------------------------------------------------------

TAXA <- TAXA / rowSums(TAXA)

#------------------------------------------------------------------------------
# Final quality control
#------------------------------------------------------------------------------

stopifnot(sum(is.na(TAXA)) == 0)

#------------------------------------------------------------------------------
# Prepare metadata
#------------------------------------------------------------------------------

final_data <- final_data[
  !duplicated(final_data$Patient),
]

metadata_filtered <- final_data[
  final_data$Patient %in% rownames(TAXA),
]

metadata_filtered <- metadata_filtered[
  match(
    rownames(TAXA),
    metadata_filtered$Patient
  ),
]

stopifnot(
  all(metadata_filtered$Patient == rownames(TAXA))
)

#------------------------------------------------------------------------------
# Summary
#------------------------------------------------------------------------------

cat("Taxonomic matrix successfully generated.\n")
cat("Samples:", nrow(TAXA), "\n")
cat("Genera :", ncol(TAXA), "\n")

#------------------------------------------------------------------------------
# Optional: save filtered metadata
#------------------------------------------------------------------------------

# write.csv(
#   metadata_filtered,
#   "C:/Users/utente/Desktop/RICERCA/IBS Mintea/datasets/PRJNA812699/metadata_filtered_new.csv",
#   row.names = FALSE
# )
# =============================================================================
# STEP 3 - Metatranscriptomic data preparation
# =============================================================================

#------------------------------------------------------------------------------
# Import metatranscriptomic dataset
#------------------------------------------------------------------------------

metagenes <- read.csv(
  "C:/Users/utente/Desktop/RICERCA/IBS Mintea/datasets/PRJNA812699/data_request/Transcripts_KO_data.csv",
  row.names = 1
)

#------------------------------------------------------------------------------
# Harmonize sample identifiers
#------------------------------------------------------------------------------

colnames(metagenes) <- trimws(colnames(metagenes))

ids <- trimws(patients$Patient)

#------------------------------------------------------------------------------
# Retain only samples included in the clinical dataset
#------------------------------------------------------------------------------

ids <- intersect(ids, colnames(metagenes))

metagenes_filtered <- metagenes[, ids]

#------------------------------------------------------------------------------
# Transpose matrix
# Samples as rows, pathways as columns
#------------------------------------------------------------------------------

metagenes_filtered <- t(metagenes_filtered)

#------------------------------------------------------------------------------
# Prefix pathway features
#------------------------------------------------------------------------------

colnames(metagenes_filtered) <- paste0(
  "G__",
  colnames(metagenes_filtered)
)

metagenes_filtered <- as.matrix(metagenes_filtered)

#------------------------------------------------------------------------------
# Remove pathways with zero variance
#------------------------------------------------------------------------------

metagenes_filtered <- metagenes_filtered[
  ,
  apply(metagenes_filtered, 2, var) != 0
]

#------------------------------------------------------------------------------
# Remove low-prevalence pathways
# Present in >15% of samples
#------------------------------------------------------------------------------

metagenes_filtered <- metagenes_filtered[
  ,
  colSums(metagenes_filtered > 0) >
    (0.15 * nrow(metagenes_filtered))
]

#------------------------------------------------------------------------------
# Create METAGENES matrix
#------------------------------------------------------------------------------

METAGENES <- metagenes_filtered

#------------------------------------------------------------------------------
# Total Sum Scaling (TSS) normalization
#------------------------------------------------------------------------------

METAGENES <- METAGENES / rowSums(METAGENES)

#------------------------------------------------------------------------------
# Retain the 500 most variable pathways
#------------------------------------------------------------------------------

vars <- apply(METAGENES, 2, var)

top_n <- 500

top_idx <- order(
  vars,
  decreasing = TRUE
)[1:top_n]

METAGENES <- METAGENES[, top_idx]

#------------------------------------------------------------------------------
# Final quality control
#------------------------------------------------------------------------------

stopifnot(sum(is.na(METAGENES)) == 0)

cat("Metatranscriptomic matrix successfully generated.\n")
cat("Samples :", nrow(METAGENES), "\n")
cat("Pathways:", ncol(METAGENES), "\n")

# =============================================================================
# STEP 4 - Metabolomic data preparation and dataset alignment
# =============================================================================

#------------------------------------------------------------------------------
# Import metabolomic dataset
#------------------------------------------------------------------------------

metabolomics <- read.csv(
  "C:/Users/utente/Desktop/RICERCA/IBS Mintea/datasets/PRJNA812699/metabolomics_data.csv",
  row.names = 1
)

#------------------------------------------------------------------------------
# Harmonize sample identifiers
#------------------------------------------------------------------------------

colnames(metabolomics) <- trimws(toupper(colnames(metabolomics)))

ids <- trimws(toupper(patients$Patient))

#------------------------------------------------------------------------------
# Retain only samples included in the clinical dataset
#------------------------------------------------------------------------------

common_ids <- intersect(ids, colnames(metabolomics))

metab_filtered <- metabolomics[, common_ids]

#------------------------------------------------------------------------------
# Transpose matrix
# Samples as rows, metabolites as columns
#------------------------------------------------------------------------------

metab_filtered <- t(metab_filtered)

#------------------------------------------------------------------------------
# Prefix metabolite features
#------------------------------------------------------------------------------

colnames(metab_filtered) <- paste0(
  "M__",
  colnames(metab_filtered)
)

metab_filtered <- as.matrix(metab_filtered)

#------------------------------------------------------------------------------
# Remove metabolites with zero variance
#------------------------------------------------------------------------------

metab_filtered <- metab_filtered[
  ,
  apply(metab_filtered, 2, var) != 0
]

#------------------------------------------------------------------------------
# Remove low-prevalence metabolites
# Present in >15% of samples
#------------------------------------------------------------------------------

metab_filtered <- metab_filtered[
  ,
  colSums(metab_filtered > 0) >
    (0.15 * nrow(metab_filtered))
]

#------------------------------------------------------------------------------
# Log transformation
#------------------------------------------------------------------------------

metab_filtered <- log1p(metab_filtered)

#------------------------------------------------------------------------------
# Create METABOLOME matrix
#------------------------------------------------------------------------------

METABOLOME <- metab_filtered

#------------------------------------------------------------------------------
# Final alignment across all omics datasets
#------------------------------------------------------------------------------

common_ids <- Reduce(
  intersect,
  list(
    rownames(TAXA),
    rownames(METAGENES),
    rownames(METABOLOME),
    metadata_filtered$Patient
  )
)

TAXA <- TAXA[common_ids, ]
METAGENES <- METAGENES[common_ids, ]
METABOLOME <- METABOLOME[common_ids, ]

metadata_filtered <- metadata_filtered[
  match(common_ids, metadata_filtered$Patient),
]

#------------------------------------------------------------------------------
# Verify sample order
#------------------------------------------------------------------------------

stopifnot(
  all(rownames(TAXA) == rownames(METAGENES))
)

stopifnot(
  all(rownames(TAXA) == rownames(METABOLOME))
)

stopifnot(
  all(rownames(TAXA) == metadata_filtered$Patient)
)

#------------------------------------------------------------------------------
# Remove highly correlated features
#------------------------------------------------------------------------------

remove_high_corr <- function(data, cutoff = 0.99) {
  
  corr <- cor(
    data,
    method = "spearman"
  )
  
  to_remove <- character(0)
  
  for (i in 1:(ncol(corr) - 1)) {
    
    for (j in (i + 1):ncol(corr)) {
      
      if (abs(corr[i, j]) > cutoff) {
        
        to_remove <- c(
          to_remove,
          colnames(data)[j]
        )
        
      }
      
    }
    
  }
  
  data[
    ,
    !colnames(data) %in% unique(to_remove)
  ]
  
}

TAXA <- remove_high_corr(TAXA)

METAGENES <- remove_high_corr(METAGENES)

METABOLOME <- remove_high_corr(METABOLOME)

#------------------------------------------------------------------------------
# Create MintTea input list
#------------------------------------------------------------------------------

data_list <- list(
  T = TAXA,
  G = METAGENES,
  M = METABOLOME
)

#------------------------------------------------------------------------------
# Summary
#------------------------------------------------------------------------------

cat("Final multi-omics datasets successfully prepared.\n")
cat("Samples      :", nrow(TAXA), "\n")
cat("Taxa         :", ncol(TAXA), "\n")
cat("Pathways     :", ncol(METAGENES), "\n")
cat("Metabolites  :", ncol(METABOLOME), "\n")

# =============================================================================
# STEP 5 - Multi-omics integration with MintTea
# =============================================================================

#------------------------------------------------------------------------------
# Preserve non-scaled matrices for downstream visualization
#------------------------------------------------------------------------------

TAXA_plot <- TAXA
METAGENES_plot <- METAGENES
METABOLOME_plot <- METABOLOME

#------------------------------------------------------------------------------
# Scale omics datasets
# Required for sgCCA implemented in MintTea
#------------------------------------------------------------------------------

TAXA <- scale(TAXA)

METAGENES <- scale(METAGENES)

METABOLOME <- scale(METABOLOME)

#------------------------------------------------------------------------------
# Build feature matrix
#------------------------------------------------------------------------------

features <- as.data.frame(
  cbind(
    TAXA,
    METAGENES,
    METABOLOME
  )
)

features <- tibble::rownames_to_column(
  features,
  "Patient"
)

#------------------------------------------------------------------------------
# Prepare metadata
#------------------------------------------------------------------------------

metadata_ready <- metadata_filtered %>%
  dplyr::select(
    Patient,
    Group
  )

metadata_ready$Group <- ifelse(
  metadata_ready$Group == "IBS",
  "disease",
  "healthy"
)

#------------------------------------------------------------------------------
# Create MintTea input table
#------------------------------------------------------------------------------

final_input <- metadata_ready %>%
  left_join(
    features,
    by = "Patient"
  )

stopifnot(
  sum(is.na(final_input$Group)) == 0
)

#------------------------------------------------------------------------------
# Prepare variable names for MintTea
#------------------------------------------------------------------------------

final_input <- final_input %>%
  dplyr::rename(
    study_group = Group
  )

colnames(final_input) <- make.names(
  colnames(final_input)
)

#------------------------------------------------------------------------------
# Run MintTea
#------------------------------------------------------------------------------

results <- MintTea(
  
  final_input,
  
  view_prefixes = c(
    "T",
    "G",
    "M"
  ),
  
  study_group_column = "study_group",
  
  sample_id_column = "Patient",
  
  param_edge_thresholds =0.75,
  
  param_diablo_keepX = 15
  
)

#------------------------------------------------------------------------------
# Extract MintTea results
#------------------------------------------------------------------------------

res <- results[[1]]

#------------------------------------------------------------------------------
# Summary
#------------------------------------------------------------------------------

cat("MintTea analysis completed successfully.\n")
cat("Identified modules:", length(res), "\n")

# =============================================================================
# STEP 6 - MintTea module characterization
# =============================================================================

# =============================================================================
# Detailed MintTea module summary
# =============================================================================

#------------------------------------------------------------------------------
# Collect modules
#------------------------------------------------------------------------------

modules <- setNames(
  res,
  paste0("Module ", seq_along(res))
)

module1 <- modules[["Module 1"]]
module2 <- modules[["Module 2"]]
module3 <- modules[["Module 3"]]
module4 <- modules[["Module 4"]]
module5 <- modules[["Module 5"]]

#------------------------------------------------------------------------------
# Function to summarize one module
#------------------------------------------------------------------------------

summarize_module <- function(
    module,
    module_name,
    dataset,
    group_column = "study_group"
) {
  
  # Check that all module features are available
  missing_features <- setdiff(
    module$features,
    colnames(dataset)
  )
  
  if (length(missing_features) > 0) {
    stop(
      module_name,
      ": missing feature(s): ",
      paste(missing_features, collapse = ", ")
    )
  }
  
  # Extract module features
  module_data <- dataset[
    ,
    module$features,
    drop = FALSE
  ]
  
  # Calculate the first principal component
  pca_model <- prcomp(
    module_data,
    center = TRUE,
    scale. = TRUE
  )
  
  module_score <- pca_model$x[, 1]
  
  # ROC based on the first principal component
  roc_object <- pROC::roc(
    response = dataset[[group_column]],
    predictor = module_score,
    quiet = TRUE
  )
  
  pc1_auc <- as.numeric(
    pROC::auc(roc_object)
  )
  
  pc1_ci <- as.numeric(
    pROC::ci.auc(
      roc_object,
      conf.level = 0.95,
      method = "delong"
    )
  )
  
  # Count features by omics type
  n_taxa <- sum(
    grepl("^T__", module$features)
  )
  
  n_pathways <- sum(
    grepl("^G__", module$features)
  )
  
  n_metabolites <- sum(
    grepl("^M__", module$features)
  )
  
  # Build summary row
  data.frame(
    Module = module_name,
    
    Number_features = length(module$features),
    
    MintTea_AUROC = as.numeric(module$auroc),
    
    PC1_AUROC = pc1_auc,
    
    CI_low = pc1_ci[1],
    
    CI_high = pc1_ci[3],
    
    AUROC_95_CI = sprintf(
      "%.3f (%.3f–%.3f)",
      pc1_auc,
      pc1_ci[1],
      pc1_ci[3]
    ),
    
    Inter_view_corr = as.numeric(
      module$inter_view_corr
    ),
    
    Taxa = n_taxa,
    
    Pathways = n_pathways,
    
    Metabolites = n_metabolites,
    
    stringsAsFactors = FALSE
  )
}

#------------------------------------------------------------------------------
# Generate summary table for all modules
#------------------------------------------------------------------------------

module_summary <- dplyr::bind_rows(
  Map(
    f = function(module, module_name) {
      summarize_module(
        module = module,
        module_name = module_name,
        dataset = final_input,
        group_column = "study_group"
      )
    },
    module = modules,
    module_name = names(modules)
  )
)

#------------------------------------------------------------------------------
# Format numeric columns
#------------------------------------------------------------------------------

module_summary <- module_summary %>%
  dplyr::mutate(
    MintTea_AUROC = round(MintTea_AUROC, 3),
    PC1_AUROC = round(PC1_AUROC, 3),
    CI_low = round(CI_low, 3),
    CI_high = round(CI_high, 3),
    Inter_view_corr = round(Inter_view_corr, 3)
  )

print(module_summary)

# =============================================================================
# PCA analysis of MintTea modules
# =============================================================================

#------------------------------------------------------------------------------
# Function to compute PCA for a MintTea module
#------------------------------------------------------------------------------

compute_module_pca <- function(
    module,
    dataset,
    group_column = "study_group"
) {
  
  module_data <- dataset[
    ,
    module$features,
    drop = FALSE
  ]
  
  stopifnot(
    ncol(module_data) > 1,
    !anyNA(module_data)
  )
  
  pca_model <- prcomp(
    module_data,
    center = TRUE,
    scale. = TRUE
  )
  
  explained_variance <- summary(
    pca_model
  )$importance[2, ]
  
  pca_scores <- data.frame(
    PC1 = pca_model$x[, 1],
    PC2 = pca_model$x[, 2],
    Group = dataset[[group_column]]
  )
  
  list(
    pca = pca_model,
    scores = pca_scores,
    variance = explained_variance
  )
}

#------------------------------------------------------------------------------
# Compute PCA for selected MintTea modules
#------------------------------------------------------------------------------

module1_pca <- compute_module_pca(
  module = module1,
  dataset = final_input
)

module2_pca <- compute_module_pca(
  module = module2,
  dataset = final_input
)

module4_pca <- compute_module_pca(
  module = module4,
  dataset = final_input
)

#------------------------------------------------------------------------------
# Extract PC1 module scores for downstream analyses
#------------------------------------------------------------------------------

module1_score <- module1_pca$scores$PC1

df_module1 <- data.frame(
  score = module1_score,
  group = module1_pca$scores$Group
)

# =============================================================================
# ROC analysis - MintTea Module 1
# =============================================================================

#------------------------------------------------------------------------------
# Compute ROC curve using Module 1 PC1 score
#------------------------------------------------------------------------------

roc_module1 <- pROC::roc(
  response = final_input$study_group,
  predictor = module1_score,
  levels = c("healthy", "disease"),
  direction = "auto",
  quiet = TRUE
)

#------------------------------------------------------------------------------
# AUROC and 95% confidence interval
#------------------------------------------------------------------------------

auc_module1 <- as.numeric(
  pROC::auc(roc_module1)
)

ci_module1 <- as.numeric(
  pROC::ci.auc(
    roc_module1,
    conf.level = 0.95,
    method = "delong"
  )
)

cat(
  "\nModule 1 AUROC:",
  round(auc_module1, 3),
  "\n95% CI:",
  paste0(
    round(ci_module1[1], 3),
    "–",
    round(ci_module1[3], 3)
  ),
  "\n"
)

# =============================================================================
# ROC analysis - MintTea Module 1
# =============================================================================

#------------------------------------------------------------------------------
# Compute ROC curve using Module 1 PC1 score
#------------------------------------------------------------------------------

roc_module1 <- pROC::roc(
  response = df_module1$group,
  predictor = df_module1$score,
  quiet = TRUE
)

#------------------------------------------------------------------------------
# AUROC and 95% confidence interval
#------------------------------------------------------------------------------

auc_module1 <- as.numeric(
  pROC::auc(roc_module1)
)

ci_module1 <- as.numeric(
  pROC::ci.auc(
    roc_module1,
    conf.level = 0.95,
    method = "delong"
  )
)

cat(
  "\nModule 1 AUROC:",
  round(auc_module1, 3),
  "\n95% CI:",
  paste0(
    round(ci_module1[1], 3),
    "–",
    round(ci_module1[3], 3)
  ),
  "\n"
)

#------------------------------------------------------------------------------
# Plot ROC curve
#------------------------------------------------------------------------------

plot(
  roc_module1,
  col = "black",
  lwd = 2,
  main = ""
)

#------------------------------------------------------------------------------
# Function to plot module PCA
#------------------------------------------------------------------------------

plot_module_pca <- function(
    pca_results,
    title
) {
  
  ggplot(
    pca_results$scores,
    aes(
      x = PC1,
      y = PC2,
      colour = Group
    )
  ) +
    geom_point(
      size = 3,
      alpha = 0.8
    ) +
    stat_ellipse(
      aes(fill = Group),
      geom = "polygon",
      alpha = 0.20,
      colour = NA
    ) +
    theme_classic(base_size = 14) +
    labs(
      title = title,
      x = paste0(
        "PC1 (",
        round(
          100 * pca_results$variance[1],
          1
        ),
        "% variance)"
      ),
      y = paste0(
        "PC2 (",
        round(
          100 * pca_results$variance[2],
          1
        ),
        "% variance)"
      )
    ) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold"
      )
    )
}

#------------------------------------------------------------------------------
# Generate PCA plots
#------------------------------------------------------------------------------

plot_module_pca(
  module1_pca,
  "PCA - MintTea Module 1"
)

plot_module_pca(
  module2_pca,
  "PCA - MintTea Module 2"
)

plot_module_pca(
  module4_pca,
  "PCA - MintTea Module 4"
)

##############################################################################
# Function: Plot MintTea network
##############################################################################

plot_module_network <- function(module,
                                cutoff = 0.6,
                                layout = "fr") {
  
  #----------------------------------------------------------
  # Extract edges
  #----------------------------------------------------------
  
  edges <- module$module_edges %>%
    dplyr::filter(edge_weight >= cutoff)
  
  #----------------------------------------------------------
  # Build graph
  #----------------------------------------------------------
  
  g <- graph_from_data_frame(
    edges,
    directed = FALSE
  )
  
  #----------------------------------------------------------
  # Node annotation
  #----------------------------------------------------------
  
  V(g)$type <- dplyr::case_when(
    grepl("^T__", V(g)$name) ~ "Taxa",
    grepl("^G__", V(g)$name) ~ "Pathway",
    grepl("^M__", V(g)$name) ~ "Metabolite",
    TRUE ~ "Unknown"
  )
  
  #----------------------------------------------------------
  # Network statistics
  #----------------------------------------------------------
  
  V(g)$degree <- degree(g)
  
  V(g)$betweenness <- betweenness(g)
  
  V(g)$size <- scales::rescale(
    V(g)$degree,
    to = c(17,21)
  )
  
  #----------------------------------------------------------
  # Plot
  #----------------------------------------------------------
  
  p <- ggraph(
    g,
    layout = layout
  ) +
    
    geom_edge_link(
      aes(width = edge_weight),
      alpha = .35,
      colour = "grey70"
    ) +
    
    scale_edge_width(
      range = c(.4,3)
    ) +
    
    geom_node_point(
      aes(
        size = size,
        colour = type
      )
    ) +
    
    geom_node_text(
      aes(
        label = ifelse(
          degree >= 1,
          name,
          ""
        )
      ),
      repel = TRUE,
      size = 2.5
    ) +
    
    scale_colour_manual(
      values = c(
        Taxa = "#4DA3FF",
        Pathway = "#4CAF50",
        Metabolite = "#FF6B6B"
      )
    ) +
    
    # più spazio attorno alla rete
    scale_x_continuous(
      expand = expansion(mult = c(0.08, 0.18))
    ) +
    
    scale_y_continuous(
      expand = expansion(mult = 0.10)
    ) +
    
    # sfondo bianco
    theme_void() +
    
    theme(
      panel.background = element_rect(
        fill = "white",
        colour = NA
      ),
      plot.background = element_rect(
        fill = "white",
        colour = NA
      ),
      legend.background = element_rect(
        fill = "white",
        colour = NA
      ),
      legend.box.background = element_rect(
        fill = "white",
        colour = NA
      ),
      legend.position = "right",
      plot.margin = margin(
        t = 15,
        r = 30,
        b = 15,
        l = 15
      )
    ) +
    
    labs(
      size = "size",
      edge_width = "edge_weight",
      colour = "type"
    )
  
  return(
    list(
      graph = g,
      plot = p
    )
  )
}

network_module1 <- plot_module_network(
  module1
)

network_module1$plot

# =============================================================================
# Permutation ROC analysis
# =============================================================================

run_auc_permutation <- function(
    score,
    group,
    n_perm = 1000,
    seed = 123
) {
  
  set.seed(seed)
  
  # ---------------------------------------------------------------------------
  # Observed ROC curve
  # ---------------------------------------------------------------------------
  
  roc_real <- pROC::roc(
    response = group,
    predictor = score,
    quiet = TRUE
  )
  
  real_auc <- as.numeric(
    pROC::auc(roc_real)
  )
  
  # ---------------------------------------------------------------------------
  # Permuted ROC curves
  # ---------------------------------------------------------------------------
  
  auc_perm <- numeric(n_perm)
  
  roc_perm <- vector(
    mode = "list",
    length = n_perm
  )
  
  for (i in seq_len(n_perm)) {
    
    shuffled_group <- sample(group)
    
    roc_tmp <- pROC::roc(
      response = shuffled_group,
      predictor = score,
      
      # Keep the same class order and direction as the observed ROC
      levels = roc_real$levels,
      direction = roc_real$direction,
      
      quiet = TRUE
    )
    
    roc_perm[[i]] <- roc_tmp
    
    auc_perm[i] <- as.numeric(
      pROC::auc(roc_tmp)
    )
  }
  
  # ---------------------------------------------------------------------------
  # Empirical permutation p-value
  # ---------------------------------------------------------------------------
  
  permutation_p <- (
    sum(auc_perm >= real_auc) + 1
  ) / (
    n_perm + 1
  )
  
  # ---------------------------------------------------------------------------
  # Return results
  # ---------------------------------------------------------------------------
  
  list(
    real_auc = real_auc,
    permutation_auc = auc_perm,
    permutation_rocs = roc_perm,
    permutation_p = permutation_p,
    ci = pROC::ci.auc(roc_real),
    threshold = pROC::coords(
      roc_real,
      x = "best"
    ),
    roc = roc_real
  )
}


# =============================================================================
# Run 1000 permutations
# =============================================================================

perm_module1 <- run_auc_permutation(
  score = module1_score,
  group = final_input$study_group,
  n_perm = 1000,
  seed = 123
)

# =============================================================================
# Plot observed and permuted ROC curves
# =============================================================================

plot_permutation_roc <- function(result) {
  
  # Create an empty ROC plotting area
  plot(
    1,
    type = "n",
    xlim = c(1, 0),
    ylim = c(0, 1),
    xlab = "Specificity",
    ylab = "Sensitivity",
    main = paste0(
      "Observed ROC and permutation ROC curves\n",
      "AUC = ",
      round(result$real_auc, 2),
      "; 95% CI ",
      round(result$ci[1], 2),
      "–",
      round(result$ci[3], 2),
      "; permutation p = ",
      format.pval(
        result$permutation_p,
        digits = 2
      )
    )
  )
  
  # Add the 1000 permuted ROC curves in grey
  for (roc_i in result$permutation_rocs) {
    
    pROC::plot.roc(
      roc_i,
      add = TRUE,
      col = grDevices::adjustcolor(
        "grey60",
        alpha.f = 0.10
      ),
      lwd = 1
    )
  }
  
  # Add the no-discrimination diagonal
  abline(
    a = 1,
    b = -1,
    lty = 2,
    col = "grey30"
  )
  
  # Add the observed ROC curve in red
  pROC::plot.roc(
    result$roc,
    add = TRUE,
    col = "red",
    lwd = 3
  )
  
  legend(
    "bottomright",
    legend = c(
      paste0(
        "Observed ROC: AUC = ",
        round(result$real_auc, 2)
      ),
      paste0(
        length(result$permutation_rocs),
        " permuted ROC curves"
      )
    ),
    col = c(
      "red",
      "grey60"
    ),
    lwd = c(
      3,
      1
    ),
    bty = "n"
  )
}


# =============================================================================
# Display plot
# =============================================================================

plot_permutation_roc(
  perm_module1
)

# =============================================================================
# STEP 7 - LASSO feature selection
# =============================================================================

#------------------------------------------------------------------------------
# Prepare predictor matrix
#------------------------------------------------------------------------------

features_mod1 <- module1$features

X <- final_input[, features_mod1, drop = FALSE]

Y <- ifelse(
  final_input$study_group == "disease",
  1,
  0
)

X <- as.matrix(X)

stopifnot(
  
  nrow(X) == length(Y),
  
  !anyNA(X),
  
  all(colnames(X) %in% features_mod1)
  
)

cat(
  "\nPredictors:", ncol(X),
  "\nSamples:", nrow(X), "\n"
)

##############################################################################
# Function: repeated LASSO cross-validation
##############################################################################

run_lasso_cv <- function(X,
                         Y,
                         seed = 123,
                         k = 5,
                         repeats = 10){
  
  set.seed(seed)
  
  folds <- createMultiFolds(
    Y,
    k = k,
    times = repeats
  )
  
  auc_values <- numeric(length(folds))
  
  selected_features <- vector(
    "list",
    length(folds)
  )
  
  i <- 1
  
  for(train_idx in folds){
    
    test_idx <- setdiff(
      seq_len(length(Y)),
      train_idx
    )
    
    model <- cv.glmnet(
      
      x = X[train_idx, ],
      
      y = Y[train_idx],
      
      family = "binomial",
      
      alpha = 1,
      
      type.measure = "auc",
      
      nfolds = 5
      
    )
    
    coef_model <- coef(
      model,
      s = "lambda.1se"
    )
    
    selected <- rownames(coef_model)[
      coef_model[,1] != 0
    ]
    
    selected <- setdiff(
      selected,
      "(Intercept)"
    )
    
    selected_features[[i]] <- selected
    
    pred <- predict(
      
      model,
      
      newx = X[test_idx, ],
      
      s = "lambda.1se",
      
      type = "response"
      
    )
    
    auc_values[i] <- as.numeric(
      
      auc(
        
        roc(
          
          Y[test_idx],
          
          as.vector(pred),
          
          quiet = TRUE
          
        )
        
      )
      
    )
    
    i <- i + 1
    
  }
  
  list(
    
    auc = auc_values,
    
    features = selected_features
    
  )
  
}
lasso_cv <- run_lasso_cv(
  X,
  Y
)

mean_auc <- mean(
  lasso_cv$auc
)

sd_auc <- sd(
  lasso_cv$auc
)

feature_frequency <-
  
  sort(
    
    table(
      
      unlist(
        lasso_cv$features
      )
      
    ),
    
    decreasing = TRUE
    
  )

cat(
  
  "\nMean AUC:",
  
  round(mean_auc,3),
  
  "\nSD:",
  
  round(sd_auc,3),
  
  "\n"
  
)

feature_frequency


##############################################################################
# Final LASSO model
##############################################################################

final_lasso <- cv.glmnet(
  
  x = X,
  
  y = Y,
  
  family = "binomial",
  
  alpha = 1,
  
  type.measure = "auc",
  
  nfolds = 5
  
)

plot(final_lasso)

coef_final <- coef(
  final_lasso,
  s = "lambda.1se"
)

selected_final <-
  
  data.frame(
    
    feature = rownames(coef_final),
    
    coefficient = as.numeric(coef_final)
    
  ) |>
  
  filter(
    
    coefficient != 0,
    
    feature != "(Intercept)"
    
  ) |>
  
  arrange(
    
    desc(abs(coefficient))
    
  )

selected_final

##############################################################################
# Function: evaluate LASSO model
##############################################################################

evaluate_lasso <- function(model,
                           X,
                           Y){
  
  pred <- predict(
    
    model,
    
    newx = X,
    
    s = "lambda.1se",
    
    type = "response"
    
  )
  
  roc_obj <- roc(
    
    Y,
    
    as.vector(pred),
    
    quiet = TRUE
    
  )
  
  list(
    
    roc = roc_obj,
    
    auc = auc(roc_obj),
    
    ci = ci.auc(roc_obj),
    
    threshold = coords(
      roc_obj,
      "best"
    )
    
  )
  
}

lasso_perf <- evaluate_lasso(
  
  final_lasso,
  
  X,
  
  Y
  
)

plot(lasso_perf$roc)

lasso_perf$auc

lasso_perf$ci

lasso_perf$threshold

##############################################################################
# Function: single-feature ROC analysis with 95% CI
##############################################################################

compute_single_feature_auc <- function(
    features,
    dataset,
    outcome
) {
  
  dplyr::bind_rows(
    
    lapply(
      
      features,
      
      function(f) {
        
        roc_obj <- pROC::roc(
          response = outcome,
          predictor = dataset[[f]],
          quiet = TRUE
        )
        
        auc_value <- as.numeric(
          pROC::auc(roc_obj)
        )
        
        auc_ci <- as.numeric(
          pROC::ci.auc(
            roc_obj,
            conf.level = 0.95,
            method = "delong"
          )
        )
        
        data.frame(
          feature = f,
          AUC = auc_value,
          CI_low = auc_ci[1],
          CI_high = auc_ci[3],
          AUC_95_CI = sprintf(
            "%.3f (%.3f–%.3f)",
            auc_value,
            auc_ci[1],
            auc_ci[3]
          ),
          AUC_abs = max(
            auc_value,
            1 - auc_value
          ),
          stringsAsFactors = FALSE
        )
      }
    )
  ) |>
    
    dplyr::arrange(
      dplyr::desc(AUC_abs)
    )
}

single_auc <- compute_single_feature_auc(
  
  selected_final$feature,
  
  final_input,
  
  final_input$study_group
  
)

single_auc

# =============================================================================
# STEP 8 - IBS reference parameters
# =============================================================================

#------------------------------------------------------------------------------
# Extract LASSO signature from IBS cohort
#------------------------------------------------------------------------------

IBS_signature <- cbind(
  
  TAXA_plot[, c(
    "T__Agathobacter",
    "T__Colidextribacter",
    "T__Monoglobus"
  )],
  
  METABOLOME_plot[, c(
    "M__salicylate",
    "M__xylose",
    "M__acesulfame",
    "M__2-keto-3-deoxy-gluconate"
  )]
  
)

IBS_signature <- as.data.frame(IBS_signature)

IBS_signature[] <- lapply(
  IBS_signature,
  function(x) as.numeric(as.character(x))
)

#------------------------------------------------------------------------------
# Compute reference scaling parameters
#------------------------------------------------------------------------------

IBS_means <- colMeans(
  IBS_signature,
  na.rm = TRUE
)

IBS_sds <- apply(
  IBS_signature,
  2,
  sd,
  na.rm = TRUE
)

summary_stats <- data.frame(
  Feature = names(IBS_means),
  Mean = IBS_means,
  SD = IBS_sds,
  row.names = NULL
)

print(summary_stats)

# =============================================================================
# IBDMDB DATASET PREPARATION
# =============================================================================

#------------------------------------------------------------------------------
# Load metagenomic taxonomic profiles
#------------------------------------------------------------------------------

tax_profiles <- read.delim(
  "C:/Users/utente/Desktop/RICERCA/IBS Mintea/datasets/IBDMDB/taxonomic_profiles_3.tsv.gz",
  sep = "\t",
  header = FALSE,
  comment.char = "#",
  stringsAsFactors = FALSE
)

sample_names <- as.character(tax_profiles[1, -1])

tax_profiles <- tax_profiles[-1, ]

colnames(tax_profiles) <- c("Taxa", sample_names)

#------------------------------------------------------------------------------
# Samples with both metagenomic and metabolomic data
#------------------------------------------------------------------------------

wanted_codes <- c(
  
  "CSM5FZ3N_P",
  "CSM5FZ3T_P",
  "CSM5MCU4_P",
  "CSM5MCWK_P",
  "CSM5MCXF_P",
  "CSM67U9V_P",
  "CSM67U9X_P",
  "CSM79HJI_P",
  "ESM5MEDZ_P",
  "ESM5MEB9_P",
  "HSM5MD7Z_P",
  "HSM5MD3L_P",
  "HSM6XRQB_P",
  "HSM7J4LP",
  "MSM5LLDI",
  "MSM5LLIC_P",
  "MSM5LLIS_P",
  "MSM5LLFG_P",
  "MSM9VZM4_P",
  "PSM6XBQM_P",
  "PSM6XBRK_P",
  "CSM67UH7",
  "HSM5MD8A_P",
  "HSM5MD82_P",
  "HSM5MD8H_P",
  "HSM5MD5Z_P",
  "HSM6XRQC_P",
  "HSM67VDX_P",
  "HSM67VDR_P",
  "HSMA33NY",
  "CSM6J2H9_P",
  "MSM6J2JF_P",
  "MSM6J2JH_P",
  "MSM6J2PO",
  "MSM6J2Q1",
  "MSM79H8D",
  "MSM79H94_P",
  "MSM9VZFJ_P",
  "MSM9VZMM",
  "CSM5FZ4A_P",
  "CSM5MCTZ_P",
  "CSM5MCXB_P",
  "CSM5MCYU_P",
  "CSM67U9T_P",
  "CSM79HQR_P",
  "HSM5MD87_P",
  "HSM6XRQE_P",
  "HSM7J4JT_P",
  "HSM7J4JV_P",
  "HSMA33NW",
  "MSM79HD6_P",
  "PSM6XBSE_P"
  
)

#------------------------------------------------------------------------------
# Match sample IDs
#------------------------------------------------------------------------------

matched_cols <- unlist(
  lapply(
    wanted_codes,
    function(x)
      grep(
        paste0("^", x),
        colnames(tax_profiles),
        value = TRUE
      )
  )
)

#------------------------------------------------------------------------------
# Filter taxonomic profiles
#------------------------------------------------------------------------------

tax_profiles_selected <- tax_profiles[, c("Taxa", matched_cols)]

dim(tax_profiles_selected)
# 932 × 53

length(matched_cols)
# 52 matched samples

setdiff(wanted_codes, matched_cols)
# Missing samples (if any)

# =============================================================================
# Filter metagenomic data at genus level
# =============================================================================

genus_df <- tax_profiles_selected[
  grepl("\\|g__[^|]+$", tax_profiles_selected$Taxa),
]

genus_df$Genus <- sub(
  ".*g__",
  "",
  genus_df$Taxa
)

rownames(genus_df) <- genus_df$Genus

genus_mat <- genus_df[
  ,
  !colnames(genus_df) %in% c("Taxa", "Genus")
]

genus_mat[] <- lapply(
  genus_mat,
  as.numeric
)

rownames(genus_mat) <- rownames(genus_df)

dim(genus_mat)
# 187 × 52

# =============================================================================
# Transpose genus matrix
# =============================================================================

genus_mat_t <- as.data.frame(
  t(genus_mat)
)

dim(genus_mat_t)
# 52 × 187

# =============================================================================
# Remove zero-variance genera
# =============================================================================

genus_mat_t <- genus_mat_t[
  ,
  apply(genus_mat_t, 2, var) != 0
]

# =============================================================================
# Convert relative abundances (%) to fractions
# =============================================================================

genus_mat_tss <- genus_mat_t / 100

rowSums(genus_mat_tss)[1:5]

# =============================================================================
# Extract IBS signature taxa
# =============================================================================

selected_genus_IBD <- genus_mat_t[
  ,
  c(
    "Monoglobus",
    "Lachnospiraceae_unclassified",
    "Oscillibacter"
  )
]

dim(selected_genus_IBD)
# 52 × 3

# =============================================================================
# Load metagenomic metadata
# =============================================================================

metadata_metagenomics <- read.delim(
  "C:/Users/utente/Desktop/RICERCA/IBS Mintea/datasets/IBDMDB/metadata_metagenomics.csv",
  sep = ";",
  stringsAsFactors = FALSE
)

metadata_metagenomics <- metadata_metagenomics[
  ,
  c(
    "SampleID",
    "diagnosis",
    "classification"
  )
]

metadata_metagenomics$SampleID <- trimws(
  metadata_metagenomics$SampleID
)

# =============================================================================
# Merge metadata and taxonomic profiles
# =============================================================================

selected_genus_IBD$SampleID <- rownames(selected_genus_IBD)

selected_genus_IBD$SampleID <- gsub(
  "_profile$",
  "",
  selected_genus_IBD$SampleID
)

selected_genus_IBD$SampleID <- trimws(
  selected_genus_IBD$SampleID
)


IBD_final <- dplyr::inner_join(
  metadata_metagenomics,
  selected_genus_IBD,
  by = "SampleID"
)

dim(IBD_final)
# 52 × 6

# =============================================================================
# Fecal metabolomics
# =============================================================================

library(reticulate)

py_require("biom-format")


py_config()

biom <- reticulate::import("biom")

table_biom <- biom$load_table(
  "C:/Users/utente/Desktop/RICERCA/IBS Mintea/datasets/IBDMDB/HMP2_metabolomics_w_metadata.biom"
)

metabolome_mat <- py_to_r(
  table_biom$to_dataframe(dense = TRUE)
)

dim(metabolome_mat)
# 81867 × 546

# =============================================================================
# Metabolite annotation
# =============================================================================

obs_md <- table_biom$metadata(axis = "observation")

meta_df <- data.frame(
  feature = rownames(metabolome_mat),
  
  metabolite = sapply(
    obs_md,
    function(x) ifelse(is.null(x$Metabolite), NA, x$Metabolite)
  ),
  
  hmdb = sapply(
    obs_md,
    function(x)
      ifelse(
        is.null(x$`HMDB (*Representative ID)`),
        NA,
        x$`HMDB (*Representative ID)`
      )
  ),
  
  mz = sapply(
    obs_md,
    function(x)
      ifelse(
        is.null(x$`m/z`),
        NA,
        x$`m/z`
      )
  ),
  
  rt = sapply(
    obs_md,
    function(x)
      ifelse(
        is.null(x$RT),
        NA,
        x$RT
      )
  ),
  
  method = sapply(
    obs_md,
    function(x)
      ifelse(
        is.null(x$Method),
        NA,
        x$Method
      )
  ),
  
  stringsAsFactors = FALSE
)

# =============================================================================
# Extract metabolites of interest
# =============================================================================

met_interest <- data.frame(
  
  patient_id = colnames(metabolome_mat),
  
  salicylate = as.numeric(
    metabolome_mat["C18n_QI66", ]
  ),
  
  acesulfame = as.numeric(
    metabolome_mat["HILn_QI21", ]
  ),
  
  xylose = as.numeric(
    metabolome_mat["HILn_QI132", ]
  ),
  
  putative_KDG = as.numeric(
    metabolome_mat["HILn_QI26540", ]
  )
  
)

# =============================================================================
# Half-minimum imputation
# =============================================================================

for(i in 2:ncol(met_interest)){
  
  min_val <- min(
    met_interest[met_interest[, i] > 0, i],
    na.rm = TRUE
  )
  
  met_interest[
    is.na(met_interest[, i]),
    i
  ] <- min_val / 2
  
}

# =============================================================================
# Log transformation
# =============================================================================

met_interest[, 2:ncol(met_interest)] <- log1p(
  met_interest[, 2:ncol(met_interest)]
)

# =============================================================================
# Keep patients with both metabolomics and metagenomics
# =============================================================================

patients_interest_clean <- sub(
  "_P$",
  "",
  wanted_codes
)

met_interest_filtered <- subset(
  met_interest,
  patient_id %in% patients_interest_clean
)

dim(met_interest_filtered)
# 51 × 5

setdiff(
  patients_interest_clean,
  met_interest$patient_id
)
# Missing sample:
# CSM6J2H9

# =============================================================================
# FINAL MULTI-OMICS DATASET (IBD)
# =============================================================================

#------------------------------------------------------------------------------
# Harmonize sample identifiers
#------------------------------------------------------------------------------

colnames(met_interest_filtered)[1] <- "SampleID"

IBD_final$SampleID <- sub(
  "_P$",
  "",
  IBD_final$SampleID
)

#------------------------------------------------------------------------------
# Merge metagenomic and metabolomic data
#------------------------------------------------------------------------------

final_multiomics <- inner_join(
  IBD_final,
  met_interest_filtered,
  by = "SampleID"
)

final_multiomics <- final_multiomics[
  ,
  !colnames(final_multiomics) %in% "patient_id"
]

dim(final_multiomics)
# 51 × 10

#------------------------------------------------------------------------------
# Extract IBS signature features
#------------------------------------------------------------------------------

IBD_signature <- final_multiomics[
  ,
  c(
    "Monoglobus",
    "Lachnospiraceae_unclassified",
    "Oscillibacter",
    "salicylate",
    "acesulfame",
    "xylose",
    "putative_KDG"
  )
]

colnames(IBD_signature) <- c(
  "T__Monoglobus",
  "T__Agathobacter",
  "T__Colidextribacter",
  "M__salicylate",
  "M__acesulfame",
  "M__xylose",
  "M__2-keto-3-deoxy-gluconate"
)

#------------------------------------------------------------------------------
# Convert to numeric
#------------------------------------------------------------------------------

IBD_signature[] <- lapply(
  IBD_signature,
  function(x) as.numeric(as.character(x))
)

#------------------------------------------------------------------------------
# Standardize using IBS reference parameters
#------------------------------------------------------------------------------

IBD_signature_scaled <- scale(
  IBD_signature,
  center = IBS_means,
  scale = IBS_sds
)

IBD_signature_scaled <- as.data.frame(
  IBD_signature_scaled
)

#------------------------------------------------------------------------------
# Cap extreme z-scores
#------------------------------------------------------------------------------

IBD_signature_scaled[
  IBD_signature_scaled > 5
] <- 5

IBD_signature_scaled[
  IBD_signature_scaled < -5
] <- -5

head(IBD_signature_scaled)

# =============================================================================
# Cross-disease transferability analysis
# =============================================================================

#------------------------------------------------------------------------------
# IBS-LASSO signature coefficients
#------------------------------------------------------------------------------

IBS_coefficients <- c(
  
  "M__salicylate" = -0.427130347,
  "T__Agathobacter" = -0.135490129,
  "M__xylose" = -0.128460994,
  "T__Colidextribacter" = -0.093303406,
  "M__acesulfame" = -0.012623396,
  "T__Monoglobus" = -0.012082215,
  "M__2.keto.3.deoxy.gluconate" = -0.008418684
)

#------------------------------------------------------------------------------
# Harmonize feature names
#------------------------------------------------------------------------------

colnames(IBD_signature_scaled) <- gsub(
  "M__2-keto-3-deoxy-gluconate",
  "M__2.keto.3.deoxy.gluconate",
  colnames(IBD_signature_scaled)
)

setdiff(
  names(IBS_coefficients),
  colnames(IBD_signature_scaled)
)

#------------------------------------------------------------------------------
# Calculate IBS-like score
#------------------------------------------------------------------------------

IBS_score_IBD <- as.matrix(
  
  IBD_signature_scaled[
    ,
    names(IBS_coefficients)
  ]
  
) %*% IBS_coefficients

IBS_score_IBD <- as.numeric(
  IBS_score_IBD
)


#------------------------------------------------------------------------------
# Add score to final dataset
#------------------------------------------------------------------------------

final_multiomics$IBS_score <- IBS_score_IBD

#------------------------------------------------------------------------------
# Group comparison
#------------------------------------------------------------------------------

aggregate(
  IBS_score ~ classification,
  data = final_multiomics,
  mean
)

aggregate(
  IBS_score ~ classification,
  data = final_multiomics,
  median
)


wilcox.test(
  IBS_score ~ classification,
  data = final_multiomics
)

#------------------------------------------------------------------------------
# Visualization
#------------------------------------------------------------------------------

ggplot(
  
  final_multiomics,
  
  aes(
    x = classification,
    y = IBS_score,
    fill = classification
  )
  
) +
  
  geom_boxplot(
    outlier.shape = NA
  ) +
  
  geom_jitter(
    width = 0.15,
    alpha = 0.5
  ) +
  
  theme_classic(base_size = 14) +
  
  labs(
    title = "Cross-disease transferability of the IBS signature",
    x = "",
    y = "IBS-like score"
  )


# =============================================================================
# ROC analysis - IBD validation cohort
# =============================================================================

roc_obj <- roc(
  
  final_multiomics$classification,
  
  final_multiomics$IBS_score
)

auc(roc_obj)
pROC::ci.auc(roc_obj, conf.level = 0.95, method = "delong")

plot(
  roc_obj,
  main = "ROC curve - IBS score applied to IBD cohort"
)

# =============================================================================
# Calculate IBS score in the discovery cohort
# =============================================================================

IBS_panel <- data.frame(
  
  M__salicylate =
    METABOLOME_plot[, "M__salicylate"],
  
  T__Agathobacter =
    TAXA_plot[, "T__Agathobacter"],
  
  M__xylose =
    METABOLOME_plot[, "M__xylose"],
  
  T__Colidextribacter =
    TAXA_plot[, "T__Colidextribacter"],
  
  M__acesulfame =
    METABOLOME_plot[, "M__acesulfame"],
  
  T__Monoglobus =
    TAXA_plot[, "T__Monoglobus"],
  
  `M__2-keto-3-deoxy-gluconate` =
    METABOLOME_plot[, "M__2-keto-3-deoxy-gluconate"]
)

# Harmonize metabolite name

colnames(IBS_panel)[
  colnames(IBS_panel) == "M__2.keto.3.deoxy.gluconate"
] <- "M__2-keto-3-deoxy-gluconate"


# Force the same feature order used in the external validation

IBS_panel <- IBS_panel[
  ,
  colnames(IBS_signature)
]


# Safety check

stopifnot(
  identical(
    colnames(IBS_panel),
    colnames(IBD_signature)
  )
)

# Convert to numeric

IBS_panel <- data.frame(
  lapply(
    IBS_panel,
    as.numeric
  )
)

# Standardization using IBS reference parameters

IBS_panel_scaled <- scale(
  IBS_panel,
  center = IBS_means,
  scale = IBS_sds
)

IBS_panel_scaled <- as.data.frame(
  IBS_panel_scaled
)

# Cap extreme values

IBS_panel_scaled[
  IBS_panel_scaled > 5
] <- 5

IBS_panel_scaled[
  IBS_panel_scaled < -5
] <- -5

# Calculate IBS score

IBS_score_original <- as.matrix(
  
  IBS_panel_scaled[
    ,
    names(IBS_coefficients)
  ]
  
) %*% IBS_coefficients

IBS_score_original <- as.numeric(
  IBS_score_original
)

metadata_filtered$IBS_score <- IBS_score_original

# =============================================================================
# IBS discovery cohort
# =============================================================================

aggregate(
  IBS_score ~ Group,
  data = metadata_filtered,
  mean
)

aggregate(
  IBS_score ~ Group,
  data = metadata_filtered,
  median
)


wilcox.test(
  IBS_score ~ Group,
  data = metadata_filtered
)

ggplot(
  metadata_filtered,
  aes(
    x = Group,
    y = IBS_score,
    fill = Group
  )
) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.5) +
  theme_classic(base_size = 14) +
  labs(
    title = "IBS score in IBS discovery cohort",
    x = "",
    y = "IBS score"
  )

ggplot(
  metadata_filtered,
  aes(
    x = Group,
    y = IBS_score,
    fill = Group
  )
) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.5) +
  scale_fill_discrete(
    labels = c("HEALTHY", "IBS")
  ) +
  scale_x_discrete(
    labels = c("HEALTHY", "IBS")
  ) +
  theme_classic(base_size = 14)

summary(metadata_filtered$IBS_score)

aggregate(
  IBS_score ~ Group,
  data = metadata_filtered,
  summary
)

apply(
  IBS_panel_scaled,
  2,
  function(x) sum(x == 5 | x == -5)
)

apply(
  IBS_panel_scaled,
  2,
  function(x) sum(x == 5 | x == -5)
)


# =============================================================================
# ROC analysis - IBS discovery cohort
# =============================================================================

metadata_filtered$IBS_binary <- ifelse(
  metadata_filtered$Group == "IBS",
  1,
  0
)

roc_IBS <- roc(
  metadata_filtered$IBS_binary,
  metadata_filtered$IBS_score
)

auc(roc_IBS)
pROC::ci.auc(roc_IBS, conf.level = 0.95, method = "delong")

plot(
  roc_IBS,
  main = "ROC curve - IBS discovery cohort"
)


# =============================================================================
# Dietary subgroup analysis
# =============================================================================

metadata_filtered$Diet_group <- NA

metadata_filtered$Diet_group[
  metadata_filtered$Diet_Pattern %in%
    c(
      "Modified_American",
      "Standard_American",
      "Mediterranean"
    )
] <- "Standard"

metadata_filtered$Diet_group[
  metadata_filtered$Diet_Pattern %in%
    c(
      "Lactose_free",
      "Gluten_free",
      "FODMAPS"
    )
] <- "Restrictive"

metadata_filtered$Diet_group[
  metadata_filtered$Diet_Pattern %in%
    c(
      "Vegetarian",
      "Vegan",
      "Pescetarian",
      "Paleo"
    )
] <- "Other"

kruskal.test(
  IBS_score ~ Diet_group,
  data = metadata_filtered
)


model_logistic_crude <- glm(
  IBS_binary ~ IBS_score,
  family = binomial,
  data = metadata_filtered
)

summary(model_logistic_crude)

exp(
  cbind(
    OR = coef(model_logistic_crude),
    confint(model_logistic_crude)
  )
)
# =============================================================================
# Multivariable models
# =============================================================================


metadata_filtered$Diet_group <- relevel(
  factor(metadata_filtered$Diet_group),
  ref = "Standard"
)

model_logistic <- glm(
  IBS_binary ~ IBS_score + Age + BMI + Diet_group,
  family = binomial,
  data = metadata_filtered
)

summary(model_logistic)

exp(
  cbind(
    OR = coef(model_logistic),
    confint(model_logistic)
  )
)

# =============================================================================
# Multicollinearity assessment
# =============================================================================

vif(model_logistic)

# =============================================================================
# LASSO Bootstrap 
# =============================================================================

set.seed(123)

n_boot <- 1000

selected_features <- list()

for(i in 1:n_boot){
  
  # bootstrap sample
  idx <- sample(
    1:nrow(X),
    size = nrow(X),
    replace = TRUE
  )
  
  X_boot <- X[idx, ]
  Y_boot <- Y[idx]
  
  # LASSO
  model <- cv.glmnet(
    x = X_boot,
    y = Y_boot,
    family = "binomial",
    alpha = 1,
    type.measure = "auc",
    nfolds = 5,
    maxit = 1e6
  )
  
  coef_model <- coef(
    model,
    s = "lambda.1se"
  )
  
  selected <- rownames(coef_model)[
    coef_model[,1] != 0
  ]
  
  selected <- selected[
    selected != "(Intercept)"
  ]
  
  selected_features[[i]] <- selected
}

all_features <- unlist(selected_features)

freq_table <- data.frame(
  feature = names(table(all_features)),
  frequency = as.numeric(table(all_features))
)

freq_table$selection_frequency <- 
  100 * freq_table$frequency / n_boot

freq_table <- freq_table[
  order(freq_table$selection_frequency,
        decreasing = TRUE),
]

freq_table

coef_matrix <- matrix(
  0,
  nrow = n_boot,
  ncol = ncol(X)
)

colnames(coef_matrix) <- colnames(X)

for(i in 1:n_boot){
  
  idx <- sample(
    1:nrow(X),
    nrow(X),
    replace = TRUE
  )
  
  X_boot <- X[idx, ]
  Y_boot <- Y[idx]
  
  model <- cv.glmnet(
    x = X_boot,
    y = Y_boot,
    family = "binomial",
    alpha = 1,
    type.measure = "auc",
    nfolds = 5,
    maxit = 1e6
  )
  
  coefs <- coef(
    model,
    s = "lambda.1se"
  )
  
  coefs <- as.matrix(coefs)
  
  coef_matrix[i,
              rownames(coefs)[-1]
  ] <- coefs[-1,1]
}

selection_freq <- colMeans(
  coef_matrix != 0
) * 100

mean_coef <- colMeans(
  coef_matrix
)

results <- data.frame(
  feature = names(selection_freq),
  selection_frequency = selection_freq,
  mean_coefficient = mean_coef
)

results <- results[
  order(results$selection_frequency,
        decreasing = TRUE),
]

results

sessionInfo()
