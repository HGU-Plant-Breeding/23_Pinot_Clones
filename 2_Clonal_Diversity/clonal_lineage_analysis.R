#!/usr/bin/env Rscript

# ==============================================================================
# SCRIPT: clonal_lineage_analysis.R
# DESCRIPTION: Reconstructs clonal lineages from Genetic (SNP, SV) and 
#              Epigenetic (CG, CHG, CHH) data.
#
# METHODOLOGY:
#   1. Calculates a Haploid Genomic Relationship Matrix (GRM) adapted for 
#      somatic variants (scaling factor: sum(p(1-p))).
#   2. Derives Distance matrices (D = 1 - G).
#   3. Constructs Dendrograms (UPGMA) and Tanglegrams (Fig 6).
#   4. Generates Heatmaps for all layers ordered by the Genetic tree (Fig S11).
#
# AUTHOR: Paolo Callipo
# DATE: 2025
# ==============================================================================

# --- LIBRARIES ---
suppressPackageStartupMessages({
  library(ape)
  library(vegan)      # For Mantel test
  library(pheatmap)   # For Heatmaps
  library(dendextend) # For Tanglegrams
  library(RColorBrewer)
  
  # Optional: For custom fonts
  if(require(showtext) & require(sysfonts)) {
    font_add_google("Roboto", "roboto")
    showtext_auto()
    FONT_FAMILY <- "roboto"
  } else {
    FONT_FAMILY <- "sans"
  }
})

# --- CONFIGURATION ---
WORKING_DIR <- "."  # Change this if running interactively
INPUT_DIR   <- "Distance_Matrix" # Subdirectory containing input tables

# Input Filenames
FILE_SNP <- file.path(INPUT_DIR, "GT_snp.tsv")
FILE_CG  <- file.path(INPUT_DIR, "CG_SMP_final.txt")
FILE_SV  <- file.path(INPUT_DIR, "GT_sv.tsv")
FILE_CHG <- file.path(INPUT_DIR, "GT_dmc_CHG.tsv")
FILE_CHH <- file.path(INPUT_DIR, "GT_dmc_CHH.tsv")
FILE_GROUPS <- file.path(INPUT_DIR, "clonal_groups.txt")

# Set Working Directory
setwd(WORKING_DIR)

# ==============================================================================
# SECTION 1: CUSTOM FUNCTIONS
# ==============================================================================

#' Compute Haploid Genomic Relationship Matrix (Fast)
#' 
#' Calculates GRM using the formula: G = ZZ' / sum(p(1-p))
#' Handles missing data via centered-mean imputation.
compute_haploid_grm_fast <- function(mat){
  mat <- as.matrix(mat)
  
  # 1. Allele Frequency
  p <- colMeans(mat, na.rm = TRUE)
  
  # Filter monomorphic
  keep <- which(p > 0 & p < 1)
  X <- mat[, keep, drop=FALSE]
  p <- p[keep]
  
  # 2. Center (Z = X - p)
  Z <- sweep(X, 2, p, FUN="-")
  
  # 3. Impute Missing (with centered mean logic)
  impute_vals <- round(p) - p
  na_idx <- which(is.na(Z), arr.ind = TRUE)
  if(nrow(na_idx) > 0){
    Z[na_idx] <- impute_vals[na_idx[, 2]]
  }
  
  # 4. Compute GRM
  G_num <- tcrossprod(Z)
  denom <- sum(p * (1 - p))
  
  return(G_num / denom)
}

#' Load and Process Input Matrices
#' 
#' Reads TSV files and formats them for GRM calculation.
#' Modes: 'genotype' (0/0, 0/1 format) or 'numeric' (0, 1 format)
load_and_process <- function(filepath, mode="numeric") {
  if(!file.exists(filepath)) {
    warning(paste("File not found:", filepath))
    return(NULL)
  }
  
  message(paste("Loading:", basename(filepath)))
  
  if(mode == "genotype") {
    # Read Genotype strings (e.g. "0/0", "0/1")
    raw <- read.table(filepath, header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
    
    # ID creation: CHR_POS
    ids <- paste(raw[,1], raw[,2], sep="_")
    if(any(duplicated(ids))) ids <- make.unique(ids)
    rownames(raw) <- ids
    
    # Drop metadata cols
    geno <- raw[, -c(1,2)]
    
    # Convert to 0/1
    mat <- apply(geno, 2, function(x) {
      x <- as.character(x)
      x[x == "0/0"] <- 0
      x[x %in% c("0/1", "1/0", "1/1")] <- 1
      x[x == "./."] <- NA
      as.numeric(x)
    })
    rownames(mat) <- rownames(geno)
    
  } else {
    # Read Pre-formatted Numeric/Matrix data
    # Assumes Row 1 = ID or First Col = ID
    raw <- read.table(filepath, header=TRUE, row.names=1, check.names=FALSE, na.strings=".")
    mat <- as.matrix(raw)
  }
  
  # Transpose for GRM (Rows=Samples, Cols=Markers)
  return(t(mat))
}

# ==============================================================================
# SECTION 2: DATA LOADING & ANALYSIS (CORE)
# ==============================================================================

# 1. Load Primary Datasets (SNP & CG)
mat_snp <- load_and_process(FILE_SNP, mode="genotype")
mat_cg  <- load_and_process(FILE_CG, mode="numeric")

# 2. Compute GRMs
message("Computing Haploid GRMs...")
G_snp <- compute_haploid_grm_fast(mat_snp)
G_cg  <- compute_haploid_grm_fast(mat_cg)

# 3. Compute Distances (D = 1 - G)
D_snp <- 1 - G_snp
D_cg  <- 1 - G_cg

# 4. Hierarchical Clustering
tree_snp <- hclust(as.dist(D_snp), method="average")
tree_cg  <- hclust(as.dist(D_cg), method="average")

# 5. Statistical Correlation (Mantel Test)
message("Performing Mantel Test (SNP vs CG)...")
mantel_res <- mantel(as.dist(D_snp), as.dist(D_cg), method="pearson", permutations=9999)
print(mantel_res)

# ==============================================================================
# SECTION 3: VISUALIZATION - TANGLEGRAM (FIGURE 6)
# ==============================================================================
message("Generating Figure 6 (Tanglegram)...")

# 1. Load Metadata & Colors
if(file.exists(FILE_GROUPS)) {
  group_data <- read.table(FILE_GROUPS, header=FALSE, sep="\t", col.names=c("Clone", "Group"))
  group_map <- setNames(group_data$Group, group_data$Clone)
  
  # Define Palette (Customize as needed)
  unique_groups <- unique(group_data$Group)
  n_groups <- length(unique_groups)
  group_colors <- setNames(brewer.pal(max(3, n_groups), "Set1")[1:n_groups], unique_groups)
  
} else {
  warning("Clonal groups file not found. Using default colors.")
  group_map <- setNames(rep("Unknown", nrow(D_snp)), rownames(D_snp))
  group_colors <- c("Unknown"="grey50")
}

# 2. Prepare Dendrograms
dend1 <- as.dendrogram(tree_snp)
dend2 <- as.dendrogram(tree_cg)

# 3. Untangle (Step2Side) for best alignment
dend_list <- dendlist(dend1, dend2) %>% untangle(method = "step2side")
dend1 <- dend_list[[1]]
dend2 <- dend_list[[2]]

# 4. Apply Colors to Branches/Labels
align_colors <- function(labs) {
  grps <- group_map[labs]
  grps[is.na(grps)] <- "Unknown"
  group_colors[grps]
}

# Left Tree (SNP)
cols1 <- align_colors(labels(dend1))
labels_colors(dend1) <- cols1
dend1 <- color_branches(dend1, col = cols1) %>% set("branches_lwd", 2)

# Right Tree (CG)
cols2 <- align_colors(labels(dend2))
labels_colors(dend2) <- cols2
dend2 <- color_branches(dend2, col = cols2) %>% set("branches_lwd", 2)

# 5. Plot Output
pdf("Figure6_Tanglegram.pdf", width=10, height=8)
par(mar=c(3,3,3,3), family=FONT_FAMILY)

tanglegram(dend1, dend2,
           main_left = "Genetic (SNP)", 
           main_right = "Epigenetic (CG-DMC)",
           color_lines = cols1, # Color connecting lines by left tree groups
           lab.cex = 0.8, 
           edge.lwd = 2, 
           margin_inner = 7,
           columns_width = c(5, 1, 5),
           fast = FALSE, 
           axes = FALSE)

# Add Legend
legend("top", legend = names(group_colors), fill = group_colors, 
       horiz = TRUE, bty = "n", cex = 0.8, inset = -0.05, xpd=TRUE)

dev.off()

# ==============================================================================
# SECTION 4: SUPPLEMENTARY HEATMAPS (ORDERED BY GENETICS)
# ==============================================================================
message("Generating Supplementary Heatmaps...")

# Define Master Order (based on SNP tree)
snp_order <- tree_snp$order
clones_ordered <- rownames(mat_snp)[snp_order]

# Helper function for consistent heatmaps
plot_ordered_heatmap <- function(mat, title, filename) {
  if(is.null(mat)) return()
  
  # Calculate Distance
  G <- compute_haploid_grm_fast(mat)
  D <- 1 - G
  diag(D) <- 0 # Clean diagonal for visual
  
  # Verify all clones exist
  common_clones <- intersect(clones_ordered, rownames(D))
  
  # Plot
  pheatmap(D[common_clones, common_clones], # Force Order
           cluster_rows = FALSE, 
           cluster_cols = FALSE,
           main = title,
           color = colorRampPalette(c("#2c7fb8","#ffffbf","#d7301f"))(100),
           border_color = NA,
           fontsize = 10,
           filename = filename,
           width = 8, height = 8)
}

# 1. SNP & CG (For reference)
plot_ordered_heatmap(mat_snp, "Genetic Distance (SNP)", "Heatmap_SNP.pdf")
plot_ordered_heatmap(mat_cg, "Epigenetic Distance (CG)", "Heatmap_CG.pdf")

# 2. Structural Variants (SV)
mat_sv <- load_and_process(FILE_SV, mode="genotype")
plot_ordered_heatmap(mat_sv, "Structural Variant Distance (SV)", "FigureS11_Heatmap_SV.pdf")

# 3. CHG Methylation
mat_chg <- load_and_process(FILE_CHG, mode="numeric")
plot_ordered_heatmap(mat_chg, "Epigenetic Distance (CHG)", "FigureS11_Heatmap_CHG.pdf")

# 4. CHH Methylation
mat_chh <- load_and_process(FILE_CHH, mode="numeric")
plot_ordered_heatmap(mat_chh, "Epigenetic Distance (CHH)", "FigureS11_Heatmap_CHH.pdf")

message("Analysis Pipeline Completed Successfully.")
