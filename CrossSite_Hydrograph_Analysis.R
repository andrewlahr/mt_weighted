# =============================================================================
# Cross-Site Analysis of Weighted Hydrographs and β(t) Curves
#
# This script loads per-site FPCA results from WorkingFPCA.R and:
#   1. Assembles β(t) curves and weighted hydrographs across all sites
#   2. Computes pairwise correlation/distance matrices
#   3. Clusters sites by β(t) similarity (hierarchical + optimal k)
#   4. Extracts summary metrics from each β(t) curve
#   5. Runs PCA/NMDS on those metrics
#   6. Produces publication-ready figures
#
# REQUIREMENTS:
#   - Run WorkingFPCA.R first for all sites so that .RData files exist
#     in output/models/
#   - Optionally provide a site_covariates.csv with columns:
#       Stream.Section, drainage_area_km2, elevation_m, latitude, ...
#     to overlay on ordination plots
# =============================================================================

library(tidyverse)
library(patchwork)
library(vegan)       # for NMDS, adonis2
library(cluster)     # for silhouette
library(dendextend)  # for coloured dendrograms
library(corrplot)    # for correlation matrix visualisation
library(ggrepel)     # for non-overlapping labels in the PCA biplot
# (mgcv no longer needed here — smoothing of β(t) happens inside WorkingFPCA.R
#  and the smoothed curve is loaded directly from the saved plot_df object)

# =============================================================================
# 0. USER CONFIGURATION
# =============================================================================

# Path to the directory containing per-site .RData files
model_dir <- "output/models"

# Optional: path to site-level covariates CSV
# Set to NULL if not available
covariate_file <- NULL   # e.g., "data/site_covariates.csv"

# Output directory for cross-site figures
out_dir <- "output/plots/cross_site"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# Site inclusion thresholds
# =============================================================================
# Sites that fail any of these thresholds are EXCLUDED from clustering, PCA,
# NMDS, and all derived figures, but they are still tracked and reported in
# `excluded_sites_summary.csv` and printed to the console so you can see who
# got dropped and why. Set a threshold to NA to disable that check.

min_years        <- 15    # minimum number of years in the model
min_n_sig_fpcs   <- 1     # minimum number of significant FPCs in the reduced model
min_r2_full      <- 0.20  # minimum R² of the full FPC regression
# Set min_r2_full <- NA to skip this check.

# DOY helpers (copied from WorkingFPCA.R for consistency)
doy_to_date <- function(doy) {
  format(as.Date(doy - 1, origin = "2001-01-01"), "%b %d")
}
key_doys   <- c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 365)
key_labels <- doy_to_date(key_doys)

# =============================================================================
# 1. LOAD ALL SITE RESULTS
# =============================================================================

rdata_files <- list.files(model_dir, pattern = "functional_regression_results\\.RData$",
                          full.names = TRUE)

cat(sprintf("Found %d site result files.\n", length(rdata_files)))

# --- 1a. First pass: load the per-site .RData into a structured list --------
# Each .RData file written by WorkingFPCA.R now contains BOTH the full-model
# objects (mod_lm, fpca_out, beta_df, wh_df, beta_boot, beta_perm_mat,
# sim_threshold, n_fpc, fpc_scores, fpc_funs) AND the reduced-model objects
# (mod_lm_reduced, sig_fpc_nums, sig_fpc_cols, plot_df, fpc_funs_reduced,
# var_explained, r2_per_fpc, p_threshold, t_threshold). plot_df already
# contains the smoothed beta(t) and significance flags.
#
# Note: mod_lm and mod_lm_reduced are now `constrained_lm` objects (not lm)
# because WorkingFPCA.R fits them with a one-sided box constraint on S_z.
# They support coef(), summary(), fitted(), and residuals() the same way.

# Selection thresholds — kept here only as fallback values for any older
# .RData files that don't have p_threshold / t_threshold saved inside them.
# Updated to match WorkingFPCA.R defaults (1.6 t, 0.05 p).
default_p_threshold <- 0.05
default_t_threshold <- 1.6

site_results <- lapply(rdata_files, function(f) {
  env <- new.env()
  load(f, envir = env)
  
  # Extract site name from filename
  site_name <- gsub("^LL_|_functional_regression_results\\.RData$", "",
                    basename(f))
  
  cat(sprintf("  Loading: %s\n", site_name))
  
  # --- Pull the reduced-model objects directly (saved by updated WorkingFPCA.R)
  # Fall back gracefully if any are missing (older RData files).
  has_plot_df       <- exists("plot_df", envir = env)
  has_sig_fpc_nums  <- exists("sig_fpc_nums", envir = env)
  
  if (!(has_plot_df && has_sig_fpc_nums)) {
    warning(sprintf(
      "%s: .RData missing reduced-model objects (plot_df, sig_fpc_nums). Skipping. Re-run WorkingFPCA.R for this site.",
      site_name))
    return(NULL)
  }
  
  plot_df       <- env$plot_df
  sig_fpc_nums  <- env$sig_fpc_nums
  beta_smooth   <- plot_df$beta_smooth
  sim_sig_pos   <- plot_df$sim_sig_pos
  sim_sig_neg   <- plot_df$sim_sig_neg
  beta_t_reduced <- plot_df$beta_raw
  
  # Weighted hydrograph (365-length vector)
  wh_vec <- env$wh_df$weighted_anomaly
  
  # Full-model beta(t) from beta_df (365-length, from the FULL model)
  beta_full <- env$beta_df$beta
  
  # Variance explained per significant FPC (read from the saved object)
  var_expl <- env$var_explained[sig_fpc_nums] * 100
  
  # Model fits — both are now constrained_lm objects
  r2_full         <- summary(env$mod_lm)$r.squared
  r2_reduced      <- summary(env$mod_lm_reduced)$r.squared
  
  # Number of years used in the model. fit_constrained_lm stores this as $n;
  # fall back to the residual length if loading from an older lm-based RData.
  n_years <- if (!is.null(env$mod_lm$n)) env$mod_lm$n
  else length(residuals(env$mod_lm))
  
  # Spawner covariate info (S_z constrained <= 0 in the model)
  full_coefs <- summary(env$mod_lm)$coefficients
  if ("S_z" %in% rownames(full_coefs)) {
    S_estimate <- unname(full_coefs["S_z", "Estimate"])
    S_t        <- unname(full_coefs["S_z", "t value"])
    S_p        <- unname(full_coefs["S_z", "Pr(>|t|)"])
    # S is "retained" if it appears in the reduced model's coefficient table
    reduced_coefs <- summary(env$mod_lm_reduced)$coefficients
    S_retained <- "S_z" %in% rownames(reduced_coefs)
  } else {
    S_estimate <- NA_real_
    S_t        <- NA_real_
    S_p        <- NA_real_
    S_retained <- FALSE
  }
  
  list(
    site          = site_name,
    n_years       = n_years,
    years_vec     = if (!is.null(env$years_vec)) env$years_vec else NA,
    fpc_scores    = env$fpc_scores,
    wh_vec        = wh_vec,
    beta_full     = beta_full,
    beta_reduced  = beta_t_reduced,
    beta_smooth   = beta_smooth,
    sim_sig_pos   = sim_sig_pos,
    sim_sig_neg   = sim_sig_neg,
    sig_fpc_nums  = sig_fpc_nums,
    n_sig         = length(sig_fpc_nums),
    var_explained = var_expl,
    r2_full       = r2_full,
    r2_reduced    = r2_reduced,
    sim_threshold = env$sim_threshold,
    S_estimate    = S_estimate,
    S_t           = S_t,
    S_p           = S_p,
    S_retained    = S_retained
  )
})

# Drop sites that failed to load
site_results <- Filter(Negate(is.null), site_results)
names(site_results) <- sapply(site_results, `[[`, "site")
n_loaded <- length(site_results)
cat(sprintf("\nLoaded %d sites successfully.\n", n_loaded))

# =============================================================================
# 1b. APPLY INCLUSION THRESHOLDS
# =============================================================================
# Decide which sites pass and which get excluded, with explicit per-site reasons.
# Excluded sites are tracked in `excluded_sites_summary.csv` and a summary table
# is printed to the console. Kept sites carry on to clustering / PCA / NMDS.

# Build a one-row-per-site bookkeeping table with all the decision inputs
inclusion_audit <- data.frame(
  site         = sapply(site_results, `[[`, "site"),
  n_years      = sapply(site_results, `[[`, "n_years"),
  n_sig_fpcs   = sapply(site_results, `[[`, "n_sig"),
  r2_full      = sapply(site_results, `[[`, "r2_full"),
  stringsAsFactors = FALSE
)

# Apply each threshold (NA threshold = check disabled)
fail_years    <- if (is.na(min_years))      rep(FALSE, nrow(inclusion_audit)) else inclusion_audit$n_years < min_years
fail_n_fpcs   <- if (is.na(min_n_sig_fpcs)) rep(FALSE, nrow(inclusion_audit)) else inclusion_audit$n_sig_fpcs < min_n_sig_fpcs
fail_r2       <- if (is.na(min_r2_full))    rep(FALSE, nrow(inclusion_audit)) else inclusion_audit$r2_full < min_r2_full

inclusion_audit$included <- !(fail_years | fail_n_fpcs | fail_r2)
inclusion_audit$reason   <- mapply(function(fy, ff, fr) {
  reasons <- c()
  if (fy) reasons <- c(reasons, sprintf("n_years < %d", min_years))
  if (ff) reasons <- c(reasons, sprintf("n_sig_fpcs < %d", min_n_sig_fpcs))
  if (fr) reasons <- c(reasons, sprintf("r2_full < %.2f", min_r2_full))
  if (length(reasons) == 0) "" else paste(reasons, collapse = "; ")
}, fail_years, fail_n_fpcs, fail_r2)

n_kept     <- sum(inclusion_audit$included)
n_excluded <- nrow(inclusion_audit) - n_kept

cat("\n=====================================================\n")
cat("INCLUSION FILTER\n")
cat("=====================================================\n")
cat(sprintf("Thresholds:  min_years = %s,  min_n_sig_fpcs = %s,  min_r2_full = %s\n",
            ifelse(is.na(min_years),      "off", as.character(min_years)),
            ifelse(is.na(min_n_sig_fpcs), "off", as.character(min_n_sig_fpcs)),
            ifelse(is.na(min_r2_full),    "off", sprintf("%.2f", min_r2_full))))
cat(sprintf("Loaded:      %d sites\n", n_loaded))
cat(sprintf("Included:    %d sites\n", n_kept))
cat(sprintf("Excluded:    %d sites\n", n_excluded))

if (n_excluded > 0) {
  cat("\nExcluded sites:\n")
  excluded_tbl <- inclusion_audit[!inclusion_audit$included,
                                  c("site", "n_years", "n_sig_fpcs",
                                    "r2_full", "reason")]
  excluded_tbl$r2_full <- round(excluded_tbl$r2_full, 3)
  print(excluded_tbl, row.names = FALSE)
}
cat("=====================================================\n")

# Save full audit table (every site, kept + excluded) so it can be inspected later
write.csv(inclusion_audit,
          file.path(out_dir, "inclusion_audit.csv"),
          row.names = FALSE)
cat(sprintf("Saved full audit: %s\n", file.path(out_dir, "inclusion_audit.csv")))

# Subset site_results to the included sites for the rest of the analysis
kept_sites   <- inclusion_audit$site[inclusion_audit$included]
site_results <- site_results[kept_sites]
n_sites      <- length(site_results)

if (n_sites < 3) {
  stop("Fewer than 3 sites passed the inclusion filter (",
       n_sites, " kept). Loosen the thresholds (min_years, min_n_sig_fpcs, ",
       "min_r2_full) at the top of the script.")
}

# =============================================================================
# 2. ASSEMBLE MATRICES
# =============================================================================

# β(t) matrix: rows = sites, cols = DOY 1:365
beta_mat <- do.call(rbind, lapply(site_results, `[[`, "beta_smooth"))
rownames(beta_mat) <- names(site_results)

# Weighted hydrograph matrix
wh_mat <- do.call(rbind, lapply(site_results, `[[`, "wh_vec"))
rownames(wh_mat) <- names(site_results)

cat(sprintf("Beta matrix: %d sites × %d DOY\n", nrow(beta_mat), ncol(beta_mat)))

# =============================================================================
# 3. PAIRWISE CORRELATION & DISTANCE
# =============================================================================

# --- 3a. Correlation matrix of β(t) curves ----------------------------------
cor_beta <- cor(t(beta_mat))

# --- 3b. L² functional distance: sqrt( ∫[β_i(t) - β_j(t)]² dt ) -----------
l2_dist <- as.dist(
  outer(1:n_sites, 1:n_sites, Vectorize(function(i, j) {
    sqrt(sum((beta_mat[i, ] - beta_mat[j, ])^2))
  }))
)
attr(l2_dist, "Labels") <- names(site_results)

# --- 3c. Also compute for weighted hydrographs ------------------------------
cor_wh <- cor(t(wh_mat))

# =============================================================================
# 4. HIERARCHICAL CLUSTERING
# =============================================================================

hc <- hclust(l2_dist, method = "ward.D2")

# Optimal number of clusters via silhouette width (k = 2 to 8)
max_k <- min(5, n_sites - 1)
sil_widths <- sapply(4:max_k, function(k) {
  cl <- cutree(hc, k = k)
  mean(silhouette(cl, l2_dist)[, "sil_width"])
})
k_opt <- which.max(sil_widths) + 3   # +3 because we started at k=4
cat(sprintf("Optimal k = %d (mean silhouette = %.3f)\n", k_opt, max(sil_widths)))

clusters <- cutree(hc, k = k_opt)

# Store cluster assignments
cluster_df <- data.frame(
  site    = names(clusters),
  cluster = factor(clusters)
)

# =============================================================================
# 5. EXTRACT SUMMARY METRICS FROM EACH β(t)
# =============================================================================

metric_df <- bind_rows(lapply(site_results, function(sr) {
  beta <- sr$beta_smooth
  sig_pos_doys <- which(sr$sim_sig_pos)
  sig_neg_doys <- which(sr$sim_sig_neg)
  
  # Identify contiguous significant windows
  pos_window_start <- if (length(sig_pos_doys) > 0) min(sig_pos_doys) else NA
  pos_window_end   <- if (length(sig_pos_doys) > 0) max(sig_pos_doys) else NA
  neg_window_start <- if (length(sig_neg_doys) > 0) min(sig_neg_doys) else NA
  neg_window_end   <- if (length(sig_neg_doys) > 0) max(sig_neg_doys) else NA
  
  data.frame(
    site              = sr$site,
    doy_peak_pos      = which.max(beta),
    doy_peak_neg      = which.min(beta),
    peak_pos_value    = max(beta),
    peak_neg_value    = min(beta),
    integrated_pos    = sum(pmax(beta, 0)),
    integrated_neg    = sum(pmin(beta, 0)),
    ratio_pos_neg     = sum(pmax(beta, 0)) / abs(sum(pmin(beta, 0))),
    n_sig_pos_days    = length(sig_pos_doys),
    n_sig_neg_days    = length(sig_neg_doys),
    centroid_pos      = if (length(sig_pos_doys) > 0) mean(sig_pos_doys) else NA,
    centroid_neg      = if (length(sig_neg_doys) > 0) mean(sig_neg_doys) else NA,
    pos_window_start  = pos_window_start,
    pos_window_end    = pos_window_end,
    pos_window_width  = ifelse(!is.na(pos_window_start), pos_window_end - pos_window_start + 1, 0),
    neg_window_start  = neg_window_start,
    neg_window_end    = neg_window_end,
    neg_window_width  = ifelse(!is.na(neg_window_start), neg_window_end - neg_window_start + 1, 0),
    n_sig_fpcs        = sr$n_sig,
    r2_full           = sr$r2_full,
    r2_reduced        = sr$r2_reduced,
    S_estimate        = sr$S_estimate,
    S_t               = sr$S_t,
    S_p               = sr$S_p,
    S_retained        = sr$S_retained,
    stringsAsFactors  = FALSE
  )
})) %>%
  left_join(cluster_df, by = "site")

cat("\n--- Summary Metrics ---\n")
print(metric_df %>% dplyr::select(site, cluster, doy_peak_pos, doy_peak_neg,
                                  n_sig_pos_days, n_sig_neg_days,
                                  r2_full, S_retained, S_estimate))

# =============================================================================
# 6. PCA ON SUMMARY METRICS
# =============================================================================

# Select numeric metrics for PCA (drop identifiers, cluster, and any NAs)
pca_vars <- c("doy_peak_pos", "doy_peak_neg", "peak_pos_value", "peak_neg_value",
              "integrated_pos", "integrated_neg", "ratio_pos_neg",
              "n_sig_pos_days", "n_sig_neg_days",
              "centroid_pos", "centroid_neg",
              "pos_window_width", "neg_window_width",
              "n_sig_fpcs", "r2_full", "r2_reduced", "S_estimate")

# Handle NAs: replace centroid NAs with 0 (no significant window)
pca_data <- metric_df %>%
  dplyr::select(all_of(pca_vars)) %>%
  mutate(across(everything(), ~replace_na(.x, 0)))

# Drop constant or zero-variance columns (can't rescale to unit variance)
col_vars   <- apply(pca_data, 2, var, na.rm = TRUE)
zero_var   <- names(col_vars[col_vars == 0 | is.na(col_vars)])
if (length(zero_var) > 0) {
  cat(sprintf("Dropping %d zero-variance column(s) from PCA: %s\n",
              length(zero_var), paste(zero_var, collapse = ", ")))
  pca_data <- pca_data %>% dplyr::select(-all_of(zero_var))
}

pca_out <- prcomp(pca_data, scale. = TRUE, center = TRUE)

pca_summary <- summary(pca_out)
cat("\nPCA variance explained:\n")
print(round(pca_summary$importance[, 1:min(5, ncol(pca_summary$importance))], 3))

n_pcs <- min(3, ncol(pca_out$x))
pca_scores <- as.data.frame(pca_out$x[, 1:n_pcs])
pca_scores$site    <- metric_df$site
pca_scores$cluster <- metric_df$cluster

# =============================================================================
# 7. NMDS ON L² DISTANCE
# =============================================================================

set.seed(42)
nmds_out <- metaMDS(l2_dist, k = 2, trymax = 100, trace = 0)
cat(sprintf("\nNMDS stress: %.4f\n", nmds_out$stress))

nmds_scores <- as.data.frame(scores(nmds_out, display = "sites"))
nmds_scores$site    <- names(site_results)
nmds_scores$cluster <- cluster_df$cluster

# =============================================================================
# 8. FIGURES
# =============================================================================

# Colour palette for clusters
# Colour palette for clusters (scales to k_opt)
clust_cols <- colorRampPalette(
  c("#4575b4", "#74add1", "#abd9e9", "#fee090",
    "#fdae61", "#f46d43", "#d73027", "#a50026")
)(k_opt)

# --- 8a. Correlation matrix heatmap ------------------------------------------
png(file.path(out_dir, "beta_correlation_matrix.png"),
    width = 10, height = 10, units = "in", res = 300)
corrplot(cor_beta, method = "color", type = "lower",
         order = "hclust", hclust.method = "ward.D2",
         tl.cex = 0.6, tl.col = "black",
         col = colorRampPalette(c("#d73027", "white", "#4575b4"))(200),
         title = "Pairwise Correlation of β(t) Curves Across Sites",
         mar = c(0, 0, 2, 0))
dev.off()
cat("Saved: beta_correlation_matrix.png\n")

# --- 8b. Dendrogram with cluster colouring -----------------------------------
dend <- as.dendrogram(hc) %>%
  set("branches_k_color", k = k_opt, value = clust_cols) %>%
  set("labels_cex", 0.6)

png(file.path(out_dir, "site_dendrogram.png"),
    width = 14, height = 7, units = "in", res = 300)
par(mar = c(8, 4, 3, 1))
plot(dend,
     main = sprintf("Hierarchical Clustering of Sites by β(t) Similarity (Ward's D², k = %d)", k_opt),
     ylab = "L² Distance", horiz = FALSE)
dev.off()
cat("Saved: site_dendrogram.png\n")

# --- 8c. All β(t) curves, coloured by cluster --------------------------------
beta_long <- data.frame(
  site    = rep(rownames(beta_mat), each = 365),
  DOY     = rep(1:365, n_sites),
  beta    = as.vector(t(beta_mat))
) %>%
  left_join(cluster_df, by = "site")

p_all_beta <- ggplot(beta_long, aes(x = DOY, y = beta, group = site, color = cluster)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_line(alpha = 0.6, linewidth = 0.5) +
  scale_color_manual(values = clust_cols, name = "Cluster") +
  scale_x_continuous(breaks = key_doys, labels = key_labels, expand = c(0.01, 0)) +
  labs(
    title    = "Smoothed β(t) Curves Across All Sites",
    subtitle = sprintf("Coloured by hierarchical cluster assignment (k = %d)", k_opt),
    x = "Day of Year", y = expression(beta(t))
  ) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title  = element_text(face = "bold"))

# --- 8d. Cluster-mean β(t) with individual curves in background --------------
cluster_means <- beta_long %>%
  group_by(cluster, DOY) %>%
  summarise(mean_beta = mean(beta), .groups = "drop")

p_cluster_means <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_line(data = beta_long,
            aes(x = DOY, y = beta, group = site, color = cluster),
            alpha = 0.15, linewidth = 0.4) +
  geom_line(data = cluster_means,
            aes(x = DOY, y = mean_beta, color = cluster),
            linewidth = 1.5) +
  scale_color_manual(values = clust_cols, name = "Cluster") +
  scale_x_continuous(breaks = key_doys, labels = key_labels, expand = c(0.01, 0)) +
  facet_wrap(~ cluster, ncol = 2) +
  labs(
    title    = "Cluster-Mean β(t) Curves with Individual Sites",
    subtitle = "Bold line = cluster mean  |  Faint lines = individual sites",
    x = "Day of Year", y = expression(beta(t))
  ) +
  theme_classic(base_size = 12) +
  theme(axis.text.x      = element_text(angle = 45, hjust = 1),
        plot.title        = element_text(face = "bold"),
        strip.text        = element_text(face = "bold"),
        strip.background  = element_rect(fill = "grey95", color = NA),
        legend.position   = "none")

fig_beta <- p_all_beta / p_cluster_means +
  plot_layout(heights = c(1, 1.2)) +
  plot_annotation(
    title = "Cross-Site Comparison of Flow–Recruitment Relationships",
    theme = theme(plot.title = element_text(size = 14, face = "bold"))
  )

ggsave(file.path(out_dir, "all_beta_curves_clustered.png"),
       fig_beta, width = 14, height = 12, dpi = 300)
cat("Saved: all_beta_curves_clustered.png\n")

# --- 8e. PCA biplot -----------------------------------------------------------
loadings_df <- as.data.frame(pca_out$rotation[, 1:2]) %>%
  rownames_to_column("variable") %>%
  mutate(
    # Scale loadings so arrows reach toward the edges of the score cloud
    PC1 = PC1 * max(abs(pca_scores$PC1)) * 0.85,
    PC2 = PC2 * max(abs(pca_scores$PC2)) * 0.85
  )

p_pca <- ggplot(pca_scores, aes(x = PC1, y = PC2)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey80") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey80") +
  # Loading arrows
  geom_segment(data = loadings_df,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
               color = "#7a4a00", linewidth = 0.6, alpha = 0.85) +
  # Site points
  geom_point(aes(fill = cluster), size = 4, shape = 21, stroke = 0.5) +
  # Site name labels (repelled so they don't pile up on a tight cluster)
  ggrepel::geom_text_repel(
    aes(label = site),
    size               = 2.6,
    color              = "grey20",
    bg.color           = "white",
    bg.r               = 0.15,
    max.overlaps       = Inf,
    box.padding        = 0.3,
    point.padding      = 0.2,
    segment.color      = "grey70",
    segment.size       = 0.25,
    seed               = 1
  ) +
  # Loading variable labels (repelled, boxed, drawn LAST so they sit on top)
  ggrepel::geom_label_repel(
    data               = loadings_df,
    aes(x = PC1, y = PC2, label = variable),
    size               = 3.2,
    fontface           = "bold",
    color              = "#7a4a00",
    fill               = scales::alpha("white", 0.9),
    label.size         = 0.3,
    label.r            = unit(0.12, "lines"),
    label.padding      = unit(0.18, "lines"),
    box.padding        = 0.6,
    point.padding      = 0.3,
    nudge_x            = 0,
    nudge_y            = 0,
    force              = 3,
    force_pull         = 0.5,
    max.overlaps       = Inf,
    segment.color      = "#7a4a00",
    segment.size       = 0.4,
    segment.alpha      = 0.6,
    segment.curvature  = -0.1,
    min.segment.length = 0,
    seed               = 2
  ) +
  scale_fill_manual(values = clust_cols, name = "Cluster") +
  labs(
    title    = "PCA of β(t) Summary Metrics",
    subtitle = sprintf("PC1: %.1f%% var.  |  PC2: %.1f%% var.",
                       pca_summary$importance[2, 1] * 100,
                       pca_summary$importance[2, 2] * 100),
    x = sprintf("PC1 (%.1f%%)", pca_summary$importance[2, 1] * 100),
    y = sprintf("PC2 (%.1f%%)", pca_summary$importance[2, 2] * 100)
  ) +
  # Add a little breathing room around the data so repelled labels don't clip
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 12) +
  theme(plot.title    = element_text(face = "bold"),
        plot.margin   = margin(10, 24, 10, 10))

ggsave(file.path(out_dir, "pca_biplot.png"),
       p_pca, width = 11, height = 9, dpi = 300)
cat("Saved: pca_biplot.png\n")

# --- 8f. NMDS ordination ------------------------------------------------------
p_nmds <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(fill = cluster), size = 4, shape = 21, stroke = 0.5) +
  ggrepel::geom_text_repel(
    aes(label = site),
    size               = 2.7,
    color              = "grey20",
    bg.color           = "white",
    bg.r               = 0.15,
    max.overlaps       = Inf,
    box.padding        = 0.3,
    point.padding      = 0.2,
    segment.color      = "grey70",
    segment.size       = 0.25,
    seed               = 1
  ) +
  scale_fill_manual(values = clust_cols, name = "Cluster") +
  labs(
    title    = "NMDS Ordination of Sites by β(t) L² Distance",
    subtitle = sprintf("Stress = %.4f  |  k = 2 dimensions", nmds_out$stress),
    x = "NMDS1", y = "NMDS2"
  ) +
  coord_equal(clip = "off") +
  theme_classic(base_size = 12) +
  theme(plot.title  = element_text(face = "bold"),
        plot.margin = margin(10, 24, 10, 10))

ggsave(file.path(out_dir, "nmds_ordination.png"),
       p_nmds, width = 10, height = 8, dpi = 300)
cat("Saved: nmds_ordination.png\n")

# --- 8g. Timing summary: peak DOY by cluster ---------------------------------
timing_df <- metric_df %>%
  dplyr::select(site, cluster, doy_peak_pos, doy_peak_neg,
                centroid_pos, centroid_neg) %>%
  pivot_longer(cols = c(doy_peak_pos, doy_peak_neg, centroid_pos, centroid_neg),
               names_to = "metric", values_to = "DOY") %>%
  filter(!is.na(DOY)) %>%
  mutate(
    type = ifelse(grepl("pos", metric), "Positive effect", "Negative effect"),
    metric_label = case_when(
      metric == "doy_peak_pos"  ~ "Peak positive β(t)",
      metric == "doy_peak_neg"  ~ "Peak negative β(t)",
      metric == "centroid_pos"  ~ "Centroid of sig. positive window",
      metric == "centroid_neg"  ~ "Centroid of sig. negative window"
    )
  )

p_timing <- ggplot(timing_df, aes(x = DOY, y = site, color = cluster)) +
  geom_point(size = 3, alpha = 0.8) +
  facet_wrap(~ metric_label, ncol = 2) +
  scale_color_manual(values = clust_cols, name = "Cluster") +
  scale_x_continuous(breaks = key_doys, labels = key_labels, expand = c(0.02, 0)) +
  labs(
    title    = "Timing of Flow Effects on Recruitment Across Sites",
    subtitle = "Peak and centroid DOY of significant β(t) windows",
    x = "Day of Year", y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(axis.text.x      = element_text(angle = 45, hjust = 1),
        axis.text.y      = element_text(size = 7),
        plot.title        = element_text(face = "bold"),
        strip.text        = element_text(face = "bold"),
        strip.background  = element_rect(fill = "grey95", color = NA))

ggsave(file.path(out_dir, "timing_summary.png"),
       p_timing, width = 14, height = max(6, n_sites * 0.3), dpi = 300)
cat("Saved: timing_summary.png\n")

# --- 8h. FPC selection heatmap: which FPCs matter at which sites? -------------
# Build a matrix of which FPC numbers were selected at each site
max_fpc_num <- max(sapply(site_results, function(sr) max(sr$sig_fpc_nums)))
fpc_selection_mat <- matrix(0, nrow = n_sites, ncol = max_fpc_num,
                            dimnames = list(names(site_results),
                                            paste0("FPC", 1:max_fpc_num)))
for (sr in site_results) {
  fpc_selection_mat[sr$site, sr$sig_fpc_nums] <- 1
}

# Count how many sites selected each FPC
fpc_freq <- colSums(fpc_selection_mat)
cat("\nFPC selection frequency across sites:\n")
print(fpc_freq[fpc_freq > 0])

fpc_heat_df <- as.data.frame(fpc_selection_mat) %>%
  rownames_to_column("site") %>%
  left_join(cluster_df, by = "site") %>%
  arrange(cluster, site) %>%
  pivot_longer(cols = starts_with("FPC"), names_to = "FPC", values_to = "selected") %>%
  mutate(FPC = factor(FPC, levels = paste0("FPC", 1:max_fpc_num)))

# Only show FPCs that were selected at least once
used_fpcs <- names(fpc_freq[fpc_freq > 0])
fpc_heat_df <- fpc_heat_df %>% filter(FPC %in% used_fpcs)

p_fpc_heat <- ggplot(fpc_heat_df, aes(x = FPC, y = fct_rev(site), fill = factor(selected))) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_manual(values = c("0" = "grey95", "1" = "#4575b4"),
                    labels = c("Not selected", "Selected"),
                    name = NULL) +
  labs(
    title    = "FPC Selection Across Sites",
    subtitle = "Which flow modes (FPCs) were statistically significant at each site?",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(axis.text.y  = element_text(size = 7),
        plot.title   = element_text(face = "bold"),
        panel.grid   = element_blank())

ggsave(file.path(out_dir, "fpc_selection_heatmap.png"),
       p_fpc_heat, width = 10, height = max(5, n_sites * 0.35), dpi = 300)
cat("Saved: fpc_selection_heatmap.png\n")

# --- 8i. Spawner (S) effect strength across sites ---------------------------
# WorkingFPCA.R now fits the spawner covariate S directly with its slope
# constrained to be <= 0. This panel shows, for each site, how strong that
# density-dependent depression of recruits-per-spawner ended up being.
S_df <- metric_df %>%
  dplyr::select(site, cluster, S_estimate, S_p, S_retained) %>%
  mutate(
    S_estimate = ifelse(is.na(S_estimate), 0, S_estimate),
    sig_label  = case_when(
      is.na(S_p)         ~ "n/a",
      S_p < 0.001        ~ "***",
      S_p < 0.01         ~ "**",
      S_p < 0.05         ~ "*",
      S_p < 0.10         ~ ".",
      TRUE               ~ ""
    ),
    fill_grp   = ifelse(S_retained, "Retained (p < 0.05, b < 0)", "Not retained")
  ) %>%
  arrange(S_estimate) %>%
  mutate(site = factor(site, levels = site))

p_spawner <- ggplot(S_df, aes(x = S_estimate, y = site, fill = fill_grp)) +
  geom_vline(xintercept = 0, color = "grey40", linewidth = 0.4) +
  geom_col(width = 0.75, alpha = 0.9) +
  geom_text(aes(label = sig_label),
            hjust = -0.2, size = 3, color = "grey25") +
  scale_fill_manual(
    values = c("Retained (p < 0.05, b < 0)" = "#1f3864",
               "Not retained"               = "grey75"),
    name = NULL
  ) +
  labs(
    title    = "Density-Dependent Spawner Effect (S, constrained ≤ 0)",
    subtitle = "Constrained slope on z-scored S in the full model. Stars indicate significance level.",
    x = "Standardized coefficient on S (negative = stronger density dependence)",
    y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(plot.title       = element_text(face = "bold"),
        axis.text.y      = element_text(size = 7),
        legend.position  = "bottom")

ggsave(file.path(out_dir, "spawner_effect_strength.png"),
       p_spawner, width = 10, height = max(5, n_sites * 0.32), dpi = 300)
cat("Saved: spawner_effect_strength.png\n")

# --- 8j. Cross-site year-by-FPC heatmap -------------------------------------
# Builds a sites × years matrix of FPC1 scores (z-scored within each site so
# scores are comparable across systems with different score variances). FPC1
# is the dominant flow-variability mode at most sites, so this view answers
# "which years had unusual flow regionally?" Years where a site's |z-score|
# exceeds 2 are marked with a bold black border.
#
# Note: each site's FPC1 captures different aspects of its hydrograph, so this
# heatmap reads more as "in which years was each site's dominant flow mode
# unusual?" rather than "which years had the same regional anomaly." Still
# very useful for spotting basin-wide drought / flood years.

# Choose which FPC to map. FPC1 is the most-comparable across sites, but you
# can change `which_fpc` to look at FPC2, FPC3, etc.
which_fpc <- 1L

# Pull FPC scores by site, z-score within site, identify extreme cells
xs_yr_df <- bind_rows(lapply(site_results, function(sr) {
  if (is.null(sr$fpc_scores) || length(sr$years_vec) < 2) return(NULL)
  if (which_fpc > ncol(sr$fpc_scores))                    return(NULL)
  scores <- sr$fpc_scores[, which_fpc]
  data.frame(
    site  = sr$site,
    Year  = sr$years_vec,
    score = as.numeric(scores)
  )
})) %>%
  group_by(site) %>%
  mutate(score_z = as.numeric(scale(score))) %>%
  ungroup()

# Order sites by hierarchical cluster so similar sites sit next to each other
xs_yr_df <- xs_yr_df %>%
  left_join(cluster_df, by = "site") %>%
  arrange(cluster, site) %>%
  mutate(site = factor(site, levels = unique(site)))

# Identify cells with |z| > 2 for the highlight overlay
extreme_cells <- xs_yr_df %>% filter(abs(score_z) > 2)

# Identify the years that hit > 2 SD across many sites (regional anomalies)
year_anomaly <- xs_yr_df %>%
  group_by(Year) %>%
  summarise(n_extreme = sum(abs(score_z) > 2, na.rm = TRUE),
            .groups = "drop") %>%
  arrange(desc(n_extreme))

cat(sprintf("\nCross-site FPC%d anomaly counts (years with most extreme z-scores across sites):\n",
            which_fpc))
print(head(year_anomaly, 10), row.names = FALSE)

score_lim_xs <- max(abs(xs_yr_df$score), na.rm = TRUE) * 1.05

p_xs_year <- ggplot(xs_yr_df, aes(x = factor(Year), y = fct_rev(site))) +
  geom_tile(aes(fill = score), color = "white", linewidth = 0.3) +
  geom_tile(data = extreme_cells,
            aes(x = factor(Year), y = fct_rev(site)),
            fill = NA, color = "black", linewidth = 0.7,
            inherit.aes = FALSE) +
  scale_fill_gradient2(
    low = "#b2182b", mid = "#f7f7f7", high = "#2166ac",
    midpoint = 0, limits = c(-score_lim_xs, score_lim_xs),
    name = sprintf("FPC%d\nscore", which_fpc)
  ) +
  labs(
    title    = sprintf("Cross-site Year × FPC%d Score Heatmap", which_fpc),
    subtitle = paste0(
      "Each cell shows that site's FPC", which_fpc, " score in that year.  ",
      "Black borders = |z-score| > 2 (extreme year for that site).  ",
      "Sites grouped by cluster on the y-axis."
    ),
    x = "Year", y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(plot.title       = element_text(face = "bold"),
        plot.subtitle    = element_text(size = 9, color = "grey40"),
        axis.text.x      = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y      = element_text(size = 7),
        axis.ticks       = element_blank(),
        panel.grid       = element_blank(),
        legend.position  = "right")

ggsave(file.path(out_dir, "year_fpc_heatmap_crosssite.png"),
       p_xs_year,
       width  = max(11, 2 + length(unique(xs_yr_df$Year)) * 0.45),
       height = max(5, n_sites * 0.32),
       dpi    = 300, limitsize = FALSE)
cat("Saved: year_fpc_heatmap_crosssite.png\n")

# Save a CSV of the anomaly counts so users can sort and look up specific years
write.csv(year_anomaly,
          file.path(out_dir, "year_anomaly_counts.csv"),
          row.names = FALSE)
cat("Saved: year_anomaly_counts.csv\n")

# =============================================================================
# 9. OPTIONAL: TEST COVARIATES AGAINST CLUSTERS / ORDINATION
# =============================================================================

if (!is.null(covariate_file) && file.exists(covariate_file)) {
  
  cat("\n--- Loading site covariates ---\n")
  covariates <- read.csv(covariate_file) %>%
    mutate(site = gsub("\\.", "", Stream.Section))   # match naming convention
  
  metric_cov <- metric_df %>%
    left_join(covariates, by = "site")
  
  # PERMANOVA: do clusters differ by covariates?
  cov_cols <- intersect(names(covariates),
                        c("drainage_area_km2", "elevation_m", "latitude",
                          "mean_annual_flow", "stream_order", "gradient"))
  
  if (length(cov_cols) > 0) {
    cat("\nPERMANOVA: testing covariate effects on β(t) distance matrix\n")
    
    # Build covariate matrix aligned with distance matrix
    cov_aligned <- metric_cov %>%
      dplyr::select(site, all_of(cov_cols)) %>%
      filter(site %in% attr(l2_dist, "Labels")) %>%
      arrange(match(site, attr(l2_dist, "Labels"))) %>%
      dplyr::select(-site)
    
    for (cv in cov_cols) {
      if (all(!is.na(cov_aligned[[cv]]))) {
        form <- as.formula(paste("l2_dist ~", cv))
        perm_result <- adonis2(form, data = cov_aligned, permutations = 999)
        cat(sprintf("  %s: R² = %.3f, p = %.3f\n",
                    cv, perm_result$R2[1], perm_result$`Pr(>F)`[1]))
      }
    }
    
    # Overlay covariates on NMDS (envfit)
    env_fit <- envfit(nmds_out, cov_aligned, permutations = 999, na.rm = TRUE)
    cat("\nenvfit results:\n")
    print(env_fit)
    
    # Plot NMDS with envfit vectors
    envfit_df <- as.data.frame(scores(env_fit, display = "vectors")) %>%
      rownames_to_column("variable") %>%
      mutate(
        NMDS1 = NMDS1 * max(abs(nmds_scores$NMDS1)) * 0.8,
        NMDS2 = NMDS2 * max(abs(nmds_scores$NMDS2)) * 0.8
      )
    
    p_nmds_env <- p_nmds +
      geom_segment(data = envfit_df,
                   aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
                   inherit.aes = FALSE,
                   arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
                   color = "darkgreen", linewidth = 0.7) +
      ggrepel::geom_label_repel(
        data               = envfit_df,
        aes(x = NMDS1, y = NMDS2, label = variable),
        inherit.aes        = FALSE,
        size               = 3.4,
        fontface           = "bold",
        color              = "darkgreen",
        fill               = scales::alpha("white", 0.9),
        label.size         = 0.3,
        label.r            = unit(0.12, "lines"),
        label.padding      = unit(0.18, "lines"),
        box.padding        = 0.6,
        force              = 3,
        max.overlaps       = Inf,
        segment.color      = "darkgreen",
        segment.size       = 0.4,
        segment.alpha      = 0.6,
        min.segment.length = 0,
        seed               = 3
      ) +
      labs(title = "NMDS with Environmental Vectors (envfit)")
    
    ggsave(file.path(out_dir, "nmds_with_covariates.png"),
           p_nmds_env, width = 10, height = 8, dpi = 300)
    cat("Saved: nmds_with_covariates.png\n")
  }
  
} else {
  cat("\nNo covariate file provided. Skipping covariate analysis.\n")
  cat("To add covariates, create a CSV with columns:\n")
  cat("  Stream.Section, drainage_area_km2, elevation_m, latitude, ...\n")
  cat("and set covariate_file at the top of this script.\n")
}

# =============================================================================
# 10. EXPORT SUMMARY TABLE
# =============================================================================

write.csv(metric_df,
          file.path(out_dir, "cross_site_summary_metrics.csv"),
          row.names = FALSE)
cat("\nSaved: cross_site_summary_metrics.csv\n")

# =============================================================================
# FINAL SUMMARY
# =============================================================================
cat("\n============================================================\n")
cat("CROSS-SITE ANALYSIS COMPLETE\n")
cat("============================================================\n")
cat(sprintf("Sites loaded:          %d\n", n_loaded))
cat(sprintf("Sites included:        %d\n", n_sites))
cat(sprintf("Sites excluded:        %d  (see inclusion_audit.csv for reasons)\n",
            n_excluded))
cat(sprintf("Optimal clusters:      %d (silhouette = %.3f)\n", k_opt, max(sil_widths)))
cat(sprintf("NMDS stress:           %.4f\n", nmds_out$stress))
cat(sprintf("PCA (PC1 + PC2):       %.1f%% of metric variance\n",
            sum(pca_summary$importance[2, 1:2]) * 100))
cat(sprintf("\nCluster membership:\n"))
print(table(clusters))
cat(sprintf("\nAll outputs saved to: %s/\n", out_dir))
cat("============================================================\n")