library(tidyverse)
library(patchwork)
LL_Sections<-read.csv("data/LL_all_df.csv")%>%pull(Stream.Section)%>%unique()

# skipSites<-c('Beaverhead.FishAndGame"')
lapply(LL_Sections,function(SS){
  # SS<-"Beaverhead.Hildreth"
  print(SS)
  
# =============================================================================
# Flow–Recruitment Functional Regression
# FPCA-based scalar-on-function regression with auto-selected reduced model
#
# Key change: the reduced model no longer hardcodes the top 2 FPCs.
# Instead, ALL statistically significant FPCs (p < 0.05 OR |t| > 1.6) from
# the full model are retained, and sign-flipping is handled automatically so
# that all retained coefficients are positive (= flow pattern benefits
# recruitment). Panel C labels show both flow variance and recruitment
# variance per FPC. Thresholds are set at the top of the selection block.
# =============================================================================
  # read.csv('weightedHydrosDatIN.csv')%>%group_by(Stream.Section)%>%summarize(n_year=length(unique(Year)))%>%filter(n_year<20)#ggplot(aes(x=Stream.Section,y=n_year))+
  #   geom_bar(stat='identity',position=position_dodge())
# --- 1. Load and inspect data -----------------------------------------------
  
  Stream<-strsplit(SS,split = '\\.')[[1]][1]
  Section<-strsplit(SS,split = '\\.')[[1]][2]
  StreamSection<-paste0(strsplit(SS,split = '\\.')[[1]][1],strsplit(SS,split = '\\.')[[1]][2])
  if(SS=="BigHole.Melrose"){
    FishDat<-read.csv(here::here('../../LL/JAGS_PVA/ModelOutput/csvs_quadratic/',paste0(Stream,Section,"allParams.csv")))%>%
      filter(SummerLag==2 & WinterLag==2)
  }else{
    FishDat <- read.csv(paste0('../../LL/JAGS_PVA/ModelOutput/csvs_quadratic/',StreamSection,"allParams.csv"),header=T)%>%filter(topMod=='yes')%>%filter(grepl("SUMMERQ",model)|grepl("Global",model))
    
  }

    head(FishDat)
    # fishYears<-FishDat%>%filter(!is.na(Year))%>%pull(Year)%>%unique()
    # flowYears<-c((min(fishYears)-3):(max(fishYears)))
    # FishDat %>% pull(name) %>% unique()
    
    N2 <- FishDat %>% filter(name == 'Estimated N2') %>% dplyr::select(Year,N50) %>% dplyr::rename(N50_R = N50)
    N34 <- FishDat %>% filter(name %in% c('Estimated N4','Estimated N3')) %>% mutate(Year = Year + 3) %>% group_by(Year) %>% 
      dplyr::summarize(N50 = sum(N50)) %>% mutate(name = 'NAD') %>% dplyr::select(Year,N50) %>% dplyr::rename(N50_AD = N50)
    
    RpS <- N2 %>% mutate(Year = as.numeric(Year)) %>% left_join(N34,by='Year') %>% mutate(R=N50_R,S=N50_AD,RpS = N50_R/N50_AD,lRpS = log(RpS)) %>% na.omit()
    
    
    recruit_df <- RpS %>% mutate(Year = Year - 2,recruit_lag=2) %>% mutate(RelWeight = lRpS%>%scale(),sdR=sd(lRpS,na.rm=T))
    
    # SS<-'BigHole.Melrose'
    if(SS=='Ruby.SilverSprings'){
      SS2<-'Ruby.Vigilante'
    }else{
      SS2<-SS
    }
df<-read.csv(paste0('imputed_output/',SS2,'_imputed.csv'))%>%left_join(recruit_df)
# df<-read.csv('weightedHydrosDatIN.csv')%>%filter(Stream.Section==SS & recruit_lag==2)
head(df)
dat_long <- df %>%
  filter(DOY!=366)%>% 
  # dplyr::select(-ch, -meanQ) %>%
  left_join(
    df %>%
      filter(DOY != 366) %>%
      group_by(DOY) %>%
      dplyr::summarize(
        ch    = mean(log(Flow_imputed), na.rm = TRUE),
        meanQ = mean(Flow_imputed,      na.rm = TRUE)
      )
  ) %>%
  mutate(flow_std = (log(Flow_imputed) - ch) %>% scale())%>%filter(!is.na(lRpS))
# write.csv(dat_long,paste0('data/LL_',SS,'_weightedhydro_dat.csv'))
cat("Dimensions:", nrow(dat_long), "rows,", ncol(dat_long), "cols\n")
cat("Years:", paste(sort(unique(dat_long$Year)), collapse = ", "), "\n")
cat("DOY range:", min(dat_long$DOY), "to", max(dat_long$DOY), "\n")
cat("flow_std range:", round(range(dat_long$flow_std, na.rm = TRUE), 3), "\n")
cat("lRpS range:", round(range(dat_long$lRpS, na.rm = TRUE), 3), "\n")

# --- 2. Handle leap years ---------------------------------------------------
# Standardize all years to DOY 1-365 by dropping DOY 366

dat_long <- dat_long %>%
  filter(DOY <= 365)

cat("\nAfter removing DOY 366:\n")
print(table(dat_long$Year))
YearsIN<-names(table(dat_long$Year)[which(table(dat_long$Year)==365)])
dat_long <- dat_long %>%
  filter(Year%in%YearsIN)
# --- 2b. Ricker residualization ---------------------------------------------
# Remove density dependence from log(R/S) before flow regression
# Model: log(R/S) = alpha - beta*S  (linearized Ricker)

sr_yearly <- dat_long %>%
  group_by(Year) %>%
  dplyr::summarise(
    lRpS = first(lRpS),
    S    = first(S),
    .groups = "drop"
  ) %>%
  filter(!is.na(S), S > 0)

mod_ricker  <- lm(lRpS ~ S, data = sr_yearly)
ricker_r2   <- summary(mod_ricker)$r.squared
ricker_beta <- coef(mod_ricker)[2]

cat(sprintf("Ricker regression: beta = %.5f, R² = %.3f\n", ricker_beta, ricker_r2))
cat(sprintf("Variance in log(R/S) explained by spawner density: %.1f%%\n", ricker_r2 * 100))

sr_yearly$lRpS_resid <- resid(mod_ricker)

dat_long <- dat_long %>%
  left_join(dplyr::select(sr_yearly, Year, lRpS_resid), by = "Year") %>%
  mutate(lRpS_raw = lRpS,
         lRpS     = lRpS)
# 
# cat(sprintf("lRpS now = Ricker residuals (SD: %.3f vs raw SD: %.3f)\n",
#             sd(sr_yearly$lRpS_resid), sd(sr_yearly$lRpS)))

# --- 3. Build functional predictor matrix X ---------------------------------
# Rows = years, Columns = DOY 1:365

X_long <- dat_long %>%
  dplyr::select(Year, DOY, flow_std) %>%
  arrange(Year, DOY)

completeness <- X_long %>%
  group_by(Year) %>%
  dplyr::summarise(n_doy = n(), n_missing = sum(is.na(flow_std)))
print(completeness)

X_wide <- X_long %>%
  pivot_wider(names_from = DOY, values_from = flow_std) %>%
  arrange(Year)

X_mat <- as.matrix(X_wide[, -1])
rownames(X_mat) <- X_wide$Year
colnames(X_mat) <- 1:365

cat("\nFunctional matrix X dimensions:", dim(X_mat), "\n")

# --- 4. Build scalar response vector y --------------------------------------

y_df <- dat_long %>%
  group_by(Year) %>%
  dplyr::summarise(lRpS = first(lRpS),
                   S = first(S)) %>%
  arrange(Year)

y_vec    <- y_df$lRpS
years_vec <- y_df$Year
s_vec <- y_df$S
cat("Response vector length:", length(y_vec), "\n")
cat("lRpS values:\n")
print(data.frame(Year = years_vec, lRpS = round(y_vec, 3)))

stopifnot(all(as.integer(rownames(X_mat)) == years_vec))
cat("X matrix rows and y vector are aligned: OK\n")

# --- 5. Compute weighted hydrograph (informal approach) ---------------------

lRpS_z             <- scale(y_vec)[, 1]
weighted_hydrograph <- colSums(X_mat * lRpS_z) / nrow(X_mat)

wh_df <- data.frame(
  DOY              = 1:365,
  weighted_anomaly = weighted_hydrograph
)

# --- 6. Fit scalar-on-function regression (FPCA approach) -------------------

library(fda)

doy_grid     <- 1:365
nbasis       <- 40
bspline_basis <- create.bspline.basis(rangeval = c(1, 365), nbasis = nbasis)

flow_fd <- smooth.basis(
  argvals  = doy_grid,
  y        = t(X_mat),
  fdParobj = bspline_basis
)$fd

cat("Flow curves smoothed onto B-spline basis\n")

fpca_out     <- pca.fd(flow_fd, nharm = 20)
var_explained <- fpca_out$varprop * 100
cum_var       <- cumsum(var_explained)

cat("\nVariance explained by each FPC:\n")
print(round(data.frame(FPC = 1:20, Pct = var_explained, Cumulative = cum_var), 1))

n_fpc <- max(5, which(cum_var >= 90)[1])
n_fpc <- min(n_fpc, 10)
cat("\nRetaining", n_fpc, "FPCs (", round(cum_var[n_fpc], 1), "% variance)\n")

fpc_scores <- fpca_out$scores[, 1:n_fpc]
colnames(fpc_scores) <- paste0("FPC", 1:n_fpc)

score_df  <- as.data.frame(fpc_scores)
score_df$y <- y_vec
score_df$S <- s_vec
mod_lm <- lm(y ~ ., data = score_df) ###add spawners in here y ~ ssb + fcpa_scores constrain slope on ssb to be negative 

cat("\n--- Full FPC Regression Summary ---\n")
print(summary(mod_lm))

fpc_funs <- eval.fd(doy_grid, fpca_out$harmonics[1:n_fpc])
b_coefs  <- coef(mod_lm)[-1]
b_coefs<-b_coefs[which(names(b_coefs)!='S')]
beta_t   <- as.vector(fpc_funs %*% b_coefs)

# Bootstrap confidence bands (full model)
set.seed(42)
n_boot    <- 1000
beta_boot <- matrix(NA, nrow = n_boot, ncol = 365)

cat("\nRunning", n_boot, "bootstrap iterations...\n")
for (i in 1:n_boot) {
  idx        <- sample(nrow(score_df), replace = TRUE)
  boot_df    <- score_df[idx, ]
  boot_lm    <- lm(y ~ ., data = boot_df)
  boot_coefs <- coef(boot_lm)[-1]
  boot_coefs<-boot_coefs[which(names(boot_coefs)!='S')]
  
  beta_boot[i, ] <- as.vector(fpc_funs %*% boot_coefs)
}
cat("Bootstrap complete\n")

beta_df <- data.frame(
  DOY      = doy_grid,
  beta     = beta_t,
  lower_95 = apply(beta_boot, 2, quantile, 0.025,na.rm=T),
  upper_95 = apply(beta_boot, 2, quantile, 0.975,na.rm=T),
  lower_80 = apply(beta_boot, 2, quantile, 0.100,na.rm=T),
  upper_80 = apply(beta_boot, 2, quantile, 0.900,na.rm=T)
) %>%
  mutate(
    sig_pos = lower_95 > 0,
    sig_neg = upper_95 < 0
  )

cat("\nDOY with significant positive beta(t) (95% bootstrap CI > 0):\n")
print(beta_df$DOY[beta_df$sig_pos])
cat("\nDOY with significant negative beta(t) (95% bootstrap CI < 0):\n")
print(beta_df$DOY[beta_df$sig_neg])

# --- 7. Permutation test (full model) ---------------------------------------

n_perm        <- 999
beta_perm_mat <- matrix(NA, nrow = n_perm, ncol = 365)

cat("\nRunning", n_perm, "permutations...\n")
set.seed(123)
for (i in 1:n_perm) {
  if (i %% 200 == 0) cat("  Permutation", i, "of", n_perm, "\n")
  perm_df       <- score_df
  perm_df$y     <- sample(score_df$y)
  perm_lm       <- lm(y ~ ., data = perm_df)
  perm_coefs    <- coef(perm_lm)[-1]
  perm_coefs<-perm_coefs[which(names(perm_coefs)!='S')]
  
  beta_perm_mat[i, ] <- as.vector(fpc_funs %*% perm_coefs)
}
cat("Permutation test complete\n")

max_dev       <- apply(abs(beta_perm_mat), 1, max)
sim_threshold <- quantile(max_dev, 0.95)

beta_df <- beta_df %>%
  mutate(
    perm_lower  = apply(beta_perm_mat, 2, quantile, 0.025),
    perm_upper  = apply(beta_perm_mat, 2, quantile, 0.975),
    sim_sig_pos = beta >  sim_threshold,
    sim_sig_neg = beta < -sim_threshold
  )

cat(sprintf("\nSimultaneous significance threshold: %.5f\n", sim_threshold))
cat("\nDOY with simultaneous significant positive effect:\n")
print(beta_df$DOY[beta_df$sim_sig_pos])
cat("\nDOY with simultaneous significant negative effect:\n")
print(beta_df$DOY[beta_df$sim_sig_neg])

# --- 8. Calendar date helpers ------------------------------------------------

doy_to_date <- function(doy) {
  format(as.Date(doy - 1, origin = "2001-01-01"), "%b %d")
}

key_doys   <- c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 365)
key_labels <- doy_to_date(key_doys)

# --- 9. Full-model diagnostic plots -----------------------------------------

p_wh <- ggplot(wh_df, aes(x = DOY, y = weighted_anomaly)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_ribbon(aes(ymin = pmin(weighted_anomaly, 0), ymax = 0),
              fill = "#d73027", alpha = 0.3) +
  geom_ribbon(aes(ymin = 0, ymax = pmax(weighted_anomaly, 0)),
              fill = "#4575b4", alpha = 0.3) +
  geom_line(color = "black", linewidth = 0.7) +
  scale_x_continuous(breaks = key_doys, labels = key_labels, expand = c(0.01, 0)) +
  labs(title = "A: Weighted Hydrograph (informal)", x = NULL,
       y = "Weighted Flow Anomaly") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_beta <- ggplot(beta_df, aes(x = DOY)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_ribbon(aes(ymin = -sim_threshold, ymax = sim_threshold),
              fill = "grey85", alpha = 0.6) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95),
              fill = "#74add1", alpha = 0.4) +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80),
              fill = "#4575b4", alpha = 0.4) +
  geom_line(aes(y = beta), color = "black", linewidth = 0.9) +
  geom_point(data = filter(beta_df, sim_sig_pos), aes(y = beta),
             color = "#4575b4", size = 0.8) +
  geom_point(data = filter(beta_df, sim_sig_neg), aes(y = beta),
             color = "#d73027", size = 0.8) +
  scale_x_continuous(breaks = key_doys, labels = key_labels, expand = c(0.01, 0)) +
  labs(
    title    = expression("B: Functional Regression "*beta(t)*" with Confidence Bands"),
    subtitle = "Grey band = simultaneous permutation null envelope (95%)\nBlue band = 95% CI from model; darker = 80% CI",
    x = NULL, y = expression(beta(t))
  ) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

wh_scaled   <- scale(wh_df$weighted_anomaly)[, 1]
beta_scaled <- scale(beta_df$beta)[, 1]

overlay_df <- data.frame(
  DOY         = 1:365,
  wh_scaled   = wh_scaled,
  beta_scaled = beta_scaled
)

p_overlay <- ggplot(overlay_df, aes(x = DOY)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_line(aes(y = wh_scaled,   color = "Weighted hydrograph"), linewidth = 0.8) +
  geom_line(aes(y = beta_scaled, color = "Functional regression β(t)"),
            linewidth = 0.8, linetype = "longdash") +
  scale_color_manual(
    values = c("Weighted hydrograph"        = "#d73027",
               "Functional regression β(t)" = "#4575b4")
  ) +
  scale_x_continuous(breaks = key_doys, labels = key_labels, expand = c(0.01, 0)) +
  labs(title = "C: Comparison (both z-scored)", x = "Day of Year",
       y = "Standardized Effect", color = NULL) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

fig_combined <- p_wh / p_beta / p_overlay +
  patchwork::plot_annotation(
    title    = "Flow–Recruitment Relationship: Single Site Analysis",
    subtitle = paste0("n = ", length(y_vec),
                      " years | Functional predictor: daily flow residuals (DOY 1–365)"),
    theme = theme(
      plot.title    = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11, color = "grey40")
    )
  )

ggsave(
  paste0("output/plots/LL_", SS, "_functional_regression_results.png"),
  fig_combined, width = 10, height = 12, dpi = 300
)
cat("\nFigure saved: functional_regression_results.png\n")

# --- 10. Summary statistics (full model) ------------------------------------

pos_sig_doys   <- beta_df$DOY[beta_df$sim_sig_pos]
neg_sig_doys   <- beta_df$DOY[beta_df$sim_sig_neg]
integrated_pos <- sum(pmax(beta_df$beta, 0))
integrated_neg <- sum(pmin(beta_df$beta, 0))
mod_summary    <- summary(mod_lm)
r2_adj         <- mod_summary$adj.r.squared
r2             <- mod_summary$r.squared

cat("\n============================================================\n")
cat("RESULTS SUMMARY (Full model)\n")
cat("============================================================\n")
cat(sprintf("FPCs retained:               %d (%.1f%% variance)\n", n_fpc, cum_var[n_fpc]))
cat(sprintf("R-squared (FPC regression):  %.3f\n", r2))
cat(sprintf("Adjusted R-squared:          %.3f\n", r2_adj))
cat(sprintf("Integrated positive beta:    %.4f\n", integrated_pos))
cat(sprintf("Integrated negative beta:    %.4f\n", integrated_neg))
cat(sprintf("DOY of peak positive beta:   %d (%s)\n",
            beta_df$DOY[which.max(beta_df$beta)],
            doy_to_date(beta_df$DOY[which.max(beta_df$beta)])))
cat(sprintf("DOY of peak negative beta:   %d (%s)\n",
            beta_df$DOY[which.min(beta_df$beta)],
            doy_to_date(beta_df$DOY[which.min(beta_df$beta)])))
if (length(pos_sig_doys) > 0) {
  cat(sprintf("Sig. positive window:        DOY %d-%d (%s to %s)\n",
              min(pos_sig_doys), max(pos_sig_doys),
              doy_to_date(min(pos_sig_doys)), doy_to_date(max(pos_sig_doys))))
} else {
  cat("Sig. positive window:        none\n")
}
if (length(neg_sig_doys) > 0) {
  cat(sprintf("Sig. negative window:        DOY %d-%d (%s to %s)\n",
              min(neg_sig_doys), max(neg_sig_doys),
              doy_to_date(min(neg_sig_doys)), doy_to_date(max(neg_sig_doys))))
} else {
  cat("Sig. negative window:        none\n")
}
cat("============================================================\n")

save(
  mod_lm, fpca_out, beta_df, wh_df,
  beta_boot, beta_perm_mat, sim_threshold,
  n_fpc, fpc_scores, fpc_funs,
  file = paste0("output/models/LL_", SS, "_functional_regression_results.RData")
)
cat("\nR objects saved: functional_regression_results.RData (full model only; reduced model appended later)\n")


# =============================================================================
# AUTO-SELECTED REDUCED MODEL
#
# Rather than hardcoding FPC1 + FPC4, the two FPCs most strongly associated
# with log(R/S) are identified from the full model by |t-statistic|.
# Any selected FPC with a negative coefficient is sign-flipped so that:
#   Positive coefficient = this flow pattern benefits recruitment.
# This makes beta(t) and the composite score consistently interpretable
# across sites regardless of which FPCs happen to dominate.
# =============================================================================
######EDIT: USE ALL STATISTICALLY SIGNIFICANT FPCs ######
# --- Step 1: Identify significant FPCs (p < 0.05 OR |t| > 1.6) from full model

# Significance thresholds (adjust here as needed)
p_threshold <- 0.05
t_threshold <- 1.6

full_model_summary <- summary(mod_lm)
t_stats            <- full_model_summary$coefficients[-1, "t value"]   # drop intercept
t_stats <- t_stats[which(names(t_stats)!='S')]
p_vals             <- full_model_summary$coefficients[-1, "Pr(>|t|)"]
p_vals <- p_vals[which(names(p_vals)!='S')]

# Select FPCs meeting EITHER criterion, ordered by |t-statistic|
sig_positions <- which(p_vals < p_threshold | abs(t_stats) > t_threshold)
sig_positions <- sig_positions[order(abs(t_stats[sig_positions]), decreasing = TRUE)]
sig_fpc_nums  <- seq_len(n_fpc)[sig_positions]                # actual FPC numbers

# Fallback: if none pass either criterion, take the single strongest
if (length(sig_fpc_nums) == 0) {
  cat(sprintf("  WARNING: No FPCs met p < %.2f or |t| > %.1f; retaining the single strongest.\n",
              p_threshold, t_threshold))
  sig_positions <- order(abs(t_stats), decreasing = TRUE)[1]
  sig_fpc_nums  <- seq_len(n_fpc)[sig_positions]
}

n_sig <- length(sig_fpc_nums)

cat("\n=====================================================\n")
cat("AUTO-SELECTED REDUCED MODEL (all significant FPCs)\n")
cat(sprintf("Criteria: p < %.2f  OR  |t| > %.1f\n", p_threshold, t_threshold))
cat("=====================================================\n")
cat(sprintf("%d FPC(s) retained:\n", n_sig))
for (k in sig_fpc_nums) {
  row_idx <- k + 1   # +1 because row 1 = intercept
  this_p <- full_model_summary$coefficients[row_idx, "Pr(>|t|)"]
  this_t <- full_model_summary$coefficients[row_idx, "t value"]
  flag   <- ifelse(this_p < p_threshold & abs(this_t) > t_threshold, "p + t",
            ifelse(this_p < p_threshold, "p only", "t only"))
  cat(sprintf("  FPC%d: t = %+.3f,  p = %.4f  [%s]\n", k, this_t, this_p, flag))
}

# --- Step 2: Auto sign-flip any FPC with a negative regression coefficient --
# We want: positive score on retained FPC  =>  better recruitment
# This is a labelling convention that makes the biology readable.

flip_flags <- coef(mod_lm)[-1][sig_positions] < 0   # TRUE = needs flipping

for (i in seq_along(sig_fpc_nums)) {
  k <- sig_fpc_nums[i]
  if (flip_flags[i]) {
    cat(sprintf("  Flipping FPC%d sign (original coefficient was negative)\n", k))
    fpca_out$harmonics$coefs[, k] <- -fpca_out$harmonics$coefs[, k]
    fpc_scores[, k]               <- -fpc_scores[, k]
  } else {
    cat(sprintf("  FPC%d sign retained (coefficient already positive)\n", k))
  }
}

# --- Step 3: Re-evaluate harmonic functions AFTER any sign flips ------------
fpc_funs <- eval.fd(doy_grid, fpca_out$harmonics[1:n_fpc])

# --- Step 4: Build reduced score dataframe ----------------------------------
sig_fpc_cols <- paste0("FPC", sig_fpc_nums)

# For backward compatibility, set fpcA/fpcB to the top two (or only one)
fpcA_num <- sig_fpc_nums[1]
fpcA_col <- sig_fpc_cols[1]
if (n_sig >= 2) {
  fpcB_num <- sig_fpc_nums[2]
  fpcB_col <- sig_fpc_cols[2]
}

score_df_reduced <- data.frame(y = y_vec)
for (i in seq_along(sig_fpc_nums)) {
  score_df_reduced[[sig_fpc_cols[i]]] <- fpc_scores[, sig_fpc_nums[i]]
}

# --- Step 5: Fit reduced model ----------------------------------------------
fpc_formula    <- as.formula(paste("y ~", paste(sig_fpc_cols, collapse = " + ")))
mod_lm_reduced <- lm(fpc_formula, data = score_df_reduced)

cat("\n--- Reduced Model Summary ---\n")
print(summary(mod_lm_reduced))
cat("\nCoefficients (both positive after sign correction):\n")
print(round(coef(mod_lm_reduced), 4))

# --- Step 6: Back-transform to beta(t) --------------------------------------
fpc_funs_reduced <- sapply(sig_fpc_nums, function(k) {
  as.vector(eval.fd(doy_grid, fpca_out$harmonics[k]))
})
# Ensure it's a matrix even with a single FPC
if (!is.matrix(fpc_funs_reduced)) {
  fpc_funs_reduced <- matrix(fpc_funs_reduced, ncol = 1)
}

b_coefs_reduced <- coef(mod_lm_reduced)[-1]
beta_t_reduced  <- as.vector(fpc_funs_reduced %*% b_coefs_reduced)

# Quick diagnostic plot
plot(doy_grid, beta_t_reduced, type = "l", lwd = 2,
     xlab = "Day of Year", ylab = "β(t)",
     main = sprintf("Reduced model β(t): %s (auto-selected)", paste(sig_fpc_cols, collapse = " + ")))
abline(h = 0, lty = 2, col = "grey50")


# =============================================================================
# Smoothed β(t) and full figure for reduced model
# =============================================================================

library(tidyverse)
library(patchwork)
library(mgcv)

x_scale <- scale_x_continuous(breaks = key_doys, labels = key_labels,
                              expand = c(0.01, 0))
base_theme <- theme_classic(base_size = 12) +
  theme(
    axis.text.x   = element_text(angle = 45, hjust = 1),
    plot.title    = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 10, color = "grey40")
  )

# --- Smooth beta(t) via penalized spline ------------------------------------
beta_raw_df   <- data.frame(doy = doy_grid, beta = beta_t_reduced)
smooth_mod    <- gam(beta ~ s(doy, bs = "ps", k = 20), data = beta_raw_df)
beta_smoothed <- fitted(smooth_mod)

# --- Bootstrap confidence bands on smoothed beta(t) ------------------------
set.seed(42)
n_boot           <- 1000
beta_boot_smooth <- matrix(NA, nrow = n_boot, ncol = 365)

cat("\nBootstrapping smoothed confidence bands...\n")
for (i in 1:n_boot) {
  idx        <- sample(nrow(score_df_reduced), replace = TRUE)
  boot_df    <- score_df_reduced[idx, ]
  boot_lm    <- lm(fpc_formula, data = boot_df)
  boot_coefs <- coef(boot_lm)[-1]
  boot_beta  <- as.vector(fpc_funs_reduced %*% boot_coefs)
  boot_gam   <- gam(beta ~ s(doy, bs = "ps", k = 20),
                    data = data.frame(doy = doy_grid, beta = boot_beta))
  beta_boot_smooth[i, ] <- fitted(boot_gam)
}
cat("Bootstrap complete.\n")

# --- Permutation null envelope on smoothed beta(t) --------------------------
set.seed(123)
n_perm           <- 999
beta_perm_smooth <- matrix(NA, nrow = n_perm, ncol = 365)

cat("Running permutation test...\n")
for (i in 1:n_perm) {
  if (i %% 200 == 0) cat("  Permutation", i, "of", n_perm, "\n")
  perm_df    <- score_df_reduced
  perm_df$y  <- sample(score_df_reduced$y)
  perm_lm    <- lm(fpc_formula, data = perm_df)
  perm_coefs <- coef(perm_lm)[-1]
  perm_beta  <- as.vector(fpc_funs_reduced %*% perm_coefs)
  perm_gam   <- gam(beta ~ s(doy, bs = "ps", k = 20),
                    data = data.frame(doy = doy_grid, beta = perm_beta))
  beta_perm_smooth[i, ] <- fitted(perm_gam)
}
cat("Permutation test complete.\n")

max_dev_smooth    <- apply(abs(beta_perm_smooth), 1, max)
sim_thresh_smooth <- quantile(max_dev_smooth, 0.95)

# --- Assemble plotting dataframe --------------------------------------------
plot_df <- data.frame(
  DOY         = doy_grid,
  beta_raw    = beta_t_reduced,
  beta_smooth = beta_smoothed,
  lower_95    = apply(beta_boot_smooth, 2, quantile, 0.025),
  upper_95    = apply(beta_boot_smooth, 2, quantile, 0.975),
  lower_80    = apply(beta_boot_smooth, 2, quantile, 0.100),
  upper_80    = apply(beta_boot_smooth, 2, quantile, 0.900)
) %>%
  mutate(
    sim_sig_pos = beta_smooth >  sim_thresh_smooth,
    sim_sig_neg = beta_smooth < -sim_thresh_smooth
  )

cat(sprintf("\nSimultaneous significance threshold (smoothed): %.5f\n", sim_thresh_smooth))
cat("Significant positive DOY:", plot_df$DOY[plot_df$sim_sig_pos], "\n")
cat("Significant negative DOY:", plot_df$DOY[plot_df$sim_sig_neg], "\n")

# =============================================================================
# PANEL A: Weighted hydrograph
# =============================================================================
p_wh <- ggplot(wh_df, aes(x = DOY, y = weighted_anomaly)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_ribbon(aes(ymin = pmin(weighted_anomaly, 0), ymax = 0),
              fill = "#d73027", alpha = 0.3) +
  geom_ribbon(aes(ymin = 0, ymax = pmax(weighted_anomaly, 0)),
              fill = "#4575b4", alpha = 0.3) +
  geom_line(color = "black", linewidth = 0.7) +
  x_scale +
  labs(
    title    = "A: Weighted Hydrograph (informal approach)",
    subtitle = "Blue = above-average flow linked to better recruitment  |  Red = linked to worse recruitment",
    x = NULL, y = "Weighted Flow Anomaly"
  ) +
  base_theme

# =============================================================================
# PANEL B: Smoothed beta(t) with uncertainty layers
# =============================================================================
p_beta <- ggplot(plot_df, aes(x = DOY)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_ribbon(aes(ymin = -sim_thresh_smooth, ymax = sim_thresh_smooth),
              fill = "grey85", alpha = 0.7) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95),
              fill = "#74add1", alpha = 0.4) +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80),
              fill = "#4575b4", alpha = 0.4) +
  geom_line(aes(y = beta_raw), color = "grey65", linewidth = 0.4) +
  geom_line(aes(y = beta_smooth), color = "black", linewidth = 1.1) +
  geom_rug(data = filter(plot_df, sim_sig_pos), aes(x = DOY),
           color = "#4575b4", sides = "b", linewidth = 0.6) +
  geom_rug(data = filter(plot_df, sim_sig_neg), aes(x = DOY),
           color = "#d73027", sides = "b", linewidth = 0.6) +
  x_scale +
  labs(
    title    = bquote("B: Smoothed "*beta(t)*" — Reduced Model ("*.(paste(sig_fpc_cols, collapse = " + "))*")"),
    subtitle = paste0(
      "Grey band = simultaneous permutation null (95%)  |  Blue bands = bootstrap 95% / 80% CI\n",
      "Faint line = raw back-transform  |  Coloured rug = significant windows"
    ),
    x = NULL, y = expression(beta(t))
  ) +
  base_theme

# =============================================================================
# PANEL C: All significant FPC harmonics (with flow + recruitment variance)
# =============================================================================
# Compute individual R² for each significant FPC against recruitment
r2_per_fpc <- sapply(seq_along(sig_fpc_nums), function(i) {
  tmp_df <- data.frame(y = y_vec, x = fpc_scores[, sig_fpc_nums[i]])
  summary(lm(y ~ x, data = tmp_df))$r.squared
})

# --- Append reduced-model objects to the saved .RData -----------------------
# This must come after sig_fpc_nums, sig_fpc_cols, plot_df, fpc_funs_reduced,
# r2_per_fpc, mod_lm_reduced, p_threshold, t_threshold all exist.
save(
  mod_lm, fpca_out, beta_df, wh_df,
  beta_boot, beta_perm_mat, sim_threshold,
  n_fpc, fpc_scores, fpc_funs,
  mod_lm_reduced, sig_fpc_nums, sig_fpc_cols, plot_df,
  fpc_funs_reduced, var_explained, r2_per_fpc,
  p_threshold, t_threshold,
  file = paste0("output/models/LL_", SS, "_functional_regression_results.RData")
)
cat("Reduced-model objects appended to .RData\n")

fpc_plot_df <- bind_rows(lapply(seq_along(sig_fpc_nums), function(i) {
  k <- sig_fpc_nums[i]
  vals <- as.vector(eval.fd(doy_grid, fpca_out$harmonics[k]))
  data.frame(
    DOY       = doy_grid,
    val       = vals,
    component = sprintf("%s  (%.1f%% flow var. | %.1f%% recruit var.)",
                        sig_fpc_cols[i], var_explained[k], r2_per_fpc[i] * 100)
  )
}))

# Determine layout: up to 3 columns
ncol_fpc <- min(n_sig, 3)

p_fpc <- ggplot(fpc_plot_df, aes(x = DOY, y = val)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_ribbon(aes(ymin = pmin(val, 0), ymax = 0),   fill = "#d73027", alpha = 0.25) +
  geom_ribbon(aes(ymin = 0, ymax = pmax(val, 0)),    fill = "#4575b4", alpha = 0.25) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ component, ncol = ncol_fpc, scales = "free_y") +
  x_scale +
  labs(
    title    = sprintf("C: Significant Flow Modes Predicting Recruitment (n = %d)", n_sig),
    subtitle = "All components have positive regression coefficients — years scoring high recruit better",
    x = NULL, y = "Harmonic loading"
  ) +
  base_theme +
  theme(strip.text       = element_text(size = 9, face = "bold"),
        strip.background = element_rect(fill = "grey95", color = NA))

# =============================================================================
# PANEL D: Observed vs predicted log(R/S)
# =============================================================================
pred_df <- data.frame(
  Year       = y_df$Year,
  observed   = y_vec,
  predicted  = fitted(mod_lm_reduced),
  fpcA_score = fpc_scores[, fpcA_num]
)

r2_lab <- sprintf("R² = %.2f  |  Adj. R² = %.2f  |  p = %.3f",
                  summary(mod_lm_reduced)$r.squared,
                  summary(mod_lm_reduced)$adj.r.squared,
                  pf(summary(mod_lm_reduced)$fstatistic[1],
                     summary(mod_lm_reduced)$fstatistic[2],
                     summary(mod_lm_reduced)$fstatistic[3],
                     lower.tail = FALSE))

p_scatter <- ggplot(pred_df, aes(x = predicted, y = observed)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(aes(color = fpcA_score), size = 3.5, alpha = 0.85) +
  geom_text(aes(label = Year), size = 2.4, vjust = -0.9, hjust = 0.5, color = "grey30") +
  scale_color_gradient2(
    low = "#d73027", mid = "grey85", high = "#4575b4", midpoint = 0,
    name = sprintf("%s score\n(flood pulse\nmagnitude)", fpcA_col) #####need to relabel
  ) +
  labs(
    title    = "D: Observed vs Predicted log(R/S)",
    subtitle = r2_lab,
    x        = "Predicted log(R/S)",
    y        = "Observed log(R/S)"
  ) +
  base_theme +
  theme(legend.position = "right", legend.title = element_text(size = 9))

# =============================================================================
# Combine into final figure
# =============================================================================
fig_final <- (p_wh | p_scatter) / p_beta / p_fpc +
  plot_layout(heights = c(1, 1.2, max(1, ceiling(n_sig / 3)))) +
  plot_annotation(
    title    = "Flow–Recruitment Relationship: FPCA Functional Regression",
    subtitle = paste0(
      SS, " Brown trout | Single site | n = ", length(y_vec), " years  |  ",
      sprintf("Reduced model: %s (auto-selected, p < %.2f or |t| > %.1f)", 
              paste(sig_fpc_cols, collapse = " + "), p_threshold, t_threshold)
    ),
    theme = theme(
      plot.title    = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10, color = "grey40")
    )
  )

ggsave(
  paste0("output/plots/LL_", SS, "_reduced_model_results.png"),
  fig_final, width = 14, height = 12 + max(0, (ceiling(n_sig / 3) - 1) * 4), dpi = 300
)
cat("\nFigure saved: reduced_model_results.png\n")

# =============================================================================
# Build model summary table page (for multi-page PDF)
# =============================================================================
library(gridExtra)
library(grid)

# Extract full model coefficient table
full_summary       <- summary(mod_lm)
coef_tab           <- as.data.frame(full_summary$coefficients)
coef_tab           <- cbind(Term = rownames(coef_tab), coef_tab)
rownames(coef_tab) <- NULL
# Round numeric columns
coef_tab[, 2:5] <- lapply(coef_tab[, 2:5], function(x) {
  ifelse(abs(x) < 1e-4 & x != 0, format(x, scientific = TRUE, digits = 3),
         formatC(x, format = "f", digits = 4))
})
# Add significance stars
p_numeric <- full_summary$coefficients[, "Pr(>|t|)"]
coef_tab$Sig <- cut(p_numeric,
                    breaks = c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
                    labels = c("***", "**", "*", ".", ""))

# Mark which FPCs were retained in the reduced model
retained_fpcs <- paste0("FPC", sig_fpc_nums)
coef_tab$Retained <- ifelse(coef_tab$Term %in% retained_fpcs, "✓", "")

# Model fit statistics
fit_stats <- data.frame(
  Statistic = c("Residual SE", "Multiple R²", "Adjusted R²",
                "F-statistic", "Numerator df", "Denominator df", "Model p-value",
                "n observations", "n FPCs in full model", "n FPCs retained (reduced)"),
  Value = c(
    formatC(full_summary$sigma, format = "f", digits = 4),
    formatC(full_summary$r.squared, format = "f", digits = 4),
    formatC(full_summary$adj.r.squared, format = "f", digits = 4),
    formatC(full_summary$fstatistic[1], format = "f", digits = 3),
    as.character(full_summary$fstatistic[2]),
    as.character(full_summary$fstatistic[3]),
    formatC(pf(full_summary$fstatistic[1],
               full_summary$fstatistic[2],
               full_summary$fstatistic[3], lower.tail = FALSE),
            format = "f", digits = 4),
    as.character(length(y_vec)),
    as.character(n_fpc),
    as.character(n_sig)
  )
)

# Build grobs
table_theme <- ttheme_minimal(
  core    = list(fg_params = list(cex = 0.75),
                 bg_params = list(fill = c("grey97", "white"), col = NA)),
  colhead = list(fg_params = list(cex = 0.85, fontface = "bold"),
                 bg_params = list(fill = "grey85", col = NA))
)

coef_grob  <- tableGrob(coef_tab, rows = NULL, theme = table_theme)
stats_grob <- tableGrob(fit_stats, rows = NULL, theme = table_theme)

# Header text
header_text <- textGrob(
  paste0(SS, " — Full Model Summary (mod_lm)"),
  gp = gpar(fontsize = 16, fontface = "bold")
)
subheader_text <- textGrob(
  sprintf("Selection criteria: p < %.2f OR |t| > %.1f  |  ✓ = retained in reduced model",
          p_threshold, t_threshold),
  gp = gpar(fontsize = 11, col = "grey40")
)
coef_label <- textGrob("Coefficients (full FPC regression)",
                       gp = gpar(fontsize = 12, fontface = "bold"), hjust = 0, x = 0.02)
stats_label <- textGrob("Model fit statistics",
                        gp = gpar(fontsize = 12, fontface = "bold"), hjust = 0, x = 0.02)
sig_legend <- textGrob(
  "Signif. codes:  *** p<0.001    ** p<0.01    * p<0.05    . p<0.1",
  gp = gpar(fontsize = 9, col = "grey40", fontface = "italic")
)

# Page 2: arrange grobs
fig_summary <- arrangeGrob(
  header_text,
  subheader_text,
  coef_label,
  coef_grob,
  sig_legend,
  stats_label,
  stats_grob,
  ncol    = 1,
  heights = unit.c(
    unit(1.5, "lines"),
    unit(1.2, "lines"),
    unit(1.5, "lines"),
    unit(0.4 + 0.35 * (nrow(coef_tab) + 1), "inches"),
    unit(1.0, "lines"),
    unit(1.5, "lines"),
    unit(0.4 + 0.35 * (nrow(fit_stats) + 1), "inches")
  )
)

# =============================================================================
# Build supplemental page: ALL FPC harmonics from the full model
# =============================================================================
# Show the harmonic loading shape for every FPC included in mod_lm,
# whether or not it was significant. Significant FPCs (retained in the
# reduced model) get a blue strip; non-significant ones get a grey strip.

# Per-FPC stats from the full model summary
full_summary_for_supp <- summary(mod_lm)
all_t_stats <- full_summary_for_supp$coefficients[-1, "t value"]
all_p_vals  <- full_summary_for_supp$coefficients[-1, "Pr(>|t|)"]

# Per-FPC R² against recruitment (one at a time)
all_r2_per_fpc <- sapply(seq_len(n_fpc), function(k) {
  summary(lm(y_vec ~ fpc_scores[, k]))$r.squared
})

# Build long-format data frame of every harmonic
all_fpc_df <- bind_rows(lapply(seq_len(n_fpc), function(k) {
  vals    <- as.vector(eval.fd(doy_grid, fpca_out$harmonics[k]))
  is_sig  <- k %in% sig_fpc_nums
  marker  <- ifelse(is_sig, " *", "")
  data.frame(
    DOY        = doy_grid,
    val        = vals,
    fpc_num    = k,
    is_sig     = is_sig,
    component  = sprintf("FPC%d%s  (%.1f%% flow var. | %.1f%% recruit var. | t = %+.2f, p = %.3f)",
                         k, marker, var_explained[k], all_r2_per_fpc[k] * 100,
                         all_t_stats[k], all_p_vals[k])
  )
}))
# Preserve FPC order in facets
all_fpc_df$component <- factor(all_fpc_df$component,
                               levels = unique(all_fpc_df$component))

# Layout: 3 columns
ncol_supp <- 3

p_all_fpc <- ggplot(all_fpc_df, aes(x = DOY, y = val)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_ribbon(aes(ymin = pmin(val, 0), ymax = 0), fill = "#d73027", alpha = 0.25) +
  geom_ribbon(aes(ymin = 0, ymax = pmax(val, 0)), fill = "#4575b4", alpha = 0.25) +
  geom_line(aes(color = is_sig), linewidth = 0.8) +
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "grey55"),
                     guide = "none") +
  facet_wrap(~ component, ncol = ncol_supp, scales = "free_y") +
  x_scale +
  labs(
    title    = sprintf("Supplemental: All %d FPC Harmonics from the Full Model", n_fpc),
    subtitle = sprintf(
      "Black lines = significant (retained in reduced model, marked with *)  |  Grey lines = not significant\nSelection criteria: p < %.2f OR |t| > %.1f",
      p_threshold, t_threshold),
    x = "Day of Year", y = "Harmonic loading"
  ) +
  base_theme +
  theme(strip.text       = element_text(size = 8, face = "bold"),
        strip.background = element_rect(fill = "grey95", color = NA),
        plot.title       = element_text(face = "bold", size = 13))

# =============================================================================
# Export per-page artifacts: PNGs for the website + per-page PDFs for the
# multi-page PDF (each page sized to its own content to eliminate whitespace gaps)
# =============================================================================

# Output directories
docs_dir   <- "docs"                              # GitHub Pages root
assets_dir <- file.path(docs_dir, "assets", "sites", SS)
dir.create(assets_dir, recursive = TRUE, showWarnings = FALSE)

# --- Page heights, each sized to its own content ----------------------------
page1_w <- 14
page1_h <- 12 + max(0, (ceiling(n_sig / 3) - 1) * 4)

page2_w <- 14
# Estimate summary page height from row count (intercept + n_fpc rows + fit-stat rows + headers)
n_coef_rows <- nrow(coef_tab)
n_stat_rows <- nrow(fit_stats)
page2_h     <- 1.5 + 0.30 * n_coef_rows + 0.30 * n_stat_rows + 2.5  # inches

page3_w <- 14
n_rows_supp <- ceiling(n_fpc / ncol_supp)
page3_h     <- max(6, 2.5 + n_rows_supp * 1.8)

# --- 1. PNGs for the website (one per page) ---------------------------------
ggsave(file.path(assets_dir, "page1_panels.png"),
       fig_final, width = page1_w, height = page1_h, dpi = 150, limitsize = FALSE)

png(file.path(assets_dir, "page2_summary.png"),
    width = page2_w, height = page2_h, units = "in", res = 150)
grid.draw(fig_summary)
dev.off()

ggsave(file.path(assets_dir, "page3_all_fpcs.png"),
       p_all_fpc, width = page3_w, height = page3_h, dpi = 150, limitsize = FALSE)

cat(sprintf("PNG pages saved to: %s\n", assets_dir))

# --- 2. Multi-page PDF, each page sized independently -----------------------
# Strategy: write each page to its own temp PDF, then concatenate with qpdf::pdf_combine.
# This is the only reliable way to get heterogeneous page sizes in a single PDF
# from base R, since pdf() forces all pages to share the device dimensions.
pdf_path <- paste0("output/plots/LL_", SS, "_reduced_model_results.pdf")

tmp_p1 <- tempfile(fileext = ".pdf")
tmp_p2 <- tempfile(fileext = ".pdf")
tmp_p3 <- tempfile(fileext = ".pdf")

# Page 1
pdf(tmp_p1, width = page1_w, height = page1_h)
print(fig_final)
dev.off()

# Page 2
pdf(tmp_p2, width = page2_w, height = page2_h)
grid.draw(fig_summary)
dev.off()

# Page 3
pdf(tmp_p3, width = page3_w, height = page3_h)
print(p_all_fpc)
dev.off()

# Merge into one PDF with no whitespace gaps
if (requireNamespace("qpdf", quietly = TRUE)) {
  qpdf::pdf_combine(c(tmp_p1, tmp_p2, tmp_p3), output = pdf_path)
  file.remove(tmp_p1, tmp_p2, tmp_p3)
  cat(sprintf("Multi-page PDF (gap-free) saved: %s\n", pdf_path))
} else {
  warning("Package 'qpdf' not installed; falling back to single-device PDF with potential whitespace.\n",
          "  Run install.packages('qpdf') to enable gap-free multi-page output.")
  pdf(pdf_path, width = page1_w, height = max(page1_h, page2_h, page3_h))
  print(fig_final)
  grid.newpage(); grid.draw(fig_summary)
  grid.newpage(); print(p_all_fpc)
  dev.off()
}

# --- 3. Site metadata as JSON (consumed by the website) ---------------------
# Encode the model summary so the website can show it as an HTML table
site_json <- list(
  site            = SS,
  n_years         = length(y_vec),
  n_fpc_full      = n_fpc,
  n_fpc_retained  = n_sig,
  retained_fpcs   = sig_fpc_cols,
  p_threshold     = p_threshold,
  t_threshold     = t_threshold,
  r_squared       = unname(round(summary(mod_lm)$r.squared, 4)),
  adj_r_squared   = unname(round(summary(mod_lm)$adj.r.squared, 4)),
  reduced_r2      = unname(round(summary(mod_lm_reduced)$r.squared, 4)),
  reduced_adj_r2  = unname(round(summary(mod_lm_reduced)$adj.r.squared, 4)),
  pages           = list(
    panels      = "page1_panels.png",
    summary     = "page2_summary.png",
    all_fpcs    = "page3_all_fpcs.png"
  )
)

if (requireNamespace("jsonlite", quietly = TRUE)) {
  jsonlite::write_json(site_json,
                       file.path(assets_dir, "metadata.json"),
                       auto_unbox = TRUE, pretty = TRUE)
} else {
  warning("Package 'jsonlite' not installed; skipping metadata.json.\n",
          "  Run install.packages('jsonlite') to enable.")
}

# =============================================================================
# Partial R² decomposition
# =============================================================================
r2_individual <- sapply(sig_fpc_cols, function(col) {
  summary(lm(as.formula(paste("y ~", col)), data = score_df_reduced))$r.squared
})
r2_both <- summary(mod_lm_reduced)$r.squared
r2_adj  <- summary(mod_lm_reduced)$adj.r.squared

cat("=====================================================\n")
cat("Variance in log(R/S) explained:\n")
cat("=====================================================\n")
for (i in seq_along(sig_fpc_cols)) {
  cat(sprintf("%s alone:               %.1f%%\n", sig_fpc_cols[i], r2_individual[i] * 100))
}
cat(sprintf("All combined (R²):     %.1f%%\n",  r2_both  * 100))
cat(sprintf("All combined (adj R²): %.1f%%\n", r2_adj   * 100))
cat(sprintf("Sum of individual R²s (confirms orthogonality): %.1f%%\n",
            sum(r2_individual) * 100))
cat("=====================================================\n")
cat("Of the variance explained by the model:\n")
for (i in seq_along(sig_fpc_cols)) {
  cat(sprintf("  Attributable to %s:  %.1f%%\n", sig_fpc_cols[i], r2_individual[i] / r2_both * 100))
}
cat("=====================================================\n")

# --- Individual scatter plots (base R) --------------------------------------
par(mfrow = c(1, n_sig))

for (i in seq_along(sig_fpc_nums)) {
  k <- sig_fpc_nums[i]
  plot(fpc_scores[, k], y_vec,
       xlab = sprintf("%s score", sig_fpc_cols[i]),
       ylab = "log(R/S)",
       main = sprintf("%s (rank %d)", sig_fpc_cols[i], i),
       pch = 19, col = "#4575b4")
  abline(lm(y_vec ~ fpc_scores[, k]), col = "#d73027", lwd = 2)
  text(fpc_scores[, k], y_vec, labels = years_vec, cex = 0.7, pos = 3, col = "grey40")
}

# =============================================================================
# FPC score time series
# =============================================================================
fpc_ts_df <- data.frame(Year = years_vec, logRS = y_vec)
for (i in seq_along(sig_fpc_nums)) {
  fpc_ts_df[[sig_fpc_cols[i]]] <- fpc_scores[, sig_fpc_nums[i]]
}

fpc_long <- fpc_ts_df %>%
  pivot_longer(
    cols      = all_of(sig_fpc_cols),
    names_to  = "component",
    values_to = "score"
  ) %>%
  mutate(
    component = {
      labels <- setNames(
        sapply(seq_along(sig_fpc_nums), function(i) {
          sprintf("%s — flow mode %d (%.1f%% flow var. | %.1f%% recruit var.)",
                  sig_fpc_cols[i], i, var_explained[sig_fpc_nums[i]], r2_per_fpc[i] * 100)
        }),
        sig_fpc_cols
      )
      labels[component]
    }
  )

p_ts <- ggplot(fpc_long, aes(x = Year, y = score)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_col(aes(fill = score > 0), width = 0.7, alpha = 0.8) +
  geom_line(
    data = fpc_long %>%
      group_by(component) %>%
      mutate(scaled_logRS = scale(logRS)[, 1] * (max(abs(score)) * 0.8)) %>%
      ungroup(),
    aes(y = scaled_logRS),
    color = "black", linewidth = 0.7
  ) +
  scale_fill_manual(
    values = c("TRUE" = "#4575b4", "FALSE" = "#d73027"),
    labels = c("TRUE" = "Positive", "FALSE" = "Negative"),
    name   = "Score direction"
  ) +
  scale_x_continuous(breaks = seq(min(years_vec), max(years_vec), by = 3)) +
  facet_wrap(~ component, ncol = 1, scales = "free_y") +
  labs(
    title    = "Time Series of Significant FPC Scores",
    subtitle = "Bars = FPC score by year  |  Black line = z-scored log(R/S) for reference",
    x = "Year", y = "FPC Score"
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.text       = element_text(size = 10, face = "bold"),
    strip.background = element_rect(fill = "grey95", color = NA),
    legend.position  = "bottom",
    axis.text.x      = element_text(angle = 45, hjust = 1),
    plot.title       = element_text(face = "bold")
  )

print(p_ts)

ggsave(
  paste0("output/plots/LL_", SS, "_fpc_timeseries.png"),
  p_ts, width = 12, height = 3 + n_sig * 2.5, dpi = 300
)
cat("Saved: fpc_timeseries.png\n")

# =============================================================================
# Composite flow score
# =============================================================================
std_formula <- as.formula(paste(
  "y ~", paste(sprintf("scale(%s)", sig_fpc_cols), collapse = " + ")
))
b_std <- coef(lm(std_formula, data = score_df_reduced))[-1]
names(b_std) <- sig_fpc_cols

cat("\nStandardized coefficients:\n")
for (col in sig_fpc_cols) {
  cat(sprintf("  %s: %.3f\n", col, b_std[col]))
}

# Standardize all significant FPC scores
fpc_std_list <- lapply(sig_fpc_nums, function(k) scale(fpc_scores[, k])[, 1])
names(fpc_std_list) <- sig_fpc_cols

composite_score <- Reduce("+", mapply(function(b, s) b * s, b_std, fpc_std_list, SIMPLIFY = FALSE))
composite_z     <- scale(composite_score)[, 1]

composite_df <- data.frame(
  Year      = years_vec,
  composite = as.vector(composite_z),
  logRS_z   = scale(y_vec)[, 1]
)
# Add individual z-scored FPC columns
for (col in sig_fpc_cols) {
  composite_df[[paste0(col, "_z")]] <- fpc_std_list[[col]]
}
composite_df <- composite_df %>%
  mutate(flow_quality = ifelse(composite > 0, "Flow favourable", "Flow unfavourable"))

r_composite <- cor(composite_df$composite, composite_df$logRS_z)
cat(sprintf("\nCorrelation of composite flow score with log(R/S): r = %.3f\n", r_composite))

p_composite <- ggplot(composite_df, aes(x = Year, y = composite)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_col(aes(fill = flow_quality), width = 0.7, alpha = 0.85) +
  geom_line(aes(y = logRS_z), color = "black", linewidth = 0.8) +
  geom_point(aes(y = logRS_z), color = "black", size = 2,
             shape = 21, fill = "white", stroke = 0.8) +
  scale_fill_manual(
    values = c("Flow favourable" = "#4575b4", "Flow unfavourable" = "#d73027"),
    name   = NULL
  ) +
  scale_x_continuous(breaks = seq(min(years_vec), max(years_vec), by = 2)) +
  annotate("text",
           x     = min(years_vec) + 1,
           y     = max(c(composite_df$composite, composite_df$logRS_z)) * 0.92,
           label = sprintf("r = %.2f", r_composite),
           size  = 4, color = "grey30", hjust = 0) +
  labs(
    title    = "Composite Flow Score vs Observed Recruitment",
    subtitle = paste0(
      "Bars = composite hydrograph score (",
      paste(sprintf("%s × %.2f", sig_fpc_cols, b_std[sig_fpc_cols]), collapse = " + "),
      ", all z-scored)\n",
      "Line = observed log(R/S)  |  Both standardized to z-scores"
    ),
    x = "Year", y = "Standardized Score"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "bottom",
    axis.text.x     = element_text(angle = 45, hjust = 1),
    plot.title      = element_text(face = "bold", size = 13),
    plot.subtitle   = element_text(size = 10, color = "grey40")
  )

# --- Stacked bar decomposition ----------------------------------------------
# Generate a color palette for N components
fpc_colors <- colorRampPalette(c("#4575b4", "#74add1", "#abd9e9", "#fee090"))(n_sig)

component_df <- bind_rows(lapply(seq_along(sig_fpc_cols), function(i) {
  data.frame(
    Year         = years_vec,
    component    = sprintf("%s: flow mode %d", sig_fpc_cols[i], i),
    contribution = b_std[sig_fpc_cols[i]] * fpc_std_list[[sig_fpc_cols[i]]]
  )
}))

color_map <- setNames(fpc_colors, 
  sapply(seq_along(sig_fpc_cols), function(i) sprintf("%s: flow mode %d", sig_fpc_cols[i], i)))

p_stacked <- ggplot(component_df, aes(x = Year, y = contribution, fill = component)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_col(width = 0.7, alpha = 0.85, position = "stack") +
  geom_line(data = composite_df, aes(x = Year, y = logRS_z),
            inherit.aes = FALSE, color = "black", linewidth = 0.8) +
  geom_point(data = composite_df, aes(x = Year, y = logRS_z),
             inherit.aes = FALSE, color = "black", size = 2,
             shape = 21, fill = "white", stroke = 0.8) +
  scale_fill_manual(values = color_map, name = "Flow component") +
  scale_x_continuous(breaks = seq(min(years_vec), max(years_vec), by = 2)) +
  labs(
    title    = "Decomposed Flow Contributions to Recruitment by Year",
    subtitle = "Stacked bars show relative contribution of each flow mode  |  Line = observed log(R/S)",
    x = "Year", y = "Weighted Flow Score (z-scored)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position  = "bottom",
    axis.text.x      = element_text(angle = 45, hjust = 1),
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(size = 10, color = "grey40")
  )

fig_ts <- p_composite / p_stacked +
  plot_annotation(
    title    = "Hydrograph Quality Score: Flow Conditions Supporting Recruitment",
    subtitle = paste0(
      SS, " Brown trout | Single site | n = ", length(years_vec), " years  |  ",
      sprintf("Composite score derived from FPCA functional regression (%s)",
              paste(sig_fpc_cols, collapse = " + "))
    ),
    theme = theme(
      plot.title    = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10, color = "grey40")
    )
  )

print(fig_ts)

ggsave(
  paste0("output/plots/LL_", SS, "_composite_flow_score_timeseries.png"),
  fig_ts, width = 12, height = 9, dpi = 300
)
cat("Saved: composite_flow_score_timeseries.png\n")

# --- Final environment summary ----------------------------------------------
cat("\nObjects available for further analysis:\n")
cat("  plot_df           — smoothed beta(t) with all CI columns\n")
cat("  pred_df           — observed vs predicted per year\n")
cat("  sim_thresh_smooth — simultaneous significance threshold\n")
cat("  composite_df      — composite flow score per year\n")
cat(sprintf("  sig_fpc_nums — auto-selected FPC indices (%s)\n",
            paste(sig_fpc_nums, collapse = ", ")))
print(fig_final)

return(list(fig_final,fig_ts))
})

# =============================================================================
# Build site index for the GitHub Pages website
# =============================================================================
# Scans docs/assets/sites/ for completed sites and writes docs/sites_index.json,
# which is consumed by the dropdown menu on index.html.

if (requireNamespace("jsonlite", quietly = TRUE)) {
  docs_dir       <- "docs"
  sites_root_dir <- file.path(docs_dir, "assets", "sites")

  if (dir.exists(sites_root_dir)) {
    site_dirs <- list.dirs(sites_root_dir, recursive = FALSE)
    site_index <- lapply(site_dirs, function(d) {
      meta_path <- file.path(d, "metadata.json")
      if (file.exists(meta_path)) {
        meta <- jsonlite::read_json(meta_path)
        list(
          id             = basename(d),
          site           = meta$site,
          n_years        = meta$n_years,
          n_fpc_retained = meta$n_fpc_retained,
          r_squared      = meta$r_squared
        )
      } else NULL
    })
    site_index <- Filter(Negate(is.null), site_index)
    # Sort alphabetically by site ID
    site_index <- site_index[order(sapply(site_index, function(x) x$id))]

    jsonlite::write_json(
      list(sites = site_index, generated = format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      file.path(docs_dir, "sites_index.json"),
      auto_unbox = TRUE, pretty = TRUE
    )
    cat(sprintf("\nSite index written: %s (%d sites)\n",
                file.path(docs_dir, "sites_index.json"), length(site_index)))
  }
}
