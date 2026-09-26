# Step A — Setup and Load Data

library(tidyverse)
library(FactoMineR)
library(factoextra)
library(mixOmics)
library(ggrepel)
library(openxlsx)
library(knitr)
library(kableExtra)

# Set working directory
setwd("C:/Users/WiNDOWS/Desktop/Metabolomics-analysis/Positive/Supernatant")

# Load peak matrix
raw_sup <- read.csv("BNC5_supernatant_pos_peakmatrix.csv",
                    row.names   = 1,
                    check.names = FALSE)

# Inspect
cat("Dimensions (features x samples):", dim(raw_sup), "\n")
cat("Column names:\n")
print(colnames(raw_sup))
cat("Total zeros:", sum(raw_sup == 0, na.rm = TRUE), "\n")
cat("Total NAs:  ", sum(is.na(raw_sup)), "\n")
print(raw_sup[1:3, 1:6])



# Step B — Zero Conversion and Sample Metadata

# Convert zeros to NA
data_sup_na <- raw_sup
data_sup_na[data_sup_na == 0] <- NA

cat("NAs after zero conversion:", sum(is.na(data_sup_na)), "\n")

# Define sample metadata
sample_info_sup <- data.frame(
  sample = colnames(data_sup_na),
  group  = c(rep("BNC5", 6),
             rep("BNC5RSw", 6),
             rep("RSw", 6)),
  lot    = c(rep("Lot1", 3), rep("Lot2", 3),
             rep("Lot1", 3), rep("Lot2", 3),
             rep("Lot1", 3), rep("Lot2", 3)),
  stringsAsFactors = FALSE
)

cat("\nSample metadata:\n")
print(sample_info_sup)


# Step C — Critical Batch Check

cat("\nBatch balance check:\n")
print(table(sample_info_sup$group, sample_info_sup$lot))


# Step D — Missing Value Filter

# Keep features detected in ≥50% of samples in at least one group
keep_feature_sup <- function(feature_row, sample_info) {
  groups <- unique(sample_info$group)
  for (grp in groups) {
    cols       <- sample_info$sample[sample_info$group == grp]
    n_detected <- sum(!is.na(feature_row[cols]))
    n_total    <- length(cols)
    if (n_detected >= ceiling(n_total * 0.5)) return(TRUE)
  }
  return(FALSE)
}

cat("Applying 50% per-group presence filter...\n")
keep_sup <- apply(data_sup_na, 1, 
                  keep_feature_sup, 
                  sample_info = sample_info_sup)
data_sup_filtered <- data_sup_na[keep_sup, ]

cat("Features before filter:", nrow(data_sup_na), "\n")
cat("Features after filter: ", nrow(data_sup_filtered), "\n")
cat("Features removed:      ", 
    nrow(data_sup_na) - nrow(data_sup_filtered), "\n")



# Step E — Imputation

data_sup_imputed <- data_sup_filtered

for (feat in rownames(data_sup_imputed)) {
  row_vals <- as.numeric(data_sup_imputed[feat, ])
  if (any(is.na(row_vals))) {
    min_val  <- min(row_vals, na.rm = TRUE)
    row_vals[is.na(row_vals)] <- min_val / 2
    data_sup_imputed[feat, ] <- row_vals
  }
}

cat("NAs after imputation:", sum(is.na(data_sup_imputed)), "\n")
# Must be 0 — do not proceed if any NAs remain


# Step F — TIC Normalisation, Log2, Autoscaling

# TIC normalisation
tic_sup    <- colSums(data_sup_imputed)
cat("TIC values per sample:\n")
print(round(tic_sup))

data_sup_tic <- sweep(data_sup_imputed, 2, tic_sup, "/")
cat("TIC after normalisation (all should be 1.0):\n")
print(round(colSums(data_sup_tic), 4))

# Log2 transformation
data_sup_log2 <- log2(data_sup_tic)
cat("Log2 range:", round(min(data_sup_log2), 2), 
    "to", round(max(data_sup_log2), 2), "\n")

# Autoscaling
data_sup_scaled <- t(scale(t(data_sup_log2)))
cat("Mean of first feature after scaling:", 
    round(mean(data_sup_scaled[1,]), 6), "\n")
cat("SD of first feature after scaling:  ", 
    round(sd(data_sup_scaled[1,]), 6), "\n")

# Rebuild full data frames with metadata
sup_log2   <- cbind(sample_info_sup[, c("sample","group","lot")],
                    as.data.frame(t(data_sup_log2)))
sup_scaled <- cbind(sample_info_sup[, c("sample","group","lot")],
                    as.data.frame(t(data_sup_scaled)))

# Transpose for PCA and PLS-DA
X_sup_log2   <- t(data_sup_log2)
X_sup_scaled <- t(data_sup_scaled)

# Remove zero variance features
zero_var_sup       <- apply(X_sup_scaled, 2, var, na.rm=TRUE) == 0
X_sup_clean        <- X_sup_scaled[, !zero_var_sup]
cat("Zero variance features removed:", sum(zero_var_sup), "\n")
cat("Final matrix dimensions:", dim(X_sup_clean), "\n")

groups_sup <- sample_info_sup$group
lots_sup   <- sample_info_sup$lot



# Step G — Before/After Normalisation Boxplot

# Convert to long format for proper boxplot across all samples
raw_long <- data_sup_imputed %>%
  as.data.frame() %>%
  rownames_to_column("Feature") %>%
  pivot_longer(-Feature, names_to = "Sample", values_to = "Intensity") %>%
  left_join(sample_info_sup, by = c("Sample" = "sample"))

log2_long <- data_sup_log2 %>%
  as.data.frame() %>%
  rownames_to_column("Feature") %>%
  pivot_longer(-Feature, names_to = "Sample", values_to = "Intensity") %>%
  left_join(sample_info_sup, by = c("Sample" = "sample"))

scaled_long <- data_sup_scaled %>%
  as.data.frame() %>%
  rownames_to_column("Feature") %>%
  pivot_longer(-Feature, names_to = "Sample", values_to = "Intensity") %>%
  left_join(sample_info_sup, by = c("Sample" = "sample"))

# Plot sample-level boxplots — one box per sample
p_raw <- ggplot(raw_long,
                aes(x = Sample, y = Intensity, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  scale_fill_manual(values = c("BNC5"="#2ca02c",
                               "BNC5RSw"="#1f77b4",
                               "RSw"="#d62728")) +
  labs(title = "Raw Data", y = "Intensity", x = "") +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7),
        legend.position = "none")

p_log2 <- ggplot(log2_long,
                 aes(x = Sample, y = Intensity, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  scale_fill_manual(values = c("BNC5"="#2ca02c",
                               "BNC5RSw"="#1f77b4",
                               "RSw"="#d62728")) +
  labs(title = "Log2 Transformed", y = "Log2 Intensity", x = "") +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7),
        legend.position = "none")

p_scaled <- ggplot(scaled_long,
                   aes(x = Sample, y = Intensity, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  scale_fill_manual(values = c("BNC5"="#2ca02c",
                               "BNC5RSw"="#1f77b4",
                               "RSw"="#d62728")) +
  labs(title = "Autoscaled", y = "Z-Score", x = "") +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7),
        legend.position = "none")

library(patchwork)
p_raw + p_log2 + p_scaled

ggsave("Supernatant_Normalisation_Check.pdf",
       width = 14, height = 5, dpi = 300)
cat("Normalisation check plot saved\n")



# Step H — PCA Quality Control

pca_sup <- PCA(X_sup_clean, graph = FALSE, scale.unit = FALSE)

var_exp_sup <- round(pca_sup$eig[, 2], 1)
cat("PC1:", var_exp_sup[1], "%\n")
cat("PC2:", var_exp_sup[2], "%\n")
cat("Combined:", var_exp_sup[1] + var_exp_sup[2], "%\n")

# PCA by Lot — batch check
p_lot_sup <- fviz_pca_ind(pca_sup,
                          geom.ind      = "point",
                          col.ind       = lots_sup,
                          palette       = c("#ff7f0e", "#9467bd"),
                          addEllipses   = TRUE,
                          ellipse.level = 0.95,
                          legend.title  = "Lot",
                          title = "Supernatant PCA — Lot Effect Check") +
  theme_bw(base_size = 12)

print(p_lot_sup)
ggsave("Supernatant_PCA_Lot.pdf", plot = p_lot_sup,
       width = 8, height = 6, dpi = 300)

# PCA by Group — biological separation
p_group_sup <- fviz_pca_ind(pca_sup,
                            geom.ind      = "point",
                            col.ind       = groups_sup,
                            palette       = c("#2ca02c","#1f77b4","#d62728"),
                            addEllipses   = TRUE,
                            ellipse.level = 0.95,
                            legend.title  = "Group",
                            title = paste0("Supernatant PCA — PC1: ",
                                           var_exp_sup[1], "%, PC2: ",
                                           var_exp_sup[2], "%")) +
  theme_bw(base_size = 12)

print(p_group_sup)
ggsave("Supernatant_PCA_Group.pdf", plot = p_group_sup,
       width = 8, height = 6, dpi = 300)

cat("PCA plots saved\n")



# Step I — ANOVA Across All Features

cat("Running ANOVA across", ncol(X_sup_clean), "features...\n")
cat("This may take 2-5 minutes...\n")

feat_names_sup  <- colnames(X_sup_clean)
groups_sup_fac  <- as.factor(groups_sup)

anova_sup <- data.frame(
  Feature    = feat_names_sup,
  ANOVA_p    = NA_real_,
  ANOVA_p_BH = NA_real_
)

for (i in seq_along(feat_names_sup)) {
  feat  <- feat_names_sup[i]
  model <- aov(X_sup_clean[, feat] ~ groups_sup_fac)
  anova_sup$ANOVA_p[i] <-
    summary(model)[[1]]["groups_sup_fac", "Pr(>F)"]
}

anova_sup$ANOVA_p_BH <- p.adjust(anova_sup$ANOVA_p, method = "BH")
anova_sup$Significant <- anova_sup$ANOVA_p_BH < 0.05
anova_sup <- anova_sup[order(anova_sup$ANOVA_p_BH), ]

cat("ANOVA complete\n")
cat("Total features tested:", nrow(anova_sup), "\n")
cat("Significant (BH < 0.05):", sum(anova_sup$Significant), "\n")
cat("\nTop 10 most significant features:\n")
print(head(anova_sup, 10))



# Step J — Pairwise Comparisons

# Extract log2 matrices per group for significant features only
sig_feats_sup <- anova_sup$Feature[anova_sup$Significant]
cat("Running pairwise comparisons on", 
    length(sig_feats_sup), "features...\n")

log2_BNC5_sup    <- X_sup_log2[groups_sup == "BNC5",
                               sig_feats_sup, drop = FALSE]
log2_BNC5RSw_sup <- X_sup_log2[groups_sup == "BNC5RSw",
                               sig_feats_sup, drop = FALSE]
log2_RSw_sup     <- X_sup_log2[groups_sup == "RSw",
                               sig_feats_sup, drop = FALSE]

# Reusable comparison function
run_comp_sup <- function(g1_data, g2_data, feat_names) {
  results <- data.frame(
    Feature = feat_names,
    mean_g1 = colMeans(g1_data, na.rm = TRUE),
    mean_g2 = colMeans(g2_data, na.rm = TRUE),
    log2FC  = colMeans(g2_data, na.rm=TRUE) -
      colMeans(g1_data, na.rm=TRUE),
    p_value = NA_real_
  )
  for (feat in feat_names) {
    results$p_value[results$Feature == feat] <-
      t.test(g1_data[, feat], g2_data[, feat])$p.value
  }
  results$p_BH        <- p.adjust(results$p_value, method = "BH")
  results$Significant <- results$p_BH < 0.05
  results$raw_FC      <- 2^results$log2FC
  results$neg_log10p  <- -log10(results$p_value)
  return(results[order(results$p_BH), ])
}

# Three comparisons
comp_sup_treatment  <- run_comp_sup(log2_BNC5_sup,
                                    log2_BNC5RSw_sup,
                                    sig_feats_sup)
comp_sup_background <- run_comp_sup(log2_RSw_sup,
                                    log2_BNC5RSw_sup,
                                    sig_feats_sup)
comp_sup_sanity     <- run_comp_sup(log2_BNC5_sup,
                                    log2_RSw_sup,
                                    sig_feats_sup)

cat("\n── BNC5RSw vs BNC5 (primary) ──\n")
cat("Significant:", sum(comp_sup_treatment$Significant), "\n")
cat("Elevated (log2FC > 0.58):",
    sum(comp_sup_treatment$Significant &
          comp_sup_treatment$log2FC > 0.58), "\n")
cat("Reduced (log2FC < -0.58):",
    sum(comp_sup_treatment$Significant &
          comp_sup_treatment$log2FC < -0.58), "\n")

cat("\n── BNC5RSw vs RSw (background check) ──\n")
cat("Significant:", sum(comp_sup_background$Significant), "\n")

cat("\n── BNC5 vs RSw (sanity check) ──\n")
cat("Significant:", sum(comp_sup_sanity$Significant), "\n")



# Step K — RSw Background Filter

# Step 1: Features elevated in BNC5RSw vs BNC5
elevated_sup <- comp_sup_treatment %>%
  filter(Significant == TRUE & log2FC > 0.58)

# Step 2: Features also in RSw background
rsw_bg_sup <- comp_sup_background %>%
  filter(Significant == TRUE & log2FC > 0) %>%
  pull(Feature)

# Step 3: True BNC5 secreted response
true_response_sup <- elevated_sup %>%
  filter(!Feature %in% rsw_bg_sup)

# Reduced features
reduced_sup <- comp_sup_treatment %>%
  filter(Significant == TRUE & log2FC < -0.58)

cat("── Supernatant RSw Background Filter ──\n")
cat("Elevated in BNC5RSw vs BNC5:          ",
    nrow(elevated_sup), "\n")
cat("Of which from RSw background:          ",
    sum(elevated_sup$Feature %in% rsw_bg_sup), "\n")
cat("TRUE BNC5 secreted response features: ",
    nrow(true_response_sup), "\n")
cat("Reduced in BNC5RSw vs BNC5:           ",
    nrow(reduced_sup), "\n")

cat("\nTop 10 true secreted response features:\n")
print(head(true_response_sup[,
                             c("Feature","log2FC","raw_FC","p_BH")], 10))



# Step L — PLS-DA and VIP Scores

true_feats_sup <- true_response_sup$Feature

X_plsda_sup <- X_sup_clean[, true_feats_sup, drop = FALSE]
Y_plsda_sup <- as.factor(groups_sup)

cat("PLS-DA matrix:", dim(X_plsda_sup), "\n")

plsda_sup <- plsda(X_plsda_sup, Y_plsda_sup, ncomp = 2)

# Scores plot
plotIndiv(plsda_sup,
          comp      = c(1, 2),
          group     = Y_plsda_sup,
          ind.names = FALSE,
          ellipse   = TRUE,
          legend    = TRUE,
          title     = "PLS-DA: True BNC5 Secreted Response",
          col       = c("#2ca02c", "#1f77b4", "#d62728"))

# VIP scores
vip_sup <- vip(plsda_sup)

vip_df_sup <- data.frame(
  Feature   = rownames(vip_sup),
  VIP_Comp1 = vip_sup[, 1]
)
vip_df_sup <- vip_df_sup[order(-vip_df_sup$VIP_Comp1), ]

cat("VIP > 1.0:", sum(vip_df_sup$VIP_Comp1 > 1.0), "\n")
cat("VIP > 1.5:", sum(vip_df_sup$VIP_Comp1 > 1.5), "\n")
cat("\nTop 15 features by VIP:\n")
print(head(vip_df_sup, 15))

# VIP bar plot
ggplot(head(vip_df_sup, 30),
       aes(x = reorder(Feature, VIP_Comp1),
           y = VIP_Comp1,
           fill = VIP_Comp1 > 1.0)) +
  geom_col() + coord_flip() +
  geom_hline(yintercept = 1.0, linetype="dashed",
             colour="red", linewidth=0.8) +
  scale_fill_manual(values = c("grey70","#1f77b4"),
                    labels = c("VIP < 1","VIP > 1")) +
  labs(title = "PLS-DA VIP Scores — Top 30 Secreted Response Features",
       x = "Feature (m/z_RT)", y = "VIP Score", fill = "Importance") +
  theme_bw(base_size = 10) +
  theme(axis.text.y = element_text(size = 7))

ggsave("Supernatant_VIP_Scores.pdf",
       width = 10, height = 8, dpi = 300)
cat("VIP plot saved\n")



# Step M — Volcano Plot

volcano_sup <- comp_sup_treatment %>%
  mutate(
    Direction = case_when(
      Significant & log2FC >  0.58 &
        !Feature %in% rsw_bg_sup ~ "True BNC5 Secreted",
      Significant & log2FC >  0.58 &
        Feature %in% rsw_bg_sup ~ "RSw Background",
      Significant & log2FC < -0.58 ~ "Reduced in BNC5RSw",
      TRUE ~ "Not Significant"
    ),
    Label = ifelse(Feature %in%
                     head(true_response_sup$Feature, 15),
                   Feature, "")
  )

cat("Volcano categories:\n")
print(table(volcano_sup$Direction))

ggplot(volcano_sup,
       aes(x = log2FC, y = neg_log10p,
           colour = Direction, label = Label)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_text_repel(size = 2.8, max.overlaps = 20,
                  colour = "black", box.padding = 0.3) +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed", colour = "grey40") +
  geom_vline(xintercept = c(-0.58, 0.58),
             linetype = "dashed", colour = "grey40") +
  scale_colour_manual(values = c(
    "True BNC5 Secreted" = "#1f77b4",
    "RSw Background"     = "#ff7f0e",
    "Reduced in BNC5RSw" = "#d62728",
    "Not Significant"    = "grey70"
  )) +
  labs(
    title    = "Volcano Plot: BNC5RSw vs BNC5 — Supernatant",
    subtitle = paste0("Blue = True BNC5 secreted response (n=",
                      nrow(true_response_sup),
                      ") | Orange = RSw background removed"),
    x        = "Log2 Fold Change (BNC5RSw vs BNC5)",
    y        = "-log10(BH-adjusted p-value)",
    colour   = "Category"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

ggsave("Supernatant_Volcano.pdf",
       width = 10, height = 7, dpi = 300)
cat("Volcano plot saved\n")



# Step N — Master Excel Output with EIC Values

# Extract raw EIC areas for true response features
eic_sup_raw        <- data_sup_imputed[true_feats_sup, ]
eic_sup_df         <- as.data.frame(round(eic_sup_raw, 0))
eic_sup_df$Feature <- rownames(eic_sup_df)

# Build master table
master_sup <- true_response_sup %>%
  dplyr::select(Feature, log2FC, raw_FC, p_BH) %>%
  left_join(
    vip_df_sup %>% rename(VIP_Score = VIP_Comp1),
    by = "Feature"
  ) %>%
  mutate(
    mz          = as.numeric(sub("M(.+)T.+", "\\1", Feature)),
    RT_min      = round(as.numeric(
      sub("M.+T(.+)", "\\1", Feature)) / 60, 2),
    Fold_Change = round(raw_FC, 2),
    Log2FC      = round(log2FC, 3),
    VIP_Score   = round(VIP_Score, 3),
    p_BH        = signif(p_BH, 4)
  ) %>%
  dplyr::select(Feature, mz, RT_min, Log2FC,
                Fold_Change, p_BH, VIP_Score) %>%
  left_join(eic_sup_df, by = "Feature") %>%
  arrange(desc(VIP_Score))

cat("Master table dimensions:", dim(master_sup), "\n")
cat("\nTop 10 features by VIP score:\n")
print(head(master_sup[, 1:7], 10))

# Build multi-sheet Excel workbook
wb_sup <- createWorkbook()

addWorksheet(wb_sup, "True_Secreted_Response_EIC")
writeData(wb_sup, "True_Secreted_Response_EIC", master_sup)

addWorksheet(wb_sup, "All_ANOVA_Significant")
writeData(wb_sup, "All_ANOVA_Significant",
          anova_sup[anova_sup$Significant, ])

addWorksheet(wb_sup, "BNC5RSw_vs_BNC5_Primary")
writeData(wb_sup, "BNC5RSw_vs_BNC5_Primary", comp_sup_treatment)

addWorksheet(wb_sup, "RSw_Background_Removed")
writeData(wb_sup, "RSw_Background_Removed",
          elevated_sup %>% filter(Feature %in% rsw_bg_sup))

addWorksheet(wb_sup, "Reduced_in_BNC5RSw")
writeData(wb_sup, "Reduced_in_BNC5RSw", reduced_sup)

addWorksheet(wb_sup, "Sample_Metadata")
writeData(wb_sup, "Sample_Metadata", sample_info_sup)

saveWorkbook(wb_sup,
             "BNC5_Supernatant_Positive_Complete_Results.xlsx",
             overwrite = TRUE)
cat("Excel workbook saved with 6 sheets\n")



# Step O — Save Complete Analysis Session

save(raw_sup, data_sup_na, data_sup_filtered,
     data_sup_imputed, data_sup_tic, data_sup_log2,
     data_sup_scaled, X_sup_log2, X_sup_scaled,
     X_sup_clean, sample_info_sup, groups_sup,
     lots_sup, groups_sup_fac,
     pca_sup, var_exp_sup,
     anova_sup, sig_feats_sup,
     comp_sup_treatment, comp_sup_background,
     comp_sup_sanity, elevated_sup,
     rsw_bg_sup, true_response_sup, true_feats_sup,
     reduced_sup, plsda_sup, vip_df_sup,
     master_sup, volcano_sup,
     file = "BNC5_supernatant_pos_analysis_session.RData")

cat("Complete analysis session saved\n")
cat("Reload: load('BNC5_supernatant_pos_analysis_session.RData')\n")
