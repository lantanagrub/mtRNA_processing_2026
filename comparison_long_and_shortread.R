# =============================================================================
# Cross-Platform DGE Comparison — Nanopore dRNA-seq (oarfish) vs Illumina (Salmon)
# C. elegans aging timecourse, paired samples, each age vs Day 1
# Reads the per-comparison CSVs written by both DGE scripts.
#
# Manuscript figures produced by this script:
#   platform_DEG_counts.pdf     Supp. Fig. 1A
#   platform_LFC_scatter.pdf    Supp. Fig. 1B
# =============================================================================

library(tidyverse)

# =============================================================================
# 0. PATHS AND PARAMETERS — edit this block
# =============================================================================

nanopore_dir <- "~/aging_rnaseq/nanopore_dge"   # out_dir of the oarfish script
illumina_dir <- "~/aging_rnaseq/illumina_dge"   # out_dir of the Salmon script
out_dir      <- "~/aging_rnaseq/platform_comparison"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ages            <- c(1, 2, 3, 4, 5, 7, 10, 15)
comparison_ages <- ages[ages != 1]
comparison_lvls <- paste0("Day", comparison_ages)

# DEG definition (same padj and LFC cutoffs used for DE calls in both DGE scripts)
padj_thresh <- 0.05
lfc_thresh  <- 1

base_size <- 7.5

# Opaque, pre-blended colors (no alpha)
platform_colors <- c("Nanopore" = "#1b9e77", "Illumina" = "#d95f02")

# -----------------------------------------------------------------------------
# LABELS — every user-visible string in this script is defined here.
# -----------------------------------------------------------------------------

# Overlap classes written to platform_paired_results.csv
overlap_levels <- c(
  "Shared",
  "Nanopore only",
  "Nanopore only, not tested in Illumina",
  "Illumina only",
  "Illumina only, not tested in Nanopore"
)

# DE classes for the LFC scatter
scatter_levels <- c("Neither", "Nanopore only", "Illumina only", "Shared")
scatter_colors <- setNames(
  c("grey80", "#1b9e77", "#d95f02", "#7b3294"),
  scatter_levels
)

any_age_label <- "Any age"

lab <- list(
  counts_title   = paste0("DE genes vs Day 1 by platform  (|LFC| > ",
                          lfc_thresh, ", FDR < ", padj_thresh, ")"),
  counts_y       = "Number of genes (up above zero, down below)",
  lfc_title      = "log2 fold change vs Day 1: Nanopore vs Illumina",
  lfc_x          = "Illumina log2 fold change",
  lfc_y          = "Nanopore log2 fold change",
  lfc_color      = "DE on"
)

theme_cmp <- theme_bw(base_size = base_size) +
  theme(
    strip.background = element_rect(fill = "grey92"),
    legend.position  = "bottom"
  )

# =============================================================================
# 1. LOAD PER-COMPARISON RESULTS FROM BOTH PLATFORMS
# =============================================================================

read_platform <- function(dir, platform) {
  bind_rows(lapply(comparison_ages, function(a) {
    f <- file.path(dir, paste0("DGE_Day", a, "_vs_Day1.csv"))
    if (!file.exists(f)) stop("File not found: ", f)
    read_csv(f, show_col_types = FALSE) %>%
      mutate(comparison = paste0("Day", a))
  })) %>%
    mutate(
      platform   = platform,
      comparison = factor(comparison, levels = comparison_lvls),
      # "tested" = has an adjusted p-value (passed low-count filter,
      # independent filtering, and outlier checks)
      tested     = !is.na(padj),
      direction  = case_when(
        tested & padj < padj_thresh & log2FoldChange >  lfc_thresh ~ "Up",
        tested & padj < padj_thresh & log2FoldChange < -lfc_thresh ~ "Down",
        tested ~ "NS",
        TRUE ~ NA_character_
      ),
      deg = direction %in% c("Up", "Down")
    )
}

res_np <- read_platform(nanopore_dir, "Nanopore")
res_il <- read_platform(illumina_dir, "Illumina")
res_both <- bind_rows(res_np, res_il) %>%
  mutate(platform = factor(platform, levels = names(platform_colors)))

# =============================================================================
# 2. PAIRED GENE x COMPARISON TABLE
# =============================================================================

pair_cols <- c("gene_id", "display_name", "comparison", "baseMean",
               "log2FoldChange", "lfcSE", "padj", "tested", "direction", "deg")

paired <- full_join(
  res_np %>% select(all_of(pair_cols)),
  res_il %>% select(all_of(pair_cols)),
  by     = c("gene_id", "comparison"),
  suffix = c("_np", "_il")
) %>%
  mutate(
    display_name = coalesce(display_name_np, display_name_il, gene_id),
    across(c(tested_np, tested_il, deg_np, deg_il), ~ coalesce(.x, FALSE))
  ) %>%
  select(-display_name_np, -display_name_il) %>%
  relocate(gene_id, display_name, comparison)

classify_overlap <- function(df) {
  df %>%
    mutate(overlap = case_when(
      deg_np & deg_il    ~ overlap_levels[1],
      deg_np & tested_il ~ overlap_levels[2],
      deg_np             ~ overlap_levels[3],
      deg_il & tested_np ~ overlap_levels[4],
      deg_il             ~ overlap_levels[5],
      TRUE               ~ NA_character_
    ))
}

paired <- paired %>%
  classify_overlap() %>%
  mutate(concordant_direction = if_else(deg_np & deg_il,
                                        direction_np == direction_il,
                                        NA))

write_csv(paired, file.path(out_dir, "platform_paired_results.csv"))

# =============================================================================
# 3. SUMMARY TABLES
# =============================================================================

# --- DEG counts per platform (data for Supp. 1A) ---
deg_counts <- res_both %>%
  group_by(platform, comparison) %>%
  summarise(
    tested = sum(tested),
    up     = sum(direction == "Up",   na.rm = TRUE),
    down   = sum(direction == "Down", na.rm = TRUE),
    total  = up + down,
    .groups = "drop"
  )

# Genes DE at any age on each platform (CSV only, not plotted)
deg_counts_any <- res_both %>%
  group_by(platform, gene_id) %>%
  summarise(deg = any(deg), .groups = "drop") %>%
  group_by(platform) %>%
  summarise(total = sum(deg), .groups = "drop") %>%
  mutate(comparison = any_age_label)

write_csv(bind_rows(deg_counts %>% mutate(comparison = as.character(comparison)),
                    deg_counts_any),
          file.path(out_dir, "platform_DEG_counts.csv"))

# --- Concordance metrics per comparison ---
# lfc_spearman_tested_both supplies the rho labels on the LFC scatter (Supp. 1B)
safe_spearman <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(NA_real_)
  cor(x[ok], y[ok], method = "spearman")
}

metrics <- paired %>%
  group_by(comparison) %>%
  summarise(
    n_tested_both        = sum(tested_np & tested_il),
    n_deg_np             = sum(deg_np),
    n_deg_il             = sum(deg_il),
    n_shared             = sum(deg_np & deg_il),
    n_shared_same_dir    = sum(concordant_direction, na.rm = TRUE),
    n_shared_opposite    = sum(!concordant_direction, na.rm = TRUE),
    jaccard_all          = n_shared / sum(deg_np | deg_il),
    jaccard_tested_both  = n_shared / sum((deg_np | deg_il) & tested_np & tested_il),
    lfc_spearman_tested_both = safe_spearman(
      log2FoldChange_np[tested_np & tested_il],
      log2FoldChange_il[tested_np & tested_il]),
    lfc_spearman_union_deg   = safe_spearman(
      log2FoldChange_np[(deg_np | deg_il) & tested_np & tested_il],
      log2FoldChange_il[(deg_np | deg_il) & tested_np & tested_il]),
    .groups = "drop"
  )

print(metrics)
write_csv(metrics, file.path(out_dir, "platform_concordance_metrics.csv"))

# =============================================================================
# 4. VISUALIZATIONS
# =============================================================================

# ---- 4a. DEG counts per platform, up above zero, down below (Supp. 1A) -----

counts_long <- deg_counts %>%
  transmute(platform, comparison, Up = up, Down = -down) %>%
  pivot_longer(c(Up, Down), names_to = "direction", values_to = "count")

dodge <- position_dodge(width = 0.8)

p_counts <- ggplot(counts_long,
                   aes(x = comparison, y = count, fill = platform, group = platform)) +
  geom_col(position = dodge, width = 0.75) +
  geom_text(
    aes(label = abs(count), vjust = if_else(count >= 0, -0.3, 1.3)),
    position = dodge, size = base_size / .pt
  ) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +
  scale_fill_manual(values = platform_colors) +
  scale_y_continuous(expand = expansion(mult = 0.12)) +
  labs(title = lab$counts_title, x = NULL, y = lab$counts_y, fill = NULL) +
  theme_cmp

ggsave(file.path(out_dir, "platform_DEG_counts.pdf"), p_counts, width = 6, height = 3.5)

# ---- 4b. LFC scatter, genes tested on both platforms (Supp. 1B) ------------

lfc_df <- paired %>%
  filter(tested_np, tested_il) %>%
  mutate(
    de_class = case_when(
      deg_np & deg_il ~ "Shared",
      deg_np          ~ "Nanopore only",
      deg_il          ~ "Illumina only",
      TRUE            ~ "Neither"
    ),
    de_class = factor(de_class, levels = scatter_levels)
  ) %>%
  arrange(de_class)   # draw Neither first, Shared last

rho_labels <- metrics %>%
  transmute(comparison,
            label = paste0("rho == ", sprintf("%.2f", lfc_spearman_tested_both)))

p_lfc <- ggplot(lfc_df, aes(x = log2FoldChange_il, y = log2FoldChange_np)) +
  geom_hline(yintercept = c(-lfc_thresh, lfc_thresh),
             linetype = "dashed", color = "grey60", linewidth = 0.25) +
  geom_vline(xintercept = c(-lfc_thresh, lfc_thresh),
             linetype = "dashed", color = "grey60", linewidth = 0.25) +
  geom_abline(slope = 1, intercept = 0, color = "grey30", linewidth = 0.3) +
  geom_point(aes(color = de_class), size = 0.5) +
  geom_text(data = rho_labels, aes(x = -Inf, y = Inf, label = label),
            parse = TRUE, hjust = -0.15, vjust = 1.4,
            size = base_size / .pt, inherit.aes = FALSE) +
  scale_color_manual(values = scatter_colors, breaks = rev(scatter_levels)) +
  facet_wrap(~ comparison, ncol = 4) +
  coord_equal() +
  labs(title = lab$lfc_title, x = lab$lfc_x, y = lab$lfc_y, color = lab$lfc_color) +
  theme_cmp +
  guides(color = guide_legend(override.aes = list(size = 2)))

ggsave(file.path(out_dir, "platform_LFC_scatter.pdf"), p_lfc, width = 6, height = 4)

cat("Platform comparison outputs written to", out_dir, "\n")
