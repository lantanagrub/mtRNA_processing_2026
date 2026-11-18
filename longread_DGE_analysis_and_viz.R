# =============================================================================
# Oarfish dRNA-seq DGE Pipeline — Gene-Level, Nanopore
# C. elegans aging timecourse: Day 1, 2, 3, 4, 5, 7, 10, 15
# Comparisons: each age vs Day 1
#
# Manuscript figures produced by this script:
#   PCA_gene_level.pdf          Supp. Fig. 2A
#   volcano_gene_level.pdf      Supp. Fig. 2C
#   GOI_expression.pdf          Fig. 2B
#   trajectory_AllGroups.pdf    Fig. 2D
#   heatmap_AllComplexes.pdf    Fig. 2E
#
# The per-comparison DGE CSVs (section 7) are also read by platform_comparison.R.
# =============================================================================

library(DESeq2)
library(tidyverse)
library(rtracklayer)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)

# =============================================================================
# 0. PATHS AND PARAMETERS — edit this block
# =============================================================================

setwd("~/aging_rnaseq")

data_dir        <- "~/aging_rnaseq/oarfish/"
gtf_path        <- "~/aging_rnaseq/reference/c_elegans.PRJNA13758.WS279.canonical_geneset.gtf.gz"
gene_names_path <- "~/aging_rnaseq/reference/c_elegans.PRJNA13758.WS279.geneIDs.txt.gz"
out_dir         <- "~/aging_rnaseq/nanopore_dge"                  # output directory for all files
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Set to TRUE to restrict DGE to protein-coding genes only
coding_only <- TRUE

ages <- c(1, 2, 3, 4, 5, 7, 10, 15)
reps <- c(1, 2, 3)

# Low-count filter: gene must have >= this many counts in >= this many samples
min_counts  <- 10
min_samples <- 3

# LFC threshold for calling significance in summaries and plots
lfc_thresh <- 1

# ---- Font sizes (pt) --------------------------------------------------------
# base_font sets the ggplot theme. theme_bw draws axis text, facet strip text
# and legend text at 0.8 x base_font, so base_font = 8 gives 6.4 pt for those.
base_font  <- 8
font_small <- 6   # volcano gene labels

# ---- Figure sizes (inches), keyed by output file name -----------------------
# Any file not listed here falls back to fig_size_default.
# heatmap_AllComplexes.pdf sets its own size in its section below.
# The expression-plot size depends on how many genes are present, so adjust
# it after a first run.
fig_size_default <- c(width = 6, height = 4)

fig_sizes <- list(
  "volcano_gene_level.pdf"        = c(width = 6,    height = 3),
  "PCA_gene_level.pdf"            = c(width = 2.35, height = 2.6),
  "GOI_expression.pdf"            = c(width = 4.5,  height = 2.5),
  "trajectory_AllGroups.pdf"      = c(width = 3,    height = 2.5)
)

fig_size <- function(fname) {
  if (fname %in% names(fig_sizes)) fig_sizes[[fname]] else fig_size_default
}

save_gg <- function(p, fname) {
  sz <- fig_size(fname)
  ggsave(file.path(out_dir, fname), p,
         width = sz[["width"]], height = sz[["height"]], units = "in")
}

# =============================================================================
# 1. BUILD TX2GENE FROM GTF
# =============================================================================

cat("Loading GTF...\n")
gtf <- import(gtf_path)

tx2gene_base <- as.data.frame(gtf) %>%
  filter(type == "transcript") %>%
  transmute(
    tname   = transcript_id,
    gene_id = gene_id
  ) %>%
  filter(!is.na(tname), !is.na(gene_id)) %>%
  distinct()

cat("Transcripts in GTF:", nrow(tx2gene_base), "\n")
cat("Unique genes in GTF:", n_distinct(tx2gene_base$gene_id), "\n")

# =============================================================================
# 2. LOAD WORMBASE GENE NAME / BIOTYPE TABLE
# =============================================================================
# Columns: species_id, gene_id, locus_name, sequence_name, status, biotype
# biotype options include: protein_coding_gene, gene, pseudogene

cat("Loading gene name table...\n")
gene_names <- read_delim(
  gene_names_path,
  delim      = ",",
  col_names  = c("species_id", "gene_id", "locus_name", "sequence_name", "status", "biotype"),
  show_col_types = FALSE
) %>%
  select(gene_id, locus_name, sequence_name, biotype) %>%
  mutate(
    # Prefer locus name (e.g. nduo-6), fall back to sequence name (e.g. 2RSSE.1),
    # then WBGene ID as last resort
    display_name = case_when(
      !is.na(locus_name)    & locus_name    != "" ~ locus_name,
      !is.na(sequence_name) & sequence_name != "" ~ sequence_name,
      TRUE ~ gene_id
    )
  ) %>%
  distinct()

cat("Biotype breakdown:\n")
print(table(gene_names$biotype, useNA = "ifany"))

# =============================================================================
# 3. JOIN TX2GENE WITH GENE NAMES
# =============================================================================

tx2gene <- tx2gene_base %>%
  left_join(gene_names, by = "gene_id") %>%
  mutate(display_name = coalesce(display_name, gene_id))

# Lookup table used for adding display names to results later
gene_id_to_name <- tx2gene %>%
  select(gene_id, display_name, biotype) %>%
  distinct()

# Optionally restrict to protein-coding genes
if (coding_only) {
  cat("Restricting to protein_coding_gene biotype.\n")
  tx2gene_use <- tx2gene %>% filter(biotype == "protein_coding_gene")
} else {
  cat("Retaining all biotypes.\n")
  tx2gene_use <- tx2gene
}

cat("Transcripts used for aggregation:", nrow(tx2gene_use), "\n")
cat("Genes used for aggregation:", n_distinct(tx2gene_use$gene_id), "\n")

# =============================================================================
# 4. LOAD AND MERGE OARFISH COUNTS
# =============================================================================

cat("Loading oarfish count files...\n")

coldata <- expand.grid(age = ages, rep = reps) %>%
  mutate(
    filename = paste0("day_", age, "_rep_", rep, "_nanopore_counts.quant"),
    sname    = paste0("Day", age, "_rep", rep),
    age_fct  = factor(paste0("Day", age), levels = paste0("Day", ages))
  ) %>%
  arrange(age, rep)
rownames(coldata) <- coldata$sname

count_list <- lapply(seq_len(nrow(coldata)), function(i) {
  f <- file.path(data_dir, coldata$filename[i])
  s <- coldata$sname[i]
  if (!file.exists(f)) stop("File not found: ", f)
  read_tsv(f, show_col_types = FALSE) %>%
    select(tname, num_reads) %>%
    rename(!!s := num_reads)
})

counts_tx <- reduce(count_list, full_join, by = "tname") %>%
  column_to_rownames("tname") %>%
  replace(is.na(.), 0)

cat("Transcripts loaded:", nrow(counts_tx), "\n")

# Coverage check
quant_tx <- rownames(counts_tx)
mapped   <- quant_tx[quant_tx %in% tx2gene_use$tname]
unmapped <- quant_tx[!quant_tx %in% tx2gene_use$tname]

cat("Quantified transcripts:", length(quant_tx), "\n")
cat("Mapped to gene:", length(mapped), "\n")
cat("Unmapped (dropped):", length(unmapped), "\n")
if (length(unmapped) > 0 && length(unmapped) <= 20) {
  cat("Unmapped IDs:\n"); print(unmapped)
} else if (length(unmapped) > 20) {
  cat("First 20 unmapped IDs:\n"); print(head(unmapped, 20))
}

# =============================================================================
# 5. AGGREGATE TRANSCRIPT COUNTS TO GENE LEVEL
# =============================================================================

counts_gene <- counts_tx %>%
  rownames_to_column("tname") %>%
  inner_join(tx2gene_use %>% select(tname, gene_id), by = "tname") %>%
  select(-tname) %>%
  group_by(gene_id) %>%
  summarise(across(everything(), sum), .groups = "drop") %>%
  column_to_rownames("gene_id")

# Round fractional EM counts to integers for DESeq2
counts_gene <- round(counts_gene)

cat("Genes after aggregation:", nrow(counts_gene), "\n")

stopifnot(all(colnames(counts_gene) == rownames(coldata)))

# =============================================================================
# 6. DESeq2
# =============================================================================

cat("Running DESeq2...\n")

dds <- DESeqDataSetFromMatrix(
  countData = counts_gene,
  colData   = coldata,
  design    = ~ age_fct
)

dds$age_fct <- relevel(dds$age_fct, ref = "Day1")

keep <- rowSums(counts(dds) >= min_counts) >= min_samples
dds  <- dds[keep, ]
cat("Genes after low-count filter:", sum(keep), "\n")

dds <- DESeq(dds)

# =============================================================================
# 7. EXTRACT PAIRWISE RESULTS AND SAVE
# =============================================================================

comparison_ages <- ages[ages != 1]

results_list <- lapply(comparison_ages, function(a) {
  label <- paste0("Day", a)
  res <- results(
    dds,
    contrast      = c("age_fct", label, "Day1"),
    alpha         = 0.05,
    pAdjustMethod = "BH"
  )
  as.data.frame(res) %>%
    rownames_to_column("gene_id") %>%
    mutate(comparison = label)
})
names(results_list) <- paste0("Day", comparison_ages, "_vs_Day1")

# Join display names onto all results
results_list <- lapply(results_list, function(df) {
  df %>%
    left_join(gene_id_to_name %>% select(gene_id, display_name), by = "gene_id") %>%
    mutate(display_name = coalesce(display_name, gene_id)) %>%
    select(gene_id, display_name, everything())
})

# Combined long-format results (used internally for plots)
results_all <- bind_rows(results_list) %>%
  mutate(comparison = factor(comparison, levels = paste0("Day", comparison_ages)))

# --- Per-comparison CSVs (wide within each file: gene_id, display_name + stats) ---
# Also read by platform_comparison.R
for (nm in names(results_list)) {
  results_list[[nm]] %>%
    select(gene_id, display_name,
           baseMean, log2FoldChange, lfcSE, stat, pvalue, padj) %>%
    write_csv(file.path(out_dir, paste0("DGE_", nm, ".csv")))
}

# --- All-comparisons CSV (one row per gene, comparisons as column groups) ---
results_wide <- results_all %>%
  select(gene_id, display_name, comparison,
         baseMean, log2FoldChange, pvalue, padj) %>%
  group_by(gene_id) %>%
  mutate(baseMean = first(baseMean)) %>%
  ungroup() %>%
  pivot_wider(
    names_from  = comparison,
    values_from = c(log2FoldChange, pvalue, padj)
  ) %>%
  # Build column order explicitly: for each day, all three stats together
  select(
    gene_id, display_name, baseMean,
    unlist(lapply(paste0("Day", comparison_ages), function(d) {
      c(paste0("log2FoldChange_", d),
        paste0("pvalue_", d),
        paste0("padj_", d))
    }))
  )

write_csv(results_wide,
          file.path(out_dir, "DGE_gene_level_all_comparisons.csv"))

# =============================================================================
# 8. SUMMARY TABLE
# =============================================================================
# sig_total counts padj < 0.05 with no LFC cutoff; sig_up / sig_down also
# require |LFC| > lfc_thresh (the same calls as the volcano plot).

sig_summary <- results_all %>%
  filter(!is.na(padj)) %>%
  group_by(comparison) %>%
  summarise(
    total_tested = n(),
    sig_total    = sum(padj < 0.05, na.rm = TRUE),
    sig_up       = sum(padj < 0.05 & log2FoldChange >  lfc_thresh, na.rm = TRUE),
    sig_down     = sum(padj < 0.05 & log2FoldChange < -lfc_thresh, na.rm = TRUE),
    .groups = "drop"
  )

print(sig_summary)
write_csv(sig_summary, file.path(out_dir, "DGE_gene_summary.csv"))

# =============================================================================
# 9. GENOME-WIDE PLOTS — VOLCANO (Supp. 2C) AND PCA (Supp. 2A)
# =============================================================================

age_palette <- colorRampPalette(c("#f5c842", "#800026"))(length(ages))
names(age_palette) <- paste0("Day", ages)

de_colors <- c("Up" = "#d73027", "Down" = "#4575b4", "NS" = "grey70")

# ---- 9a. Volcano plots (faceted) ------------------------------------------

# Annotate direction for plot data
volcano_data <- results_all %>%
  filter(!is.na(padj), !is.na(log2FoldChange)) %>%
  mutate(
    direction = case_when(
      padj < 0.05 & log2FoldChange >  lfc_thresh ~ "Up",
      padj < 0.05 & log2FoldChange < -lfc_thresh ~ "Down",
      TRUE ~ "NS"
    )
  )

# Top 5 up + top 5 down per comparison (lowest padj) for volcano labels
top_labels <- volcano_data %>%
  filter(direction != "NS") %>%
  group_by(comparison, direction) %>%
  slice_min(padj, n = 5) %>%
  ungroup()

p_volcano <- ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = direction), alpha = 0.5, size = 0.7) +
  geom_text_repel(
    data          = top_labels,
    aes(label     = display_name),
    size          = font_small / .pt,
    max.overlaps  = 12,
    segment.color = "grey50"
  ) +
  scale_color_manual(values = de_colors) +
  geom_vline(xintercept = c(-lfc_thresh, lfc_thresh),
             linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed", color = "grey40", linewidth = 0.4) +
  facet_wrap(~ comparison, ncol = 6) +
  labs(
    title = "Differential Gene Expression vs Day 1",
    x     = expression(log[2]~"Fold Change"),
    y     = expression(-log[10]~"(adj. p-value)"),
    color = NULL
  ) +
  theme_bw(base_size = base_font) +
  theme(
    strip.background = element_rect(fill = "grey92"),
    legend.position  = "bottom"
  )

save_gg(p_volcano, "volcano_gene_level.pdf")

# ---- 9b. PCA ---------------------------------------------------------------

# VST (not blind to design); also used for all plots in section 10
vsd <- vst(dds, blind = FALSE)

pca_data <- plotPCA(vsd, intgroup = "age_fct", returnData = TRUE)
pct_var  <- round(100 * attr(pca_data, "percentVar"))

# Sample-name labels omitted; age is shown by color and the legend
p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = age_fct)) +
  geom_point(size = 3.5) +
  scale_color_manual(values = age_palette, name = "Age") +
  labs(
    title = "PCA — VST Gene Expression",
    x     = paste0("PC1: ", pct_var[1], "% variance"),
    y     = paste0("PC2: ", pct_var[2], "% variance")
  ) +
  theme_bw(base_size = base_font)

save_gg(p_pca, "PCA_gene_level.pdf")


# =============================================================================
# 10. GENE-SET PLOTS — mtDNA-ENCODED GENES, ETC COMPLEXES, HOUSEKEEPING
# =============================================================================
# BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
library(circlize)

# VST matrix used throughout this section
vsd_mat <- assay(vsd)

# Heatmap color scale (z-score, -3 to 3)
shared_colors <- colorRamp2(seq(-3, 3, length.out = 11),
                            colorRampPalette(rev(brewer.pal(11, "RdBu")))(11))

# Column order: samples arranged by age then rep
col_order <- coldata %>% arrange(age, rep) %>% pull(sname)

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

# Expression plot — ids_present sets which genes are plotted; display_order
# (display names) sets facet order.
#   relative_to_day1 = TRUE plots each gene's VST minus its own Day 1 mean.
#   shared_y = TRUE uses one y-axis scale for all panels, so axis labels are
#   drawn only on the outer panels of the grid.
#   point_size sets both replicate points and mean points.
#   axis_title_pt, when set, sets axis titles and facet strip text.
#   other_pt, when set, sets the plot title and tick labels (and strip text,
#   if axis_title_pt is not set).
#   An asterisk marks ages where the gene has padj < 0.05 vs Day 1
#   (results_subset$sig).
make_expression_plot <- function(ids_present, display_order,
                                 results_subset, title, filename,
                                 relative_to_day1 = FALSE,
                                 shared_y         = relative_to_day1,
                                 ncol             = 6,
                                 x_breaks         = ages,
                                 y_label          = expression(Delta*"VST Expression (~"*log[2]*" scale)"),
                                 point_size       = 3,
                                 axis_title_pt    = NULL,
                                 other_pt         = NULL) {
  expr_df <- vsd_mat[ids_present, ] %>%
    as.data.frame() %>%
    rownames_to_column("gene_id") %>%
    pivot_longer(-gene_id, names_to = "sname", values_to = "vst") %>%
    left_join(coldata %>% select(sname, age, age_fct), by = "sname") %>%
    left_join(gene_id_to_name %>% select(gene_id, display_name), by = "gene_id") %>%
    mutate(
      display_name = coalesce(display_name, gene_id),
      # Apply specified display order as factor for facet ordering
      display_name = factor(display_name, levels = display_order)
    )

  # Express each gene relative to its own Day 1 mean (delta VST)
  if (relative_to_day1) {
    expr_df <- expr_df %>%
      group_by(gene_id) %>%
      mutate(vst = vst - mean(vst[age == 1], na.rm = TRUE)) %>%
      ungroup()
  }

  sig_df <- results_subset %>%
    filter(sig) %>%
    mutate(display_name = factor(coalesce(display_name, gene_id),
                                 levels = display_order))

  p <- ggplot(expr_df, aes(x = age, y = vst)) +
    {
      if (relative_to_day1)
        geom_hline(yintercept = 0, linetype = "dashed",
                   color = "grey60", linewidth = 0.3)
    } +
    stat_summary(aes(group = display_name), fun = mean,
                 geom = "line", color = "grey40", linewidth = 0.6) +
    geom_jitter(aes(color = age_fct), width = 0.15, size = point_size, alpha = 0.8) +
    stat_summary(fun = mean, geom = "point", size = point_size,
                 shape = 21, fill = "white", color = "black") +
    geom_text(
      data        = sig_df,
      aes(x = age_num, y = Inf, label = "*"),
      vjust       = 1.5, size = 5, color = "black",
      inherit.aes = FALSE
    ) +
    scale_color_manual(values = age_palette, name = "Age") +
    scale_x_continuous(breaks = x_breaks) +
    facet_wrap(~ display_name, ncol = ncol,
               scales = if (shared_y) "fixed" else "free_y") +
    labs(title = title, x = "Age (days)", y = y_label) +
    theme_bw(base_size = base_font) +
    theme(
      strip.background = element_rect(fill = "grey92"),
      legend.position  = "none",
      panel.spacing    = unit(1.5, "mm")
    )

  # Optional explicit font sizes (pt)
  if (!is.null(other_pt)) {
    p <- p + theme(
      plot.title  = element_text(size = other_pt),
      axis.text   = element_text(size = other_pt),
      strip.text  = element_text(size = axis_title_pt %||% other_pt)
    )
  }
  if (!is.null(axis_title_pt)) {
    p <- p + theme(axis.title = element_text(size = axis_title_pt))
  }

  save_gg(p, filename)
  cat("Expression plot written:", filename, "\n")
}

# Wide-format DE summary CSV — one row per gene, days as columns
make_sig_summary <- function(results_subset, filename) {
  results_subset %>%
    mutate(display_name = coalesce(display_name, gene_id)) %>%
    select(gene_id, display_name, comparison, log2FoldChange, padj, sig) %>%
    pivot_wider(
      names_from  = comparison,
      values_from = c(log2FoldChange, padj, sig)
    ) %>%
    select(
      gene_id, display_name,
      unlist(lapply(paste0("Day", comparison_ages), function(d) {
        c(paste0("log2FoldChange_", d),
          paste0("padj_", d),
          paste0("sig_", d))
      }))
    ) %>%
    arrange(display_name) %>%
    write_csv(file.path(out_dir, filename))
  cat("Summary CSV written:", filename, "\n")
}

# Resolve gene name/ID set to gene IDs, warn on misses
resolve_ids <- function(gene_set) {
  ids <- gene_id_to_name %>%
    filter(display_name %in% gene_set | gene_id %in% gene_set) %>%
    pull(gene_id) %>%
    unique()
  not_found <- gene_set[
    !gene_set %in% c(gene_id_to_name$display_name, gene_id_to_name$gene_id)
  ]
  if (length(not_found) > 0) cat("Not found in annotation:", not_found, "\n")
  ids
}

# Scale and clip VST matrix rows, replace rownames with display names
prepare_scaled_mat <- function(ids, col_order) {
  mat    <- vsd_mat[ids, col_order, drop = FALSE]
  scaled <- t(scale(t(mat)))
  scaled <- pmin(pmax(scaled, -3), 3)
  rownames(scaled) <- coalesce(
    gene_id_to_name$display_name[match(rownames(scaled), gene_id_to_name$gene_id)],
    rownames(scaled)
  )
  scaled
}

# =============================================================================
# DEFINE GENE SETS
# =============================================================================

# mtDNA-encoded genes
genes_of_interest <- c("ndfl-4",  "nduo-6", "nduo-2", "nduo-3", "nduo-4",
                       "nduo-5",  "atp-6",  "nduo-1", "ctb-1",  "ctc-2",
                       "ctc-3",   "ctc-1")

# Housekeeping genes (used only in the trajectory plot)
housekeeping_genes <- c("rps-23",  "rps-27", "rps-16", "rps-26", "rps-4",
                        "rps-2",   "rps-17", "rpl-24.1", "rpl-15", "rpl-35",
                        "rpl-36",  "rpl-27", "rpl-33")

# Display order for mtDNA-encoded genes (GOI expression plot facets)
goi_order <- c("nduo-6", "ndfl-4", "nduo-1", "atp-6",  "nduo-2",
               "ctb-1",  "ctc-3",  "nduo-4", "ctc-1",  "ctc-2",
               "nduo-3", "nduo-5")

# Resolve to gene IDs
goi_ids <- resolve_ids(genes_of_interest)
hk_ids  <- resolve_ids(housekeeping_genes)

goi_ids_present <- goi_ids[goi_ids %in% rownames(vsd_mat)]
hk_ids_present  <- hk_ids[hk_ids  %in% rownames(vsd_mat)]

# GOI IDs in goi_order
goi_order_ids <- gene_id_to_name %>%
  filter(display_name %in% goi_order) %>%
  mutate(display_name = factor(display_name, levels = goi_order)) %>%
  arrange(display_name) %>%
  pull(gene_id) %>%
  .[. %in% rownames(vsd_mat)]

# Display names in goi_order for expression plot facet ordering
goi_display_order <- coalesce(
  gene_id_to_name$display_name[match(goi_order_ids, gene_id_to_name$gene_id)],
  goi_order_ids
)

# GOI results subset. sig = padj < 0.05 with no LFC cutoff (unlike the
# volcano, which also requires |LFC| > lfc_thresh).
goi_results <- results_all %>%
  filter(gene_id %in% goi_ids_present) %>%
  mutate(
    sig     = !is.na(padj) & padj < 0.05,
    age_num = as.integer(str_remove(as.character(comparison), "Day"))
  )

# =============================================================================
# mtDNA-ENCODED GENES — EXPRESSION PLOT (Fig. 2B) AND CSV
# =============================================================================

cat("\n--- Genes of Interest ---\n")

make_sig_summary(goi_results, "GOI_DE_summary.csv")

# Delta VST relative to each gene's Day 1 mean, shared y-axis, 2 x 6 grid
make_expression_plot(
  ids_present      = goi_order_ids,
  display_order    = goi_display_order,
  results_subset   = goi_results,
  title            = NULL,
  filename         = "GOI_expression.pdf",
  relative_to_day1 = TRUE,
  ncol             = 6,
  x_breaks         = c(1, 5, 10, 15),
  y_label          = expression(Delta*"VST Expression (~"*log[2]*" scale)"),
  point_size       = 1,
  axis_title_pt    = 8,
  other_pt         = 6
)


# =============================================================================
# MITOCHONDRIAL COMPLEX GENE SETS (nuclear-encoded)
# =============================================================================

complex_gene_sets <- list(
  "Complex_I" = c(
    "WBGene00003831", "WBGene00009992", "WBGene00021562", "WBGene00001520",
    "WBGene00006463", "WBGene00020417", "WBGene00003061", "WBGene00021839",
    "WBGene00009051", "WBGene00012376", "WBGene00020636", "WBGene00009294",
    "WBGene00013406", "WBGene00007880", "WBGene00013308", "WBGene00009740",
    "WBGene00021849", "WBGene00021800", "WBGene00019401", "WBGene00022380",
    "WBGene00016393", "WBGene00013237", "WBGene00018151", "WBGene00009712",
    "WBGene00007684", "WBGene00012166", "WBGene00016118", "WBGene00014086",
    "WBGene00008414", "WBGene00013094", "WBGene00015810", "WBGene00010326",
    "WBGene00018361", "WBGene00022169", "WBGene00015501", "WBGene00022075",
    "WBGene00007192"
  ),
  "Complex_II" = c(
    "sdha-1", "sdha-2", "sdhb-1", "mev-1", "sdhd-1"
  ),
  "Complex_III" = c(
    "ucr-1", "ucr-2.1", "ucr-2.2", "ucr-2.3", "isp-1", "cyc-1",
    "T02H6.11", "F45H10.2", "R07E4.3", "T27E9.2", "ucr-11"
  ),
  "Complex_IV" = c(
    "cox-4", "cox-5A", "cox-5B", "cox-6A", "cox-6B", "cox-6C", "cox-7C"
  ),
  "Complex_V" = c(
    "atp-1", "atp-2", "Y69A2AR.18", "F58F12.1", "hpo-18", "R05D3.6",
    "ZC262.5", "Y82E9BR.3", "asb-1", "asb-2", "atp-5", "atp-4", "atp-3",
    "R04F11.2", "R53.4", "asg-1", "asg-2", "mai-1", "mai-2", "T26E3.7"
  )
)

# =============================================================================
# RESOLVE IDS AND BUILD RESULTS FOR EACH COMPLEX
# =============================================================================

complex_data <- lapply(names(complex_gene_sets), function(cx) {
  gene_set <- complex_gene_sets[[cx]]

  ids <- resolve_ids(gene_set)
  ids_present <- ids[ids %in% rownames(vsd_mat)]

  if (length(ids_present) == 0) {
    cat("No genes found in VST matrix for", cx, "— skipping.\n")
    return(NULL)
  }

  results_sub <- results_all %>%
    filter(gene_id %in% ids_present) %>%
    mutate(
      sig     = !is.na(padj) & padj < 0.05,
      age_num = as.integer(str_remove(as.character(comparison), "Day"))
    )

  list(
    ids_present = ids_present,
    results_sub = results_sub
  )
})
names(complex_data) <- names(complex_gene_sets)

# Remove any NULLs
complex_data <- Filter(Negate(is.null), complex_data)

# =============================================================================
# ALL-COMPLEX DE SUMMARY CSV
# =============================================================================
# Nuclear-encoded complex genes only; the mtDNA-encoded genes are in
# GOI_DE_summary.csv.

all_complex_results <- bind_rows(lapply(names(complex_data), function(cx) {
  complex_data[[cx]]$results_sub %>% mutate(complex = cx)
})) %>%
  distinct(gene_id, comparison, .keep_all = TRUE)

make_sig_summary(all_complex_results, "AllComplexes_DE_summary.csv")

# =============================================================================
# COMBINED HEATMAP — mtDNA-encoded genes + nuclear Complexes I–V (Fig. 2E)
# =============================================================================
# Everything for heatmap_AllComplexes.pdf is in this section. It uses objects
# made earlier in the script: vsd_mat, coldata, col_order, goi_order_ids,
# complex_data, prepare_scaled_mat(), shared_colors, out_dir.

cat("\n--- All Complexes Combined Heatmap ---\n")

# --- Settings ---------------------------------------------------------------
all_heat_file   <- "heatmap_AllComplexes.pdf"
all_heat_width  <- 3.6   # in
all_heat_height <- 2.5   # in
all_heat_font   <- 7     # pt, column (day) labels
all_heat_rowlab <- 7     # pt, row group labels (mtDNA / Complex)
all_heat_legend <- 7     # pt, z-score legend
all_heat_lwd    <- 0.5   # line width for borders and dendrogram
all_heat_gap    <- 0.5   # mm, gap between row groups and between days
all_heat_pad    <- 1     # mm, padding around the whole figure

# Row group labels as displayed ("\n" = line break)
split_labels <- c(
  "mtDNA encoded" = "mtDNA\nencoded",
  "Complex I"     = "Complex I",
  "Complex II"    = "Complex II",
  "Complex III"   = "Complex III",
  "Complex IV"    = "Complex IV",
  "Complex V"     = "Complex V"
)

# --- Data: z-scored VST, one matrix per group, stacked ----------------------
mtdna_scaled <- prepare_scaled_mat(goi_order_ids, col_order)

complex_scaled_list <- lapply(names(complex_data), function(cx) {
  prepare_scaled_mat(complex_data[[cx]]$ids_present, col_order)
})
names(complex_scaled_list) <- names(complex_data)

all_scaled <- do.call(rbind, c(list(mtdna_scaled), complex_scaled_list))

# --- Row groups -------------------------------------------------------------
# Rows are clustered within each group, and the groups themselves are ordered
# by clustering (ComplexHeatmap default), not by the factor levels below.
row_split_raw <- c(
  rep("mtDNA encoded", nrow(mtdna_scaled)),
  unlist(lapply(names(complex_data), function(cx) {
    rep(gsub("_", " ", cx), nrow(complex_scaled_list[[cx]]))
  }))
)

row_split_vec <- factor(
  unname(split_labels[row_split_raw]),
  levels = unname(rev(split_labels))
) %>% droplevels()

# --- Columns: age numbers only, centered over each age group ----------------
ordered_coldata <- coldata %>% arrange(age, rep)
day_labels <- ifelse(ordered_coldata$rep == 2,
                     as.character(ordered_coldata$age), "")
# Labels sit at the bottom of their annotation row, and the annotation row is
# only as tall as the text, so they sit close to the heatmap.
day_gp  <- gpar(fontsize = all_heat_font)
col_ann <- HeatmapAnnotation(
  Day = anno_text(day_labels, rot = 0, just = "bottom",
                  location = unit(0, "npc"), gp = day_gp,
                  height = max_text_height(day_labels, gp = day_gp)),
  show_annotation_name = FALSE
)
col_split <- ordered_coldata$age_fct

# --- Heatmap ----------------------------------------------------------------
# The dotted vertical line in the dendrogram marks where the group-level tree
# joins each group's own tree. It is controlled by a global option, so turn it
# off here and restore the previous setting after drawing.
old_parent_line <- ht_opt$show_parent_dend_line
ht_opt$show_parent_dend_line <- FALSE

# Gap between the day labels and the heatmap body (set explicitly; 1 mm is
# also the ComplexHeatmap default)
old_col_anno_pad <- ht_opt$COLUMN_ANNO_PADDING
ht_opt$COLUMN_ANNO_PADDING <- unit(1, "mm")

# Row group labels: ComplexHeatmap right-aligns horizontal row titles. The
# built-in titles are drawn invisible (col = NA) so they still reserve space,
# and centered labels are drawn on top of them after the heatmap is drawn.
row_lab_gp <- gpar(fontsize = all_heat_rowlab, fontface = "plain")

ht_all <- Heatmap(
  all_scaled,
  name              = "z-score",
  col               = shared_colors,
  top_annotation    = col_ann,
  show_column_names = FALSE,
  show_row_names    = FALSE,          # no individual gene labels
  cluster_columns   = FALSE,
  cluster_rows      = TRUE,
  row_split         = row_split_vec,  # groups ordered by clustering
  row_dend_gp       = gpar(lwd = all_heat_lwd),
  row_title_gp      = gpar(fontsize = all_heat_rowlab, col = NA),  # see above
  row_title_rot     = 0,              # horizontal group labels
  row_gap           = unit(all_heat_gap, "mm"),
  column_split      = col_split,
  column_title      = NULL,           # no title, no per-day split titles
  column_gap        = unit(all_heat_gap, "mm"),
  border            = TRUE,
  border_gp         = gpar(lwd = all_heat_lwd),
  heatmap_legend_param = list(
    title         = "z-score",
    title_gp      = gpar(fontsize = all_heat_legend, fontface = "plain"),
    labels_gp     = gpar(fontsize = all_heat_legend),
    at            = c(-3, -1.5, 0, 1.5, 3),
    legend_height = unit(30, "mm")
  )
)

pdf(file.path(out_dir, all_heat_file),
    width = all_heat_width, height = all_heat_height)
ht_drawn <- draw(ht_all, padding = unit(rep(all_heat_pad, 4), "mm"))

# Centered group labels, one per row group, in the order the groups are drawn
slice_names <- names(row_order(ht_drawn))
if (is.null(slice_names)) stop("Row group names not found; labels not drawn.")
for (k in seq_along(slice_names)) {
  decorate_row_title("z-score", {
    grid.text(slice_names[k], x = unit(0.5, "npc"), y = unit(0.5, "npc"),
              just = "center", gp = row_lab_gp)
  }, slice = k)
}
dev.off()

ht_opt$show_parent_dend_line <- old_parent_line
ht_opt$COLUMN_ANNO_PADDING   <- old_col_anno_pad

cat("Combined heatmap written:", all_heat_file, "\n")


# =============================================================================
# EXPRESSION TRAJECTORY RELATIVE TO DAY 1 — ALL GENE GROUPS (Fig. 2D)
# =============================================================================

group_ids <- c(
  list("mtDNA\nencoded" = goi_order_ids),
  list("House-\nkeeping" = hk_ids_present),
  # Names taken from complex_data so labels stay correct if a complex was skipped
  setNames(
    lapply(complex_data, `[[`, "ids_present"),
    gsub("_", " ", names(complex_data))
  )
)

group_labels <- c("mtDNA\nencoded",
                  "House-\nkeeping",
                  paste0("Complex ", c("I", "II", "III", "IV", "V"))
                  )

group_palette <- c(
  "mtDNA\nencoded"  = "#e41a1c",
  "House-\nkeeping" = "#999999",
  "Complex I"       = "#377eb8",
  "Complex II"      = "#ff7f00",
  "Complex III"     = "#4daf4a",
  "Complex IV"      = "#984ea3",
  "Complex V"       = "#a65628"
)

expr_long <- bind_rows(lapply(names(group_ids), function(grp) {
  ids <- group_ids[[grp]]
  ids <- ids[ids %in% rownames(vsd_mat)]
  if (length(ids) == 0) return(NULL)

  vsd_mat[ids, col_order, drop = FALSE] %>%
    as.data.frame() %>%
    rownames_to_column("gene_id") %>%
    pivot_longer(-gene_id, names_to = "sname", values_to = "vst") %>%
    left_join(coldata %>% select(sname, age, age_fct), by = "sname") %>%
    mutate(group = grp)
})) %>%
  # Subtract each gene's Day 1 mean (delta VST)
  group_by(gene_id) %>%
  mutate(
    day1_mean = mean(vst[age == 1], na.rm = TRUE),
    vst       = vst - day1_mean
  ) %>%
  ungroup() %>%
  filter(is.finite(vst))

# Per-gene mean at each age
gene_means <- expr_long %>%
  group_by(group, gene_id, age) %>%
  summarise(gene_mean = mean(vst, na.rm = TRUE), .groups = "drop") %>%
  mutate(group = factor(group, levels = group_labels))

# Group mean and SE across genes
group_summary <- gene_means %>%
  group_by(group, age) %>%
  summarise(
    group_mean = mean(gene_mean, na.rm = TRUE),
    group_se   = sd(gene_mean, na.rm = TRUE) / sqrt(n()),
    .groups    = "drop"
  ) %>%
  mutate(group = factor(group, levels = group_labels))

p_traj_all <- ggplot() +
  # Reference line at 0 = no change from Day 1 (values are VST minus Day 1 mean)
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey40", linewidth = 0.4) +
  # Individual gene lines
  geom_line(
    data      = gene_means,
    aes(x     = age, y = gene_mean,
        group = interaction(gene_id, group),
        color = group),
    alpha     = 0.15,
    linewidth = 0.3
  ) +
  # SE ribbon
  geom_ribbon(
    data  = group_summary,
    aes(x    = age,
        ymin = group_mean - group_se,
        ymax = group_mean + group_se,
        fill = group),
    alpha = 0.2
  ) +
  # Group mean line
  geom_line(
    data      = group_summary,
    aes(x     = age, y = group_mean, color = group),
    linewidth = 0.8
  ) +
  scale_color_manual(values = group_palette, name = NULL) +
  scale_fill_manual(values  = group_palette, name = NULL) +
  scale_x_continuous(breaks = ages) +
  labs(
    x     = "Age (days)",
    y     = expression(Delta*"VST Expression (~"*log[2]*" scale)")
  ) +
  theme_bw(base_size = base_font) +
  theme(
    legend.position        = c(0.5, 0.02),   # centered horizontally, near bottom
    legend.justification   = c(0.5, 0),       # anchor point at bottom center of legend
    legend.direction       = "horizontal",
    legend.background      = element_rect(fill = alpha("white", 0), color = NA),
    legend.margin          = margin(2, 2, 2, 2),
    legend.key.width       = unit(2, "mm")
  )

save_gg(p_traj_all, "trajectory_AllGroups.pdf")

cat("Trajectory plot written: trajectory_AllGroups.pdf\n")
