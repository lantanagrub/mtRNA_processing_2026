# =============================================================================
# Salmon RNA-seq DGE Pipeline — Gene-Level, Illumina
# C. elegans aging timecourse: Day 1, 2, 3, 4, 5, 7, 10, 15
# Comparisons: each age vs Day 1
#
# Data processing and DESeq2 (sections 1–8) mirror the Nanopore oarfish
# pipeline; differences are marked "ILLUMINA:".
#
# Manuscript figure produced by this script:
#   PCA_gene_level.pdf          Supp. Fig. 2B
#
# The per-comparison DGE CSVs (section 7) are also read by platform_comparison.R.
# =============================================================================

library(DESeq2)
library(tidyverse)
library(rtracklayer)
library(ggplot2)
library(ggrepel)

# =============================================================================
# 0. PATHS AND PARAMETERS — edit this block
# =============================================================================

setwd("~/aging_rnaseq")

data_dir        <- "~/aging_rnaseq/salmon"
gtf_path        <- "~/aging_rnaseq/reference/c_elegans.PRJNA13758.WS279.canonical_geneset.gtf.gz"
gene_names_path <- "~/aging_rnaseq/reference/c_elegans.PRJNA13758.WS279.geneIDs.txt.gz"

# ILLUMINA: separate output directory so file names do not collide with the
# Nanopore outputs (which are written to ~/aging_rnaseq/nanopore_dge).
out_dir <- "~/aging_rnaseq/illumina_dge"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ILLUMINA: path to each sample's quant.sf, relative to data_dir.
# Edit to match your Salmon output layout.
quant_path <- function(age, rep) {
  file.path(paste0("illumina_day", age, "_rep", rep, "_quant"), "quant.sf")
}

# ILLUMINA: optional regex stripped from Salmon transcript names so they match
# GTF transcript_id (e.g. "^transcript:"). NULL = no change. The coverage check
# in section 4 will show whether this is needed.
tx_id_strip_regex <- NULL

# ILLUMINA: aggregation method.
#   FALSE = sum NumReads per gene and round (identical to the oarfish pipeline;
#           keeps the two platforms methodologically parallel)
#   TRUE  = tximport with average-transcript-length offsets (standard for
#           Salmon; corrects for isoform-usage-driven changes in effective
#           gene length across samples)
use_tximport <- FALSE

# Set to TRUE to restrict DGE to protein-coding genes only
coding_only <- TRUE

ages <- c(1, 2, 3, 4, 5, 7, 10, 15)
reps <- c(1, 2, 3)

# Low-count filter: gene must have >= this many counts in >= this many samples
min_counts  <- 10
min_samples <- 3

# LFC threshold for calling significance in the summary table
lfc_thresh <- 1

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
# 4. LOAD AND MERGE SALMON COUNTS
# =============================================================================

cat("Loading Salmon quant.sf files...\n")

clean_tx_id <- function(x) {
  if (!is.null(tx_id_strip_regex)) x <- str_remove(x, tx_id_strip_regex)
  x
}

coldata <- expand.grid(age = ages, rep = reps) %>%
  mutate(
    filename = map2_chr(age, rep, quant_path),
    sname    = paste0("Day", age, "_rep", rep),
    age_fct  = factor(paste0("Day", age), levels = paste0("Day", ages))
  ) %>%
  arrange(age, rep)
rownames(coldata) <- coldata$sname

# ILLUMINA: quant.sf columns are Name, Length, EffectiveLength, TPM, NumReads
count_list <- lapply(seq_len(nrow(coldata)), function(i) {
  f <- file.path(data_dir, coldata$filename[i])
  s <- coldata$sname[i]
  if (!file.exists(f)) stop("File not found: ", f)
  read_tsv(f, show_col_types = FALSE) %>%
    transmute(tname = clean_tx_id(Name), !!s := NumReads)
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

if (!use_tximport) {

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

} else {

  # ILLUMINA: tximport route. Read at transcript level, clean IDs, then
  # summarize to gene level so the same tx2gene_use filtering applies.
  library(tximport)

  files <- setNames(file.path(data_dir, coldata$filename), coldata$sname)
  txi_tx <- tximport(files, type = "salmon", txOut = TRUE)
  for (m in c("abundance", "counts", "length")) {
    rownames(txi_tx[[m]]) <- clean_tx_id(rownames(txi_tx[[m]]))
  }
  txi <- summarizeToGene(txi_tx, tx2gene = tx2gene_use %>% select(tname, gene_id))

  cat("Genes after aggregation:", nrow(txi$counts), "\n")

  stopifnot(all(colnames(txi$counts) == rownames(coldata)))
}

# =============================================================================
# 6. DESeq2
# =============================================================================

cat("Running DESeq2...\n")

if (!use_tximport) {
  dds <- DESeqDataSetFromMatrix(
    countData = counts_gene,
    colData   = coldata,
    design    = ~ age_fct
  )
} else {
  dds <- DESeqDataSetFromTximport(
    txi,
    colData = coldata,
    design  = ~ age_fct
  )
}

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

# Combined long-format results
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
# require |LFC| > lfc_thresh.

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
# 9. PCA (Supp. 2B)
# =============================================================================

age_palette <- colorRampPalette(c("#f5c842", "#800026"))(length(ages))
names(age_palette) <- paste0("Day", ages)

# VST (not blind to design)
vsd <- vst(dds, blind = FALSE)

pca_data <- plotPCA(vsd, intgroup = "age_fct", returnData = TRUE)
pct_var  <- round(100 * attr(pca_data, "percentVar"))

# ILLUMINA: points are labeled with sample names (the Nanopore PCA is not)
p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = age_fct, label = name)) +
  geom_point(size = 3.5) +
  geom_text_repel(size = 3, max.overlaps = 20) +
  scale_color_manual(values = age_palette, name = "Age") +
  labs(
    title = "PCA — VST Gene Expression",
    x     = paste0("PC1: ", pct_var[1], "% variance"),
    y     = paste0("PC2: ", pct_var[2], "% variance")
  ) +
  theme_bw(base_size = 12)

ggsave(file.path(out_dir, "PCA_gene_level.pdf"), p_pca, width = 5, height = 4.5)
