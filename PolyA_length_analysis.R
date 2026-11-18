# =============================================================================
# poly(A) Tail Length Visualization Script
# C. elegans aging dataset — nuclear vs. mitochondrial comparison
# =============================================================================
# Edit the CONFIG block below, then Source the script (Ctrl+Shift+S / Cmd+Shift+S)
# or run it section by section interactively.
#
# Manuscript figures produced by this script (PNG):
#   01_polya_distribution_mt_vs_nu                       Fig. 4A
#   02_polya_boxplot_mt_vs_nu                            Fig. 4B
#   03_timecourse_mean_polya                             Fig. 4C
#   04_mtdna_per_gene_trajectories                       Fig. 4D
#   05_polya_per_day_boxplot                             Supp. Fig. 4A
#   06_correlation_polya_readcount_mt_all_days_pooled    Fig. 4E
#
# Tables:
#   polya_summary_by_origin_day.csv           data behind plot 03
#   polya_wilcox_mt_vs_nu_by_day.csv          per-day mt vs nu tests (plot 05)
#   polya_readcount_spearman_mt_by_day.csv    per-day Spearman, mtDNA (plot 06)
# =============================================================================

# ── Font: Liberation Sans, registered with systemfonts for the ragg device ───
if (!requireNamespace("fontLiberation", quietly = TRUE)) install.packages("fontLiberation")
if (!requireNamespace("ragg",           quietly = TRUE)) install.packages("ragg")

font_dir <- system.file(package = "fontLiberation")
find_ttf <- function(f) list.files(font_dir, pattern = f, recursive = TRUE, full.names = TRUE)

systemfonts::register_font(
  name       = "Liberation Sans",
  plain      = find_ttf("^LiberationSans-Regular\\.ttf$"),
  bold       = find_ttf("^LiberationSans-Bold\\.ttf$"),
  italic     = find_ttf("^LiberationSans-Italic\\.ttf$"),
  bolditalic = find_ttf("^LiberationSans-BoldItalic\\.ttf$")
)

systemfonts::registry_fonts()   # check: Liberation Sans should be listed


# ── CONFIG — edit these paths before running ──────────────────────────────────
opt <- list(
  excel  = "~/aging_rnaseq/06_all_days_polya_summary.xlsx",
  gtf    = "~/aging_rnaseq/reference/c_elegans.PRJNA13758.WS279.canonical_geneset.gtf.gz",
  outdir = "~/aging_rnaseq/polya_plots_2",
  dpi    = 600
)

file.exists(path.expand(opt$excel))
file.exists(path.expand(opt$gtf))

# ─────────────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(purrr)
  library(ragg)
})

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# ── Shared plot theme ─────────────────────────────────────────────────────────
theme_polya <- function() {
  theme_classic(base_size = 6) +
    theme(
      plot.title       = element_text(size = 8, hjust = 0.5),
      axis.title       = element_text(size = 8),
      axis.text        = element_text(size = 7, color = "black"),
      axis.ticks       = element_line(color = "black"),
      legend.title     = element_text(size = 8),
      legend.text      = element_text(size = 8),
      legend.key.size  = unit(0.3, "cm"),
      plot.margin      = margin(3, 3, 3, 3),
      strip.text       = element_text(size = 8),
      text = element_text(family = "Liberation Sans")
    )
}

MT_COLOR <- "#D55E00"   # orange-red — mtDNA
NU_COLOR <- "#0072B2"   # blue       — nuDNA

# ── Save helper (PNG via ragg, so the registered font is used) ────────────────
save_plot <- function(p, name, width, height) {
  path <- file.path(opt$outdir, paste0(name, ".png"))
  ragg::agg_png(path, width = width, height = height, units = "in",
                res = opt$dpi, background = "white")
  print(p)
  dev.off()
  message("  Saved: ", path)
}

# =============================================================================
# 1. PARSE GTF — extract gene_id → biotype and MtDNA gene IDs
# =============================================================================
message("\n[1/5] Parsing GTF for gene biotypes and MtDNA gene IDs...")

parse_gtf_genes <- function(gtf_path) {
  gtf_path <- path.expand(gtf_path)
  con <- if (grepl("\\.gz$", gtf_path)) gzcon(file(gtf_path, "rb")) else file(gtf_path, "r")
  on.exit(close(con))

  result <- list()

  repeat {
    lines <- readLines(con, n = 50000, warn = FALSE)
    if (length(lines) == 0) break

    # Keep only gene-level rows
    gene_lines <- lines[grepl("\tgene\t", lines) & !grepl("^#", lines)]
    if (length(gene_lines) == 0) next

    # GTF columns: seqname, source, feature, start, end, score, strand, frame, attributes
    fields <- strsplit(gene_lines, "\t")

    # Extract chromosome (col 1)
    chroms <- sapply(fields, `[[`, 1)

    # Extract gene_id
    ids <- regmatches(gene_lines,
                      regexpr('gene_id "([^"]+)"', gene_lines))
    ids <- sub('gene_id "', "", ids)
    ids <- sub('"$', "", ids)

    # Extract biotype (gene_biotype field in WormBase GTFs)
    biotypes <- ifelse(
      grepl('gene_biotype "', gene_lines),
      sub('.*gene_biotype "([^"]+)".*', "\\1", gene_lines),
      ifelse(
        grepl('biotype "', gene_lines),
        sub('.*biotype "([^"]+)".*', "\\1", gene_lines),
        NA_character_
      )
    )

    # Flag MtDNA chromosome
    is_mt <- chroms == "MtDNA"

    result[[length(result) + 1]] <- data.frame(
      gene_id  = ids,
      chrom    = chroms,
      biotype  = biotypes,
      is_mt    = is_mt,
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, result)
}

gtf_genes   <- parse_gtf_genes(opt$gtf)
mt_gene_ids <- gtf_genes$gene_id[gtf_genes$is_mt]

message("  Found ", length(mt_gene_ids), " MtDNA gene IDs in GTF.")
message("  GTF biotype table: ", nrow(gtf_genes), " genes.")
message("  Biotypes present: ", paste(sort(unique(gtf_genes$biotype)), collapse = ", "))

# Build lookup: gene_id → gtf_biotype
biotype_lookup <- gtf_genes %>%
  select(gene_id, gtf_biotype = biotype) %>%
  distinct()

# =============================================================================
# 2. LOAD & TIDY EXCEL DATA
# =============================================================================
message("\n[2/5] Loading Excel sheets...")

day_sheets <- excel_sheets(path.expand(opt$excel))
day_sheets <- day_sheets[grepl("^day", day_sheets, ignore.case = TRUE)]

read_day_sheet <- function(sheet) {
  df <- read_xlsx(path.expand(opt$excel), sheet = sheet,
                  col_names = TRUE, skip = 1)
  df <- df %>%
    rename(
      wb_id        = "WB_ID",
      isoform_name = "isoform name",
      read_count   = "read count ( 3 reps combined)",
      median_polya = "median poly(A)",
      min_polya    = "min poly(A)",
      max_polya    = "max poly(A)",
      day          = "day",
      gene_name    = "gene name",
      biotype      = "biotype"
    )
  df <- df %>%
    mutate(across(c(read_count, median_polya, min_polya, max_polya),
                  as.numeric))
  df
}

raw <- map_dfr(day_sheets, read_day_sheet)

# Drop unannotated isoforms (non-WBGene rows)
n_before <- nrow(raw)
raw <- raw %>% filter(grepl("^WBGene", wb_id))
n_after  <- nrow(raw)
message("  Dropped ", n_before - n_after, " unannotated isoforms; ",
        n_after, " rows retained.")

# Parse day number and origin; add GTF biotype
raw <- raw %>%
  mutate(
    day_num = as.integer(gsub("[^0-9]", "", day)),
    is_mt   = wb_id %in% mt_gene_ids,
    origin  = if_else(is_mt, "mtDNA", "nuDNA")
  ) %>%
  left_join(biotype_lookup, by = c("wb_id" = "gene_id"))
# gtf_biotype (not the Excel biotype column) is used for all filtering below

message("  GTF biotype join complete.")
message("  GTF biotype breakdown (all rows):")
print(table(raw$gtf_biotype, useNA = "ifany"))

# =============================================================================
# 3. FILTERED VIEWS
#    - tx_df     : all rows (transcript-level, no collapsing)
#    - pc_df     : protein_coding only (both origins) — plots 01–03, 05, 06
#    - mt_pc_df  : mtDNA protein_coding only          — plot 04
# =============================================================================
message("\n[3/5] Building filtered data views...")

tx_df    <- raw   # no collapsing; each row = one transcript observation

pc_df    <- tx_df %>% filter(gtf_biotype == "protein_coding")
mt_pc_df <- pc_df %>% filter(origin == "mtDNA")

days_ordered <- sort(unique(tx_df$day_num))

message("  protein_coding rows (both origins): ", nrow(pc_df))
message("  mtDNA protein_coding rows: ", nrow(mt_pc_df))

# =============================================================================
# 4. PLOTS
# =============================================================================
message("\n[4/5] Generating plots...")

# ── Plot 01: poly(A) length distribution — mt vs nu (protein_coding) ─────────
# All days pooled: each transcript contributes one value per day it was detected.
message("  Plot 01: Distribution comparison mt vs nu (protein_coding)...")

p1 <- pc_df %>%
  ggplot(aes(x = median_polya, fill = origin, colour = origin)) +
  geom_density(alpha = 0.45, linewidth = 0.4) +
  scale_fill_manual(values   = c(mtDNA = MT_COLOR, nuDNA = NU_COLOR), name = NULL) +
  scale_colour_manual(values = c(mtDNA = MT_COLOR, nuDNA = NU_COLOR), name = NULL) +
  scale_x_continuous(limits = c(0, NA)) +
  labs(
    x     = "Median poly(A) (nt)",
    y     = "Density"
  ) +
  theme_polya() +
  theme(
    legend.position        = "inside",
    legend.position.inside = c(0.90, 0.90),
    legend.justification   = c(1, 1),
    legend.background      = element_blank(),
    plot.margin            = margin(4, 8, 4, 4)
  )

save_plot(p1, "01_polya_distribution_mt_vs_nu", width = 2.24, height = 2.05)

# ── Plot 02: Boxplot mt vs nu — protein_coding, Day 1 ────────────────────────
message("  Plot 02: Boxplot mt vs nu (protein_coding, Day 1)...")

p2_df <- pc_df %>% filter(day_num == 1)

# Two-sided Wilcoxon rank-sum test. With these sample sizes and ties, R uses
# the normal approximation with continuity correction, not an exact test.
# Full result is printed below for the figure legend; the plot shows a
# bracket with the p-value only.
p2_test <- wilcox.test(median_polya ~ origin, data = p2_df)
p2_n    <- table(p2_df$origin)
p2_med  <- tapply(p2_df$median_polya, p2_df$origin, median)
p2_avg  <- tapply(p2_df$median_polya, p2_df$origin, mean)
p2_min  <- tapply(p2_df$median_polya, p2_df$origin, min)
p2_max  <- tapply(p2_df$median_polya, p2_df$origin, max)

message("    W = ", round(p2_test$statistic, 1),
        ", p = ", format.pval(p2_test$p.value, digits = 3),
        "; n(mtDNA) = ", p2_n[["mtDNA"]], ", n(nuDNA) = ", p2_n[["nuDNA"]],
        "; median mtDNA = ", round(p2_med[["mtDNA"]], 1), " nt",
        ", nuDNA = ", round(p2_med[["nuDNA"]], 1), " nt")
message("    mean mtDNA = ", round(p2_avg[["mtDNA"]], 1),
        ", nuDNA = ", round(p2_avg[["nuDNA"]], 1), " nt",
        "; range mtDNA = ", p2_min[["mtDNA"]], "–", p2_max[["mtDNA"]],
        ", nuDNA = ", p2_min[["nuDNA"]], "–", p2_max[["nuDNA"]], " nt")

p2 <- p2_df %>%
  ggplot(aes(x = origin, y = median_polya,
             fill = origin, colour = origin)) +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.4,
               linewidth = 0.4, alpha = 0.6, width = 0.5) +
  ggpubr::geom_bracket(
    xmin       = "mtDNA",
    xmax       = "nuDNA",
    y.position = max(p2_df$median_polya, na.rm = TRUE) * 1.04,
    label = if (p2_test$p.value < 0.001) "p < 0.001" else paste0("p = ", signif(p2_test$p.value, 2)),
    label.size = 2,
    tip.length = 0.02,
    size       = 0.3,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(values   = c(mtDNA = MT_COLOR, nuDNA = NU_COLOR)) +
  scale_colour_manual(values = c(mtDNA = MT_COLOR, nuDNA = NU_COLOR)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
  guides(fill = "none", colour = "none") +
  labs(
    title = "Day 1",
    x     = NULL,
    y     = "Median poly(A) (nt)"
  ) +
  theme_polya()

save_plot(p2, "02_polya_boxplot_mt_vs_nu", width = 1.5, height = 2)


# ── Plot 03: Time-course — mean poly(A) ± SD, mt vs nu (protein_coding) ──────
# Mean and SD across transcripts of each transcript's median poly(A).
# The same values are written to polya_summary_by_origin_day.csv (section 5).
message("  Plot 03: Time-course mean poly(A) (protein_coding)...")

timecourse <- pc_df %>%
  group_by(origin, day_num) %>%
  summarise(
    mean_pa = mean(median_polya, na.rm = TRUE),
    sd_pa   = sd(median_polya,   na.rm = TRUE),
    .groups = "drop"
  )

p3 <- timecourse %>%
  ggplot(aes(x = day_num, y = mean_pa,
             colour = origin, fill = origin)) +
  geom_ribbon(aes(ymin = mean_pa - sd_pa, ymax = mean_pa + sd_pa),
              alpha = 0.2, colour = NA) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 0.8) +
  scale_x_continuous(breaks = days_ordered) +
  scale_colour_manual(values = c(mtDNA = MT_COLOR, nuDNA = NU_COLOR), name = NULL) +
  scale_fill_manual(values   = c(mtDNA = MT_COLOR, nuDNA = NU_COLOR), name = NULL) +
  labs(
    x     = "Day",
    y     = "Mean median poly(A) ± SD (nt)"
  ) +
  theme_polya() +
  theme(
    legend.position  = "top",
    legend.direction = "horizontal",
    legend.margin    = margin(0, 0, 0, 0),
    plot.margin      = margin(3, 3, 3, 3)
  )

save_plot(p3, "03_timecourse_mean_polya", width = 2.7, height = 2)

# ── Plot 04: Per-gene trajectories — mtDNA protein_coding ────────────────────
message("  Plot 04: Per-gene mtDNA trajectories (protein_coding)...")

p4 <- mt_pc_df %>%
  ggplot(aes(x = day_num, y = median_polya,
             colour = gene_name, group = gene_name)) +
  geom_line(linewidth = 0.4, alpha = 0.8) +
  geom_point(size = 0.6) +
  scale_x_continuous(breaks = days_ordered) +
  labs(
    x     = "Day",
    y     = "Median poly(A) (nt)",
    colour = "Gene"
  ) +
  theme_polya() +
  theme(legend.key.size      = unit(0.25, "cm"),
        legend.text          = element_text(size = 7),
        legend.key.spacing.y = unit(1.5, "pt"))

save_plot(p4, "04_mtdna_per_gene_trajectories", width = 3.25, height = 2.25)

# ── Plot 05: Per-day boxplot — mt vs nu (protein_coding) ─────────────────────
message("  Plot 05: Per-day boxplots mt vs nu (protein_coding)...")

p5 <- pc_df %>%
  ggplot(aes(x = factor(day_num), y = median_polya,
             fill = origin, colour = origin)) +
  geom_boxplot(outlier.size = 0.2, outlier.alpha = 0.3,
               linewidth = 0.3, alpha = 0.55,
               position = position_dodge(width = 0.7), width = 0.55) +
  scale_fill_manual(values   = c(mtDNA = MT_COLOR, nuDNA = NU_COLOR), name = NULL) +
  scale_colour_manual(values = c(mtDNA = MT_COLOR, nuDNA = NU_COLOR), name = NULL) +
  labs(
    title = "poly(A) by day: mt vs nu\n(protein-coding)",
    x     = "Day",
    y     = "Median poly(A) (nt)"
  ) +
  theme_polya()

save_plot(p5, "05_polya_per_day_boxplot", width = 4, height = 2.2)

# ── Statistics for plot 05: mt vs nu, tested separately within each day ──────
# Wilcoxon rank-sum (normal approximation with continuity correction), BH
# correction across days. Not drawn on the plot.
polya_wilcox_by_day <- pc_df %>%
  filter(!is.na(median_polya)) %>%
  group_by(day_num) %>%
  group_modify(~ {
    mt <- .x$median_polya[.x$origin == "mtDNA"]
    nu <- .x$median_polya[.x$origin == "nuDNA"]

    if (length(mt) < 3 || length(nu) < 3) {
      return(tibble(n_mt = length(mt), n_nu = length(nu),
                    median_mt = NA_real_, median_nu = NA_real_,
                    diff_median = NA_real_, W = NA_real_,
                    p = NA_real_, rank_biserial = NA_real_))
    }

    wt <- wilcox.test(mt, nu, exact = FALSE, correct = TRUE)
    U  <- unname(wt$statistic)

    tibble(
      n_mt          = length(mt),
      n_nu          = length(nu),
      median_mt     = median(mt),
      median_nu     = median(nu),
      diff_median   = median(mt) - median(nu),
      W             = U,
      p             = wt$p.value,
      # rank-biserial correlation: −1 to 1, positive = mtDNA longer
      rank_biserial = 2 * U / (length(mt) * length(nu)) - 1
    )
  }) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p, method = "BH")) %>%
  arrange(day_num)

print(as.data.frame(polya_wilcox_by_day))

write.csv(polya_wilcox_by_day,
          file.path(opt$outdir, "polya_wilcox_mt_vs_nu_by_day.csv"),
          row.names = FALSE)

# ── Plot 06: mtDNA poly(A) vs read count — all days layered, pooled fit ──────
# All timepoints overlaid on one set of axes, points colored by age, with a
# single linear fit across all days.
message("  Plot 06: mtDNA poly(A) vs read count, all days (protein_coding)...")

COR_MIN_READS <- 50   # isoforms below this combined read count are dropped

# Yellow-to-red age palette, matching the DESeq2 script
age_palette_polya <- colorRampPalette(c("#f5c842", "#800026"))(length(days_ordered))
names(age_palette_polya) <- as.character(days_ordered)

cor_all_mt <- mt_pc_df %>%
  filter(!is.na(median_polya), !is.na(read_count),
         read_count >= COR_MIN_READS) %>%
  mutate(day_f = factor(day_num, levels = days_ordered))

cat("  Isoforms per day at >=", COR_MIN_READS, "reads:\n")
print(table(cor_all_mt$day_num))

# Spearman within each day (valid per-day statistics), BH across days
spearman_mt_by_day <- cor_all_mt %>%
  group_by(day_num) %>%
  group_modify(~ {
    if (nrow(.x) < 5) {
      return(tibble(n = nrow(.x), rho = NA_real_, p = NA_real_))
    }
    ct <- cor.test(log10(.x$read_count + 1), .x$median_polya,
                   method = "spearman", exact = FALSE)
    tibble(n = nrow(.x), rho = unname(ct$estimate), p = ct$p.value)
  }) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p, method = "BH"))

print(as.data.frame(spearman_mt_by_day))
write.csv(spearman_mt_by_day,
          file.path(opt$outdir, "polya_readcount_spearman_mt_by_day.csv"),
          row.names = FALSE)

# Spearman on the pooled points, printed for reference only. It is not a valid
# test: each gene contributes up to one point per day, so points are not
# independent. Per-day statistics are in spearman_mt_by_day above.
p6_cor <- cor.test(log10(cor_all_mt$read_count + 1), cor_all_mt$median_polya,
                   method = "spearman", exact = FALSE)
cat("  Pooled Spearman rho = ", round(p6_cor$estimate, 3),
    ", p = ", format.pval(p6_cor$p.value, digits = 3),
    ", n = ", nrow(cor_all_mt), "\n", sep = "")

p6 <- ggplot(cor_all_mt, aes(x = log10(read_count + 1), y = median_polya)) +
  geom_point(aes(colour = day_f), size = 0.8, alpha = 0.75) +
  geom_smooth(method = "lm", formula = y ~ x,
              colour = "grey30", linewidth = 0.4, se = TRUE) +
  scale_colour_manual(values = age_palette_polya, name = "Day") +
  labs(
    x     = expression(log[10]~"(read count + 1)"),
    y     = "Median poly(A) (nt)"
  ) +
  theme_polya() +
  theme(legend.key.size = unit(0.25, "cm"))

save_plot(p6, "06_correlation_polya_readcount_mt_all_days_pooled",
          width = 3.1, height = 2.25)


# =============================================================================
# 5. SUMMARY TABLE (data behind plot 03)
# =============================================================================
message("\n[5/5] Writing summary table...")

summary_tbl <- pc_df %>%
  group_by(origin, day_num) %>%
  summarise(
    n_transcripts  = n(),
    mean_polya     = round(mean(median_polya,   na.rm = TRUE), 2),
    median_polya   = round(median(median_polya, na.rm = TRUE), 2),
    sd_polya       = round(sd(median_polya,     na.rm = TRUE), 2),
    total_reads    = sum(read_count, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(summary_tbl,
          file.path(opt$outdir, "polya_summary_by_origin_day.csv"),
          row.names = FALSE, quote = FALSE)

message("\nDone. All outputs written to: ", opt$outdir)
message("Files generated:")
message("  Plots 01–06 (png, ", opt$dpi, " dpi)")
message("  polya_summary_by_origin_day.csv")
message("  polya_wilcox_mt_vs_nu_by_day.csv")
message("  polya_readcount_spearman_mt_by_day.csv")

