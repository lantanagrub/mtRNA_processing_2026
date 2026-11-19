#!/usr/bin/env python3
"""
five_prime_extension.py
=======================
Quantify failure to cleave transcript 5' ends, per gene and per sample.

Direct RNA sequences 3'->5', so a read's 3' end (reference_end) is precise and
its 5' end (reference_start) is potentially subject to truncation. To avoid 
interpreting truncated reads as correctly trimmed reads, need to filter out 
truncated reads and retain full-length reads.

    1. Assign each read to a gene by its 3' end (the reliable end).
    2. Require the read to cover >= min_coverage of that gene. This removes
       reads that lost their 5' end to truncation, which would otherwise be
       miscounted as correctly trimmed.
    3. Among surviving reads, classify the 5' end:
         TRIMMED    starts at the gene's 5' boundary (+/- tolerance)
         EXTENSION  starts upstream of it -> 5' cleavage failed

    extension_fraction = EXTENSION / (EXTENSION + TRIMMED)

Manuscript figures (with --out .../5_prime_extension_mRNA, --biotypes mrna):
  <out>_per_sample.png         Fig. 5A
  <out>_dumbbell.png           Fig. 5B
  <out>_extension_counts.png   Fig. 5C
  <out>_read_integrity.png     Supp. Fig. 5A-5D

Usage
-----
    python five_prime_extension.py \
        --features $MTRNA_ROOT/ref/mtdna_features.tsv \
        --bam-dir  $MTRNA_ROOT/elac2/bam/mito \
        --samples  $MTRNA_SRC/config/samples/elac2.tsv \
        --out      $MTRNA_ROOT/elac2/results/processing/5_prime_extension_mRNA \
        [--per-sample-size 1.5 2.5] [--dumbbell-size 2.3 2.8] \
        [--extension-counts-size 2.7 2.8] [--read-integrity-size 6.0 7.7]

Figures are PNG at 600 dpi in Liberation Sans, saved at exactly the size given.
Axis and legend text is set in the LABELS block near the top of the file.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import font_manager
from matplotlib.legend_handler import HandlerTuple
from matplotlib.ticker import FuncFormatter, PercentFormatter

import pysam


GENOTYPE_PALETTE = ["#37474F", "#C62828", "#1565C0", "#2E7D32"]

# Bar segment for reads WITHOUT a 5' extension in the extension-count plot.
# Opaque, shared by every genotype; the extension segment takes the genotype
# colour.
NO_EXTENSION_COLOR = "#CFD8DC"

# First family found is used for every plot.
FONT_FAMILY = ["Liberation Sans", "Arial", "DejaVu Sans"]

# Figure text
LABELS = {
    # y axis of the per-sample (summary) plot, x axis of the dumbbell plot
    "extension_pct":   "mt-mRNA reads with 5' extension",
    # extension-count plot
    "counts_x":        "Reads (summed across replicates)",
    "with_extension":  "5' extension",
    "no_extension":    "No 5' extension",
}

# Default figure sizes in inches (width, height). Each is overridable on the
# command line with the matching --*-size W H flag. The saved PNG is exactly
# this size: layout is constrained, not cropped with bbox_inches="tight".
SIZES = {
    "dumbbell":          (2.3, 2.8),
    "extension_counts":  (2.7, 2.8),   # same height as the dumbbell, so rows align
    "per_sample":        (1.5, 2.5),
    "read_integrity":    (6.0, 7.7),
}

# Which features enter the figures and the aggregated tables. tRNAs are excluded
# by default: a tRNA's own 5' end is the RNase P site, and mature tRNAs are
# largely absent from direct-RNA libraries because they carry no polyA tail for
# adapter ligation, so tRNA rows measure a different molecular population.
BIOTYPE_SETS = {
    "mrna":      ["protein_coding"],
    "mrna+rrna": ["protein_coding", "rRNA"],
    "all":       ["protein_coding", "rRNA", "non_coding", "tRNA"],
}


# ── Feature context ──────────────────────────────────────────────────────────

def build_context(features):
    """Order features and record what sits immediately upstream of each.

    Coordinates stay 1-based inclusive, as produced by mtdna_features.py.
    The genome is circular, so the first feature's upstream neighbour is the
    last one.
    """
    f = features.sort_values("start").reset_index(drop=True)
    n = len(f)

    up_label, up_biotype, up_is_trna = [], [], []
    for i in range(n):
        j = i - 1 if i > 0 else n - 1
        up_label.append(f.loc[j, "label"])
        up_biotype.append(f.loc[j, "biotype"])
        up_is_trna.append(f.loc[j, "biotype"] == "tRNA")

    f["upstream_label"] = up_label
    f["upstream_biotype"] = up_biotype
    f["upstream_is_trna"] = up_is_trna
    return f


# ── Read classification ──────────────────────────────────────────────────────

def classify_reads(bam_path, contig, feats, tol, min_cov, min_alignment):
    """Assign reads to genes by 3' end, then classify the 5' end.

    Returns per-(gene, class) counts plus the extension length distribution.
    """
    if not Path(str(bam_path) + ".bai").exists():
        sys.exit(f"[ERROR] BAM index missing. Run: samtools index {bam_path}")

    ends = feats["end"].to_numpy()
    starts = feats["start"].to_numpy()
    lengths = (ends - starts + 1).astype(float)
    labels = feats["label"].to_numpy()

    # Three nested populations per gene, each a denominator for a different
    # question:
    #   aligned       any read aligned to >= min_alignment nt of the gene, however
    #                 its ends fall. The broadest view of reads at this locus.
    #   anchored      of those, reads whose 3' end matched this gene's 3'
    #                 boundary within tolerance
    #   trimmed +     of those, reads that also cleared the full-length filter
    #   extension     - the reads the extension fraction is computed from
    #
    # A read aligned to several genes counts toward `aligned` for each of
    # them, so that column sums to more than the library size. It is a per-gene
    # denominator, not a partition of the reads.
    counts = {lab: dict(trimmed=0, extension=0, anchored=0, below_coverage=0,
                        aligned=0)
              for lab in labels}
    extensions = {lab: [] for lab in labels}
    n_total = n_no3p = n_short = 0

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam.fetch(contig):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.reference_end is None:
                continue
            n_total += 1

            # pysam: reference_start is 0-based, reference_end is 0-based
            # exclusive. Converting to 1-based inclusive to match the GTF:
            r_end = read.reference_end
            r_start = read.reference_start + 1

            # 0. Every gene this read covers, before any filtering.
            aligned_all = (np.minimum(r_end, ends)
                           - np.maximum(r_start, starts) + 1)
            for lab_al in labels[aligned_all >= min_alignment]:
                counts[lab_al]["aligned"] += 1

            # 1. Anchor on the 3' end. Where gene ends sit close together,
            #    several may be within tolerance; the closest wins.
            diffs = np.abs(ends - r_end)
            k = int(np.argmin(diffs))
            if diffs[k] > tol:
                n_no3p += 1
                continue

            lab = labels[k]
            counts[lab]["anchored"] += 1

            # 2. Full-length requirement, so 5'-truncated reads are dropped
            #    rather than being counted as trimmed.
            aligned_len = min(r_end, ends[k]) - max(r_start, starts[k]) + 1
            if aligned_len / lengths[k] < min_cov:
                n_short += 1
                counts[lab]["below_coverage"] += 1
                continue

            # 3. Positive extension = read begins upstream of the 5' boundary.
            oh = int(starts[k] - r_start)
            if oh > tol:
                counts[lab]["extension"] += 1
                extensions[lab].append(oh)
            else:
                counts[lab]["trimmed"] += 1

    return counts, extensions, dict(total=n_total, no_3p_match=n_no3p,
                                   below_coverage=n_short)


def collect(feats, bam_dir, suffix, contig, samples, tol, min_cov,
            min_reads, min_alignment):
    records = []
    sample_stats = []
    for s in samples.itertuples():
        path = Path(bam_dir) / f"{s.sample_id}{suffix}"
        if not path.exists():
            sys.exit(f"[ERROR] BAM not found: {path}")

        counts, extensions, stats = classify_reads(
            path, contig, feats, tol, min_cov, min_alignment)

        kept = sum(c["trimmed"] + c["extension"] for c in counts.values())
        print(f"[INFO] {s.sample_id} ({s.genotype}): {stats['total']:,} reads "
              f"-> {kept:,} assigned "
              f"({stats['no_3p_match']:,} no 3' match, "
              f"{stats['below_coverage']:,} below {min_cov:.0%} coverage)")

        sample_stats.append(dict(
            sample_id=s.sample_id,
            genotype=s.genotype,
            genotype_label=getattr(s, "genotype_label", s.genotype),
            n_reads_total=stats["total"],
            n_no_3p_match=stats["no_3p_match"],
            n_below_coverage=stats["below_coverage"],
            n_qualified=kept,
        ))

        for f in feats.itertuples():
            c = counts[f.label]
            tot = c["trimmed"] + c["extension"]
            oh = extensions[f.label]
            records.append(dict(
                sample_id=s.sample_id,
                genotype=s.genotype,
                genotype_label=getattr(s, "genotype_label", s.genotype),
                label=f.label,
                gene_id=f.gene_id,
                biotype=f.biotype,
                start=f.start, end=f.end,
                upstream_label=f.upstream_label,
                upstream_biotype=f.upstream_biotype,
                upstream_is_trna=f.upstream_is_trna,
                n_aligned_to_gene=c["aligned"],
                n_anchored=c["anchored"],
                n_below_coverage=c["below_coverage"],
                n_trimmed=c["trimmed"],
                n_extension=c["extension"],
                n_assigned=tot,
                # Of ALL reads aligned to this gene, the share whose 3' end
                # sits at the gene's annotated 3' boundary. Low values mean
                # transcript 3' ends are not where the annotation puts them.
                intact_3_of_total=(c["aligned"]
                                   and c["anchored"] / c["aligned"]),
                # Of reads with an intact 3' end, the share that also reach far
                # enough back to be full-length - i.e. did not lose their 5'
                # end to truncation. This is a library-quality measure.
                intact_5_of_intact_3=(c["anchored"] and tot / c["anchored"]),
                # Of ALL reads aligned to this gene, the share with both ends
                # intact. The product of the two above, and the only reads the
                # extension fraction can be computed from.
                intact_5_and_3_of_total=(c["aligned"]
                                         and tot / c["aligned"]),
                extension_fraction=(c["extension"] / tot
                                   if tot >= min_reads else np.nan),
                median_extension_bp=(float(np.median(oh)) if oh else np.nan),
                max_extension_bp=(int(np.max(oh)) if oh else np.nan),
            ))
    return pd.DataFrame(records), pd.DataFrame(sample_stats)


# ── Subsampling ──────────────────────────────────────────────────────────────

def report_depth(samples, bam_dir, suffix, contig):
    print("\n[INFO] Mito depth per sample:")
    depths = {}
    for s in samples.itertuples():
        path = Path(bam_dir) / f"{s.sample_id}{suffix}"
        with pysam.AlignmentFile(str(path), "rb") as bam:
            n = sum(1 for r in bam.fetch(contig)
                    if not (r.is_unmapped or r.is_secondary or r.is_supplementary))
        depths[s.sample_id] = n
        print(f"         {s.sample_id} ({s.genotype}): {n:>9,}")
    lo, hi = min(depths.values()), max(depths.values())
    print(f"         range: {lo:,} - {hi:,} ({hi / lo:.1f}x)")
    print("         Rates below are per-gene proportions, so depth affects "
          "precision,\n         not expectation. Genes below --min-reads are "
          "reported as NA.\n")
    return depths


# ── Plot ─────────────────────────────────────────────────────────────────────

RC = {
    "font.size": 7, "axes.labelsize": 7.5, "axes.titlesize": 8,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 6.5,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "pdf.fonttype": 42, "ps.fonttype": 42,
}


def setup_fonts():
    """Use the first installed family in FONT_FAMILY for all figure text."""
    found = None
    for fam in FONT_FAMILY:
        try:
            font_manager.findfont(font_manager.FontProperties(family=fam),
                                  fallback_to_default=False)
            found = fam
            break
        except ValueError:
            continue
    found = found or "DejaVu Sans"
    if found != FONT_FAMILY[0]:
        print(f"[WARN] Font '{FONT_FAMILY[0]}' not found; using '{found}'.")
    matplotlib.rcParams.update({"font.family": "sans-serif",
                                "font.sans-serif": [found]})


def save_figure(fig, out, suffix, formats, dpi):
    for fmt in formats:
        path = f"{out}_{suffix}.{fmt}"
        fig.savefig(path, dpi=dpi, facecolor="white")
        print(f"[INFO] Figure -> {path}")
    plt.close(fig)


def dumbbell_table(df, biotypes):
    """Rows and gene order shared by the dumbbell and extension-count plots.

    Returns (keep, glabels, means), or None if there is nothing to plot.
    `means` is indexed by gene in plotting order: y = 0 at the bottom, largest
    genotype difference at the top. Both plots draw from this, so their rows
    line up when the figures are placed side by side.
    """
    keep = (df[df["biotype"].isin(biotypes)]
            .dropna(subset=["extension_fraction"]))
    if keep.empty:
        return None
    glabels = list(dict.fromkeys(keep["genotype_label"]))
    means = (keep.groupby(["label", "genotype_label"])["extension_fraction"]
                 .mean().unstack())
    # Genes measured in only one genotype have no difference to sort on.
    means = means.dropna(how="any")
    if means.empty or len(glabels) < 2:
        return None
    means["_diff"] = means[glabels[-1]] - means[glabels[0]]
    means = means.sort_values("_diff")          # largest ends up at the top
    # Replicates of genes dropped above must not reach the plots either.
    keep = keep[keep["label"].isin(means.index)]
    return keep, glabels, means


def plot_dumbbell(df, gcolor, width, height, biotypes, show_reps=True):
    """Dumbbell plot: each gene one row, genotype means joined by a line.

    Genes are ordered by the DIFFERENCE between genotype means, largest at top.

    Individual replicates are drawn as small open markers beside each mean.

    Returns (fig, [ax]) unsaved; main() aligns its rows with the
    extension-count figure before saving.
    """
    tbl = dumbbell_table(df, biotypes)
    if tbl is None:
        print("[WARN] Need two genotypes with shared genes passing "
              "--min-reads; skipping dumbbell plot.")
        return
    keep, glabels, means = tbl
    labels = means.index.tolist()
    y = np.arange(len(labels))

    fig, ax = plt.subplots(figsize=(width, height), layout="constrained")
    mean_size = 26

    # Connector first, so markers sit on top of it.
    ax.hlines(y, means[glabels[0]], means[glabels[-1]],
              color="#9E9E9E", linewidth=0.8, zorder=1)

    markers = ["o", "s", "^", "D"]
    for i, g in enumerate(glabels):
        ax.scatter(means[g], y, s=mean_size, marker=markers[i % len(markers)],
                   color=gcolor[g], label=g, zorder=5 - i,
                   edgecolor="none")

        if show_reps:
            sub = keep[keep.genotype_label == g]
            for yi, lab in enumerate(labels):
                v = sub.loc[sub.label == lab, "extension_fraction"].to_numpy()
                if not len(v):
                    continue
                # Offset above or below the mean marker by genotype, so the
                # two genotypes' replicate points stay separate.
                off = 0.22 if i % 2 else -0.22
                ax.scatter(v, np.full(len(v), yi + off), s=5,
                           facecolor="white", edgecolor=gcolor[g],
                           linewidth=0.5, zorder=2)

    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xlabel(LABELS["extension_pct"])
    # Blank title: reserves the same top margin the extension-count plot uses
    # for its genotype names, so rows line up when the two sit side by side.
    ax.set_title(" ", fontsize=matplotlib.rcParams["axes.labelsize"])
    # Legend before the draw below, so the constrained layout has already made
    # room for it when the axes width is measured for the padding.
    fig.legend(*ax.get_legend_handles_labels(), loc="outside lower center",
               ncol=len(glabels), frameon=False,
               handlelength=1.0, columnspacing=1.4)

    # Percentage ticks at 10-point steps, with the axis rounded UP to the next
    # step above the largest value plotted - replicate markers included, since
    # those sit outside the means. Matplotlib's default locator stops at the
    # last tick below the data, which clipped the top of the range.
    hi = max(float(keep["extension_fraction"].max()),
             float(means[glabels].to_numpy().max()))
    upper = np.ceil(hi * 10) / 10


    fig.canvas.draw()
    ax_w_in = ax.get_position().width * fig.get_size_inches()[0]
    radius_in = (np.sqrt(mean_size) / 2) / 72.0        # points -> inches
    pad = (radius_in / ax_w_in) * upper * 1.05         # 5% margin
    ax.set_xlim(-pad, upper + pad)
    ax.set_xticks(np.arange(0, upper + 1e-9, 0.1))
    ax.xaxis.set_major_formatter(
        FuncFormatter(lambda v, _: f"{v * 100:.0f}%"))

    ax.set_ylim(-0.8, len(labels) - 0.2)
    for xt in ax.get_xticks():
        if abs(xt) < 1e-9:
            continue
        ax.axvline(xt, color="#E0E0E0", linewidth=0.5, zorder=0)
    ax.set_axisbelow(True)

    # Left spine at the padded axis limit, so it sits just clear of 0% and
    # markers at zero are not clipped by it.
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines["left"].set_visible(True)
    # One tick per gene on the left spine, matching the x-axis tick style.
    ax.tick_params(axis="y", length=2.5, width=0.6, direction="out")

    return fig, [ax]


def plot_extension_counts(df, gcolor, out, width, height, biotypes):
    """Companion to the dumbbell: reads with vs without a 5' extension.

    One panel per genotype, one row per gene in the dumbbell's order and with
    its y limits, so the two figures can sit side by side with rows aligned.
    Each bar is stacked: reads with an extension (genotype colour) from zero,
    then reads without one (grey). Counts are summed over the same replicates
    that enter the dumbbell - those passing --min-reads for that gene - so each
    bar is the number of reads behind the corresponding dumbbell point.
    Panels share the x axis, so bar lengths compare directly across genotypes.

    Writes <out>_extension_counts.tsv and returns (fig, axes) unsaved; main()
    aligns rows with the dumbbell, adds the shared x label, then saves.
    """
    tbl = dumbbell_table(df, biotypes)
    if tbl is None:
        print("[WARN] Need two genotypes with shared genes passing "
              "--min-reads; skipping extension-count plot.")
        return
    keep, glabels, means = tbl
    labels = means.index.tolist()
    y = np.arange(len(labels))

    sums = (keep.groupby(["label", "genotype_label"])[["n_extension",
                                                       "n_trimmed"]]
                .sum())

    # The values behind the plot, in its row order (top row first).
    tbl_out = (keep.groupby(["label", "genotype", "genotype_label"],
                            sort=False)
                   .agg(n_replicates=("sample_id", "size"),
                        samples=("sample_id", lambda s: ",".join(s)),
                        reads_with_extension=("n_extension", "sum"),
                        reads_without_extension=("n_trimmed", "sum"))
                   .reset_index())
    tbl_out["total_reads"] = (tbl_out["reads_with_extension"]
                              + tbl_out["reads_without_extension"])
    tbl_out["plot_row"] = tbl_out["label"].map(
        {lab: len(labels) - i for i, lab in enumerate(labels)})
    tbl_out["_g"] = tbl_out["genotype_label"].map(
        {g: i for i, g in enumerate(glabels)})
    tbl_out = (tbl_out.sort_values(["plot_row", "_g"])
                      .drop(columns="_g"))
    tsv = f"{out}_extension_counts.tsv"
    tbl_out.to_csv(tsv, sep="\t", index=False)
    print(f"[INFO] Extension counts -> {tsv}")


    hi = float(tbl_out["total_reads"].max()) or 1.0
    step = 25000 * max(1, int(np.ceil(hi / 25000 / 4)))
    xtop = step * np.ceil(hi / step)
    xticks = np.arange(0, xtop + 1, step)

    fig, axes = plt.subplots(1, len(glabels), figsize=(width, height),
                             sharey=True, sharex=True, layout="constrained",
                             squeeze=False)
    # Room between panels so the last tick label of one does not run into the
    # "0" of the next.
    fig.get_layout_engine().set(wspace=0.12)
    axes = axes[0]
    for ax, g in zip(axes, glabels):
        sub = sums.xs(g, level="genotype_label").reindex(labels).fillna(0)
        ext = sub["n_extension"].to_numpy(dtype=float)
        noext = sub["n_trimmed"].to_numpy(dtype=float)
        ax.barh(y, ext, height=0.7, color=gcolor[g], linewidth=0, zorder=2)
        ax.barh(y, noext, left=ext, height=0.7, color=NO_EXTENSION_COLOR,
                linewidth=0, zorder=2)
        ax.set_title(g, fontsize=matplotlib.rcParams["axes.labelsize"])
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_xticks(xticks)
        # Thousands of reads: 0, 10k, 20k, ...
        ax.xaxis.set_major_formatter(
            FuncFormatter(lambda v, _: "0" if v == 0 else f"{v / 1000:g}k"))
        ax.set_xlim(0, xtop)
        ax.tick_params(axis="y", length=2.5, width=0.6, direction="out")
    for ax in axes[1:]:
        ax.tick_params(axis="y", left=False)

    axes[0].set_yticks(y)
    axes[0].set_yticklabels(labels)
    axes[0].set_ylim(-0.8, len(labels) - 0.2)    # same as the dumbbell
    # Row alignment with the dumbbell: both figures reserve the same space
    # around the axes
    for ax in axes:
        ax.set_xlabel(" ")

    # One legend entry for the extension segment showing every genotype's
    # colour, one for the shared grey.
    ext_handle = tuple(mpatches.Patch(color=gcolor[g]) for g in glabels)
    no_handle = mpatches.Patch(color=NO_EXTENSION_COLOR)
    fig.legend([ext_handle, no_handle],
               [LABELS["with_extension"], LABELS["no_extension"]],
               handler_map={tuple: HandlerTuple(ndivide=None, pad=0)},
               loc="outside lower center", ncol=2, frameon=False,
               handlelength=1.6, columnspacing=1.4)

    return fig, list(axes)


def freeze_layout(fig):
    """Run the constrained layout once, then fix every element in place."""
    fig.canvas.draw()
    fig.set_layout_engine("none")


def place_counts_xlabel(fig, axes):
    """Draw the shared x label on the line the blank per-panel labels hold."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    inv = fig.transFigure.inverted()
    lb = axes[0].xaxis.label.get_window_extent(r).transformed(inv)
    x0 = axes[0].get_position().x0
    x1 = axes[-1].get_position().x1
    fig.text((x0 + x1) / 2, (lb.y0 + lb.y1) / 2, LABELS["counts_x"],
             ha="center", va="center",
             fontsize=matplotlib.rcParams["axes.labelsize"])


def align_rows(fig_a, axes_a, fig_b, axes_b):
    ha, hb = fig_a.get_size_inches()[1], fig_b.get_size_inches()[1]
    freeze_layout(fig_a)
    freeze_layout(fig_b)
    if abs(ha - hb) > 1e-6:
        print("[WARN] Dumbbell and extension-count figures differ in height; "
              "rows will not line up side by side.")
        return
    pa, pb = axes_a[0].get_position(), axes_b[0].get_position()
    bottom = max(pa.y0, pb.y0)
    top = min(pa.y1, pb.y1)
    for ax in list(axes_a) + list(axes_b):
        p = ax.get_position()
        ax.set_position([p.x0, bottom, p.width, top - bottom])


def plot_per_sample(sample_df, gcolor, out, formats, dpi, width, height):
    glabels = list(dict.fromkeys(sample_df["genotype_label"]))
    fig, ax = plt.subplots(figsize=(width, height), layout="constrained")

    for i, g in enumerate(glabels):
        sub = sample_df[sample_df.genotype_label == g]
        ax.bar(i, sub["extension_fraction"].mean(), 0.55,
               yerr=sub["extension_fraction"].std(), capsize=2.5,
               color=gcolor[g], error_kw=dict(lw=0.6, ecolor="#444444"),
               zorder=2)
        # Replicates drawn individually: with n=3 the reader should see the
        # actual spread rather than infer it from an error bar.
        jitter = np.linspace(-0.12, 0.12, len(sub))
        ax.scatter(i + jitter, sub["extension_fraction"],
                   s=9, color="white", edgecolor="#222222",
                   linewidth=0.5, zorder=3)

    ax.set_xticks(range(len(glabels)))
    ax.set_xticklabels(glabels)
    ax.set_ylabel(LABELS["extension_pct"])
    # Top of the axis rounded UP to the next 10% above the highest bar, error
    # bar or replicate point, so the range ends on a round value.
    stats = sample_df.groupby("genotype_label")["extension_fraction"]
    hi = max(float(sample_df["extension_fraction"].max()),
             float((stats.mean() + stats.std().fillna(0)).max()))
    ax.set_ylim(0, max(np.ceil(hi * 10 - 1e-9) / 10, 0.1))
    # Values are stored as fractions; displayed as whole percentages.
    ax.yaxis.set_major_formatter(PercentFormatter(xmax=1, decimals=0))
    ax.spines[["top", "right"]].set_visible(False)

    save_figure(fig, out, "per_sample", formats, dpi)


def plot_read_integrity(per_gene_qual, overall, gcolor, out, formats, dpi,
                        width, height, biotypes):
    """How many reads at each gene have intact ends, at three levels.

      intact_3_of_total        of all reads aligned to the gene, those whose 3'
                               end sits at its annotated 3' boundary. Direct RNA
                               sequences 3'->5', so this end is where sequencing
                               began and is accurate. A low value means real 3'
                               ends are not where the annotation puts them.

      intact_5_of_intact_3     of those, the ones that span far enough
                                to be considered full-length - they did not lose their
                               5' end to truncation. 

      intact_5_and_3_of_total  both ends intact, out of all reads aligned to
                               the gene. The product of the two above, and the
                               only reads whose 5' end can be interpreted.
    """
    keep = per_gene_qual[per_gene_qual["biotype"].isin(biotypes)]
    if keep.empty:
        print("[WARN] No genes in the selected biotypes; "
              "skipping read-integrity plot.")
        return

    glabels = list(dict.fromkeys(keep["genotype_label"]))
    labels = (keep[["label", "start"]].drop_duplicates("label")
                                      .sort_values("start")["label"].tolist())
    x = np.arange(len(labels))
    w = 0.8 / max(len(glabels), 1)

    panels = [
        ("mean_intact_3_of_total", "sd_intact_3_of_total",
         "Proportion of total\naligned reads with\nintact 3' end"),
        ("mean_intact_5_of_intact_3", "sd_intact_5_of_intact_3",
         "Proportion of 3'\nintact reads that\nhave 5' end intact"),
        ("mean_intact_5_and_3_of_total", "sd_intact_5_and_3_of_total",
         "Proportion of total\naligned reads with\nintact 3' and 5' end"),
    ]

    n_p = len(panels)
    fig = plt.figure(figsize=(width, height), layout="constrained")
    fig.get_layout_engine().set(hspace=0.01)
    gs = fig.add_gridspec(n_p + 2, 1,
                          height_ratios=[1] * n_p + [0.18, 1.1])
    axes = [fig.add_subplot(gs[i]) for i in range(n_p)]
    spacer = fig.add_subplot(gs[n_p])
    spacer.axis("off")
    axes.append(fig.add_subplot(gs[n_p + 1]))

    for j, (ax, (col, errcol, ylab)) in enumerate(zip(axes, panels)):
        for i, g in enumerate(glabels):
            sub = (keep[keep.genotype_label == g]
                   .set_index("label").reindex(labels))
            ax.bar(x + i * w - 0.4 + w / 2, sub[col], w,
                   yerr=(sub[errcol] if errcol else None), capsize=1.5,
                   color=gcolor[g], label=g,
                   error_kw=dict(lw=0.5, ecolor="#444444"))
        ax.set_xticks(x)
        ax.set_xlim(-0.7, len(labels) - 0.3)
        if j == len(panels) - 1:
            ax.set_xticklabels(labels, rotation=30, ha="right",
                               rotation_mode="anchor")
        else:
            ax.set_xticklabels([])
        ax.set_ylabel(ylab)
        ax.set_ylim(0, 1)
        ax.spines[["top", "right"]].set_visible(False)


    handles, lbls = axes[0].get_legend_handles_labels()
    spacer.legend(handles, lbls, loc="center", ncol=len(glabels),
                  frameon=False, handlelength=1.2, columnspacing=1.4)


    ax = axes[-1]
    for i, g in enumerate(glabels):
        sub = overall[overall.genotype_label == g]
        ax.bar(i, sub["intact_5_and_3_of_total"].mean(), 0.28,
               yerr=sub["intact_5_and_3_of_total"].std(), capsize=2.5,
               color=gcolor[g], error_kw=dict(lw=0.6, ecolor="#444444"),
               zorder=2)
        jitter = np.linspace(-0.07, 0.07, len(sub))
        ax.scatter(i + jitter, sub["intact_5_and_3_of_total"], s=9, color="white",
                   edgecolor="#222222", linewidth=0.5, zorder=3)
    ax.set_xticks(range(len(glabels)))
    ax.set_xticklabels(glabels)
    ax.set_xlim(-0.7, len(glabels) - 0.3)
    ax.set_ylabel("Proportion of total reads\nwith intact 3' and 5' end")
    ax.set_ylim(0, None)
    ax.spines[["top", "right"]].set_visible(False)

    save_figure(fig, out, "read_integrity", formats, dpi)


def test_genotypes(per_sample):
    """Welch's t-test on the whole-sample extension fraction.

    """
    from scipy import stats

    glabels = list(dict.fromkeys(per_sample["genotype_label"]))
    if len(glabels) < 2:
        print("[WARN] Only one genotype; no statistical test run.")
        return pd.DataFrame()

    geno_code = (per_sample.drop_duplicates("genotype_label")
                           .set_index("genotype_label")["genotype"])
    ref = glabels[0]
    f0 = per_sample.loc[per_sample.genotype_label == ref,
                        "extension_fraction"].to_numpy(dtype=float)
    rows = []
    for g in glabels[1:]:
        f1 = per_sample.loc[per_sample.genotype_label == g,
                            "extension_fraction"].to_numpy(dtype=float)
        if len(f0) < 2 or len(f1) < 2:
            print(f"[WARN] {g} vs {ref}: need at least two libraries per "
                  "genotype; test skipped.")
            continue
        n0, n1 = len(f0), len(f1)
        v0, v1 = f0.var(ddof=1) / n0, f1.var(ddof=1) / n1
        se = np.sqrt(v0 + v1)
        # Welch-Satterthwaite degrees of freedom
        df_w = (v0 + v1) ** 2 / (v0 ** 2 / (n0 - 1) + v1 ** 2 / (n1 - 1))
        diff = f1.mean() - f0.mean()
        res = stats.ttest_ind(f1, f0, equal_var=False)
        tcrit = stats.t.ppf(0.975, df_w)
        rows.append(dict(
            reference=geno_code[ref], genotype=geno_code[g],
            reference_label=ref, genotype_label=g,
            n_libraries_reference=n0, n_libraries=n1,
            mean_fraction_reference=f0.mean(), sd_reference=f0.std(ddof=1),
            mean_fraction=f1.mean(), sd=f1.std(ddof=1),
            difference=diff,
            diff_ci95_low=diff - tcrit * se,
            diff_ci95_high=diff + tcrit * se,
            welch_t=res.statistic, welch_df=df_w, welch_p=res.pvalue,
        ))
    return pd.DataFrame(rows)


def write_excel(path, tables):
    """Write every summary table to one workbook, a sheet each."""
    try:
        with pd.ExcelWriter(path, engine="openpyxl") as xl:
            for sheet, tbl in tables.items():
                tbl.to_excel(xl, sheet_name=sheet[:31], index=False)
        print(f"[INFO] Excel -> {path}")
    except ImportError:
        print("[WARN] openpyxl not installed; skipping Excel output. "
              "Install with: pip install openpyxl")


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--features", required=True)
    p.add_argument("--bam-dir", required=True)
    p.add_argument("--samples", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--bam-suffix", default="_mito.bam")
    p.add_argument("--contig", default="MtDNA")
    p.add_argument("--tolerance", type=int, default=10,
                   help="bp slack when matching a read end to a boundary "
                        "(default 10)")
    p.add_argument("--min-coverage", type=float, default=0.95,
                   help="Fraction of the gene a read must cover to be counted. "
                        "Removes 5'-truncated reads. Default 0.95.")
    p.add_argument("--min-reads", type=int, default=10,
                   help="Minimum assigned reads before a fraction is reported")
    p.add_argument("--min-alignment", type=int, default=10,
                   help="nt a read must align within a gene to count toward "
                        "that gene's aligned-read total (default 10)")
    p.add_argument("--biotypes", choices=list(BIOTYPE_SETS),
                   default="mrna",
                   help="Which features appear in the figures and aggregated "
                        "tables. 'mrna' (default) = protein-coding only; "
                        "'mrna+rrna' adds rRNA; 'all' also adds the D-loop, "
                        "NCRs and tRNAs. The raw per-gene table "
                        "always contains every feature regardless.")
    p.add_argument("--sample-scope", choices=list(BIOTYPE_SETS),
                   default=None,
                   help="Features pooled for the per-sample fraction. Defaults "
                        "to whatever --biotypes is set to.")
    # Figure sizes: width and height in inches, saved at exactly this size.
    for key in SIZES:
        flag = "--" + key.replace("_", "-") + "-size"
        w, h = SIZES[key]
        p.add_argument(flag, nargs=2, type=float, metavar=("W", "H"),
                       default=SIZES[key],
                       help=f"{key.replace('_', '-')} figure size in inches "
                            f"(default {w} {h})")
    p.add_argument("--no-replicate-points", action="store_true",
                   help="Draw genotype means only, without the individual "
                        "replicate markers")
    p.add_argument("--excel", action="store_true",
                   help="Also write all summary tables to one .xlsx workbook")
    p.add_argument("--formats", nargs="+", default=["png"],
                   help="Output formats (default png)")
    p.add_argument("--dpi", type=int, default=600)
    p.add_argument("--no-plot", action="store_true")
    args = p.parse_args()

    matplotlib.rcParams.update(RC)
    setup_fonts()

    features = pd.read_csv(args.features, sep="\t")
    samples = pd.read_csv(args.samples, sep="\t")
    for col in ("sample_id", "genotype"):
        if col not in samples.columns:
            sys.exit(f"[ERROR] Sample sheet missing column: {col}")

    # `genotype` is the analysis code (short, ASCII, safe as a model factor
    # level); `genotype_label` is what appears on figures.
    if "genotype_label" not in samples.columns:
        samples["genotype_label"] = samples["genotype"]
        print("[WARN] No genotype_label column; using genotype values as "
              "figure labels.")

    if "non_coding" not in set(features["biotype"]):
        print("[WARN] No non_coding features in the table. Rerun "
              "mtdna_features.py so the D-loop and NCRs are included.")

    biotypes = BIOTYPE_SETS[args.biotypes]
    scope_biotypes = BIOTYPE_SETS[args.sample_scope or args.biotypes]

    feats = build_context(features)
    print(f"[INFO] {len(feats)} features "
          f"({(feats.biotype == 'tRNA').sum()} tRNA, "
          f"{(feats.upstream_is_trna).sum()} genes with a tRNA-punctuated "
          f"5' boundary)")
    print(f"[INFO] tolerance={args.tolerance} bp, "
          f"min_coverage={args.min_coverage:.0%}, min_reads={args.min_reads}")
    print(f"[INFO] biotypes={args.biotypes} ({', '.join(biotypes)})")
    if set(biotypes) - {"protein_coding"} or set(scope_biotypes) - {"protein_coding"}:
        print("[WARN] Axis labels say 'mt-mRNA', but the plotted features "
              "include non-mRNA biotypes. Use --biotypes mrna (and leave "
              "--sample-scope unset) for the labels to be accurate.")

    report_depth(samples, args.bam_dir, args.bam_suffix, args.contig)

    df, sample_stats = collect(feats, args.bam_dir, args.bam_suffix,
                               args.contig, samples, args.tolerance,
                               args.min_coverage, args.min_reads,
                               args.min_alignment)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    gpath = f"{out}_per_gene.tsv"
    df.to_csv(gpath, sep="\t", index=False, float_format="%.5f")
    print(f"\n[INFO] Per-gene table -> {gpath}")

    coding = df[df["biotype"].isin(biotypes)]
    per_gene = (coding.dropna(subset=["extension_fraction"])
                      .groupby(["label", "genotype", "genotype_label"],
                               as_index=False)
                      .agg(n_rep=("extension_fraction", "size"),
                           mean_extension_fraction=("extension_fraction", "mean"),
                           median_extension_fraction=("extension_fraction", "median"),
                           sd=("extension_fraction", "std"),
                           total_assigned=("n_assigned", "sum")))
    pos = df[["label", "start"]].drop_duplicates("label")
    per_gene = (per_gene.merge(pos, on="label")
                        .sort_values(["start", "genotype"])
                        .drop(columns="start"))

    gsum = f"{out}_per_gene_summary.tsv"
    per_gene.to_csv(gsum, sep="\t", index=False, float_format="%.5f")

    print(f"\n=== Per-gene extension fraction "
          f"({args.biotypes}; mean over replicates) ===\n")
    print(per_gene.to_string(index=False, float_format=lambda v: f"{v:.4f}"))
    print(f"\n[INFO] Per-gene summary -> {gsum}")


    scope = df[df["biotype"].isin(scope_biotypes)]
    per_sample = (scope.groupby(["sample_id", "genotype", "genotype_label"],
                                as_index=False)
                       .agg(n_extension=("n_extension", "sum"),
                            n_trimmed=("n_trimmed", "sum")))
    per_sample["n_assigned"] = (per_sample["n_extension"]
                                + per_sample["n_trimmed"])
    per_sample["extension_fraction"] = (per_sample["n_extension"]
                                       / per_sample["n_assigned"])

    spath = f"{out}_per_sample.tsv"
    per_sample.to_csv(spath, sep="\t", index=False, float_format="%.5f")

    print(f"\n=== Per-sample extension fraction "
          f"({args.sample_scope or args.biotypes}) ===\n")
    print(per_sample.to_string(index=False, float_format=lambda v: f"{v:.4f}"))

    gmean = (per_sample.groupby(["genotype", "genotype_label"], as_index=False)
                       .agg(n_samples=("extension_fraction", "size"),
                            mean_extension_fraction=("extension_fraction", "mean"),
                            sd=("extension_fraction", "std"),
                            total_assigned=("n_assigned", "sum")))
    print("\n  Genotype means (of the per-sample values above):\n")
    print(gmean.to_string(index=False, float_format=lambda v: f"{v:.4f}"))

    gpath2 = f"{out}_genotype_summary.tsv"
    gmean.to_csv(gpath2, sep="\t", index=False, float_format="%.5f")
    print(f"\n[INFO] Per-sample -> {spath}")
    print(f"[INFO] Genotype summary -> {gpath2}")

    # ── Genotype test on the whole-sample fractions ─────────────────────────
    gtest = test_genotypes(per_sample)
    if not gtest.empty:
        tpath = f"{out}_genotype_test.tsv"
        gtest.to_csv(tpath, sep="\t", index=False, float_format="%.6g")
        print("\n=== Welch's t-test, whole-sample extension fraction "
              "(one value per library) ===\n")
        for r in gtest.itertuples():
            print(f"  {r.genotype_label} vs {r.reference_label} "
                  f"({r.n_libraries} vs {r.n_libraries_reference} libraries)")
            print(f"    mean {r.mean_fraction:.2%} vs "
                  f"{r.mean_fraction_reference:.2%}; difference "
                  f"{r.difference * 100:.2f} percentage points "
                  f"(95% CI {r.diff_ci95_low * 100:.2f} to "
                  f"{r.diff_ci95_high * 100:.2f})")
            print(f"    t = {r.welch_t:.2f}, df = {r.welch_df:.2f}, "
                  f"p = {r.welch_p:.3g}")
        print(f"\n[INFO] Genotype test -> {tpath}")

    # ── Read integrity: how many reads at each gene have intact ends 
    per_gene_qual = (df.groupby(["label", "biotype", "genotype",
                                 "genotype_label"], as_index=False)
                       .agg(n_rep=("intact_5_and_3_of_total", "size"),
                            mean_intact_3_of_total=("intact_3_of_total", "mean"),
                            sd_intact_3_of_total=("intact_3_of_total", "std"),
                            mean_intact_5_of_intact_3=("intact_5_of_intact_3", "mean"),
                            sd_intact_5_of_intact_3=("intact_5_of_intact_3", "std"),
                            mean_intact_5_and_3_of_total=("intact_5_and_3_of_total", "mean"),
                            sd_intact_5_and_3_of_total=("intact_5_and_3_of_total", "std"),
                            n_aligned_to_gene=("n_aligned_to_gene", "sum"),
                            n_intact_3=("n_anchored", "sum"),
                            n_intact_5_and_3=("n_assigned", "sum")))
    per_gene_qual = (per_gene_qual.merge(pos, on="label")
                                  .sort_values(["start", "genotype"]))

    qgpath = f"{out}_read_integrity_per_gene.tsv"
    per_gene_qual.drop(columns="start").to_csv(
        qgpath, sep="\t", index=False, float_format="%.5f")


    overall = sample_stats.copy()
    overall = overall.rename(columns={"n_qualified": "n_intact_5_and_3"})
    overall["intact_5_and_3_of_total"] = (overall["n_intact_5_and_3"]
                                          / overall["n_reads_total"])
    overall["pct_no_intact_3"] = (100 * overall["n_no_3p_match"]
                                  / overall["n_reads_total"])
    overall["pct_no_intact_5"] = (100 * overall["n_below_coverage"]
                                  / overall["n_reads_total"])
    qspath = f"{out}_read_integrity_per_sample.tsv"
    overall.to_csv(qspath, sep="\t", index=False, float_format="%.5f")

    qual_geno = (overall.groupby(["genotype", "genotype_label"], as_index=False)
                        .agg(n_samples=("intact_5_and_3_of_total", "size"),
                             mean_intact_5_and_3_of_total=("intact_5_and_3_of_total", "mean"),
                             sd=("intact_5_and_3_of_total", "std"),
                             total_reads=("n_reads_total", "sum"),
                             total_intact_5_and_3=("n_intact_5_and_3", "sum")))
    qgeno_path = f"{out}_read_integrity_genotype.tsv"
    qual_geno.to_csv(qgeno_path, sep="\t", index=False, float_format="%.5f")

    print("\n=== Reads with both ends intact, per library ===\n")
    print("  pct_no_intact_3  3' end not at any annotated gene boundary")
    print("  pct_no_intact_5  3' end matched, but the read was not full length\n")
    print(overall[["sample_id", "genotype", "n_reads_total",
                   "n_intact_5_and_3", "intact_5_and_3_of_total",
                   "pct_no_intact_3", "pct_no_intact_5"]]
          .to_string(index=False, float_format=lambda v: f"{v:.4f}"))
    print("\n  Genotype means:\n")
    print(qual_geno[["genotype", "n_samples", "mean_intact_5_and_3_of_total",
                     "sd", "total_reads", "total_intact_5_and_3"]]
          .to_string(index=False, float_format=lambda v: f"{v:.4f}"))
    print("\n  Denominator here is every read in the library, counted once.")
    print("  The per-gene rates use reads aligned to each gene, and a read")
    print("  spanning several features counts toward each - so those")
    print("  denominators sum to more than the library size.")
    print("  pct_no_intact_3 counts reads whose 3' end is not within tolerance")
    print("  of any gene boundary; they still align to a gene.")
    print(f"\n=== Read integrity per gene "
          f"({args.biotypes}, mean over replicates) ===\n")
    show_q = ["label", "genotype", "mean_intact_3_of_total",
              "mean_intact_5_of_intact_3", "mean_intact_5_and_3_of_total",
              "n_aligned_to_gene", "n_intact_3", "n_intact_5_and_3"]
    print(per_gene_qual[per_gene_qual["biotype"].isin(biotypes)][show_q]
          .to_string(index=False, float_format=lambda v: f"{v:.4f}"))
    print("\n  intact_3_of_total        3' end at the gene's annotated 3'")
    print("                           boundary, of all reads aligned to it")
    print("  intact_5_of_intact_3     also full-length (>= --min-coverage of")
    print("                           the gene), of those with an intact 3'")
    print("  intact_5_and_3_of_total  both ends intact, of all reads aligned")
    print("                           to the gene - the reads the extension")
    print("                           fraction is computed from")
    print("\n  A read aligned to several genes counts toward each one's total,")
    print("  so n_aligned_to_gene sums to more than the library size. These are")
    print("  per-gene rates; the library-level rate is reported separately.")

    print(f"\n[INFO] Read integrity per gene     -> {qgpath}")
    print(f"[INFO] Read integrity per sample   -> {qspath}")
    print(f"[INFO] Read integrity by genotype  -> {qgeno_path}")

    if args.excel:
        write_excel(f"{out}.xlsx", {
            "per_gene_raw": df,
            "extension_per_gene": per_gene,
            "extension_per_sample": per_sample,
            "extension_genotype": gmean,
            "extension_genotype_test": gtest,
            "integrity_per_gene": per_gene_qual.drop(columns="start"),
            "integrity_per_sample": overall,
            "integrity_genotype": qual_geno,
        })

    if not args.no_plot:
        glabels = list(dict.fromkeys(df["genotype_label"]))
        gcolor = {g: GENOTYPE_PALETTE[i % len(GENOTYPE_PALETTE)]
                  for i, g in enumerate(glabels)}
        # Dumbbell and extension counts are built together and aligned before
        # saving, so their gene rows sit at the same heights.
        db = plot_dumbbell(df, gcolor, *args.dumbbell_size, biotypes,
                           show_reps=not args.no_replicate_points)
        ec = plot_extension_counts(df, gcolor, out,
                                   *args.extension_counts_size, biotypes)
        if db and ec:
            align_rows(*db, *ec)
        else:
            for f in (db, ec):
                if f:
                    freeze_layout(f[0])
        if ec:
            place_counts_xlabel(*ec)
            save_figure(ec[0], out, "extension_counts", args.formats, args.dpi)
        if db:
            save_figure(db[0], out, "dumbbell", args.formats, args.dpi)
        plot_per_sample(per_sample, gcolor, out, args.formats, args.dpi,
                        *args.per_sample_size)
        plot_read_integrity(per_gene_qual, overall, gcolor, out, args.formats,
                            args.dpi, *args.read_integrity_size, biotypes)


if __name__ == "__main__":
    main()
