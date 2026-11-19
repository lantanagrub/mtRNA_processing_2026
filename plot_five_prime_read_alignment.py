#!/usr/bin/env python3
"""
plot_five_prime_read_map.py
===========================
Read maps of mtDNA direct-RNA transcripts, per genotype, with 5' processing
state shown per read. Uses the same read classification as
five_prime_extension.py:

  1. Assign each read to a gene by its 3' end.
  2. Discard reads covering less than --min-coverage of that gene.
  3. Classify the 5' end:
       TRIMMED    begins at the gene's 5' boundary (+/- tolerance)
       EXTENSION  begins upstream of it

Within each genotype panel, reads are grouped by gene, with extension reads
on the lowest tracks and trimmed reads above them.

Manuscript figure: five_prime_read_map.png, Fig. 5D (default settings).

Usage
-----
    python plot_five_prime_read_map.py \
        --features $MTRNA_ROOT/ref/mtdna_features.tsv \
        --bam-dir  $MTRNA_ROOT/elac2/bam/mito \
        --samples  $MTRNA_SRC/config/samples/elac2.tsv \
        --out      $MTRNA_ROOT/elac2/results/processing/five_prime_read_map

PNG at 600 dpi, exactly --width x --height inches (default 6.6 x 3).
Liberation Sans: gene labels 6 pt, all other text 7 pt.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
import matplotlib.colors as mcolors

import pysam


# ── Palette ──────────────────────────────────────────────────────────────────

BIOTYPE_COLOR = {
    "protein_coding": "#DCDCF0",
    "tRNA":           "#FFD9A8",
    "rRNA":           "#D6EDD8",
    "non_coding":     "#EAEAEA",
}
BG_COLOR         = "#F5F5F5"
READ_COLOR       = "#1A1A1A"
UNASSIGNED_COLOR = "#BDBDBD"
DIVIDER_COLOR    = "#757575"    # feature boundary lines in read panels
DIVIDER_LW       = 0.2          # points

DEFAULT_WIDTH_IN = 6.6
DEFAULT_HEIGHT_IN = 3.0

# First family found is used for all text.
FONT_FAMILY = ["Liberation Sans", "Arial", "DejaVu Sans"]
GENE_LABEL_PT = 6      # gene names in the gene track
TEXT_PT = 7            # everything else

RC = {
    "font.size": TEXT_PT,
    "axes.labelsize": TEXT_PT,
    "axes.titlesize": TEXT_PT,
    "xtick.labelsize": TEXT_PT,
    "ytick.labelsize": TEXT_PT,
    "legend.fontsize": TEXT_PT,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "xtick.major.size": 2.5,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
}


def setup_fonts():
    """Use the first installed family in FONT_FAMILY."""
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


def _blend(hex_color, alpha, bg=(1.0, 1.0, 1.0)):
    """Pre-blend a colour against white (avoids dark seams where features abut)."""
    c = np.array(mcolors.to_rgb(hex_color))
    return tuple(c * alpha + np.array(bg) * (1 - alpha))


def shade_features(ax, features, alpha=0.55):
    """Feature shading for a read panel."""
    for f in features.itertuples():
        base = BIOTYPE_COLOR.get(f.biotype, BG_COLOR)
        ax.axvspan(f.start, f.end, color=_blend(base, alpha),
                   linewidth=0, zorder=1)


# ── Read classification ──────────────────────────────────────────────────────

def classify_reads(bam_path, contig, feats, tol, min_cov, keep_unassigned):
    """Assign reads to genes by 3' end and classify the 5' end."""
    if not Path(str(bam_path) + ".bai").exists():
        sys.exit(f"[ERROR] BAM index missing. Run: samtools index {bam_path}")

    ends    = feats["end"].to_numpy()
    starts  = feats["start"].to_numpy()
    lengths = (ends - starts + 1).astype(float)
    labels  = feats["label"].to_numpy()

    rows = []
    n_total = n_no3p = n_short = 0

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam.fetch(contig):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.reference_end is None:
                continue
            n_total += 1

            # 1-based inclusive, to match the feature table
            r_end = read.reference_end
            r_start = read.reference_start + 1

            diffs = np.abs(ends - r_end)
            k = int(np.argmin(diffs))

            if diffs[k] > tol:
                n_no3p += 1
                if keep_unassigned:
                    rows.append((r_start, r_end, "unassigned", "", 0))
                continue

            aligned_len = min(r_end, ends[k]) - max(r_start, starts[k]) + 1
            if aligned_len / lengths[k] < min_cov:
                n_short += 1
                if keep_unassigned:
                    rows.append((r_start, r_end, "unassigned", "", 0))
                continue

            oh = int(starts[k] - r_start)
            if oh > tol:
                rows.append((r_start, r_end, "extension", labels[k], oh))
            else:
                rows.append((r_start, r_end, "trimmed", labels[k], 0))

    df = pd.DataFrame(rows, columns=["start", "end", "klass", "gene", "extension"])

    stats = dict(total=n_total, no_3p_match=n_no3p, below_coverage=n_short,
                 assigned=int((df["klass"] != "unassigned").sum()))
    return df, stats


def load_all(bam_dir, suffix, contig, feats, samples, tol, min_cov,
             keep_unassigned):
    frames = []
    for s in samples.itertuples():
        path = Path(bam_dir) / f"{s.sample_id}{suffix}"
        if not path.exists():
            sys.exit(f"[ERROR] BAM not found: {path}")
        df, st = classify_reads(path, contig, feats, tol, min_cov,
                                keep_unassigned)
        df["sample_id"] = s.sample_id
        df["genotype"] = s.genotype
        n_oh = int((df["klass"] == "extension").sum())
        print(f"[INFO] {s.sample_id} ({s.genotype}): {st['total']:,} reads -> "
              f"{st['assigned']:,} assigned, {n_oh:,} with 5' extension "
              f"({st['no_3p_match']:,} no 3' match, "
              f"{st['below_coverage']:,} below {min_cov:.0%} coverage)")
        frames.append(df)
    return pd.concat(frames, ignore_index=True)


# ── Depth handling ───────────────────────────────────────────────────────────

def subsample(reads, mode, seed):
    """Equalise assigned-read depth across samples."""
    counts = reads.groupby("sample_id").size()

    print("\n[INFO] Reads available per sample after filtering:")
    for sid, n in counts.items():
        print(f"         {sid}: {n:>9,}")
    print(f"         range: {counts.min():,} - {counts.max():,} "
          f"({counts.max() / counts.min():.1f}x)")

    if mode == "none":
        print("[INFO] No subsampling.\n")
        return reads, None

    target = int(counts.min()) if mode == "min" else int(mode)
    if target > counts.min():
        print(f"[WARN] Requested depth {target:,} exceeds the shallowest "
              f"sample ({counts.min():,}); it will be used in full.")

    # Per group rather than groupby().apply(), which drops sample_id in pandas 3
    parts = [g.sample(n=min(target, len(g)), random_state=seed)
             for _, g in reads.groupby("sample_id")]
    out = pd.concat(parts, ignore_index=True)
    print(f"[INFO] Subsampled to {target:,} reads per sample "
          f"({len(out):,} total)\n")
    return out, target


# ── Panels ───────────────────────────────────────────────────────────────────

def pack_tracks(df, feature_order, min_gap):
    """Pack reads onto tracks, grouped by gene, extension reads first.

    First-fit over reads ordered by (gene position, extension first, start).
    """
    rank = {lab: i for i, lab in enumerate(feature_order)}
    key = pd.DataFrame({
        "gene_rank": df["gene"].map(rank).fillna(len(rank)).astype(int),
        "not_extension": df["klass"] != "extension",   # False sorts first
        "start": df["start"],
    })
    order = key.sort_values(["gene_rank", "not_extension", "start"]).index

    starts = df["start"].to_numpy()
    ends = df["end"].to_numpy()
    pos = {ix: i for i, ix in enumerate(df.index)}

    track_end = []                      # highest end placed on each track
    track_of = np.empty(len(df), dtype=int)

    for ix in order:
        i = pos[ix]
        placed = False
        for t in range(len(track_end)):
            if starts[i] > track_end[t] + min_gap:
                track_end[t] = max(track_end[t], ends[i])
                track_of[i] = t
                placed = True
                break
        if not placed:
            track_of[i] = len(track_end)
            track_end.append(ends[i])

    return track_of, len(track_end)


def draw_read_map(ax, sub, label, feature_order, region, n_show, seed,
                  min_gap_frac=0.002):
    """One genotype panel. n_show reads, split evenly across replicates."""
    lo, hi = region
    sub = sub[(sub["end"] > lo) & (sub["start"] < hi)]

    if len(sub):
        # Remainder goes to the first replicates (e.g. 400 over 3 -> 134, 133, 133)
        groups = list(sub.groupby("sample_id"))
        reps = len(groups)
        base, extra = divmod(n_show, reps)
        parts = [g.sample(n=min(base + (j < extra), len(g)), random_state=seed)
                 for j, (_, g) in enumerate(groups)]
        sub = pd.concat(parts, ignore_index=True)

    if not len(sub):
        ax.set_yticks([])
        ax.set_ylabel(label, fontsize=TEXT_PT)
        return 0

    min_gap = (hi - lo) * min_gap_frac
    track_of, n_tracks = pack_tracks(sub, feature_order, min_gap)

    starts = sub["start"].to_numpy()
    ends = sub["end"].to_numpy()
    is_un = (sub["klass"] == "unassigned").to_numpy()

    if is_un.any():
        ax.hlines(track_of[is_un], starts[is_un], ends[is_un],
                  color=UNASSIGNED_COLOR, linewidth=0.5, zorder=2)
    keep = ~is_un
    if keep.any():
        ax.hlines(track_of[keep], starts[keep], ends[keep],
                  color=READ_COLOR, linewidth=0.6, zorder=3)

    ax.set_ylim(-1, n_tracks)
    ax.set_yticks([])
    ax.set_ylabel(f"{label}\nn={len(sub):,}", fontsize=TEXT_PT)
    ax.tick_params(axis="x", bottom=False, labelbottom=False)
    ax.spines[["top", "right", "left"]].set_visible(False)
    return n_tracks


def draw_gene_track(ax, features, region, fig, label_size=GENE_LABEL_PT):
    """Gene track with labels inside the boxes.

    Each label is tried horizontal, then split onto two lines, then rotated;
    left unlabelled if none fit. tRNAs are not labelled.
    """
    lo, hi = region
    span = max(hi - lo, 1)

    ax_w_in = ax.get_position().width * fig.get_size_inches()[0]
    in_per_bp = ax_w_in / span

    renderer = fig.canvas.get_renderer()

    def text_extent_in(txt):
        """(width, height) of the label in inches."""
        t = ax.text(0, 0, txt, fontsize=label_size,
                    **({"linespacing": 1.0} if "\n" in txt else {}))
        bb = t.get_window_extent(renderer)
        t.remove()
        return bb.width / fig.dpi, bb.height / fig.dpi

    box_h_in = ax.get_position().height * fig.get_size_inches()[1]
    cap_height_in = label_size * 0.95 / 72.0
    margin = 0.92       # label must fit inside 92% of its box

    ax.axvspan(lo, hi, color=BG_COLOR, zorder=1)
    for sp in ax.spines.values():
        sp.set_linewidth(0.5)

    for f in features.itertuples():
        color = BIOTYPE_COLOR.get(f.biotype, BG_COLOR)
        ax.add_patch(plt.Rectangle(
            (f.start, 0), f.end - f.start, 1,
            facecolor=color, edgecolor="black", linewidth=0.3, zorder=3))

        if f.biotype == "tRNA":
            continue

        box_in = (f.end - f.start) * in_per_bp
        text = str(f.label)

        def fits(txt):
            w, h = text_extent_in(txt)
            return w <= box_in * margin and h <= box_h_in * margin

        # Two-line version splits at the space nearest the middle
        two_line = None
        if " " in text:
            words = text.split()
            k = min(range(1, len(words)),
                    key=lambda i: abs(len(" ".join(words[:i]))
                                      - len(" ".join(words[i:]))))
            two_line = " ".join(words[:k]) + "\n" + " ".join(words[k:])

        if fits(text):
            rotation = 0
        elif two_line and fits(two_line):
            text, rotation = two_line, 0
        elif box_in >= cap_height_in * 0.6:
            rotation = 90
        else:
            continue

        ax.text((f.start + f.end) / 2, 0.5, text,
                ha="center", va="center", fontsize=label_size,
                rotation=rotation, zorder=4,
                **({"linespacing": 1.0} if "\n" in text else {}))

    ax.set_ylim(0, 1)
    ax.set_yticks([])


def close_track_gap(axes):
    """Gene track and panel spacing settings."""
    reads, track = axes[:-1], axes[-1]
    pos = [ax.get_position() for ax in axes]
    top, bottom = pos[0].y1, pos[-1].y0
    gap = (pos[0].y0 - pos[1].y1) if len(reads) > 1 else 0.0
    ratio = pos[-1].height / pos[0].height
    n = len(reads)
    ph = (top - bottom - gap * (n - 1)) / (n + ratio)
    y = top
    for ax, p in zip(reads, pos):
        y -= ph
        ax.set_position([p.x0, y, p.width, ph])
        y -= gap
    y += gap
    p = pos[-1]
    track.set_position([p.x0, bottom, p.width, y - bottom])


def add_boundary_dividers(ax, features, region, max_gap=20,
                          linewidth=DIVIDER_LW):
    """Dashed line at each feature boundary in a read panel.

    Features up to max_gap bp apart share one line at the gap midpoint.
    """
    lo, hi = region
    f = features.sort_values("start").reset_index(drop=True)
    xs = []
    for i in range(len(f) - 1):
        a_end, b_start = f.loc[i, "end"], f.loc[i + 1, "start"]
        if b_start - a_end - 1 <= max_gap:
            xs.append((a_end + b_start) / 2.0)
        else:
            xs += [a_end + 0.5, b_start - 0.5]
    for x in xs:
        if lo < x < hi:
            ax.axvline(x, color=DIVIDER_COLOR, linewidth=linewidth,
                       linestyle=(0, (3, 2)), zorder=2)


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--features", required=True)
    p.add_argument("--bam-dir", required=True)
    p.add_argument("--samples", required=True)
    p.add_argument("--out", required=True, help="Output prefix (no file suffix)")
    p.add_argument("--bam-suffix", default="_mito.bam")
    p.add_argument("--contig", default="MtDNA")
    p.add_argument("--genome-length", type=int, default=13794)
    p.add_argument("--region", nargs=2, type=int, default=None,
                   metavar=("START", "END"))
    p.add_argument("--tolerance", type=int, default=10,
                   help="bp slack matching a read end to a boundary (default 10)")
    p.add_argument("--min-coverage", type=float, default=0.95,
                   help="Fraction of the assigned gene a read must cover "
                        "(default 0.95)")
    p.add_argument("--subsample", default="min",
                   help="'min' (default), 'none', or an integer read count")
    p.add_argument("--n-reads", type=int, default=300,
                   help="Reads drawn per genotype panel, split evenly across "
                        "replicates (default 300)")
    p.add_argument("--show-unassigned", action="store_true",
                   help="Also draw reads that failed 3' anchoring or the "
                        "coverage filter, in grey")
    p.add_argument("--width", type=float, default=DEFAULT_WIDTH_IN,
                   help=f"Figure width in inches (default {DEFAULT_WIDTH_IN})")
    p.add_argument("--height", type=float, default=DEFAULT_HEIGHT_IN,
                   help=f"Figure height in inches (default {DEFAULT_HEIGHT_IN})")
    p.add_argument("--track-ratio", type=float, default=0.36,
                   help="Gene track height relative to one read panel "
                        "(default 0.36)")
    p.add_argument("--title", default=None, help="Optional figure title")
    p.add_argument("--no-legend", action="store_true")
    p.add_argument("--formats", nargs="+", default=["png"])
    p.add_argument("--dpi", type=int, default=600)
    p.add_argument("--seed", type=int, default=21)
    args = p.parse_args()

    matplotlib.rcParams.update(RC)
    setup_fonts()

    features = pd.read_csv(args.features, sep="\t").sort_values("start")
    samples = pd.read_csv(args.samples, sep="\t")
    for col in ("sample_id", "genotype"):
        if col not in samples.columns:
            sys.exit(f"[ERROR] Sample sheet missing column: {col}")

    # `genotype` is the analysis code; `genotype_label` appears on the figure.
    if "genotype_label" not in samples.columns:
        samples["genotype_label"] = samples["genotype"]
        print("[WARN] No genotype_label column in the sample sheet; "
              "using genotype values verbatim as figure labels.")

    if "non_coding" not in set(features["biotype"]):
        print("[WARN] No non_coding features present. Rerun mtdna_features.py "
              "so the D-loop and NCRs are included.")

    print(f"[INFO] {len(features)} features | tolerance={args.tolerance} bp, "
          f"min_coverage={args.min_coverage:.0%}\n")

    reads = load_all(args.bam_dir, args.bam_suffix, args.contig, features,
                     samples, args.tolerance, args.min_coverage,
                     args.show_unassigned)
    reads, depth = subsample(reads, args.subsample, args.seed)

    region = tuple(args.region) if args.region else (0, args.genome_length)
    feature_order = features["label"].tolist()
    genotypes = list(dict.fromkeys(samples["genotype"]))
    display = (samples.drop_duplicates("genotype")
                      .set_index("genotype")["genotype_label"].to_dict())

    n_geno = len(genotypes)
    fig, axes = plt.subplots(
        n_geno + 1, 1,
        figsize=(args.width, args.height),
        gridspec_kw={"height_ratios": [1] * n_geno + [args.track_ratio]},
        sharex=True, layout="constrained")
    fig.get_layout_engine().set(hspace=0.01, h_pad=0.02)
    fig.patch.set_facecolor("white")

    for i, g in enumerate(genotypes):
        ax = axes[i]
        shade_features(ax, features)
        n_tracks = draw_read_map(ax, reads[reads.genotype == g],
                                 display.get(g, g), feature_order, region,
                                 args.n_reads, args.seed)
        print(f"[INFO] {display.get(g, g)}: packed into {n_tracks} track(s)")
        add_boundary_dividers(ax, features, region)

    axes[-1].set_xlim(*region)
    axes[-1].xaxis.set_major_locator(MaxNLocator(nbins=7, integer=True))
    axes[-1].set_xlabel("Genome position (bp)")
    axes[-1].set_yticks([])

    if not args.no_legend:
        legend = [
            Line2D([0], [0], color=READ_COLOR, lw=1.4, label="read"),
            plt.Rectangle((0, 0), 1, 1, fc=BIOTYPE_COLOR["protein_coding"],
                          ec="#999999", lw=0.3, label="protein coding"),
            plt.Rectangle((0, 0), 1, 1, fc=BIOTYPE_COLOR["tRNA"],
                          ec="#999999", lw=0.3, label="tRNA"),
            plt.Rectangle((0, 0), 1, 1, fc=BIOTYPE_COLOR["rRNA"],
                          ec="#999999", lw=0.3, label="rRNA"),
            plt.Rectangle((0, 0), 1, 1, fc=BIOTYPE_COLOR["non_coding"],
                          ec="#999999", lw=0.3, label="non-coding"),
        ]
        if args.show_unassigned:
            legend.append(Line2D([0], [0], color=UNASSIGNED_COLOR, lw=1.4,
                                 label="unassigned"))
        fig.legend(handles=legend, loc="outside lower center",
                   ncol=len(legend), frameon=False,
                   handlelength=1.2, columnspacing=1.0, handletextpad=0.4)

    if args.title:
        fig.suptitle(args.title, fontsize=TEXT_PT, fontweight="bold")

    # Freeze layout first so gene labels are fitted to the final axes width
    fig.canvas.draw()
    fig.set_layout_engine("none")
    close_track_gap(axes)
    draw_gene_track(axes[-1], features, region, fig)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    for fmt in args.formats:
        path = f"{out}.{fmt}"
        fig.savefig(path, dpi=args.dpi, facecolor="white")
        print(f"[INFO] Figure -> {path} "
              f"({args.width:g} x {args.height:g} in, {args.dpi} dpi)")
    plt.close(fig)


if __name__ == "__main__":
    main()
