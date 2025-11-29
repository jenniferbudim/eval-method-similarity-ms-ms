#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogFormatterSciNotation


# ---------------- Helper Functions ---------------- #

def flatten_similarity_matrix(df):
    """Convert a square similarity matrix to a flat vector (upper triangle only)."""
    arr = df.to_numpy(dtype=float)
    tri_upper = arr[np.triu_indices_from(arr, k=1)]
    return tri_upper


def ensure_1d_valid(arr):
    """Ensure array is 1D, finite, and clipped to [0, 1]."""
    arr = np.asarray(arr, dtype=float)
    arr = np.clip(arr, 0.0, 1.0)
    arr = arr[~np.isnan(arr)]
    return arr


def prepare_logscale_counts(H_counts):
    """Prepare safe log normalization and handle zeros (clamped to ≥10^0)."""
    H = np.asarray(H_counts, dtype=float)
    H_for_color = H.copy()

    # Replace zeros and very small counts with 1 (10^0)
    H_for_color[H_for_color < 1] = 1.0

    # Clamp log scale between 10^0 and 10^6
    vmin = 1.0
    vmax = max(1e6, np.max(H_for_color))
    norm = LogNorm(vmin=vmin, vmax=vmax)

    # Colorbar ticks: 10^0 → 10^6
    cb_ticks = [10 ** e for e in range(0, 7)]
    return H_for_color, norm, cb_ticks


def plot_heatmap(H_counts, xedges, yedges, xbins, ybins, output_prefix, annotated=True):
    """Plot a viridis log-scale heatmap with optional log10(count) annotations."""
    H_for_color, norm, cb_ticks = prepare_logscale_counts(H_counts)

    fig, ax = plt.subplots(figsize=(10, 8))
    mesh = ax.pcolormesh(
        xedges, yedges, H_for_color.T,
        cmap="viridis", norm=norm, shading="auto"
    )

    # Colorbar with log ticks
    cbar = plt.colorbar(mesh, ax=ax, ticks=cb_ticks)
    cbar.ax.yaxis.set_major_formatter(LogFormatterSciNotation(base=10.0))
    cbar.set_label("Count (log scale)")

    # Axes labels and title
    ax.set_xticks(xbins)
    ax.set_yticks(ybins)
    ax.set_xticklabels([f"{x:.2f}" for x in xbins], rotation=45)
    ax.set_yticklabels([f"{y:.1f}" for y in ybins])

    # Annotate each bin with log10(count)
    if annotated:
        for i in range(len(xedges) - 1):
            for j in range(len(yedges) - 1):
                count = H_counts[i, j]
                # Values < 1 shown as 0.00 (dark purple)
                text = f"{np.log10(count):.2f}" if count > 0 else "0.00"
                ax.text(
                    (xedges[i] + xedges[i + 1]) / 2,
                    (yedges[j] + yedges[j + 1]) / 2,
                    text,
                    ha="center", va="center",
                    color="white", fontsize=10, fontweight="bold"
                )

    plt.tight_layout()

    # Save PNG in current directory (for Nextflow)
    out_file = f"{output_prefix}{'_annotated' if annotated else '_clean'}.png"
    plt.savefig(out_file, dpi=300)
    plt.close()
    print(f"Saved: {out_file}")
    return out_file


# ---------------- Main ---------------- #

def main():
    parser = argparse.ArgumentParser(description="Generate train/test similarity heatmaps (matrix format).")
    parser.add_argument("--train_test_similarities", required=True, help="Wide CSV matrix for train-test similarities.")
    parser.add_argument("--test_similarities", required=True, help="Wide CSV matrix for test-test similarities.")
    parser.add_argument("--test_csv", required=True, help="CSV file containing test metadata (not used).")
    args = parser.parse_args()

    # Load data
    train_test = pd.read_csv(args.train_test_similarities, index_col=0)
    test_pairs = pd.read_csv(args.test_similarities, index_col=0)

    # Flatten to 1D arrays
    train_values = ensure_1d_valid(flatten_similarity_matrix(train_test))
    test_values = ensure_1d_valid(flatten_similarity_matrix(test_pairs))

    # Fixed 13 train-test similarity bins (0.40 → 1.00, step 0.05)
    xbins = np.array([0.40, 0.45, 0.50, 0.55, 0.60, 0.65,
                      0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00])
    # Test pairwise similarity bins (0.0 → 1.0, step 0.1)
    ybins = np.arange(0.0, 1.0 + 0.1, 0.1)

    # 2D histogram (counts)
    H_counts, xedges, yedges = np.histogram2d(train_values, test_values, bins=(xbins, ybins))

    # Generate annotated and clean plots
    plot_heatmap(H_counts, xedges, yedges, xbins, ybins,
                 output_prefix="train_test_similarity_plot",
                 annotated=True)
    plot_heatmap(H_counts, xedges, yedges, xbins, ybins,
                 output_prefix="train_test_similarity_plot",
                 annotated=False)


if __name__ == "__main__":
    main()