#!/usr/bin/env python3
"""
Title:         prosift_plot_utils.py
Project:       ProSIFT (PROtein Statistical Integration and Filtering Tool)
Author:        Reina Hastings (reinahastings13@gmail.com)
Created:       2026-03-27
Last Modified: 2026-03-27
Purpose:       Shared Plotly plotting utilities used by Module 02 (prenorm_qc.py)
               and Module 03 (normalize.py). Provides color palette helpers, a
               save_plot dispatcher (PNG + HTML), and the three diagnostic plots
               that appear in both pre- and post-normalization QC: intensity box
               plots, PCA scatter, and sample-to-sample correlation heatmap.
Inputs:        (library module -- imported, not invoked directly)
Outputs:       (library module -- no output files)
Usage:
  from prosift_plot_utils import (
      make_color_map, save_plot,
      plot_intensity_boxplots, compute_and_plot_pca, compute_and_plot_correlation,
  )
"""

import logging
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

# ============================================================
# COLOR PALETTE
# ============================================================

# D3 qualitative palette -- colorblind-friendly, 10 colors.
_PALETTE = px.colors.qualitative.D3


def make_color_map(groups: list[str]) -> dict[str, str]:
    """Assign a consistent D3 color to each unique group label."""
    unique_groups = sorted(set(groups))
    return {g: _PALETTE[i % len(_PALETTE)] for i, g in enumerate(unique_groups)}


# ============================================================
# PLOT I/O
# ============================================================

def save_plot(fig: go.Figure, base_path: Path) -> None:
    """
    Save a Plotly figure as both a static PNG (via kaleido) and a
    standalone interactive HTML file with embedded plotly.js.

    base_path must NOT have an extension -- the extension is appended as a
    string rather than using Path.with_suffix, because Path.with_suffix
    replaces only the last dot-separated token in the stem (e.g.,
    "run.intensity_boxplots" would lose "intensity_boxplots").
    """
    png_path  = Path(str(base_path) + ".png")
    html_path = Path(str(base_path) + ".html")

    fig.write_image(str(png_path))
    fig.write_html(str(html_path), full_html=True, include_plotlyjs=True)
    logging.info(f"  Saved: {png_path.name}, {html_path.name}")


# ============================================================
# SHARED DIAGNOSTIC PLOTS
# ============================================================

def plot_intensity_boxplots(
    log2_df: pd.DataFrame,
    summary_df: pd.DataFrame,
    color_map: dict[str, str],
    run_id: str,
    label: str = "Pre-Normalization",
) -> go.Figure:
    """
    Box plot of log2 intensity distributions. One box per sample, samples
    ordered by group. Boxes show median, IQR, 1.5xIQR whiskers, and outlier
    points. Used by Module 02 (pre-normalization) and Module 03 (post-normalization).

    Parameters
    ----------
    log2_df : pd.DataFrame
        Log2 abundance matrix; index = protein_id, columns = sample_ids.
    summary_df : pd.DataFrame
        Per-sample summary with at minimum 'sample_id' and 'group' columns.
    color_map : dict[str, str]
        {group_label: hex_color} from make_color_map().
    run_id : str
        Run identifier used in the plot title.
    label : str
        Phase label for the plot title, e.g. "Pre-Normalization" or
        "Post-Normalization".
    """
    sample_order = (
        summary_df.sort_values(["group", "sample_id"])["sample_id"].tolist()
    )

    fig = go.Figure()
    shown_groups: set[str] = set()
    for sid in sample_order:
        grp = str(summary_df.loc[summary_df["sample_id"] == sid, "group"].iloc[0])
        show_legend = grp not in shown_groups
        shown_groups.add(grp)
        values = log2_df[sid].dropna().values
        fig.add_trace(
            go.Box(
                y=values,
                name=sid,
                marker_color=color_map[grp],
                legendgroup=grp,
                legendgrouptitle_text=grp if show_legend else None,
                showlegend=show_legend,
                boxpoints="outliers",
                hoverinfo="y+name",
            )
        )

    fig.update_layout(
        title=f"Intensity Distributions ({label}): {run_id}",
        xaxis_title="Sample",
        yaxis_title="log2 Intensity",
        height=520,
        width=800,
        showlegend=True,
    )
    return fig


def compute_and_plot_pca(
    log2_df: pd.DataFrame,
    summary_df: pd.DataFrame,
    color_map: dict[str, str],
    run_id: str,
    label: str = "Pre-Normalization",
    note: str | None = None,
) -> tuple[pd.DataFrame, go.Figure]:
    """
    PCA on complete-case log2 data (proteins with any NaN dropped).
    Centers and scales each protein (correlation-matrix PCA) so all proteins
    contribute equally regardless of absolute abundance.

    Parameters
    ----------
    log2_df : pd.DataFrame
        Log2 abundance matrix; index = protein_id, columns = sample_ids.
    summary_df : pd.DataFrame
        Per-sample summary with at minimum 'sample_id' and 'group' columns.
    color_map : dict[str, str]
        {group_label: hex_color} from make_color_map().
    run_id : str
        Run identifier used in the plot title.
    label : str
        Phase label for the plot title, e.g. "Pre-Normalization".
    note : str or None
        Optional italic annotation text displayed below the plot axis. Pass
        None (default) to omit the annotation.

    Returns
    -------
    pca_df : pd.DataFrame
        Per-sample PCA coordinates (PC1-PC3) and variance explained.
    fig : go.Figure
        PC1 vs PC2 scatter plot.
    """
    sample_ids = log2_df.columns.tolist()
    group_of = {
        row["sample_id"]: str(row["group"])
        for _, row in summary_df.iterrows()
    }

    # --- Drop proteins with any missing value (complete cases only) ---
    complete_df = log2_df.dropna(axis=0, how="any")
    n_proteins_used = len(complete_df)
    n_proteins_total = len(log2_df)
    n_dropped = n_proteins_total - n_proteins_used
    logging.info(
        f"  PCA: using {n_proteins_used} complete-case proteins "
        f"({n_dropped} dropped due to missing values)"
    )

    if n_proteins_used < 2:
        raise ValueError(
            f"Only {n_proteins_used} complete-case proteins available for PCA. "
            "Too many missing values to compute PCA."
        )

    # --- Center and scale: rows=samples, columns=proteins ---
    X = complete_df.T.values  # shape: (n_samples, n_proteins)
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    # --- Fit PCA ---
    n_components = min(len(sample_ids), n_proteins_used)
    pca = PCA(n_components=n_components)
    coords = pca.fit_transform(X_scaled)  # shape: (n_samples, n_components)

    var_exp = pca.explained_variance_ratio_
    pc1_var = float(var_exp[0]) if len(var_exp) > 0 else np.nan
    pc2_var = float(var_exp[1]) if len(var_exp) > 1 else np.nan
    pc3_var = float(var_exp[2]) if len(var_exp) > 2 else np.nan

    groups = [group_of[s] for s in sample_ids]

    pca_df = pd.DataFrame({
        "sample_id":              sample_ids,
        "group":                  groups,
        "PC1":                    coords[:, 0],
        "PC2":                    coords[:, 1] if coords.shape[1] > 1 else np.nan,
        "PC3":                    coords[:, 2] if coords.shape[1] > 2 else np.nan,
        "variance_explained_PC1": pc1_var,
        "variance_explained_PC2": pc2_var,
        "variance_explained_PC3": pc3_var,
        "n_proteins_used":        n_proteins_used,
    })

    # --- PC1 vs PC2 scatter plot ---
    pc1_label = f"PC1 ({pc1_var * 100:.1f}%)"
    pc2_label = f"PC2 ({pc2_var * 100:.1f}%)"

    fig = go.Figure()
    shown_groups: set[str] = set()
    for grp, grp_df in pca_df.groupby("group"):
        grp = str(grp)
        show_legend = grp not in shown_groups
        shown_groups.add(grp)
        fig.add_trace(
            go.Scatter(
                x=grp_df["PC1"],
                y=grp_df["PC2"],
                mode="markers+text",
                text=grp_df["sample_id"],
                textposition="top center",
                textfont=dict(size=11),
                marker=dict(color=color_map[grp], size=14,
                            line=dict(width=1, color="white")),
                name=grp,
                legendgroup=grp,
                showlegend=show_legend,
                hovertemplate=(
                    "<b>%{text}</b><br>"
                    f"PC1: %{{x:.3f}}<br>"
                    f"PC2: %{{y:.3f}}<extra></extra>"
                ),
            )
        )

    layout_kwargs: dict = dict(
        title=f"PCA ({label}): {run_id}",
        xaxis_title=pc1_label,
        yaxis_title=pc2_label,
        height=580,
        width=700,
        legend_title_text="Group",
    )
    if note is not None:
        layout_kwargs["annotations"] = [
            dict(
                text=f"<i>{note}</i>",
                xref="paper", yref="paper",
                x=0.0, y=-0.18,
                showarrow=False,
                font=dict(size=11, color="gray"),
                align="left",
            )
        ]

    fig.update_layout(**layout_kwargs)
    return pca_df, fig


def compute_and_plot_correlation(
    log2_df: pd.DataFrame,
    summary_df: pd.DataFrame,
    color_map: dict[str, str],
    run_id: str,
    label: str = "Pre-Normalization",
) -> tuple[pd.DataFrame, go.Figure]:
    """
    Pearson correlation matrix using pairwise-complete observations.
    Color scale anchored to [min_off_diagonal_corr, 1.0] to amplify subtle
    within-group vs between-group differences. Samples ordered by group.

    Parameters
    ----------
    log2_df : pd.DataFrame
        Log2 abundance matrix; index = protein_id, columns = sample_ids.
    summary_df : pd.DataFrame
        Per-sample summary with at minimum 'sample_id' and 'group' columns.
    color_map : dict[str, str]
        {group_label: hex_color} -- not used for the heatmap itself but
        accepted for API consistency.
    run_id : str
        Run identifier used in the plot title.
    label : str
        Phase label for the plot title, e.g. "Pre-Normalization".

    Returns
    -------
    corr_df : pd.DataFrame
        Square Pearson correlation matrix (sample_id index and columns).
    fig : go.Figure
        Annotated heatmap.
    """
    # Order samples: group together for clear block structure
    sample_order = (
        summary_df.sort_values(["group", "sample_id"])["sample_id"].tolist()
    )
    log2_ordered = log2_df[sample_order]

    # Pairwise-complete Pearson correlation (pandas default)
    corr_df = log2_ordered.corr(method="pearson")

    # Anchor color scale to the range of off-diagonal values
    mask = ~np.eye(len(sample_order), dtype=bool)
    min_corr = float(corr_df.values[mask].min())

    # Axis labels: "Group: SampleID"
    group_of = dict(zip(summary_df["sample_id"], summary_df["group"].astype(str)))
    axis_labels = [f"{group_of[s]}: {s}" for s in sample_order]

    # Cell annotations with 2dp correlation values
    annotations = []
    for i, row_sid in enumerate(sample_order):
        for j, col_sid in enumerate(sample_order):
            val = corr_df.loc[row_sid, col_sid]
            annotations.append(
                dict(
                    x=j, y=i,
                    text=f"{val:.2f}",
                    showarrow=False,
                    font=dict(size=10, color="black"),
                )
            )

    fig = go.Figure(
        data=go.Heatmap(
            z=corr_df.values,
            x=axis_labels,
            y=axis_labels,
            colorscale="RdBu",
            zmin=min_corr,
            zmax=1.0,
            colorbar=dict(title="Pearson r"),
            hovertemplate="%{y}<br>%{x}<br>r = %{z:.4f}<extra></extra>",
        )
    )
    fig.update_layout(
        title=f"Sample Correlation Heatmap ({label}): {run_id}",
        annotations=annotations,
        height=520,
        width=640,
        xaxis=dict(tickangle=40),
    )
    return corr_df, fig
