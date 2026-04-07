#!/usr/bin/env python3
"""
Title:         differential_abundance.py
Project:       ProSIFT (PROtein Statistical Integration and Filtering Tool)
Author:        Reina Hastings (reinahastings13@gmail.com)
Created:       2026-03-30
Last Modified: 2026-03-30
Purpose:       Module 04 DIFFERENTIAL_ABUNDANCE process. Parses and validates
               contrasts from params.yml, fits a linear model per protein using
               limma empirical Bayes via rpy2, applies DEqMS peptide-count-aware
               variance correction, and produces a per-protein results table (14
               columns), a plain-text summary, and volcano / MA diagnostic plots
               (static PNG + interactive HTML) for each contrast.
Inputs:
  --matrix       {run_id}.imputed_matrix.parquet    (Module 03 IMPUTE)
  --metadata     {run_id}.validated_metadata.parquet (Module 01 VALIDATE_INPUTS)
  --id-mapping   {run_id}.id_mapping.parquet         (Module 01 UNIPROT_MAPPING)
  --params       {run_id}_params.yml
Outputs:
  {run_id}.diff_abundance_results.parquet
  {run_id}.diff_abundance_results.csv
  {run_id}.diff_abundance_summary.txt
  {run_id}.{contrast}.volcano_plot.png/.html  (one pair per contrast)
  {run_id}.{contrast}.ma_plot.png/.html       (one pair per contrast)
Usage:
  differential_abundance.py \
    --matrix     CTXcyto_WT_vs_CTXcyto_KO.imputed_matrix.parquet \
    --metadata   CTXcyto_WT_vs_CTXcyto_KO.validated_metadata.parquet \
    --id-mapping CTXcyto_WT_vs_CTXcyto_KO.id_mapping.parquet \
    --params     CTXcyto_WT_vs_CTXcyto_KO_params.yml \
    --run-id     CTXcyto_WT_vs_CTXcyto_KO \
    --outdir     .
"""

import argparse
import datetime
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import yaml

try:
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.conversion import localconverter
    from rpy2.robjects.packages import importr
    _HAVE_RPY2 = True
except ImportError:
    _HAVE_RPY2 = False

from prosift_plot_utils import save_plot


# ============================================================
# CONSTANTS
# ============================================================

# Colors for volcano / MA plots -- fixed per direction, not group-based
_COLOR_UP = "#d62728"  # D3 red: significant, positive FC
_COLOR_DN = "#1f77b4"  # D3 blue: significant, negative FC
_COLOR_NS = "#cccccc"  # gray: not significant

_PEPTIDE_PREFIX = "peptide_count_"


# ============================================================
# ARGUMENT PARSING
# ============================================================

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="differential_abundance.py",
        description=(
            "ProSIFT Module 04 DIFFERENTIAL_ABUNDANCE: "
            "limma + DEqMS differential abundance analysis"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Example:\n"
            "  differential_abundance.py \\\n"
            "    --matrix     run.imputed_matrix.parquet \\\n"
            "    --metadata   run.validated_metadata.parquet \\\n"
            "    --id-mapping run.id_mapping.parquet \\\n"
            "    --params     run_params.yml \\\n"
            "    --run-id     run \\\n"
            "    --outdir     ."
        ),
    )
    parser.add_argument("--matrix",     required=True,
                        help="Imputed abundance matrix (Parquet)")
    parser.add_argument("--metadata",   required=True,
                        help="Validated metadata (Parquet)")
    parser.add_argument("--id-mapping", required=True, dest="id_mapping",
                        help="ID mapping table (Parquet, from UNIPROT_MAPPING)")
    parser.add_argument("--params",     required=True,
                        help="Run params.yml")
    parser.add_argument("--run-id",     required=True, dest="run_id",
                        help="Run identifier (used as output file prefix)")
    parser.add_argument("--outdir",     required=True,
                        help="Output directory")
    return parser.parse_args()


def setup_logging() -> None:
    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )


# ============================================================
# DATA LOADING AND PREPARATION
# ============================================================

def load_params(params_path: Path) -> dict:
    with open(params_path) as f:
        return yaml.safe_load(f)


def extract_abundance_and_peptide_cols(
    matrix_df: pd.DataFrame,
    params: dict,
) -> tuple[pd.DataFrame, pd.DataFrame, list[str]]:
    """
    Split the imputed matrix into abundance and peptide count DataFrames.

    Returns:
        abund_df   -- index=protein_id, columns=bare sample IDs (prefix stripped)
        pep_df     -- index=protein_id, columns=peptide_count_{sample_id} originals
                      (empty DataFrame if no peptide count columns present)
        sample_ids -- ordered list of bare sample IDs matching abund_df columns
    """
    abundance_prefix = params["input"].get("abundance_prefix", "")

    all_cols     = matrix_df.columns.tolist()
    peptide_cols = [c for c in all_cols if c.startswith(_PEPTIDE_PREFIX)]
    abund_cols   = [
        c for c in all_cols
        if c != "protein_id" and not c.startswith(_PEPTIDE_PREFIX)
    ]

    if not abund_cols:
        raise ValueError("No abundance columns found in imputed matrix.")

    sample_ids = (
        [c[len(abundance_prefix):] for c in abund_cols]
        if abundance_prefix
        else list(abund_cols)
    )

    abund_df = matrix_df.set_index("protein_id")[abund_cols].copy()
    abund_df.columns = sample_ids

    pep_df = (
        matrix_df.set_index("protein_id")[peptide_cols].copy()
        if peptide_cols
        else pd.DataFrame(index=matrix_df.set_index("protein_id").index)
    )

    return abund_df, pep_df, sample_ids


def build_group_map(metadata_df: pd.DataFrame, params: dict) -> dict[str, str]:
    """Return {sample_id: group_label} from validated metadata."""
    group_col = params["design"]["group_column"]
    if group_col not in metadata_df.columns:
        raise ValueError(
            f"design.group_column '{group_col}' not found in metadata. "
            f"Available columns: {metadata_df.columns.tolist()}"
        )
    return dict(
        zip(
            metadata_df["sample_id"].astype(str),
            metadata_df[group_col].astype(str),
        )
    )


# ============================================================
# CONTRAST PARSING AND VALIDATION
# ============================================================

def parse_and_validate_contrasts(
    params: dict,
    metadata_df: pd.DataFrame,
) -> list[tuple[str, str, str, str]]:
    """
    Parse and validate contrast strings from params.yml design.contrasts.

    Each contrast string must use the 'numerator_vs_denominator' format
    (split on first occurrence of '_vs_'). Both group names must match
    values found in the group_column of the metadata.

    Returns a list of 4-tuples:
        (contrast_user, numerator, denominator, r_contrast_str)
        e.g. ("KO_vs_WT", "KO", "WT", "KO - WT")

    Raises ValueError with a specific message on any validation failure.
    """
    group_col = params["design"]["group_column"]
    available_groups = sorted(
        metadata_df[group_col].dropna().astype(str).unique()
    )

    raw_contrasts = params.get("design", {}).get("contrasts", [])
    if not raw_contrasts:
        raise ValueError(
            "No contrasts defined in params.yml under design.contrasts. "
            "Add at least one contrast in 'numerator_vs_denominator' format "
            "(e.g., 'KO_vs_WT')."
        )

    parsed: list[tuple[str, str, str, str]] = []
    for contrast_user in raw_contrasts:
        # Split on first _vs_ to handle (unlikely) group names containing '_vs_'
        idx = contrast_user.find("_vs_")
        if idx == -1:
            raise ValueError(
                f"Contrast '{contrast_user}' does not contain '_vs_' delimiter. "
                "Use format 'numerator_vs_denominator' (e.g., 'KO_vs_WT')."
            )
        numerator   = contrast_user[:idx]
        denominator = contrast_user[idx + 4:]

        if not denominator:
            raise ValueError(
                f"Contrast '{contrast_user}': denominator is empty after "
                "splitting on '_vs_'."
            )

        for name, role in [(numerator, "numerator"), (denominator, "denominator")]:
            if name not in available_groups:
                raise ValueError(
                    f"Contrast '{contrast_user}': {role} group '{name}' not found "
                    f"in metadata column '{group_col}'. "
                    f"Available groups: {available_groups}"
                )

        r_contrast_str = f"{numerator} - {denominator}"
        parsed.append((contrast_user, numerator, denominator, r_contrast_str))

    return parsed


# ============================================================
# PEPTIDE COUNT SUMMARIZATION
# ============================================================

def summarize_peptide_counts(pep_df: pd.DataFrame) -> pd.Series:
    """
    Per-protein peptide count: minimum of nonzero values across all samples.

    Counts of 0 correspond to imputed positions (no DIA-NN detection); these
    are excluded from the minimum. Every protein is guaranteed to have at least
    one nonzero value because Module 01's detection filter removes ABSENT
    proteins.

    Returns pd.Series with index=protein_id, dtype int64.
    """
    pep_values = pep_df.values.astype(float)
    pep_values[pep_values == 0] = np.nan        # treat zeros as "no measurement"
    min_nonzero = np.nanmin(pep_values, axis=1)  # shape: (n_proteins,)

    return pd.Series(
        min_nonzero.astype(np.int64),
        index=pep_df.index,
        name="n_peptides",
    )


# ============================================================
# STATISTICAL ANALYSIS -- rpy2 bridge
# ============================================================

def _run_one_contrast_r(
    abund_df: pd.DataFrame,
    groups: list[str],
    unique_groups: list[str],
    pep_counts: "pd.Series | None",
    r_contrast_str: str,
    use_deqms: bool,
) -> tuple[pd.DataFrame, str]:
    """
    Run limma + DEqMS (or limma only) for one contrast via rpy2.

    All data is passed to R via the global environment with 'prosift_' prefix.
    R code is executed as string blocks. Results are converted back to pandas
    via pandas2ri.

    Parameters
    ----------
    abund_df      : index=protein_id, columns=sample_ids. No NaN.
    groups        : group label per sample, aligned to abund_df.columns.
    unique_groups : unique group levels (controls factor level order in design).
    pep_counts    : per-protein min-nonzero peptide count. None if limma-only.
    r_contrast_str: R-side contrast expression, e.g. "KO - WT".
    use_deqms     : True to run spectraCounteBayes after eBayes.

    Returns
    -------
    raw_df      : DataFrame with protein_id and raw R column names.
                  Columns: protein_id, logFC, AveExpr, t, P.Value, adj.P.Val, B
                  (plus sca.t, sca.P.Value, sca.adj.pval, count if DEqMS).
    method_used : "DEqMS" or "limma".
    """
    if not _HAVE_RPY2:
        raise RuntimeError(
            "rpy2 is not installed. Cannot run statistical analysis. "
            "Install rpy2 and ensure R (with limma and DEqMS) is accessible."
        )

    # --- Load R packages (fail fast with informative message) ---
    try:
        importr("limma")
        if use_deqms:
            importr("DEqMS")
    except Exception as exc:
        raise RuntimeError(
            f"Failed to load R packages: {exc}. "
            "Ensure limma and DEqMS are installed in the R library."
        ) from exc

    n_proteins = len(abund_df)
    n_samples  = len(abund_df.columns)

    # --- Pass data to R global environment ---
    # Flatten in row-major (C) order; R matrix built with byrow=TRUE below.
    mat_flat = abund_df.values.astype(float).flatten(order="C")

    ro.globalenv["prosift_mat_values"]    = ro.FloatVector(mat_flat.tolist())
    ro.globalenv["prosift_protein_ids"]   = ro.StrVector(abund_df.index.tolist())
    ro.globalenv["prosift_sample_ids"]    = ro.StrVector(list(abund_df.columns))
    ro.globalenv["prosift_groups"]        = ro.StrVector(groups)
    ro.globalenv["prosift_unique_groups"] = ro.StrVector(unique_groups)
    ro.globalenv["prosift_contrast_str"]  = ro.StrVector([r_contrast_str])

    # --- Build design matrix and fit linear model ---
    ro.r(f"""
        # Reconstruct protein x sample matrix (row-major values, byrow=TRUE)
        prosift_mat <- matrix(
            prosift_mat_values,
            nrow  = {n_proteins},
            ncol  = {n_samples},
            byrow = TRUE
        )
        rownames(prosift_mat) <- prosift_protein_ids
        colnames(prosift_mat) <- prosift_sample_ids

        # Means model: one coefficient per group, no intercept
        group_f <- factor(prosift_groups, levels = prosift_unique_groups)
        design  <- model.matrix(~ 0 + group_f)
        colnames(design) <- prosift_unique_groups

        # Fit -> contrast -> empirical Bayes moderation
        fit  <- limma::lmFit(prosift_mat, design)
        cmat <- limma::makeContrasts(
                    contrasts = prosift_contrast_str[1],
                    levels    = design
                )
        fit2 <- limma::contrasts.fit(fit, cmat)
        fit3 <- limma::eBayes(fit2)
    """)

    # --- DEqMS peptide count correction or limma-only ---
    if use_deqms and pep_counts is not None:
        # Align peptide counts to the protein order in abund_df
        pep_aligned = pep_counts.reindex(abund_df.index).fillna(1).astype(int)
        ro.globalenv["prosift_pep_counts"] = ro.IntVector(pep_aligned.tolist())

        ro.r("""
            fit3$count <- as.integer(prosift_pep_counts)
            fit4 <- DEqMS::spectraCounteBayes(fit3)
            prosift_results <- DEqMS::outputResult(fit4, coef_col = 1)
            # Add protein_id as a column so it survives the rpy2 conversion
            prosift_results$protein_id <- rownames(prosift_results)
        """)
        method_used = "DEqMS"

    else:
        ro.r("""
            # sort.by="none" preserves original protein order (rowname alignment)
            prosift_results <- limma::topTable(
                fit3,
                number  = Inf,
                sort.by = "none",
                coef    = 1
            )
            prosift_results$protein_id <- rownames(prosift_results)
        """)
        method_used = "limma"

    # --- Convert R data frame to pandas ---
    with localconverter(ro.default_converter + pandas2ri.converter):
        raw_df = ro.conversion.rpy2py(ro.globalenv["prosift_results"])
    raw_df = raw_df.reset_index(drop=True)

    return raw_df, method_used


# ============================================================
# RESULT ASSEMBLY
# ============================================================

def assemble_results(
    raw_df: pd.DataFrame,
    id_mapping_df: pd.DataFrame,
    params: dict,
    method_used: str,
    contrast_user: str,
) -> pd.DataFrame:
    """
    Map R column names to the ProSIFT output schema, add gene_symbol,
    compute significance and direction, and attach the contrast label.

    Output schema (14 columns):
        protein_id, gene_symbol, log2_fc, avg_abundance,
        limma_t, limma_pvalue, limma_adj_pvalue,
        deqms_t, deqms_pvalue, deqms_adj_pvalue,
        n_peptides, significant, direction, contrast
    """
    da_cfg      = params.get("differential_abundance", {})
    sig_cfg     = da_cfg.get("significance", {})
    fdr_thresh  = float(sig_cfg.get("fdr_threshold", 0.05))
    fc_thresh   = float(sig_cfg.get("fc_threshold", 1.0))

    # --- Rename R columns to ProSIFT schema names ---
    rename_map = {
        "logFC":       "log2_fc",
        "AveExpr":     "avg_abundance",
        "t":           "limma_t",
        "P.Value":     "limma_pvalue",
        "adj.P.Val":   "limma_adj_pvalue",
    }
    if method_used == "DEqMS":
        rename_map.update({
            "sca.t":       "deqms_t",
            "sca.P.Value": "deqms_pvalue",
            "sca.adj.pval":"deqms_adj_pvalue",
            "count":       "n_peptides",
        })

    out = raw_df.rename(columns=rename_map).copy()

    # Drop extra R columns (B statistic, etc.) not in the output schema
    keep = {"protein_id", "log2_fc", "avg_abundance",
            "limma_t", "limma_pvalue", "limma_adj_pvalue"}
    if method_used == "DEqMS":
        keep.update({"deqms_t", "deqms_pvalue", "deqms_adj_pvalue", "n_peptides"})
    out = out[[c for c in out.columns if c in keep]].copy()

    # --- Fill DEqMS columns with NA if limma-only ---
    for col in ("deqms_t", "deqms_pvalue", "deqms_adj_pvalue", "n_peptides"):
        if col not in out.columns:
            out[col] = pd.NA

    # --- Enforce schema dtypes ---
    # R returns integer counts as int32; cast to int64 per schema spec.
    if out["n_peptides"].notna().any():
        out["n_peptides"] = out["n_peptides"].astype("Int64")

    # --- Add gene symbol from ID mapping ---
    gene_col = None
    for candidate in ("gene_symbol_mouse", "gene_symbol"):
        if candidate in id_mapping_df.columns:
            gene_col = candidate
            break

    if gene_col is not None:
        sym_map = id_mapping_df.set_index('protein_id')[gene_col]
        out["gene_symbol"] = out["protein_id"].map(sym_map)
    else:
        logging.warning(
            "ID mapping table has no 'gene_symbol_mouse' or 'gene_symbol' column. "
            "gene_symbol will be null for all proteins."
        )
        out["gene_symbol"] = pd.NA

    # --- Determine primary p-value column for significance ---
    primary_adj_pval = (
        "deqms_adj_pvalue" if method_used == "DEqMS" else "limma_adj_pvalue"
    )

    # --- Significance call ---
    adj_pval = out[primary_adj_pval].astype(float)
    fc       = out["log2_fc"].astype(float)

    passes_fdr = adj_pval < fdr_thresh
    passes_fc  = (fc.abs() > fc_thresh) if fc_thresh > 0 else pd.Series(True, index=out.index)
    out["significant"] = passes_fdr & passes_fc

    # --- Direction ---
    def _direction(row: pd.Series) -> str:
        if not row["significant"]:
            return "ns"
        return "up" if row["log2_fc"] > 0 else "down"

    out["direction"] = out.apply(_direction, axis=1)

    # --- Contrast label ---
    out["contrast"] = contrast_user

    # --- Canonical column order ---
    col_order = [
        "protein_id", "gene_symbol",
        "log2_fc", "avg_abundance",
        "limma_t", "limma_pvalue", "limma_adj_pvalue",
        "deqms_t", "deqms_pvalue", "deqms_adj_pvalue",
        "n_peptides",
        "significant", "direction", "contrast",
    ]
    out = out[col_order]

    return out


# ============================================================
# SUMMARY TEXT
# ============================================================

def write_summary_txt(
    contrast_results: list[tuple[str, str, str, pd.DataFrame]],
    run_id: str,
    method_used: str,
    params: dict,
    outdir: Path,
) -> None:
    """
    Write a plain-text summary following the template in Section 4.6 of the
    Module 04 spec.

    contrast_results is a list of (contrast_user, numerator, denominator, df)
    for each contrast, in order.
    """
    now     = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    da_cfg  = params.get("differential_abundance", {})
    sig_cfg = da_cfg.get("significance", {})
    fdr_threshold = sig_cfg.get("fdr_threshold", 0.05)
    fc_threshold  = sig_cfg.get("fc_threshold", 1.0)

    primary_pval_col = "deqms_adj_pvalue" if method_used == "DEqMS" else "limma_adj_pvalue"
    primary_pval_label = "deqms_adj_pvalue" if method_used == "DEqMS" else "limma_adj_pvalue"

    lines = [
        "========================================",
        "DIFFERENTIAL ABUNDANCE SUMMARY",
        "========================================",
        "",
        f"Run:              {run_id}",
        f"Date:             {now}",
        f"Method:           {method_used}",
        "",
    ]

    for contrast_user, numerator, denominator, df in contrast_results:
        n_proteins   = len(df)
        n_numerator  = int((df["direction"] != "ns").any())  # placeholder; compute below
        group_col    = params["design"]["group_column"]

        # Count samples per group
        n_sig  = int(df["significant"].sum())
        n_up   = int((df["direction"] == "up").sum())
        n_down = int((df["direction"] == "down").sum())
        n_ns   = n_proteins - n_sig

        pct_sig  = f"{100 * n_sig / n_proteins:.1f}" if n_proteins > 0 else "0.0"
        pct_up   = f"{100 * n_up / n_proteins:.1f}"  if n_proteins > 0 else "0.0"
        pct_down = f"{100 * n_down / n_proteins:.1f}" if n_proteins > 0 else "0.0"

        lines += [
            "----------------------------------------",
            f"CONTRAST: {contrast_user}  ({numerator} - {denominator})",
            "----------------------------------------",
            "",
            "INPUT",
            "----------------------------------------",
            f"Proteins tested:          {n_proteins}",
            "",
        ]

        # Peptide count section (DEqMS only)
        if method_used == "DEqMS" and "n_peptides" in df.columns:
            pep = df["n_peptides"].dropna().astype(float)
            lines += [
                "PEPTIDE COUNTS (DEqMS)",
                "----------------------------------------",
                f"Summary method:           min of nonzero values per protein",
                f"Median count:             {pep.median():.1f}",
                f"Range:                    {int(pep.min())} - {int(pep.max())}",
                "",
            ]
        else:
            if method_used == "limma":
                lines += [
                    "NOTE: Method = limma (no peptide count correction applied).",
                    "",
                ]

        lines += [
            "SIGNIFICANCE THRESHOLDS",
            "----------------------------------------",
            f"FDR threshold:            {fdr_threshold} (BH-corrected)",
            f"log2 FC threshold:        {fc_threshold} (absolute)",
            f"Primary p-value:          {primary_pval_label}",
            "",
            "RESULTS",
            "----------------------------------------",
            f"Significant proteins:     {n_sig} / {n_proteins} ({pct_sig}%)",
            f"  Up-regulated:           {n_up} ({pct_up}%)",
            f"  Down-regulated:         {n_down} ({pct_down}%)",
            f"Not significant:          {n_ns}",
            "",
        ]

        # Top 10 by adjusted p-value
        top10 = (
            df.sort_values(primary_pval_col, ascending=True)
              .head(10)[["protein_id", "gene_symbol", "log2_fc", primary_pval_col]]
        )
        lines += [
            "TOP 10 BY SIGNIFICANCE",
            "----------------------------------------",
            f"  {'protein_id':<24} {'gene_symbol':<16} {'log2_fc':>8}  {'adj_pvalue':>12}",
        ]
        for _, row in top10.iterrows():
            gene = str(row["gene_symbol"]) if pd.notna(row["gene_symbol"]) else "NA"
            pval = row[primary_pval_col]
            pval_str = f"{pval:.3e}" if pd.notna(pval) else "NA"
            lines.append(
                f"  {str(row['protein_id']):<24} {gene:<16} "
                f"{row['log2_fc']:>8.3f}  {pval_str:>12}"
            )

        lines.append("")

    lines.append("========================================")

    path = outdir / f"{run_id}.diff_abundance_summary.txt"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    logging.info(f"  Saved: {path.name}")


# ============================================================
# DIAGNOSTIC PLOTS
# ============================================================

def plot_volcano(
    df: pd.DataFrame,
    contrast_user: str,
    run_id: str,
    params: dict,
) -> go.Figure:
    """
    Volcano plot: x = log2_fc, y = -log10(primary adj p-value).
    Points colored by direction (up/down/ns). Threshold lines drawn for
    both FDR and FC cutoffs. Hover text shows protein_id, gene_symbol,
    log2_fc, and adj_pvalue.
    """
    da_cfg      = params.get("differential_abundance", {})
    sig_cfg     = da_cfg.get("significance", {})
    fdr_thresh  = float(sig_cfg.get("fdr_threshold", 0.05))
    fc_thresh   = float(sig_cfg.get("fc_threshold", 1.0))

    # Primary adjusted p-value column
    primary_col = (
        "deqms_adj_pvalue"
        if "deqms_adj_pvalue" in df.columns and df["deqms_adj_pvalue"].notna().any()
        else "limma_adj_pvalue"
    )

    neg_log_p = -np.log10(df[primary_col].clip(lower=1e-300).astype(float))
    log2_fc   = df["log2_fc"].astype(float)

    color_map = {"up": _COLOR_UP, "down": _COLOR_DN, "ns": _COLOR_NS}
    traces: dict[str, dict] = {"up": {"x": [], "y": [], "text": []},
                                "down": {"x": [], "y": [], "text": []},
                                "ns":  {"x": [], "y": [], "text": []}}

    for _, row in df.iterrows():
        direction = str(row["direction"])
        gene      = str(row["gene_symbol"]) if pd.notna(row["gene_symbol"]) else "NA"
        pval      = row[primary_col]
        pval_str  = f"{pval:.3e}" if pd.notna(pval) else "NA"
        hover     = (
            f"{row['protein_id']}<br>"
            f"Gene: {gene}<br>"
            f"log2FC: {row['log2_fc']:.2f}<br>"
            f"adj.p: {pval_str}"
        )
        traces[direction]["x"].append(float(row["log2_fc"]))
        traces[direction]["y"].append(float(neg_log_p[row.name]))
        traces[direction]["text"].append(hover)

    label_map = {"up": "Up-regulated", "down": "Down-regulated", "ns": "Not significant"}
    fig = go.Figure()
    for direction in ("ns", "down", "up"):   # ns drawn first (background)
        d = traces[direction]
        if not d["x"]:
            continue
        opacity = 0.4 if direction == "ns" else 0.8
        fig.add_trace(
            go.Scatter(
                x=d["x"],
                y=d["y"],
                mode="markers",
                name=label_map[direction],
                marker=dict(color=color_map[direction], size=5, opacity=opacity),
                text=d["text"],
                hoverinfo="text",
            )
        )

    # Threshold lines
    y_fdr = -np.log10(fdr_thresh)
    fig.add_hline(y=y_fdr, line_dash="dash", line_color="#888888", line_width=1)
    if fc_thresh > 0:
        fig.add_vline(x=fc_thresh,  line_dash="dash", line_color="#888888", line_width=1)
        fig.add_vline(x=-fc_thresh, line_dash="dash", line_color="#888888", line_width=1)

    fig.update_layout(
        title=f"Volcano Plot: {contrast_user}  ({run_id})",
        xaxis_title="log2 Fold Change",
        yaxis_title=f"-log10({primary_col})",
        height=560,
        width=700,
        legend_title_text="Direction",
    )
    return fig


def plot_ma(
    df: pd.DataFrame,
    contrast_user: str,
    run_id: str,
    params: dict,
) -> go.Figure:
    """
    MA plot: x = avg_abundance (AveExpr), y = log2_fc.
    Same coloring scheme and threshold lines as the volcano plot.
    Hover text identical to volcano plot.
    """
    da_cfg      = params.get("differential_abundance", {})
    sig_cfg     = da_cfg.get("significance", {})
    fc_thresh   = float(sig_cfg.get("fc_threshold", 1.0))

    primary_col = (
        "deqms_adj_pvalue"
        if "deqms_adj_pvalue" in df.columns and df["deqms_adj_pvalue"].notna().any()
        else "limma_adj_pvalue"
    )

    color_map  = {"up": _COLOR_UP, "down": _COLOR_DN, "ns": _COLOR_NS}
    traces: dict[str, dict] = {"up": {"x": [], "y": [], "text": []},
                                "down": {"x": [], "y": [], "text": []},
                                "ns":  {"x": [], "y": [], "text": []}}

    for _, row in df.iterrows():
        direction = str(row["direction"])
        gene      = str(row["gene_symbol"]) if pd.notna(row["gene_symbol"]) else "NA"
        pval      = row[primary_col]
        pval_str  = f"{pval:.3e}" if pd.notna(pval) else "NA"
        hover     = (
            f"{row['protein_id']}<br>"
            f"Gene: {gene}<br>"
            f"log2FC: {row['log2_fc']:.2f}<br>"
            f"adj.p: {pval_str}"
        )
        traces[direction]["x"].append(float(row["avg_abundance"]))
        traces[direction]["y"].append(float(row["log2_fc"]))
        traces[direction]["text"].append(hover)

    label_map = {"up": "Up-regulated", "down": "Down-regulated", "ns": "Not significant"}
    fig = go.Figure()
    for direction in ("ns", "down", "up"):
        d = traces[direction]
        if not d["x"]:
            continue
        opacity = 0.4 if direction == "ns" else 0.8
        fig.add_trace(
            go.Scatter(
                x=d["x"],
                y=d["y"],
                mode="markers",
                name=label_map[direction],
                marker=dict(color=color_map[direction], size=5, opacity=opacity),
                text=d["text"],
                hoverinfo="text",
            )
        )

    if fc_thresh > 0:
        fig.add_hline(y=fc_thresh,  line_dash="dash", line_color="#888888", line_width=1)
        fig.add_hline(y=-fc_thresh, line_dash="dash", line_color="#888888", line_width=1)

    fig.update_layout(
        title=f"MA Plot: {contrast_user}  ({run_id})",
        xaxis_title="Average log2 Abundance (AveExpr)",
        yaxis_title="log2 Fold Change",
        height=520,
        width=700,
        legend_title_text="Direction",
    )
    return fig


# ============================================================
# MAIN
# ============================================================

def main() -> None:
    args = parse_args()
    setup_logging()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    run_id = args.run_id

    logging.info(f"Module 04 DIFFERENTIAL_ABUNDANCE -- {run_id}")

    # --- Load params ---
    params = load_params(Path(args.params))

    # --- Load inputs ---
    logging.info("Loading input files...")
    matrix_df  = pd.read_parquet(args.matrix)
    metadata_df = pd.read_parquet(args.metadata)
    mapping_df  = pd.read_parquet(args.id_mapping)

    # --- Split abundance and peptide count columns ---
    abund_df, pep_df, sample_ids = extract_abundance_and_peptide_cols(matrix_df, params)
    n_proteins = len(abund_df)
    n_samples  = len(sample_ids)
    logging.info(f"  {n_proteins} proteins, {n_samples} samples")

    # --- Build group map ---
    group_map = build_group_map(metadata_df, params)
    missing_samples = [s for s in sample_ids if s not in group_map]
    if missing_samples:
        raise ValueError(
            f"Samples in matrix not found in metadata: {missing_samples}"
        )

    groups        = [group_map[s] for s in sample_ids]
    unique_groups = sorted(set(groups))
    logging.info(f"  Groups: {unique_groups}")

    # --- Validate contrasts ---
    contrasts = parse_and_validate_contrasts(params, metadata_df)
    logging.info(f"  Contrasts: {[c[0] for c in contrasts]}")

    # --- Determine method and peptide count availability ---
    da_cfg      = params.get("differential_abundance", {})
    method      = da_cfg.get("method", "deqms").lower()
    has_pep     = not pep_df.empty

    if method == "deqms" and not has_pep:
        logging.warning(
            "DEqMS requested but no peptide count columns found. "
            "Falling back to limma."
        )
        use_deqms   = False
        method_used = "limma"
    elif method == "deqms":
        use_deqms   = True
        method_used = "DEqMS"
    else:
        use_deqms   = False
        method_used = "limma"

    logging.info(f"  Statistical method: {method_used}")

    # --- Summarize peptide counts (DEqMS path only) ---
    pep_counts: "pd.Series | None" = None
    if use_deqms:
        logging.info("Summarizing peptide counts (min of nonzero per protein)...")
        pep_counts = summarize_peptide_counts(pep_df)
        logging.info(
            f"  Peptide count range: {int(pep_counts.min())} - {int(pep_counts.max())}, "
            f"median: {pep_counts.median():.1f}"
        )

    # --- Run statistical analysis per contrast ---
    contrast_results: list[tuple[str, str, str, pd.DataFrame]] = []

    for contrast_user, numerator, denominator, r_contrast_str in contrasts:
        logging.info(f"Running contrast: {contrast_user}  ({r_contrast_str})...")

        raw_df, actual_method = _run_one_contrast_r(
            abund_df, groups, unique_groups,
            pep_counts, r_contrast_str, use_deqms,
        )

        result_df = assemble_results(
            raw_df, mapping_df, params, actual_method, contrast_user
        )

        n_sig  = int(result_df["significant"].sum())
        n_up   = int((result_df["direction"] == "up").sum())
        n_down = int((result_df["direction"] == "down").sum())
        logging.info(
            f"  {n_sig}/{n_proteins} significant "
            f"({n_up} up, {n_down} down)"
        )

        contrast_results.append((contrast_user, numerator, denominator, result_df))

        # --- Generate diagnostic plots ---
        logging.info(f"  Generating plots for {contrast_user}...")
        volcano_fig = plot_volcano(result_df, contrast_user, run_id, params)
        save_plot(volcano_fig, outdir / f"{run_id}.{contrast_user}.volcano_plot")

        ma_fig = plot_ma(result_df, contrast_user, run_id, params)
        save_plot(ma_fig, outdir / f"{run_id}.{contrast_user}.ma_plot")

    # --- Combine all contrast results ---
    combined_df = pd.concat(
        [df for _, _, _, df in contrast_results], ignore_index=True
    )

    # --- Write results table ---
    logging.info("Writing outputs...")
    pq_path  = outdir / f"{run_id}.diff_abundance_results.parquet"
    csv_path = outdir / f"{run_id}.diff_abundance_results.csv"
    combined_df.to_parquet(pq_path,  index=False)
    combined_df.to_csv(csv_path,     index=False)
    logging.info(f"  Saved: {pq_path.name}")
    logging.info(f"  Saved: {csv_path.name}")

    # --- Write text summary ---
    write_summary_txt(contrast_results, run_id, method_used, params, outdir)

    # --- Final log ---
    total_sig = int(combined_df["significant"].sum())
    logging.info(
        f"Module 04 DIFFERENTIAL_ABUNDANCE complete. "
        f"{len(contrasts)} contrast(s), "
        f"{total_sig}/{n_proteins} proteins significant "
        f"(method: {method_used})."
    )


if __name__ == "__main__":
    main()
