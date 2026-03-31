#!/usr/bin/env python3
"""
Title:         enrichment.py
Project:       ProSIFT (PROtein Statistical Integration and Filtering Tool)
Author:        Reina Hastings (reinahastings13@gmail.com)
Created:       2026-03-31
Last Modified: 2026-03-31
Purpose:       Module 05 ENRICHMENT process. Runs overrepresentation analysis (ORA)
               via gseapy.enrich() and preranked GSEA via gseapy.prerank() against
               local MSigDB GMT files. Operates on gene symbols from Module 04's
               differential abundance results table. Produces a unified enrichment
               results table, a protein-term mapping table, lollipop plots (Plotly,
               PNG + HTML) per contrast per library, GSEA running score plots (gseapy
               built-in, PNG) for top significant terms, and a summary text file.
Inputs:
  --results      {run_id}.diff_abundance_results.parquet  (Module 04 DIFFERENTIAL_ABUNDANCE)
  --params       {run_id}_params.yml
Outputs:
  {run_id}.enrichment_results.parquet
  {run_id}.enrichment_results.csv
  {run_id}.protein_term_mapping.parquet
  {run_id}.enrichment_summary.txt
  {run_id}.{contrast}.{library}.ora_lollipop.png/.html   (per contrast, per library)
  {run_id}.{contrast}.{library}.gsea_lollipop.png/.html  (per contrast, per library)
  {run_id}.{contrast}.{library}.gsea_running_score.{term}.png  (top N terms)
Usage:
  enrichment.py --results CTXcyto_WT_vs_CTXcyto_KO.diff_abundance_results.parquet \
                --params CTXcyto_WT_vs_CTXcyto_KO_params.yml \
                --run-id CTXcyto_WT_vs_CTXcyto_KO \
                --outdir .
"""

import argparse
import datetime
import logging
import math
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import gseapy
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import yaml

# ============================================================
# ARGUMENT PARSING
# ============================================================

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="enrichment.py",
        description="ProSIFT Module 05 ENRICHMENT: ORA and GSEA enrichment analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Example:\n"
            "  enrichment.py \\\n"
            "    --results CTXcyto_WT_vs_CTXcyto_KO.diff_abundance_results.parquet \\\n"
            "    --params CTXcyto_WT_vs_CTXcyto_KO_params.yml \\\n"
            "    --run-id CTXcyto_WT_vs_CTXcyto_KO \\\n"
            "    --outdir .\n"
        ),
    )
    parser.add_argument("--results",  required=True, help="Module 04 diff_abundance_results.parquet")
    parser.add_argument("--params",   required=True, help="Run params.yml")
    parser.add_argument("--run-id",   required=True, dest="run_id", help="Run identifier prefix for output files")
    parser.add_argument("--outdir",   required=True, help="Output directory")
    return parser.parse_args()


# ============================================================
# LOGGING
# ============================================================

def setup_logging() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        stream=sys.stdout,
    )


# ============================================================
# LIBRARY SHORT NAME MAPPING
# ============================================================

# Maps substrings found in GMT filenames to short library identifiers.
# Used for output file naming and the `library` column in results tables.
# Checked in order; first match wins. Falls back to the filename stem.
_LIBRARY_NAME_MAP: List[Tuple[str, str]] = [
    ("go.bp",           "GO_BP"),
    ("go.mf",           "GO_MF"),
    ("go.cc",           "GO_CC"),
    ("kegg_medicus",    "KEGG"),
    ("kegg_legacy",     "KEGG_LEGACY"),
    ("reactome",        "REACTOME"),
    ("hallmark",        "HALLMARK"),
    ("mh.all",          "HALLMARK"),
]

def _library_short_name(gmt_path: str) -> str:
    """Derive a short library identifier from a GMT file path."""
    stem = Path(gmt_path).stem.lower()
    for substring, short_name in _LIBRARY_NAME_MAP:
        if substring in stem:
            return short_name
    # Fallback: uppercase the filename stem, truncated
    return Path(gmt_path).stem.upper()[:20]


# ============================================================
# PARAMETER LOADING
# ============================================================

def load_params(params_path: str) -> dict:
    """Load and validate enrichment parameters from params.yml."""
    params_dir = Path(params_path).parent.resolve()

    with open(params_path) as fh:
        params = yaml.safe_load(fh)

    enr = params.get("enrichment", {})

    # Required: gene_set_libraries must be a list of paths.
    # Paths are resolved relative to the params.yml file's directory so that
    # relative paths work the same way during standalone testing (where the
    # working directory may differ from the run directory) and during Nextflow
    # execution (where the params.yml is staged into the work directory).
    libraries = enr.get("gene_set_libraries", [])
    if not libraries:
        logging.error("params.yml: enrichment.gene_set_libraries is empty or missing. "
                      "Provide at least one GMT file path.")
        sys.exit(1)
    resolved = []
    for p in libraries:
        resolved_p = Path(p) if Path(p).is_absolute() else (params_dir / p).resolve()
        if not resolved_p.exists():
            logging.error("Gene set library not found: %s (resolved from %s)", resolved_p, p)
            sys.exit(1)
        resolved.append(str(resolved_p))
    enr["gene_set_libraries"] = resolved

    # Defaults with explicit type coercion
    enr.setdefault("run_ora",           True)
    enr.setdefault("run_gsea",          True)
    enr.setdefault("background",        "detected")
    enr.setdefault("gsea_ranking",      "t_statistic")
    enr.setdefault("min_gene_set_size", 15)
    enr.setdefault("max_gene_set_size", 500)
    enr.setdefault("fdr_threshold",     0.05)
    enr.setdefault("plot_top_n",        20)
    enr.setdefault("plot_top_gsea_traces", 10)
    enr.setdefault("gsea_permutations", 1000)
    enr.setdefault("gsea_seed",         42)

    params["enrichment"] = enr
    return params


# ============================================================
# GENE SYMBOL PREPARATION
# ============================================================

def prepare_gene_symbols(
    da_df: pd.DataFrame,
    contrast: str,
) -> Tuple[pd.DataFrame, dict]:
    """
    Filter and deduplicate the Module 04 results for one contrast.

    Returns:
        contrast_df: filtered, deduplicated DataFrame for this contrast
        stats: dict of counts for summary reporting
    """
    df = da_df[da_df["contrast"] == contrast].copy()
    n_total = len(df)

    # --- Drop unmapped proteins (null gene_symbol) ---
    n_unmapped = df["gene_symbol"].isna().sum()
    df = df[df["gene_symbol"].notna()].copy()
    if n_unmapped > 0:
        logging.warning(
            "Contrast %s: dropped %d proteins with null gene_symbol", contrast, n_unmapped
        )

    # --- Determine primary adjusted p-value column ---
    # Use deqms_adj_pvalue if present and non-null, otherwise limma_adj_pvalue
    if "deqms_adj_pvalue" in df.columns and df["deqms_adj_pvalue"].notna().any():
        pval_col = "deqms_adj_pvalue"
    else:
        pval_col = "limma_adj_pvalue"

    # --- Deduplicate: one row per gene_symbol, keep lowest adj_pvalue ---
    n_before_dedup = len(df)
    df = (
        df.sort_values(pval_col, ascending=True)
          .drop_duplicates(subset="gene_symbol", keep="first")
          .reset_index(drop=True)
    )
    n_collapsed = n_before_dedup - len(df)
    if n_collapsed > 0:
        logging.info(
            "Contrast %s: collapsed %d duplicate gene symbols (kept most significant per gene)",
            contrast, n_collapsed,
        )

    stats = {
        "n_total":     n_total,
        "n_unmapped":  n_unmapped,
        "n_mapped":    n_total - n_unmapped,
        "n_collapsed": n_collapsed,
        "n_unique":    len(df),
        "pval_col":    pval_col,
    }
    return df, stats


# ============================================================
# GSEA RANKING METRIC
# ============================================================

def build_ranked_series(df: pd.DataFrame, ranking: str, pval_col: str) -> pd.Series:
    """
    Build the gene-symbol-indexed ranked Series for gseapy.prerank().

    ranking options:
      "t_statistic"   -- deqms_t if available, else limma_t
      "signed_log10p" -- sign(log2_fc) * -log10(raw_pvalue), clamped for p=0
      "log2fc"        -- log2_fc directly
    """
    if ranking == "t_statistic":
        t_col = "deqms_t" if ("deqms_t" in df.columns and df["deqms_t"].notna().any()) else "limma_t"
        rnk = df.set_index("gene_symbol")[t_col].astype(float)

    elif ranking == "signed_log10p":
        # Use deqms_pvalue or limma_pvalue (raw, not adjusted -- GSEA does its own FDR)
        raw_p_col = "deqms_pvalue" if ("deqms_pvalue" in df.columns and df["deqms_pvalue"].notna().any()) else "limma_pvalue"
        pvals = df[raw_p_col].astype(float)
        # Clamp p=0 to the smallest nonzero p-value to avoid log(0) = -inf
        min_p = pvals[pvals > 0].min()
        pvals = pvals.clip(lower=min_p)
        signs = np.sign(df["log2_fc"].astype(float))
        scores = signs * (-np.log10(pvals))
        rnk = pd.Series(scores.values, index=df["gene_symbol"])

    elif ranking == "log2fc":
        rnk = df.set_index("gene_symbol")["log2_fc"].astype(float)

    else:
        logging.error("Unknown gsea_ranking: %s. Choose t_statistic, signed_log10p, or log2fc.", ranking)
        sys.exit(1)

    return rnk.dropna()


# ============================================================
# ORA (gseapy.enrich)
# ============================================================

def run_ora(
    sig_genes: List[str],
    background_genes: List[str],
    gmt_path: str,
    library_name: str,
    contrast: str,
    params: dict,
) -> pd.DataFrame:
    """
    Run ORA for one library against one contrast's significant gene set.
    Returns a DataFrame in the unified enrichment results schema.
    Returns an empty DataFrame if no significant genes or no significant terms.
    """
    enr = params["enrichment"]
    fdr_threshold = enr["fdr_threshold"]

    if not sig_genes:
        logging.warning(
            "Contrast %s, library %s: no significant genes -- ORA skipped", contrast, library_name
        )
        return pd.DataFrame()

    try:
        result = gseapy.enrich(
            gene_list=sig_genes,
            gene_sets=gmt_path,
            background=background_genes,
            no_plot=True,
            verbose=False,
            cutoff=1.0,  # Return all terms; we filter ourselves for consistent handling
        )
    except Exception as exc:
        logging.warning("ORA failed for contrast %s, library %s: %s", contrast, library_name, exc)
        return pd.DataFrame()

    res = result.res2d
    if res is None or res.empty:
        logging.info("Contrast %s, library %s ORA: no terms returned", contrast, library_name)
        return pd.DataFrame()

    # --- Map gseapy columns to unified schema ---
    # gseapy enrich columns: Gene_set, Term, Overlap, P-value, Adjusted P-value,
    #                        Odds Ratio, Combined Score, Genes
    def _parse_overlap_size(overlap_str: str) -> int:
        """Parse '3/500' -> 3"""
        try:
            return int(str(overlap_str).split("/")[0])
        except (ValueError, AttributeError):
            return 0

    def _parse_gene_set_size(overlap_str: str) -> int:
        """Parse '3/500' -> 500"""
        try:
            return int(str(overlap_str).split("/")[1])
        except (ValueError, AttributeError, IndexError):
            return 0

    out = pd.DataFrame({
        "term_id":         res["Term"],
        "term_name":       res["Term"],
        "library":         library_name,
        "analysis_type":   "ORA",
        "contrast":        contrast,
        "pvalue":          res["P-value"].astype(float),
        "adj_pvalue":      res["Adjusted P-value"].astype(float),
        "enrichment_score": np.nan,
        "odds_ratio":      res["Odds Ratio"].astype(float),
        "combined_score":  res["Combined Score"].astype(float),
        "gene_set_size":   res["Overlap"].apply(_parse_gene_set_size),
        "overlap_size":    res["Overlap"].apply(_parse_overlap_size),
        "overlap_genes":   res["Genes"].str.replace(";", ";", regex=False),
    })

    logging.info(
        "Contrast %s, library %s ORA: %d terms tested, %d significant (FDR < %.2f)",
        contrast, library_name, len(out),
        (out["adj_pvalue"] < fdr_threshold).sum(), fdr_threshold,
    )
    return out.reset_index(drop=True)


# ============================================================
# GSEA (gseapy.prerank)
# ============================================================

def run_gsea(
    ranked_series: pd.Series,
    gmt_path: str,
    library_name: str,
    contrast: str,
    params: dict,
) -> Tuple[pd.DataFrame, Optional[object]]:
    """
    Run preranked GSEA for one library against one contrast's full ranked list.
    Returns (results_df, prerank_result_object).
    results_df is in the unified enrichment results schema.
    prerank_result_object is needed for running score plots; None on failure.
    """
    enr = params["enrichment"]
    fdr_threshold = enr["fdr_threshold"]

    if ranked_series.empty:
        logging.warning("Contrast %s, library %s: empty ranked list -- GSEA skipped", contrast, library_name)
        return pd.DataFrame(), None

    try:
        pre_res = gseapy.prerank(
            rnk=ranked_series,
            gene_sets=gmt_path,
            min_size=enr["min_gene_set_size"],
            max_size=enr["max_gene_set_size"],
            permutation_num=enr["gsea_permutations"],
            ascending=False,   # highest scores (most upregulated) first
            no_plot=True,
            verbose=False,
            seed=enr["gsea_seed"],
            threads=1,         # deterministic; parallel threads can affect permutation results
        )
    except Exception as exc:
        logging.warning("GSEA failed for contrast %s, library %s: %s", contrast, library_name, exc)
        return pd.DataFrame(), None

    res = pre_res.res2d
    if res is None or res.empty:
        logging.info("Contrast %s, library %s GSEA: no terms returned", contrast, library_name)
        return pd.DataFrame(), pre_res

    # --- Map gseapy columns to unified schema ---
    # gseapy prerank res2d columns:
    #   Name, Term, ES, NES, NOM p-val, FDR q-val, FWER p-val,
    #   Tag %, Gene %, Lead_genes
    def _leading_edge_size(lead_genes_str) -> int:
        if pd.isna(lead_genes_str) or lead_genes_str == "":
            return 0
        return len(str(lead_genes_str).split(";"))

    def _gene_set_size_from_results(term: str, pre_res_obj) -> int:
        """Extract matched gene count from the results dict."""
        try:
            matched = pre_res_obj.results[term].get("matched_genes", [])
            return len(matched) if matched else 0
        except (KeyError, TypeError):
            return 0

    out_rows = []
    for _, row in res.iterrows():
        term = row["Term"]
        lead_genes = row.get("Lead_genes", "")
        lead_size = _leading_edge_size(lead_genes)
        gene_set_size = _gene_set_size_from_results(term, pre_res)

        out_rows.append({
            "term_id":          term,
            "term_name":        term,
            "library":          library_name,
            "analysis_type":    "GSEA",
            "contrast":         contrast,
            "pvalue":           float(row["NOM p-val"]),
            "adj_pvalue":       float(row["FDR q-val"]),
            "enrichment_score": float(row["NES"]),
            "odds_ratio":       np.nan,
            "combined_score":   np.nan,
            "gene_set_size":    gene_set_size,
            "overlap_size":     lead_size,
            "overlap_genes":    str(lead_genes) if pd.notna(lead_genes) else "",
        })

    out = pd.DataFrame(out_rows)

    logging.info(
        "Contrast %s, library %s GSEA: %d terms tested, %d significant (FDR < %.2f)",
        contrast, library_name, len(out),
        (out["adj_pvalue"] < fdr_threshold).sum(), fdr_threshold,
    )
    return out.reset_index(drop=True), pre_res


# ============================================================
# PROTEIN-TERM MAPPING TABLE
# ============================================================

def build_protein_term_mapping(
    da_df: pd.DataFrame,
    gmt_paths: List[str],
    library_names: List[str],
    enrichment_results: pd.DataFrame,
    gsea_results_by_key: Dict[Tuple[str, str], object],
) -> pd.DataFrame:
    """
    Build the many-to-many protein-term mapping table.

    Includes all protein-term pairs for terms that were actually tested
    (i.e., appear in the enrichment results table). Annotates each pair
    with significance status and GSEA leading edge membership.
    """
    if enrichment_results.empty:
        return pd.DataFrame(columns=[
            "gene_symbol", "protein_id", "term_id", "term_name",
            "library", "in_significant_set", "is_leading_edge", "contrast",
        ])

    # Build protein->gene_symbol lookup (all proteins, not just significant)
    protein_gene = (
        da_df[["protein_id", "gene_symbol", "contrast", "significant"]]
        .dropna(subset=["gene_symbol"])
        .drop_duplicates()
    )

    # Build set of significant gene symbols per contrast (for in_significant_set)
    sig_genes_by_contrast: Dict[str, set] = {}
    for contrast in da_df["contrast"].unique():
        sig = da_df[(da_df["contrast"] == contrast) & (da_df["significant"] == True)]
        sig_genes_by_contrast[contrast] = set(sig["gene_symbol"].dropna())

    # Build leading edge sets per (contrast, library, term)
    leading_edge: Dict[Tuple[str, str, str], set] = {}
    for (contrast, lib_name), pre_res in gsea_results_by_key.items():
        if pre_res is None:
            continue
        for term, term_data in pre_res.results.items():
            lead_str = term_data.get("lead_genes", "")
            if lead_str:
                leading_edge[(contrast, lib_name, term)] = set(str(lead_str).split(";"))

    # Parse GMT files to get gene-to-term mapping for tested terms only
    tested_terms: set = set(enrichment_results["term_id"].unique())

    rows = []
    for gmt_path, lib_name in zip(gmt_paths, library_names):
        with open(gmt_path) as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 3:
                    continue
                term_id = parts[0]
                if term_id not in tested_terms:
                    continue
                genes_in_term = set(parts[2:])

                for contrast in da_df["contrast"].unique():
                    # Proteins in this run that are annotated to this term
                    pg = protein_gene[protein_gene["contrast"] == contrast]
                    annotated = pg[pg["gene_symbol"].isin(genes_in_term)]

                    for _, prow in annotated.iterrows():
                        gsym = prow["gene_symbol"]
                        rows.append({
                            "gene_symbol":      gsym,
                            "protein_id":       prow["protein_id"],
                            "term_id":          term_id,
                            "term_name":        term_id,
                            "library":          lib_name,
                            "in_significant_set": gsym in sig_genes_by_contrast.get(contrast, set()),
                            "is_leading_edge":  gsym in leading_edge.get((contrast, lib_name, term_id), set()),
                            "contrast":         contrast,
                        })

    if not rows:
        return pd.DataFrame(columns=[
            "gene_symbol", "protein_id", "term_id", "term_name",
            "library", "in_significant_set", "is_leading_edge", "contrast",
        ])

    mapping_df = pd.DataFrame(rows).drop_duplicates().reset_index(drop=True)
    logging.info("Protein-term mapping: %d rows across %d contrasts", len(mapping_df), da_df["contrast"].nunique())
    return mapping_df


# ============================================================
# VISUALIZATION: LOLLIPOP PLOTS
# ============================================================

# Color scales and visual constants
_LOLLIPOP_COLORSCALE = "Blues_r"     # darker = more significant
_LOLLIPOP_STEM_COLOR = "#cccccc"
_MAX_TERM_LABEL_LEN = 55             # truncate long GO term names for static PNG

def _truncate_label(s: str, maxlen: int = _MAX_TERM_LABEL_LEN) -> str:
    return s if len(s) <= maxlen else s[:maxlen - 3] + "..."

def _make_lollipop_fig(
    terms: List[str],
    x_vals: List[float],
    dot_colors: List[float],   # values mapped to color scale (adj_pvalue)
    dot_sizes: List[int],      # overlap or leading edge size
    x_label: str,
    title: str,
    hover_texts: List[str],
    colorbar_title: str,
) -> go.Figure:
    """
    Build a Plotly lollipop figure.
    Terms are displayed on the y-axis (top = most significant).
    """
    labels = [_truncate_label(t) for t in terms]

    # Stems: horizontal lines from x=0 to each dot
    shapes = []
    for i, xv in enumerate(x_vals):
        shapes.append(dict(
            type="line",
            x0=0, x1=xv,
            y0=i, y1=i,
            line=dict(color=_LOLLIPOP_STEM_COLOR, width=1.5),
        ))

    # Dots
    scatter = go.Scatter(
        x=x_vals,
        y=list(range(len(terms))),
        mode="markers",
        marker=dict(
            color=dot_colors,
            colorscale=_LOLLIPOP_COLORSCALE,
            reversescale=False,
            size=[max(6, min(20, 6 + s // 5)) for s in dot_sizes],
            colorbar=dict(title=colorbar_title, thickness=12, len=0.6),
            line=dict(width=0.5, color="#555555"),
        ),
        text=hover_texts,
        hoverinfo="text",
    )

    fig = go.Figure(data=[scatter])
    fig.update_layout(
        title=dict(text=title, font=dict(size=13)),
        xaxis=dict(title=x_label, zeroline=True, zerolinewidth=1, zerolinecolor="#999999"),
        yaxis=dict(
            tickmode="array",
            tickvals=list(range(len(labels))),
            ticktext=labels,
            autorange="reversed",
        ),
        shapes=shapes,
        plot_bgcolor="#ffffff",
        paper_bgcolor="#ffffff",
        margin=dict(l=300, r=80, t=60, b=60),
        height=max(300, 30 * len(terms) + 120),
    )
    return fig


def _save_lollipop(fig: go.Figure, base_path: str) -> None:
    """Save a lollipop figure as both PNG (kaleido) and HTML (Plotly)."""
    try:
        fig.write_image(base_path + ".png", scale=2)
    except Exception as exc:
        logging.warning("Failed to write PNG %s: %s", base_path + ".png", exc)
    try:
        fig.write_html(base_path + ".html", include_plotlyjs="cdn")
    except Exception as exc:
        logging.warning("Failed to write HTML %s: %s", base_path + ".html", exc)


def plot_ora_lollipop(
    ora_df: pd.DataFrame,
    contrast: str,
    library_name: str,
    run_id: str,
    outdir: Path,
    params: dict,
) -> None:
    """Generate ORA lollipop plot for one contrast + library."""
    enr = params["enrichment"]
    fdr_threshold = enr["fdr_threshold"]
    top_n = enr["plot_top_n"]

    sig = ora_df[ora_df["adj_pvalue"] < fdr_threshold].copy()
    if sig.empty:
        logging.info("Contrast %s, library %s ORA: no significant terms -- lollipop plot skipped", contrast, library_name)
        return

    sig = sig.nsmallest(top_n, "adj_pvalue")
    sig = sig.sort_values("adj_pvalue", ascending=False)  # best at top of plot

    # x-axis: combined score (Enrichr-style effect size)
    x_vals      = sig["combined_score"].tolist()
    dot_colors  = sig["adj_pvalue"].tolist()
    dot_sizes   = sig["overlap_size"].tolist()

    hover_texts = [
        (f"<b>{row['term_id']}</b><br>"
         f"adj p-value: {row['adj_pvalue']:.3e}<br>"
         f"Odds Ratio: {row['odds_ratio']:.2f}<br>"
         f"Overlap: {row['overlap_size']} genes<br>"
         f"Genes: {row['overlap_genes']}")
        for _, row in sig.iterrows()
    ]

    fig = _make_lollipop_fig(
        terms=sig["term_id"].tolist(),
        x_vals=x_vals,
        dot_colors=dot_colors,
        dot_sizes=dot_sizes,
        x_label="Combined Score",
        title=f"ORA: {contrast} | {library_name}",
        hover_texts=hover_texts,
        colorbar_title="-log10(adj p)",
    )

    # Add vertical line at x=0
    fig.add_vline(x=0, line_width=1, line_color="#888888")

    safe_contrast = contrast.replace("/", "_")
    base = str(outdir / f"{run_id}.{safe_contrast}.{library_name}.ora_lollipop")
    _save_lollipop(fig, base)
    logging.info("ORA lollipop saved: %s.png/.html", base)


def plot_gsea_lollipop(
    gsea_df: pd.DataFrame,
    contrast: str,
    library_name: str,
    run_id: str,
    outdir: Path,
    params: dict,
) -> None:
    """Generate GSEA lollipop plot for one contrast + library."""
    enr = params["enrichment"]
    fdr_threshold = enr["fdr_threshold"]
    top_n = enr["plot_top_n"]

    sig = gsea_df[gsea_df["adj_pvalue"] < fdr_threshold].copy()
    if sig.empty:
        logging.info("Contrast %s, library %s GSEA: no significant terms -- lollipop plot skipped", contrast, library_name)
        return

    # Order by absolute NES, show top_n
    sig = sig.reindex(sig["enrichment_score"].abs().nlargest(top_n).index)
    sig = sig.sort_values("enrichment_score", ascending=True)  # negative NES at bottom

    x_vals     = sig["enrichment_score"].tolist()
    dot_colors = sig["adj_pvalue"].tolist()
    dot_sizes  = sig["overlap_size"].tolist()

    hover_texts = [
        (f"<b>{row['term_id']}</b><br>"
         f"NES: {row['enrichment_score']:.3f}<br>"
         f"FDR: {row['adj_pvalue']:.3e}<br>"
         f"Leading edge: {row['overlap_size']} genes<br>"
         f"Genes: {row['overlap_genes'][:200]}{'...' if len(str(row['overlap_genes'])) > 200 else ''}")
        for _, row in sig.iterrows()
    ]

    fig = _make_lollipop_fig(
        terms=sig["term_id"].tolist(),
        x_vals=x_vals,
        dot_colors=dot_colors,
        dot_sizes=dot_sizes,
        x_label="Normalized Enrichment Score (NES)",
        title=f"GSEA: {contrast} | {library_name}",
        hover_texts=hover_texts,
        colorbar_title="FDR q-val",
    )

    # Reference line at NES=0
    fig.add_vline(x=0, line_width=1.5, line_color="#444444")

    safe_contrast = contrast.replace("/", "_")
    base = str(outdir / f"{run_id}.{safe_contrast}.{library_name}.gsea_lollipop")
    _save_lollipop(fig, base)
    logging.info("GSEA lollipop saved: %s.png/.html", base)


# ============================================================
# VISUALIZATION: GSEA RUNNING SCORE PLOTS
# ============================================================

def plot_gsea_running_scores(
    gsea_df: pd.DataFrame,
    pre_res: object,
    contrast: str,
    library_name: str,
    run_id: str,
    outdir: Path,
    params: dict,
) -> None:
    """
    Generate GSEA running enrichment score plots for top N significant terms.
    Uses gseapy's built-in gseaplot() (matplotlib, PNG only).
    """
    enr = params["enrichment"]
    fdr_threshold = enr["fdr_threshold"]
    top_n_traces = enr["plot_top_gsea_traces"]

    sig = gsea_df[gsea_df["adj_pvalue"] < fdr_threshold].copy()
    if sig.empty:
        return

    top_terms = sig.nsmallest(top_n_traces, "adj_pvalue")["term_id"].tolist()

    for term in top_terms:
        try:
            term_data = pre_res.results.get(term)
            if term_data is None:
                continue

            # Sanitize term name for use in a filename
            safe_term = term.replace("/", "_").replace(" ", "_").replace(":", "_")
            safe_contrast = contrast.replace("/", "_")
            ofname = str(
                outdir / f"{run_id}.{safe_contrast}.{library_name}.gsea_running_score.{safe_term}.png"
            )

            gseapy.gseaplot(
                term=term,
                hits=term_data["hits"],
                nes=term_data["nes"],
                pval=term_data["pval"],
                fdr=term_data["fdr"],
                RES=term_data["RES"],
                rank_metric=pre_res.ranking,
                ofname=ofname,
                figsize=(8, 5),
            )
            logging.info("GSEA running score plot: %s", ofname)
        except Exception as exc:
            logging.warning("Failed to generate running score plot for term %s: %s", term, exc)


# ============================================================
# SUMMARY TEXT FILE
# ============================================================

def write_summary(
    run_id: str,
    outdir: Path,
    da_df: pd.DataFrame,
    params: dict,
    gmt_paths: List[str],
    library_names: List[str],
    gene_sym_stats: Dict[str, dict],
    enrichment_results: pd.DataFrame,
    timestamp: str,
) -> None:
    enr = params["enrichment"]
    contrasts = sorted(da_df["contrast"].unique())
    lines = []

    lines += [
        "=" * 40,
        "ENRICHMENT ANALYSIS SUMMARY",
        "=" * 40,
        "",
        f"Run:              {run_id}",
        f"Date:             {timestamp}",
        "",
    ]

    # --- Input ---
    # Use stats from the first contrast as a proxy for shared protein counts
    first_stats = next(iter(gene_sym_stats.values())) if gene_sym_stats else {}
    n_total    = first_stats.get("n_total", len(da_df) // max(len(contrasts), 1))
    n_unmapped = first_stats.get("n_unmapped", 0)
    n_mapped   = first_stats.get("n_mapped", n_total - n_unmapped)
    n_unique   = first_stats.get("n_unique", n_mapped)

    lines += [
        "-" * 40,
        "INPUT",
        "-" * 40,
        "",
        f"Proteins from Module 04:      {n_total}",
        f"  With gene symbol:           {n_mapped} ({100*n_mapped//max(n_total,1)}%)",
        f"  Unmapped (dropped):         {n_unmapped}",
        f"Unique gene symbols:          {n_unique}",
        f"Background (ORA):             {n_unique} unique gene symbols (detected proteins)",
        "",
    ]

    # --- Libraries ---
    lines += [
        "-" * 40,
        "GENE SET LIBRARIES",
        "-" * 40,
        "",
    ]
    header = f"  {'Library':<15} {'GMT file':<55} {'Size filter'}"
    lines.append(header)
    for gmt, lib in zip(gmt_paths, library_names):
        lines.append(
            f"  {lib:<15} {Path(gmt).name:<55} "
            f"{enr['min_gene_set_size']}-{enr['max_gene_set_size']} genes"
        )
    lines.append("")

    # --- Parameters ---
    lines += [
        "-" * 40,
        "PARAMETERS",
        "-" * 40,
        "",
        f"ORA:              {'enabled' if enr['run_ora'] else 'disabled'}",
        f"GSEA:             {'enabled' if enr['run_gsea'] else 'disabled'}",
        f"GSEA ranking:     {enr['gsea_ranking']}",
        f"GSEA permutations:{enr['gsea_permutations']}",
        f"FDR threshold:    {enr['fdr_threshold']}",
        f"Correction:       Benjamini-Hochberg, per-library",
        "",
    ]

    # --- Per-contrast results ---
    for contrast in contrasts:
        n_sig = int(da_df[(da_df["contrast"] == contrast) & (da_df["significant"] == True)]["gene_symbol"].notna().sum())
        n_genes = gene_sym_stats.get(contrast, {}).get("n_unique", "?")

        lines += [
            "-" * 40,
            f"RESULTS: {contrast}",
            "-" * 40,
            "",
            f"Significant genes (ORA input):  {n_sig} / {n_genes}",
            "",
        ]

        if not enrichment_results.empty:
            cdf = enrichment_results[enrichment_results["contrast"] == contrast]

            if enr["run_ora"]:
                lines.append("ORA results:")
                lines.append(f"  {'Library':<12} {'Sig terms':<12} Top term")
                for lib in library_names:
                    ldf = cdf[(cdf["library"] == lib) & (cdf["analysis_type"] == "ORA")]
                    sig_count = (ldf["adj_pvalue"] < enr["fdr_threshold"]).sum()
                    if sig_count > 0 and not ldf.empty:
                        top = ldf.nsmallest(1, "adj_pvalue").iloc[0]
                        top_str = f"{_truncate_label(top['term_id'], 40)} (adj_p={top['adj_pvalue']:.2e})"
                    else:
                        top_str = "none"
                    lines.append(f"  {lib:<12} {sig_count:<12} {top_str}")
                lines.append("")

            if enr["run_gsea"]:
                lines.append("GSEA results:")
                lines.append(f"  {'Library':<12} {'Sig terms':<12} Top term (NES)")
                for lib in library_names:
                    ldf = cdf[(cdf["library"] == lib) & (cdf["analysis_type"] == "GSEA")]
                    sig_count = (ldf["adj_pvalue"] < enr["fdr_threshold"]).sum()
                    if sig_count > 0 and not ldf.empty:
                        top = ldf.reindex(ldf["enrichment_score"].abs().nlargest(1).index).iloc[0]
                        top_str = (
                            f"{_truncate_label(top['term_id'], 35)} "
                            f"(NES={top['enrichment_score']:.2f}, FDR={top['adj_pvalue']:.2e})"
                        )
                    else:
                        top_str = "none"
                    lines.append(f"  {lib:<12} {sig_count:<12} {top_str}")
                lines.append("")
        else:
            lines.append("No enrichment results produced.")
            lines.append("")

    out_path = outdir / f"{run_id}.enrichment_summary.txt"
    out_path.write_text("\n".join(lines) + "\n")
    logging.info("Summary written: %s", out_path)


# ============================================================
# MAIN
# ============================================================

def main() -> None:
    setup_logging()
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    logging.info("=" * 60)
    logging.info("ProSIFT Module 05 ENRICHMENT: %s", args.run_id)
    logging.info("=" * 60)

    # --- Load inputs ---
    params = load_params(args.params)
    enr = params["enrichment"]

    logging.info("Loading Module 04 results: %s", args.results)
    da_df = pd.read_parquet(args.results)
    contrasts = sorted(da_df["contrast"].unique())
    logging.info("Contrasts: %s | Proteins: %d", contrasts, len(da_df) // len(contrasts))

    gmt_paths    = enr["gene_set_libraries"]
    library_names = [_library_short_name(p) for p in gmt_paths]
    logging.info("Libraries: %s", library_names)

    # --- Per-contrast enrichment ---
    all_enrichment: List[pd.DataFrame] = []
    gene_sym_stats: Dict[str, dict] = {}
    gsea_results_by_key: Dict[Tuple[str, str], object] = {}  # (contrast, lib_name) -> pre_res

    for contrast in contrasts:
        logging.info("--- Contrast: %s ---", contrast)

        # Gene symbol preparation (shared across libraries for this contrast)
        contrast_df, stats = prepare_gene_symbols(da_df, contrast)
        gene_sym_stats[contrast] = stats

        sig_genes     = contrast_df[contrast_df["significant"] == True]["gene_symbol"].tolist()
        background_genes = contrast_df["gene_symbol"].tolist()

        ranked_series = build_ranked_series(
            contrast_df,
            ranking=enr["gsea_ranking"],
            pval_col=stats["pval_col"],
        )

        for gmt_path, lib_name in zip(gmt_paths, library_names):
            logging.info("Library: %s", lib_name)

            # ORA
            if enr["run_ora"]:
                ora_df = run_ora(
                    sig_genes=sig_genes,
                    background_genes=background_genes,
                    gmt_path=gmt_path,
                    library_name=lib_name,
                    contrast=contrast,
                    params=params,
                )
                if not ora_df.empty:
                    all_enrichment.append(ora_df)
                    plot_ora_lollipop(
                        ora_df=ora_df,
                        contrast=contrast,
                        library_name=lib_name,
                        run_id=args.run_id,
                        outdir=outdir,
                        params=params,
                    )

            # GSEA
            if enr["run_gsea"]:
                gsea_df, pre_res = run_gsea(
                    ranked_series=ranked_series,
                    gmt_path=gmt_path,
                    library_name=lib_name,
                    contrast=contrast,
                    params=params,
                )
                if not gsea_df.empty:
                    all_enrichment.append(gsea_df)
                    gsea_results_by_key[(contrast, lib_name)] = pre_res
                    plot_gsea_lollipop(
                        gsea_df=gsea_df,
                        contrast=contrast,
                        library_name=lib_name,
                        run_id=args.run_id,
                        outdir=outdir,
                        params=params,
                    )
                    if pre_res is not None:
                        plot_gsea_running_scores(
                            gsea_df=gsea_df,
                            pre_res=pre_res,
                            contrast=contrast,
                            library_name=lib_name,
                            run_id=args.run_id,
                            outdir=outdir,
                            params=params,
                        )

    # --- Combine enrichment results ---
    if all_enrichment:
        enrichment_results = pd.concat(all_enrichment, ignore_index=True)
        # Enforce schema dtypes
        enrichment_results["gene_set_size"] = enrichment_results["gene_set_size"].astype("Int64")
        enrichment_results["overlap_size"]  = enrichment_results["overlap_size"].astype("Int64")
    else:
        logging.warning("No enrichment results produced for any contrast or library.")
        enrichment_results = pd.DataFrame(columns=[
            "term_id", "term_name", "library", "analysis_type", "contrast",
            "pvalue", "adj_pvalue", "enrichment_score", "odds_ratio",
            "combined_score", "gene_set_size", "overlap_size", "overlap_genes",
        ])

    # --- Protein-term mapping ---
    logging.info("Building protein-term mapping table...")
    mapping_df = build_protein_term_mapping(
        da_df=da_df,
        gmt_paths=gmt_paths,
        library_names=library_names,
        enrichment_results=enrichment_results,
        gsea_results_by_key=gsea_results_by_key,
    )

    # --- Write outputs ---
    results_parquet = outdir / f"{args.run_id}.enrichment_results.parquet"
    results_csv     = outdir / f"{args.run_id}.enrichment_results.csv"
    mapping_parquet = outdir / f"{args.run_id}.protein_term_mapping.parquet"
    summary_txt     = outdir / f"{args.run_id}.enrichment_summary.txt"

    enrichment_results.to_parquet(results_parquet, index=False)
    enrichment_results.to_csv(results_csv, index=False)
    mapping_df.to_parquet(mapping_parquet, index=False)
    logging.info("Enrichment results: %s (%d rows)", results_parquet, len(enrichment_results))
    logging.info("Protein-term mapping: %s (%d rows)", mapping_parquet, len(mapping_df))

    write_summary(
        run_id=args.run_id,
        outdir=outdir,
        da_df=da_df,
        params=params,
        gmt_paths=gmt_paths,
        library_names=library_names,
        gene_sym_stats=gene_sym_stats,
        enrichment_results=enrichment_results,
        timestamp=timestamp,
    )

    logging.info("Module 05 ENRICHMENT complete.")


if __name__ == "__main__":
    main()
