# ProSIFT (PROtein Statistical Integration and Filtering Tool)

## Project context
ProSIFT is a pipeline and interactive exploration platform for discovery proteomics. It takes protein abundance data through QC, statistical analysis, functional enrichment, and multi-database annotation (UniProt, DisGeNET, DGIdb, CTD, PubMed), producing comprehensive per-protein profiles compiled into a queryable SQLite database. Sole developer: Reina Hastings, Blanco-Suárez Lab, SDSU.

## Current status
Pre-implementation. All four phases (core skeleton/QC, statistical analysis, database integration, interactive frontend) are not yet started. Architecture, module specs, and database schema are designed and documented.

## Stack
- Primary: Python 3.11+
- R via rpy2 for DEqMS and limma
- Nextflow (DSL2) for pipeline orchestration, following nf-core conventions
- Singularity/Apptainer for containerized reproducibility
- Target cluster: mesx.sdsu.edu (Torque/PBS, not SLURM)
- SQLite for results database, Parquet for intermediates
- gseapy for ORA/GSEA enrichment

## Architecture (three layers)
1. Nextflow pipeline (compute layer) - runs on HPC cluster
2. Structured results directory (data contract between pipeline and frontend)
3. Interactive frontend (technology TBD - Shiny for Python, Streamlit, Panel, or static HTML under consideration)

## Key design conventions
- Centralized UniProt ID mapping at pipeline start, cached as Parquet, used by all downstream modules
- Mouse-to-human ortholog mapping for human-centric databases (DGIdb, DisGeNET)
- All API responses cached locally with configurable lifetime (default 30 days)
- Scientific parameters live in params.yml, infrastructure config in nextflow.config
- Proteins with no human ortholog get empty DGIdb/DisGeNET results (not silently dropped)

## Development workflow
- Code locally on macOS, rsync to cluster for execution
- Git used locally for milestones only (not for shuttling changes to cluster)
- Singularity container planned but not yet built; Conda environment during development

## Working with the master documentation (ProSIFT_project_documentation.Rmd)
- Check it for existing design decisions before proposing alternatives
- Flag contradictions between conversation and what the doc records
- Date all decision register and development log entries
- If the doc is not uploaded, remind me to log updates afterward
