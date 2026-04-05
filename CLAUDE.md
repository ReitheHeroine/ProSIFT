# ProSIFT (PROtein Statistical Integration and Filtering Tool)

## Project context
ProSIFT is a pipeline and interactive exploration platform for discovery proteomics. It takes protein abundance data through QC, statistical analysis, functional enrichment, and multi-database annotation (UniProt, DisGeNET, DGIdb, CTD, PubMed), producing comprehensive per-protein profiles compiled into a queryable SQLite database. Sole developer: Reina (Rei) Hastings, Blanco-Suarez Lab, SDSU.

## Current status
See `proSIFT_handoff.Rmd` for current dev state and task queue.

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
3. Interactive frontend (technology TBD)

## Key design conventions
- All API responses cached locally with configurable lifetime (default 30 days)
- Scientific parameters live in params.yml, infrastructure config in nextflow.config

## Development workflow
- Code locally on macOS, rsync to cluster for execution

## Git policy
- **Do not commit or push to Git unless explicitly asked.** Build and verify first; Rei will say when to commit.
- Branch: `main` only (no feature branches unless discussed).

## Build and test output
- Keep build/test output quiet by default. Use `--silent`, `--quiet`, or redirect stderr unless debugging a failure.

## Efficiency
- Batch API and bash calls where possible instead of making many tiny calls. For example, create multiple files in sequence rather than one at a time with confirmation between each, and combine related shell commands.
- When reading files, read only the sections you need. Use line ranges or grep rather than dumping entire files.

## Project documentation

### Documentation files
- `proSIFT_master_doc.Rmd` - Current-state reference for how ProSIFT works. Architecture, module summaries, infrastructure, roadmap, testing, limitations, references. Check it for existing design decisions before proposing alternatives.
- `proSIFT_project_history.Rmd` - Chronological record of cross-cutting design decisions (Decision Register) and development work (Development Log). Append-only.
- `proSIFT_handoff.Rmd` - Operational coordination. Read at every session start.
- `proSIFT_session_log.Rmd` - Append-only session history. Read when catching up on recent work.
- `proSIFT_environment.Rmd` - Platform, rsync, conda, git reference. Read when working with the cluster or committing.
- `proSIFT_file_map.Rmd` - Directory tree and scripts reference. Read when creating or locating files.
- Module specs in `project_documentation/modules/` (01 through 05, plus 03b QC Report Assembly) - Per-module design and implementation details.

### Reading module documentation
Module specs follow a standard structure: Purpose (1), Inputs/Outputs (2), Workflow Diagram (3), Process Specs (4), Design Decisions (5), Parameters (6), Limitations (7), Testing Notes (8), Change History (9). Read sections based on the task type:

**Implementing a new process:** Read Sections 1-4, 6 (full spec). Skip 5, 7-9 (not yet written or not needed). Also read Handoff Section 5 (data contracts) for upstream/downstream interfaces.

**Modifying or debugging an existing process:** Read Section 4 (process specs) for the relevant process. Read Section 8 (testing notes) to understand existing test coverage. Skim Section 7 (limitations) for known edge cases.

**Adding tests or validation:** Read Sections 4 and 8.

**Cross-module work (wiring, integration):** Read Section 2 (I/O) of all involved modules. Read Handoff Section 5 (data contracts).

### General documentation rules
- Flag contradictions between the conversation and what the docs record.
- Date all decision register and development log entries. Cross-cutting decisions (affecting multiple modules or pipeline architecture) go in `proSIFT_project_history.Rmd`. Module-internal decisions go in the relevant module spec's Section 5.
- Use the @prosift-docs subagent to locate relevant sections in .Rmd documentation files before reading them directly. The agent returns file paths and line ranges.

## Session logging
During a session, log work in the handoff file's "Current Session" section (Section 6). At session end, follow the flush procedure in the handoff's session protocol: prepend the entry to `proSIFT_session_log.Rmd` (read only the first ~35 lines to get the structure, insert after the `# Active Session Log` header), then clear Section 6 back to the empty template. Always write to the session log before clearing.

## Code style
- No em dashes anywhere in the codebase or message strings (use single hyphens)
- Use single quotes (`'`) instead of double quotes (`"`) where language allows
- Include argument parsing (`argparse` for Python, `getopts`/`getopt` for Bash)
- Include usage messages so scripts are easy to use by anyone
- Comment with section headers, numbered steps for multi-step algorithms, explain the *why* for non-obvious logic
- Write detailed inline comments for complex operations or domain-specific logic
- For new standalone scripts, use the header template in `.claude/templates/script_header.txt`
- New bin/ scripts must have `chmod +x` applied immediately after creation
