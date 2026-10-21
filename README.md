# Prague Cyclotron Exp 2 DOGMA Workspace

This is the Mac-side mirror of the DOGMA/NCAL Prague Cyclotron Experiment 2 project that was previously organized from the Dell laptop.

## Migration Status

- Source on Dell: `/home/dachi/Documents/NCAL_FSD/Prague_Cyclotron_Exp_2/`
- Mac mirror: `/Users/dachi_macbookpro/Documents/PhD/NCAL/Prague2/Dogma/`
- Mirrored on 2026-06-30 with `rsync` over `ssh dell`.
- Excluded from the Mac mirror: Linux conda environments (`.conda/`), nested `.git/` metadata, `.DS_Store`, and AppleDouble `._*` sidecar files.
- No existing GitHub remote was found in the Dell project tree. GitHub CLI/auth is not ready yet on this Mac, and Dell's saved `gh` token is invalid.

## Project Layout

- `Dogma_Analysis/Analysis/`: active cleaned-data analysis chain and runners.
- `Dogma_Analysis/Reference/`: experiment context, RF phase scan notes, and colleague reference material.
- `Dogma_Analysis/Dogma_test_Data/`: raw DOGMA test input; treat as read-only.
- `Dogma_Analysis/Results/`: small mirrored result artifacts only. Full result production belongs on `/data6`.
- `Documents/Latex/`: analysis write-up and supporting figures.
- `Documents/`: detector, electronics, articles, tables, and presentation references.
- `BeamtimeScreens/` and `CAEN_DT5725S/`: setup screenshots and CAEN-specific support material.

## Local Results

All local analysis result copies belong under [Dogma_Analysis/Results/](Dogma_Analysis/Results/), grouped by run type:

- [Position Runs](<Dogma_Analysis/Results/Position Runs/>): Prague position scans, comparison collections, and derived gamma/neutron diagnostics.
  - The canonical [5 m matched RF-method analysis](<Dogma_Analysis/Results/Position Runs/NCAL_20us_Pos_5m_0000/README.md>) has numbered Steps01-10, shared upstream data, two unshifted Step04 alternatives, complete NCAL pulse TTrees, exact histogram-replay TTrees and per-step reports. Start with the [side-by-side 3P comparison](<Dogma_Analysis/Results/Position Runs/NCAL_20us_Pos_5m_0000/Comparison/NCAL_3P_ToT_comparison.png>). Previous canonical ROOT products and their copy information are preserved unchanged in the [archive](<Dogma_Analysis/Results/Position Runs/NCAL_20us_Pos_5m_0000/Archive/2026-10-19_Previous_5m_results/Folded_RF_3x/>).
- [Pulser Runs](<Dogma_Analysis/Results/Pulser Runs/>): Prague external-trigger pulser calibration.
- [Source Runs](<Dogma_Analysis/Results/Source Runs/>): Prague radioactive-source control.
- [Zebra Runs](<Dogma_Analysis/Results/Zebra Runs/>): Zebra laboratory pulser tests.

Keep acquisition names and the existing per-run step, board, and report directories inside these groups. Future local copies must use the corresponding group rather than create separate result roots. This organization is local only: remote production remains under `/data6/Dogma_analysis_by_Dachi/Results/<run_dir>/...`.

## Working Model

- Edit code and documentation on the Mac workspace.
- Run DOGMA analysis code only on `padiwa`, using the real `/data6` filesystem.
- Do not run DOGMA Python runners, C++ analyzers, compilers, notebooks, or data-processing commands locally on the Mac unless this rule is explicitly overridden for a non-DOGMA task.
- Write DOGMA result sets first to `/data6/Dogma_analysis_by_Dachi/Results/<run_dir>/...`.
- Mirror back only small artifacts by default: PNG, PDF, small text summaries/tables, and small ROOT files.
- Keep raw input data read-only unless explicitly changing the data management plan.
- Activate `conda_env` before running DOGMA Python on Linux-side machines that provide that environment.

## GitHub Sync Plan

The recommended first GitHub commit is code, project-authored notes, LaTeX source, setup screenshots, small tables, and small selected result artifacts. Raw DOGMA data, generated bulk results, Linux conda environments, and bulky third-party PDF/manual libraries remain local-only or on `/data6` unless deliberately added later, preferably with Git LFS for large binary material.

Once a private GitHub repository URL is available and authentication is configured on this Mac:

```bash
git init
git branch -M main
git add .
git commit -m "Initial Prague Cyclotron Exp 2 DOGMA mirror"
git remote add origin <private-github-repo-url>
git push -u origin main
```

Then clone or pull the same repository on Dell/padiwa-side working directories so code and documentation stay synchronized through GitHub, while heavy processing outputs remain on `/data6`.