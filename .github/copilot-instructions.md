# Repository Instructions

## Project Context

- This repository is the Mac-side mirror of the Prague Cyclotron Experiment 2 DOGMA/NCAL workspace formerly organized on Dell.
- The Dell source path was `/home/dachi/Documents/NCAL_FSD/Prague_Cyclotron_Exp_2/`.
- The heavy compute/data target is `padiwa` with result storage under `/data6/Dogma_analysis_by_Dachi/Results/`.
- Create and edit code locally on the Mac, but run DOGMA analysis code only on `padiwa`.
- Do not run DOGMA Python runners, C++ analyzers, compilers, notebooks, or data-processing commands locally on the Mac unless the user explicitly overrides this rule for a non-DOGMA task.

## Layout

- `Dogma_Analysis/Analysis/` contains active analysis scripts, C++ analyzers, and Python runners.
- `Dogma_Analysis/Archive/` and `Dogma_Analysis/Analysis/Archive/` contain older or replaced analysis code.
- `Dogma_Analysis/Reference/experiment_context.md` is the first timing/detector context note to read.
- `Dogma_Analysis/Reference/ncal_rf_phase_scan_walkthrough.md` explains the RF phase scan at a human level.
- `Documents/Latex/` contains the write-up.
- `data6/` is a local mount/reference path and must not be committed.

## Analysis Rules

- Treat `Dogma_Analysis/Dogma_test_Data/` as read-only raw/test input.
- Cleaning and validation happen before downstream analysis.
- Downstream analyses should consume cleaned data, not raw DOGMA export, unless explicitly requested.
- Keep result layout run-key based: `Results/<run_dir>/...`.
- Preserve the remote result layout on `/data6`.
- For ch0-referenced analysis, use the guarded first-valid ch0 hit in file order rather than simply the minimum-time ch0 hit.

## Execution Model

- Edit locally, inspect code locally, and commit/sync through GitHub.
- Execute DOGMA analysis commands only on `padiwa`, using the real `/data6` filesystem.
- The MacBook Pro workspace does not require `conda_env`; use it for editing, git operations, documentation, and lightweight inspection only.
- Activate `conda_env` before DOGMA Python commands when working from the Dell laptop/Linux-side workflow where that environment exists.
- Current cleaned-data entrypoints are:
  - `step01_run_dogma_raw_data_refiner.py`
  - `step02_run_cleaned_all_channel_rates.py`
  - `step03_run_cleaned_all_channel_time_tot.py`
  - `step04_run_cleaned_all_channel_rf_period.py`
  - `run_all_cleaned_analysis_chain.py`

## Git/Data Hygiene

- Do not commit `.conda/`, `data6/`, Python caches, LaTeX build products, or generated bulk result trees.
- GitHub should synchronize code, notes, LaTeX sources, reference docs, and small selected artifacts.
- Full raw data and full generated result sets stay on `padiwa`/`/data6` unless the user explicitly asks for a different storage plan.

## Worker Protocol (MiMo code worker)

This section defines the multi-model collaboration protocol established on 2026-10-19. It applies to all work in this repository.

### Roles

- **Code worker (this model)**: MiMo-V2.6-Pro-RL on Juwels Booster (BLABLADOR endpoint). Implements all requested work: code, analysis, documentation, and fixes.
- **Reviewer**: ChatGPT Sol (GitHub Copilot). Reviews worker output frequently on GitHub.
- **Orchestrator**: Claude Opus. Manages and orchestrates global project progress; its task assignments are treated like user tasks and logged the same way.
- **User**: assigns tasks and is the source of truth for priorities.

### Quality Standard

- Be professional and scientific; hold work to academic and publication-grade standards.
- Deliver reliable research with scientific quality: correct methodology, reproducible results, honest uncertainty, and no overclaiming.
- Code must be pro-level: clear structure, precise naming, minimal side effects, and defensible numerical methods.

### Task Logging (mandatory)

- Every task prompt received from the user (or orchestrator) MUST be saved verbatim, with the date/time received and all details needed to track it, under `Documents/Worker_Tasks/`.
- One file per task, named `YYYY-MM-DD_HHMMSS_<slug>.md` (local time, Europe/Prague unless stated otherwise).
- Each entry MUST contain: task ID, received timestamp, source (user/orchestrator), assigned model, reviewer, status, the verbatim prompt, assumptions made, work performed, files changed, commits created, and validation evidence.
- Update the entry when the task progresses or completes so the reviewer can follow the full trail without re-running anything.
- See [Documents/Worker_Tasks/README.md](../../Documents/Worker_Tasks/README.md) for the entry template.

### Commit Attribution (mandatory)

- Every commit made by this worker MUST state that it was done by this model, e.g. a trailer line: `Done by: MiMo-V2.6-Pro-RL on Juwels Booster (BLABLADOR)`.
- Keep the standard `Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>` trailer as well.

### Reviewer-Friendly Work Style (minimize Copilot review cost)

- Make small, focused commits: one logical change per commit, no unrelated churn or reformatting.
- Write commit messages that state the rationale, not just the diff.
- Keep diffs minimal and surgical; never mix refactors with behavior changes.
- Self-review each diff before handover and record validation evidence (commands run, results) in the task log.
- Prefer clear names, small functions, and explicit data flow so the reviewer can verify correctness by reading, not by running.
- Reference the task-log entry in commit messages when practical so the reviewer can connect code to requirements.
