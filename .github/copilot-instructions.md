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

## Multi-Model Workflow (cost-aware, quality-gated)

This repository uses the collaboration protocol agreed on 2026-10-19.
Minimize paid-model work without weakening scientific or software validation.

### Roles

- **Implementation worker**: MiMo-V2.6-Pro-RL on Juwels Booster (Blablador endpoint). Prefer this worker for suitable substantial implementation, investigation, tests, and fixes.
- **Latency fallback worker**: DeepSeek-V4.1-Flash on JUWELS Booster (Blablador endpoint). Use when MiMo is unavailable, stalls, or takes too long for the task; apply the same acceptance criteria and independent verification.
- **Reviewer/orchestrator**: the main Copilot chat. Defines scope and acceptance criteria, delegates bounded work, reviews changes, independently verifies consequential behavior, and reports to the user. Do not assume the main chat is MiMo.
- **Senior escalation**: Claude Opus for difficult architecture/scientific questions, unresolved disagreements, or repeated worker failures. Do not require a third management layer for routine tasks. If the user selects Opus as the main chat model, it can combine orchestration and review.
- **User**: assigns tasks and is the source of truth for priorities.

### Routing and Handover

- Invoke MiMo as a model-pinned subagent in the current workflow when supported; do not require separate chats, manual copy/paste, or a GitHub round trip to review local changes.
- The configured worker model identifier is `customendpoint/Blablador/80 - MiMo-V2.6-Pro-RL on Juwels Booster`. If unavailable or tool use fails, report the limitation; do not silently substitute a paid model or claim a delegation succeeded.
- The authorized fallback identifier is `customendpoint/Blablador/15 - DeepSeek-V4.1-Flash on JUWELS Booster`. Set a task-specific latency expectation in the delegation brief; if MiMo exceeds it or stalls, report the switch and hand over the remaining bounded work rather than restarting completed work.
- Before reassigning the same files, confirm the first worker has stopped writing. If cancellation is unavailable, isolate the fallback output and integrate only one result; never allow competing workers to edit the same files concurrently.
- Handle trivial edits, short lookups, and conversational replies directly when delegation would cost more than it saves. Do not offload this routing policy recursively.
- Before delegating, estimate whether orchestration, review and likely corrections will cost less than direct implementation. Delegate only when a meaningful saving is expected; treat review overhead below the estimated solo effort as a routing target, not an unmeasured guarantee.
- Give each worker a bounded objective, relevant context, allowed files, execution constraints, measurable acceptance criteria, and a stopping condition. Avoid sending the full conversation.
- Require a concise handover: changes made, files changed, exact validation commands and outcomes, assumptions, and unresolved issues. Never treat a worker's success claim alone as verification.
- Review the patch instead of redoing the implementation. Independently check critical behavior and inspect actual validation evidence. For scientific analysis, verify population accounting, guarded reference selection, phase conventions, numerical assumptions, and reproducibility.
- Allow at most two focused correction rounds before reassessing or escalating; do not allow unbounded agent discussion. Workers must not launch additional agents unless explicitly assigned to do so.
- Reassess sooner if review/correction overhead threatens the expected saving. Stop redundant delegation or finish directly instead of reimplementing a worker's task through repeated review. Never omit required validation to fit the budget.
- Keep review focused on the patch, critical assumptions and targeted checks. Reuse existing evidence and one concise task log; do not generate duplicate audits, repeated full-context reviews or per-task model benchmarks unless requested.
- Do not send credentials, confidential source, raw experimental data, or other sensitive material to cloud models. All existing local-editing, padiwa-only execution, data-preservation, and Git hygiene rules still apply.
- Model selection must use actual runtime capabilities; instructions do not change the main chat model. Do not promise guaranteed quality or cost savings. Track meaningful validation and correction overhead when evaluating the arrangement.

### Short User-Facing Execution Summary

- End task reports with one or two short lines identifying who performed the work. For direct work, say `Execution: direct (Copilot); no delegation`.
- For delegated work, name the actual worker and any fallback, report elapsed dispatch-to-handover time, worker input/output tokens when available, and corrected defects with a brief description. Distinguish worker correction rounds from reviewer fixes, and disclose unresolved failures.
- Use existing runtime telemetry; report missing timing or token counts as `not reported`, never estimates presented as measurements. Separate worker usage from parent review usage when available; do not double-count reasoning tokens or imply token counts are monetary cost.
- Keep this footer concise. Do not create a new benchmark, scan the full conversation history, or generate extra audit files just to populate it.

### Quality Standard

- Be professional and scientific; hold work to academic and publication-grade standards.
- Deliver reliable research with scientific quality: correct methodology, reproducible results, honest uncertainty, and no overclaiming.
- Code must be pro-level: clear structure, precise naming, minimal side effects, and defensible numerical methods.

### Worker Task Logging (mandatory)

- Every task prompt received by a worker from the user (or orchestrator) MUST be saved verbatim, with the date/time received and all details needed to track it, under `Documents/Worker_Tasks/`. Include the original user prompt alongside any delegated brief, and record the actual assigned model, including a fallback switch.
- One file per task, named `YYYY-MM-DD_HHMMSS_<slug>.md` (local time, Europe/Prague unless stated otherwise).
- Each entry MUST contain: task ID, received timestamp, source (user/orchestrator), assigned model, reviewer, status, the verbatim prompt, assumptions made, work performed, files changed, commits created, and validation evidence.
- Update the entry when the task progresses or completes so the reviewer can follow the full trail without re-running anything.
- See [Documents/Worker_Tasks/README.md](../../Documents/Worker_Tasks/README.md) for the entry template.

### Commit Attribution (mandatory)

- Every commit made by a worker MUST identify the model that actually performed the work, e.g. a trailer line: `Done by: MiMo-V2.6-Pro-RL on Juwels Booster (BLABLADOR)`. Use DeepSeek attribution for its work, not MiMo attribution. Do not attribute work to a model that did not perform it.
- Keep the standard `Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>` trailer as well.

### Reviewer-Friendly Work Style (minimize Copilot review cost)

- Make small, focused commits: one logical change per commit, no unrelated churn or reformatting.
- Write commit messages that state the rationale, not just the diff.
- Keep diffs minimal and surgical; never mix refactors with behavior changes.
- Self-review each diff before handover and record validation evidence (commands run, results) in the task log.
- Prefer clear names, small functions, and explicit data flow so the reviewer can verify correctness by reading, not by running.
- Reference the task-log entry in commit messages when practical so the reviewer can connect code to requirements.
