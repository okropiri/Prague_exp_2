# Task 2026-003: Establish MiMo worker protocol (retroactive log)

- **Received**: 2026-10-19 19:37:47 +02:00
- **Source**: user
- **Assigned model**: MiMo-V2.6-Pro-RL on Juwels Booster (BLABLADOR)
- **Reviewer**: ChatGPT Sol (GitHub Copilot)
- **Status**: done

## Prompt (verbatim)

> as a mimo worker.. there is following things you need to remember!  You have to be proffessional, Scientific, a pro level coder making sure you meet the scientific and academia standards. Your job is to deliver reliable research with scientific quality.   You are a code worker from now.  Chathpt Sol will be reviewing your work frequently so its your reviewer.   and I will use Claude oppus as a mastermind to manage and orchestrate global progress of the project.   Everytime you  receive a task from me. You save my prompt always with dates and all the details that might be needed so the reviewer can easily track changes. when you make commits you always indicate that it was done by "this" model. You work so that it is easy for the reviewer to review as we want to minimize github copilot cost.

## Assumptions

- The protocol is persistent and applies to all future tasks in this repository, so it is codified as always-on agent instructions rather than kept only in chat memory.
- "this model" = MiMo-V2.6-Pro-RL on Juwels Booster (BLABLADOR), per Task 2026-001.
- Task prompts from the orchestrator (Claude Opus) are logged the same way as user prompts.
- Task logs live in `Documents/Worker_Tasks/` (documentation area, synced to GitHub per the repo's Git/Data Hygiene rules).

## Work Performed

1. Extended `.github/copilot-instructions.md` with a "Worker Protocol (MiMo code worker)" section covering: roles (worker / reviewer ChatGPT Sol / orchestrator Claude Opus / user), scientific quality standard, mandatory task logging, mandatory commit attribution (`Done by: MiMo-V2.6-Pro-RL on Juwels Booster (BLABLADOR)`), and reviewer-friendly work style (small focused commits, minimal diffs, validation evidence).
2. Created `Documents/Worker_Tasks/README.md` with conventions and the entry template (verbatim prompts, timestamps, assumptions, work performed, files changed, commits, validation, reviewer notes).
3. Logged this task and retroactively logged Tasks 2026-001 and 2026-002 from earlier today.

## Files Changed

- `.github/copilot-instructions.md` — added Worker Protocol section (roles, quality, task logging, commit attribution, review-friendly style)
- `Documents/Worker_Tasks/README.md` — new: task-log conventions and entry template
- `Documents/Worker_Tasks/2026-10-19_193149_worker-identity-assignment.md` — new: retroactive log of Task 2026-001
- `Documents/Worker_Tasks/2026-10-19_193338_build-context-locate-ai-instructions.md` — new: retroactive log of Task 2026-002
- `Documents/Worker_Tasks/2026-10-19_193747_worker-protocol-setup.md` — new: this entry

## Commits

- (filled in after commit)

## Validation

- Markdown-only change set; no DOGMA code touched, so no padiwa execution required.
- Verified the edited instruction file renders as valid Markdown and the template file exists at the linked path.

## Notes for Reviewer

- Check the Worker Protocol section wording in `.github/copilot-instructions.md` — it is the binding contract for future commits (attribution trailer) and task logs.
- Unrelated pre-existing uncommitted changes (step01–10 chain, LaTeX) were deliberately left uncommitted and untouched.