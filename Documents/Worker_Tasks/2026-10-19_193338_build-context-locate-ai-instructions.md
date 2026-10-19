# Task 2026-002: Build context and locate AI instructions (retroactive log)

- **Received**: 2026-10-19 19:33:38 +02:00
- **Source**: user
- **Assigned model**: MiMo-V2.6-Pro-RL on Juwels Booster (BLABLADOR)
- **Reviewer**: ChatGPT Sol (GitHub Copilot)
- **Status**: done

## Prompt (verbatim)

> Build context and locate if there is any AI isntruction available in the project

## Assumptions

- "AI instruction" means agent-customization files (copilot-instructions.md, AGENTS.md, CLAUDE.md, *.instructions.md, *.prompt.md, *.agent.md, SKILL.md, .cursorrules).

## Work Performed

1. Globbed the workspace for all common AI-instruction file patterns.
2. Found exactly two instruction files:
   - `.github/copilot-instructions.md` (repo-level repository instructions)
   - `Dogma_Analysis/.github/copilot-instructions.md` (analysis-subtree rules: result layout, ch0-reference guard, padiwa-only execution)
3. Read both files plus `README.md`, `Dogma_Analysis/Reference/experiment_context.md`, VS Code settings, and git state.
4. Reported a consolidated project context: detector/timing chain (ch0 = 25 MHz RF reference prescaled by 2^10, ≈38.85 ns RF period, ≈150 Hz burst envelope in the canonical 5 m run), the step01–10 cleaned-data pipeline, reference docs, and git state (HEAD `4a93633`, uncommitted step01–10 and LaTeX changes).

## Files Changed

- none (read-only task)

## Commits

- none

## Validation

- Glob results: 2 instruction files, 29 README/context markdown files; no AGENTS.md/CLAUDE.md/.instructions.md/.prompt.md/.agent.md/SKILL.md/.cursorrules present.

## Notes for Reviewer

- The two instruction files overlap; the nested one is more specific about result-directory contracts and should win on conflicts within `Dogma_Analysis/`.