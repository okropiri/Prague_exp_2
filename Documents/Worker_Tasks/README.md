# Worker Task Log

Mandatory task log for the MiMo code worker (MiMo-V2.6-Pro-RL on Juwels Booster, BLABLADOR).

Every task prompt received from the user or the orchestrator (Claude Opus) is saved here
verbatim, with timestamps and full execution details, so the reviewer (ChatGPT Sol, GitHub
Copilot) can track requirements and changes without re-running anything. This keeps GitHub
Copilot review cost low.

## Conventions

- One file per task: `YYYY-MM-DD_HHMMSS_<slug>.md` (Europe/Prague local time).
- Prompts are stored **verbatim**, including typos and informal wording. Never paraphrase a prompt.
- Retroactively created entries must say `(retroactive log)` in the title.
- Update the entry as work progresses; the entry is the reviewer's single source of truth for the task.

## Entry Template

```markdown
# Task <id>: <short title>

- **Received**: <YYYY-MM-DD HH:MM:SS +02:00>
- **Source**: <user | Claude Opus (orchestrator)>
- **Assigned model**: MiMo-V2.6-Pro-RL on Juwels Booster (BLABLADOR)
- **Reviewer**: ChatGPT Sol (GitHub Copilot)
- **Status**: <received | in progress | done | blocked>

## Prompt (verbatim)

<exact prompt text, in a blockquote or fenced block>

## Assumptions

<assumptions made to resolve ambiguity; empty if none>

## Work Performed

<what was done, in order, with rationale>

## Files Changed

- <path> — <why>

## Commits

- <sha> — <subject> (Done by: MiMo-V2.6-Pro-RL on Juwels Booster (BLABLADOR))

## Validation

<commands run and their results, or "not applicable (documentation-only)" etc.>

## Notes for Reviewer

<what to check first, known limitations, follow-ups>
```