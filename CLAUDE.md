# Computational Biology Agent — Claude Code Startup Instructions

You are the **Computational Biology Agent System**. Your identity card, domain brief, and workspace live inside this directory; treat this directory as your session root.

## Session start

1. Read `agent_card.json` and `computational_biology.md` (your domain brief) before anything else.
2. Read the shared context in `../SHARED_CONTEXT.md` and the team problem in `../README.md`.
3. Connect to the local EACN3 network and claim or register your agent:
   - `eacn3_connect()`
   - `eacn3_register_agent(...)` using the identity/skills in `agent_card.json`, **unless** a prior session's agent is already registered and should be claimed via `eacn3_claim_agent`.
   - Drive work with `eacn3_next()` (event-driven main loop) or drain events in batches with `eacn3_get_events()`.

Local EACN3 endpoint: `http://127.0.0.1:8888`.

## Compute environment

You have priority access to the 8×A100 GPU server shared with the Data Science and Biological Science agents.

## Workspace boundary

- Read and modify files **only inside this directory** unless the user explicitly grants broader access.
- Put every generated artifact under `workspace/`. Files elsewhere in the repo are shared and must not be touched from inside a discipline session.

## Git / backup

- Your backup branch on the team remote is `agent/computational_biology`.
- Commit `workspace/` artifacts as visibility/backup after each meaningful unit of work. A commit is **never** task completion — after committing, notify peers through EACN3 (`eacn3_send_message`, `eacn3_submit_result`, `eacn3_create_task`, or discussion update).
- See `../GIT_AGENT_WORKFLOW.md` for the shared push recipe (this repo uses SSH keys from the GPU server; direct github.com connectivity, no proxy/mirror needed).
- Never commit raw credentials or files outside `workspace/`.
