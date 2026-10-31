# Decision 0003: shared-home credential file switching

Status: accepted for the two configured Linux work machines. This supersedes
[Decision 0002](0002-account-homes.md), retained as historical context.

The unmodified installed Codex CLI uses one home, project tree, session history,
and SQLite state. The selector moves private authentication files between the
usual `auth.json` path and inactive account slots. A private journal makes moves
recoverable after interruption without overwriting credential contents.

A shared lock remains held by each `cody` session. A switch, login, or logout
requires exclusive access and a conservative scan for other Codex processes using
the home. Existing sessions and servers are never stopped automatically. The npm
launcher closes inherited locks, so `cody` resolves the matching installed native
CLI and preserves the npm package environment before executing it directly.
There is no background monitor, prompt hook, patched binary, or custom build.

Direct `codex` invocations do not take the `cody` lock. Already-running direct
processes are checked, but a direct process starting concurrently with a switch
can race that check. Use `cody` consistently for coordinated access. Synthetic
fixtures verify native session resume, lock contention, interrupted file moves,
token refresh, and interactive resume; actual cross-account ChatGPT behavior
requires verification after sign-in.
