# Decision 0002: separate Codex account homes

Status: accepted, 2026-10-31. Source: `cody` at `b47d545`;
`tests/cody-accounts-live.py` at `e5942e0`.

The personal and business selectors use separate Codex homes, SQLite state,
file-backed login credentials, and session histories. Conflicting selectors,
environment-provided Codex credentials, and explicit storage overrides fail
rather than silently selecting a different account. Default cody behavior remains
unchanged. The account label describes the selected home, not a verified login.

On the selected work machines, symlinks share the primary configuration,
AGENTS.md, skills, rules, and plugin code cache. Credential files, session files,
SQLite state, and plugin-local secret directories are not symlinked. Shared
configuration edits affect all three homes. Real account login is performed by
its owner, and external connector authorization may be account-specific.

The two additional homes on each work machine are explicitly enrolled in the
private Mesh inventory and its existing backup policy. No hook is installed.
