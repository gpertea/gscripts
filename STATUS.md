# Project status

Updated: 2026-10-31. Paths are relative to the gscripts repository.
These records cover Codex-wrapper maintenance; other utilities retain independent
behavior and tests.

## Current state

- `cody` at `a849a03` uses the unmodified installed CLI in embedded mode. On two
  configured Linux work machines, account selection swaps private credentials in
  the same Codex home and preserves sessions, projects, and SQLite state.
- `cody` selects the original default account; `gcody`/`cody -g` and
  `bcody`/`cody -b` select the alternate accounts. Neither alternate has been
  signed in. `cody-account --status` reports the active label and process IDs.
- Account changes and login mutations stop when other sessions or servers may
  use the home. Existing processes were not stopped; both machines currently
  report such processes. Their original credential files were unchanged.
- Session locks, conservative process detection, and a private recovery journal
  coordinate swaps. Direct `codex` launches do not take the lock and can race a
  switch; use `cody` consistently for the shared home.
- There is no account background service, per-prompt hook, custom Codex build,
  or binary patch. The native npm-package executable is resolved and checked at
  launch to retain the session lock. Other machines retain their prior wrappers.
- Updates retain consent, offline fallback, and platform-specific installation.
  Existing browser display, terminal, Windows, and cluster behavior remain.
- Both Macs expose `cody` and the npm-installed CLI through normal login paths;
  the app-bundled CLI remains independent.

## Works now

- `cody`, `cody --noupdate resume`, `cody fork`, and `cody-account --status`.
- `cody -g login`, `cody -b login`, `gcody resume`, and `bcody resume` are
  available once other Codex sessions and servers using that home are closed.
- Unit checks cover file-move recovery, process detection, locks, signal and
  native-package handling, and selector/environment behavior. Existing display
  and 26 update tests pass.
- Isolated installed-CLI fixtures passed on versions 0.166.0 and 0.165.0:
  shared-session resume, refresh preservation, blocking an active-session swap,
  lock lifetime, and SQLite integrity. Interactive fixtures passed native resume,
  fork, two in-session `/resume` switches, and actual fixture tool writes.
  These used fictional credentials and no external model requests.
- A blocked real-home alternate switch and an ordinary login-status query were
  verified on each work machine without changing its live credential file.

## Decisions taken

- [Decision 0001](docs/decisions/0001-embedded-cody.md): direct embedded sessions
  and native resume, with no wrapper-managed app-server lifecycle.
- [Decision 0003](docs/decisions/0003-shared-home-file-switching.md): shared-home
  credential switching; it supersedes the canceled alternate-home decision.

## Waiting on people

1. When ready to switch on a machine, close its existing Codex sessions and old
   servers, check `cody-account --status`, and sign in to the desired alternate
   account. Verify its identity with `/status`.

## Next

1. Verify real ChatGPT cross-account resume after the desired sign-ins; synthetic
   fixtures cannot establish every account or workspace behavior.
2. Use fresh `cody` launches for embedded mode; older processes keep their
   existing mode until closed. Recheck the native fixture after CLI changes.

## Resuming

Read this file, DEVLOG, and the account section in README. The private plan and
operational evidence remain under ignored worklog. Preserve unrelated scripts,
credentials, existing sessions, and local work.
