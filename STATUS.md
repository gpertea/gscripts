# Project status

Updated: 2026-10-31. Paths are relative to the gscripts repository.
These records currently cover Codex-wrapper maintenance; other utilities retain
independent behavior and tests.

## Current state

- `cody` at `659e170` on the two work machines and `fd115d7` elsewhere launches
  interactive new/resume/fork sessions directly with
  `--no-daemon --yolo`. The display manager and resume proxy are removed.
- Updates retain consent, offline fallback, and platform-specific installation;
  server version comparisons, process inspection, and restarts are removed.
- Browser display defaults, explicit environment choices, terminal settings,
  Windows/MSYS2 behavior, and the existing cluster helper workaround remain.
- Native in-session `/resume` switches saved threads with the intended Full
  Access settings. Mesh continues passive history/process collection separately.
- The two work machines use the account-selector wrapper; other machines retain
  the previous embedded wrapper. Reviewed code is pushed.
  Windows uses a small launcher exporting its trusted optional patch-helper root.
  Native Windows checks verify normal and opt-in behavior without the old wrapper.
- Personal and business homes are configured on both work machines. Configuration,
  instructions, skills, rules, and plugin code cache are shared by symlink;
  Codex credentials, SQLite state, and sessions are independent. All three
  homes per machine are enrolled in the private Mesh inventory.
- Both Macs now expose `cody` and the same npm-installed Codex CLI version
  through their normal login-shell paths. The app-bundled CLI remains independent.

## Works now

- `cody`, `cody --noupdate resume`, `cody --noupdate resume --last`, and `cody fork`.
- `cody --update` updates the CLI and reports its version; it does not restart
  running sessions. Explicit remote/display-wrapper flags are retired.
- `bash tests/cody-env.sh` passes; `python3 tests/cody-update.py` passes 26 tests.
  `python3 tests/cody-accounts.py` passes five selector/isolation tests; the
  native `tests/cody-accounts-live.py` at `3f1ba19` passes on both work machines,
  including actual shared-instruction/skill loading and isolated session/SQLite writes.
- `python3 tests/cody-embedded-live.py` passes against installed CLI 0.166.0 with
  isolated homes and a local Responses fixture. New, explicit resume, fork, and
  two actual in-session switches verify saved identities, permissions, and tool
  writes. No account credentials or external model requests are used.

## Decisions taken

- [Decision 0001](docs/decisions/0001-embedded-cody.md): direct embedded sessions,
  native resume, and no wrapper-managed app-server lifecycle.
- [Decision 0002](docs/decisions/0002-account-homes.md): independent account
  state with explicitly shared settings and enrolled metadata.

## Waiting on people

1. Sign in to the personal and business homes on each work machine as desired;
   then check the active account with `/status`.

## Next

1. After login, verify account identity and any external connector authorization.
2. Use fresh cody launches for embedded mode; existing processes keep their current
   mode until closed. Future CLI changes require rerunning the native fixture.

## Resuming

Read this file, DEVLOG, and the embedded-session section in README. Use isolated
fixtures for regressions. Private plans and deployment evidence are under ignored
worklog. Preserve unrelated scripts and existing local work. No automatic restart
or Codex update is part of a wrapper rollout.
