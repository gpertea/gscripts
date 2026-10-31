# Project status

Updated: 2026-10-31. Paths are relative to the gscripts repository.
These records currently cover Codex-wrapper maintenance; other utilities retain
independent behavior and tests.

## Current state

- `cody` at `fd115d7` launches interactive new/resume/fork sessions directly with
  `--no-daemon --yolo`. The display manager and resume proxy are removed.
- Updates retain consent, offline fallback, and platform-specific installation;
  server version comparisons, process inspection, and restarts are removed.
- Browser display defaults, explicit environment choices, terminal settings,
  Windows/MSYS2 behavior, and the existing cluster helper workaround remain.
- Native in-session `/resume` switches saved threads with the intended Full
  Access settings. Mesh continues passive history/process collection separately.
- All seven checkouts have the same reviewed wrapper bytes; code is pushed.
  Windows uses a small launcher exporting its trusted optional patch-helper root.
  Native Windows checks verify normal and opt-in behavior without the old wrapper.
- One desktop has no terminal-visible cody/Codex commands. Its bundled CLI accepts
  embedded mode and the updated wrapper passes with an explicit executable path;
  terminal setup choice is pending. No Codex installation/update was performed.

## Works now

- `cody`, `cody --noupdate resume`, `cody --noupdate resume --last`, and `cody fork`.
- `cody --update` updates the CLI and reports its version; it does not restart
  running sessions. Explicit remote/display-wrapper flags are retired.
- `bash tests/cody-env.sh` passes; `python3 tests/cody-update.py` passes 26 tests.
- `python3 tests/cody-embedded-live.py` passes against installed CLI 0.166.0 with
  isolated homes and a local Responses fixture. New, explicit resume, fork, and
  two actual in-session switches verify saved identities, permissions, and tool
  writes. No account credentials or external model requests are used.

## Decisions taken

- [Decision 0001](docs/decisions/0001-embedded-cody.md): direct embedded sessions,
  native resume, and no wrapper-managed app-server lifecycle.

## Waiting on people

1. Choose whether to expose the app-bundled CLI, install a standalone CLI, or
   leave terminal PATH unchanged on the desktop lacking terminal commands.

## Next

1. Apply the selected terminal setup on the desktop lacking terminal commands;
   other source/launcher deployment checks are complete.
2. Use fresh cody launches for embedded mode; existing processes keep their current
   mode until closed. Future CLI changes require rerunning the native fixture.

## Resuming

Read this file, DEVLOG, and the embedded-session section in README. Use isolated
fixtures for regressions. Private plans and deployment evidence are under ignored
worklog. Preserve unrelated scripts and existing local work. No automatic restart
or Codex update is part of a wrapper rollout.
