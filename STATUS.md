# Project status

Updated: 2026-11-27. Paths are relative to the gscripts repository.
These records cover Codex-wrapper maintenance, dbPg, and Linux VNC viewers.
Wrapper verification below remains dated 2026-11-05; no wrapper code changed
in the dbPg round.

## dbPg utility

- Version 0.2, implementation 31d4b6e: legacy connection defaults remain;
  explicit options add peer/socket access, read-only mode and autocommit.
- Password selection delegates to libpq. dbQuery binds values; dbQueryHash
  returns named fields. Prepared-statement helpers remain unchanged.
- `DBPG_LIVE_TEST=1 prove tests/dbPg.t` passes all 14 checks, including an
  actual read-only rse connection. Existing pg_db.pl/get_br_demo.pl compile.
- Usage: [docs/dbPg.md](docs/dbPg.md). No pending helper implementation work.

## Current state

- The fallback `jwm-vnc-session.xml` binds Super+Shift+4 to
  `flameshot-vnc gui` (d6a10b4). The shortcut opened Flameshot in an isolated X
  session and in the active VNC desktops on gvlin and gglin.
- Linux `vnc-session-view` backgrounds its viewer after session/tunnel setup
  (dfedc9f). All `vnc-HOST` aliases use this helper, including direct mode.
  Output goes to `${XDG_STATE_HOME:-$HOME/.local/state}/vnc-session/viewer.log`.
  Setup errors and management commands remain foreground operations. Installed
  helpers on this workstation are symlinks to this checkout; no reinstall is
  needed. VNC servers, lin-config, codex-kit, and other platforms are untouched.
- `cody` uses the unmodified installed CLI in embedded mode. Shared-home
  switching was implemented at `a849a03`, retained account selection at `bdac018`,
  and status context at `57008a1`. On two
  configured Linux work machines, account selection swaps private credentials in
  the same Codex home and preserves sessions, projects, and SQLite state.
- Plain `cody` uses the selected account. `wcody`/`cody -w` selects the original
  work account; `gcody`/`cody -g` and `bcody`/`cody -b` select alternates.
  `cody-account` reports status and `cody-account w`, `g`, or `b` switches.
- `cody-account` and `--status` show PID, live working directory, verified tmux
  session/window/pane and socket, and cached Mesh project labels matched by cwd.
  Existing local schema-2 catalogs and schema-1 registries are queried read-only;
  Mesh discovery follows `MESH_STATE_DIR`, `XDG_STATE_HOME`, then the Linux default.
  Missing optional metadata never removes conservatively detected processes.
- Account changes and login mutations stop when other sessions or servers may
  use the home; they do not stop those processes.
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

- `python3 tests/vnc-viewer-background.py`: three tests pass, including six
  alias/transport combinations, viewer survival after SIGHUP, redirected I/O,
  setup failure propagation, and foreground list/status/stop. SSH and viewer
  are test doubles; SSH PID discovery is stubbed. No live GUI was opened.
- Existing VNC tunnel, host CLI, user-alias, Windows title/tunnel shell suites
  and all 30 macOS routing fixtures pass. `bash -n vnc-session-view` passes;
  ShellCheck passes with only existing SC2155/SC2015 findings excluded.
- `cody`, `cody --noupdate resume`, `cody fork`, `cody-account`, and
  `cody-account --status`.
- `cody -g login`, `cody -b login`, `wcody resume`, and `bcody resume` are
  available once other Codex sessions and servers using that home are closed.
- Unit checks cover file-move recovery, process detection, locks, signal and
  native-package handling, and selector/environment behavior. Existing display
  and 26 update tests pass.
- Status fixtures pass 15 tests, including both status forms, unmanaged homes,
  a real private tmux server, renamed sessions, shell descendants, stale pane
  hints, missing/timed-out tmux, both Mesh formats, local-node filtering,
  directory matching, discovery, read-only access, and metadata failure paths.
  Account swap/selector suites pass 23 tests; a local read-only status check
  also displayed live cwd/tmux context and cached Mesh labels.
- Previously, isolated installed-CLI fixtures passed on versions 0.166.0 and 0.165.0:
  shared-session resume, refresh preservation, blocking an active-session swap,
  lock lifetime, and SQLite integrity. Interactive fixtures passed native resume,
  fork, two in-session `/resume` switches, and actual fixture tool writes.
  These used fictional credentials and no external model requests.
- The updated installed-CLI fixture was blocked at login by conservative
  process detection and has not completed in this round.
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
Status regression command: `python3 tests/cody-account-status.py`.
VNC viewer regression command: `python3 tests/vnc-viewer-background.py`.
