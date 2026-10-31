# Development log

Newest round first. Paths are relative to the gscripts repository root.

## 2026-10-31  Isolated account homes with shared settings

Branch: `master`. Implementation: `b47d545`, documentation/native fixture:
`15a987b` and `e5942e0`. Sources: `cody`, `tests/cody-accounts.py`, and
`tests/cody-accounts-live.py` at those revisions. Private setup and verification
records remain under ignored `worklog/`.

### Done

- Added personal/business wrapper flags and symlink aliases. Explicit selection
  fences Codex home, SQLite location, and file-backed login storage; conflicting
  selectors and credential/storage overrides fail closed.
- Configured private homes on two work machines. Shared the primary settings,
  instructions, skills, rules, and plugin code cache via symlinks, without sharing
  Codex credentials, sessions, SQLite state, or plugin-local secret directories.
- Enrolled all three homes per work machine in the existing private Mesh catalog
  with explicit owner authorization. No hooks were installed.

### Measured

- Five selector tests, 26 update tests, and the environment test pass. Actual
  installed CLIs on both work machines pass the isolated account fixture: signed-out
  login checks, shared configuration writes, real fixture session and SQLite files
  in the selected home, and unchanged primary auth/history sentinels.
- All four new history sources completed without errors; monitors remained fresh
  and sync acknowledgments caught up. Actual account login remains outstanding.
- A local footer profile-label fixture works with CLI 0.166.0; CLI 0.165.0 did not
  show it. A profile label identifies configuration, not verified account identity.

### Not done

- No user login was performed; external connector authorization may need attention.
- Codex has no documented account-email footer item. Profile-label behavior was
  tested in isolated local fixtures; no persistent status-line settings were changed.

### Next

- Sign in to each desired account and verify identity using Codex's account status.

## 2026-10-31  Align Mac Codex CLI installations

Branch: `master`. Source and previous fleet deployment: `fd115d7`; this round
changed the installed CLI on one machine and did not change wrapper source.

### Done

- Installed the same global npm Codex CLI version on the second Mac using its
  existing Homebrew Node/npm installation. Authentication and session storage were
  not modified.

### Measured

- Fresh login shells on both Macs resolve `cody` from the scripts checkout and
  `codex` from Homebrew's global npm installation. Both report CLI 0.164.0 through
  `codex --version` and `cody --noupdate --version`; `--no-daemon --help` passes.

### Not done

- No real model turn or interactive resume was run in this round. Existing
  sessions retain their launch mode until reopened.

### Next

- Use a fresh cody session on either Mac for embedded operation.

## 2026-10-31  Embedded cody and native resume

Implementation: `c2e2635` and `fd115d7` on `master`. Source: `cody` and the three
current `tests/cody-*` fixtures at `fd115d7`. Approved plans and operational evidence
remain under ignored `worklog/`.

### Done

- Replaced managed display/server/proxy routing with direct embedded launches.
  Retained Full Access, native cwd/arguments, browser context, terminal choices,
  Windows installer isolation, docs opt-in, and bounded consent-based updates.
- Removed obsolete display manager, permission proxy, and their old tests.
  Update reporting no longer contacts or manages servers or running processes.
- Preserved optional Windows color-patch support through a trusted helper loaded
  only for explicit opt-in. Both update/check and launch paths now use canonical
  cody; the previous Windows wrapper is never executed.
- Updated all seven checkouts, preserving existing local work. Windows launcher
  identities were backed up; unrelated Windows repository edits remain intact.
- Added current project records and an accepted embedded-mode decision.

### Measured

- Shell syntax and environment/routing checks pass. All 26 update regression tests
  pass, including Windows doubles, failures, administrative commands, and retired
  options rejected before update work.
- Installed CLI 0.166.0 passes the loopback fixture through the actual wrapper
  at `c2e2635` (normal mode is unchanged by optional patch support at `fd115d7`):
  new, explicit resume, fork, and two in-session resume switches. Real git staging
  and outside-workspace writes confirm selected-thread permissions. The fixture
  uses 29 local requests and zero external model requests or account credentials.

### Not done

- One desktop has no terminal-visible cody/Codex commands. Its existing bundled
  CLI accepts embedded mode; native wrapper checks pass with the explicit binary
  path. Terminal exposure or standalone installation awaits a preference.
- Native Windows checks exercise installed CLI 0.160.0 help and isolated launcher
  routing, not full Windows TUI conversation switching. In-session switches were
  verified on Linux CLI 0.166.0; other installed CLI versions were not upgraded.
- Existing processes retain their launch mode; no automatic shutdown or Codex
  update was performed.

### Next

- Finish the selected desktop terminal setup and use fresh cody launches for
  embedded sessions. Preserve independent passive Mesh collection.

## 2026-10-31  Baseline reconstructed from git history

Sources: `cody`, `cody-display`, `cody-resume`, and README at `20bf213`;
related integration commit: `b00e301`. These are historical observations.

The prior implementation routed Linux display sessions through per-display
app-servers and a permission-restoration proxy for remote resume/fork. It included
server-version reporting/restart consent, terminal/browser environment defaults,
Windows standalone installation support, and protocol/native test fixtures. The
embedded round supersedes the server/proxy design while preserving platform and
update behavior.
