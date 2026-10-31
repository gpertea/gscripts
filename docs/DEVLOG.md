# Development log

Newest round first. Paths are relative to the gscripts repository root.

## 2026-10-31  Retain the selected account

Branch: `master`. Implementation: `bdac018`.

### Done

- Plain `cody` retains the active account; `wcody` and `cody -w` select the
  original work slot. `cody-account w`, `g`, or `b` switches without launching
  Codex and uses the same locks and process checks as session launches.

### Measured

- Selector and credential-swap fixtures passed 23 tests; display and argument
  routing checks passed. Fixtures cover retained selection, all shortcuts,
  explicit switches, and blocking a switch while a session holds the lock.

### Not done

- The installed-CLI loopback fixture was blocked at its initial login by
  conservative process detection; the updated native path was not verified.

### Next

- Recheck the installed-CLI fixture when process detection permits it, and
  verify real cross-account resume when practical.

## 2026-10-31  Account status and usage

Branch: `master`. Implementation: `4393dba`.

### Done

- Made `cody-account` show the selected account, available login email,
  process status, and usage when called without arguments. Authentication
  tokens are not displayed.

### Measured

- The 12 account-swap fixture tests passed, including no-argument output,
  missing setup, environment credentials, and token privacy.

### Not done

- Cross-account resume with real accounts has not been verified in this round.

### Next

- Use the account commands as needed and verify real cross-account resume when
  the alternate accounts have been signed in.

## 2026-10-31  Shared-home credential switching

Branch: `master`. Implementation: `a849a03`; the approved plan remains under
ignored `worklog/`.

### Done

- Replaced alternate-home selection with private credential-file switching in
  one Codex home. Preserved the installed CLI, project files, session history,
  and SQLite state. Added session locks, process checks, a recovery journal, and
  an account/process status command.
- Initialized private account metadata on both work machines without changing
  their existing auth files or stopping any session or server.

### Measured

- Synthetic native-CLI tests passed on versions 0.166.0 and 0.165.0 for shared
  session resume, token refresh, switch blocking, lock lifetime, and database
  integrity. Interactive fixtures passed native resume and in-session switching.
- On both live homes, a switch was blocked while existing processes remained,
  the original credential bytes were unchanged, and default login status worked.
- A local fixture measured approximately 50 ms median for coordination and
  75 ms for a switch; these are launch-time measurements, not prompt overhead.

### Not done

- Alternate account sign-in and real ChatGPT cross-account resume remain
  unverified. Direct `codex` launches cannot be fully coordinated with `cody`.

### Next

- Close relevant sessions and servers when convenient, then sign in and verify
  real cross-account resume on the desired work machines.

## 2026-10-31  Cancel alternate account homes

Branch: `master`. Documentation: `6da5320`. Sources and tests for the installed
selector remain at `659e170` and `3f1ba19`. The replacement is not deployed.

### Done

- Removed the unused alternate homes on both work machines and restored Mesh
  enrollment to the primary home only. Preserved existing sessions and logins.
- Corrected README guidance; the old alternate selectors now fail closed.

### Measured

- Both monitors are active and both producer queues are caught up.
- An isolated installed-CLI fixture resumed one session using two fictional
  credential files. A synthetic refresh updated only its selected file; its
  model turn did not complete due to a network-policy error.

### Not done

- Shared-home account selection and live-account verification are unfinished.

### Next

- Validate a maintainable per-process credential selector before deployment.

## 2026-10-31  Account-home native verification follow-up

Branch: `master`. Sources: `cody` at `659e170`, `tests/cody-accounts.py` at
`659e170`, and `tests/cody-accounts-live.py` at `3f1ba19` in this repository.
The approved plan and operational evidence remain under ignored `worklog/`.

### Done

- Hardened storage-override rejection for whitespace and quoted key variants.
- Configured the selected homes with links to primary configuration, instructions,
  skills, rules, and plugin cache. Codex credentials, session data, SQLite, and
  plugin-local secret directories remain separate.
- Enrolled both additional homes on each selected work machine in the existing
  private Mesh catalog and backup boundary after explicit authorization.

### Measured

- Five selector tests, 26 update tests, and the environment test pass. The native
  isolated fixture passes on both installed CLIs: login state and actual session/
  SQLite writes stay in the selected home, while shared instructions and skills
  reach the local fixture provider and config changes reach the shared target.
- New history sources completed without errors; monitoring and synchronization
  remained healthy. All four new account homes were signed out at verification.
- A `profile` footer item displayed a selected profile label on CLI 0.166.0 in an
  isolated TUI fixture; CLI 0.165.0 did not display it. This is a label, not proof
  of which account is authenticated.

### Not done

- No user login or external model request was made. External connector credentials
  can require separate authorization. No status-line settings were changed.

### Next

- Sign in to desired homes, verify identity with `/status`, and confirm connectors.

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
