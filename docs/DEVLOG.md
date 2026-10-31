# Development log

Newest round first. Paths are relative to the gscripts repository root.

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
