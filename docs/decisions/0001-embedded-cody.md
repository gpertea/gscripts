# Decision 0001: embedded interactive cody

Status: accepted, 2026-10-31.

Source: `cody` at `fd115d7`; tests under `tests/` at that revision.

Interactive sessions launch Codex directly with `--no-daemon --yolo`, including
resume and fork. The process inherits cwd, the Codex home, and terminal context.
Linux browser-display defaults and explicit environment choices remain available.
Administrative commands pass through without injected embedded-mode flags.

Remove display-server routing, the Unix-socket manager, permission-restoration
proxy, and wrapper-driven server version comparison/restart handling. Native
embedded in-session resume supplies the intended permissions without proxy RPCs.
Updates require consent and report only the installed CLI version. Existing
processes continue running their loaded executable until they exit.

Preserve Windows installer, terminal, and environment behavior. An existing
Windows launcher exports a trusted optional color-patch helper root; every launch
uses the canonical wrapper in embedded mode. Only explicit `--patch` loads that
helper, after the canonical updater. The previous Windows wrapper is not executed.

Mesh observes saved history and OS processes independently. It does not require
an app-server connection, cody hook, or wrapper-triggered synchronization.
