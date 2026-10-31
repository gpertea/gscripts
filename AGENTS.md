# gscripts guidance

Read [STATUS.md](STATUS.md) for current Codex-wrapper maintenance state and
[docs/DEVLOG.md](docs/DEVLOG.md) for verified rounds. The repository also contains
independent utilities; keep changes scoped to the requested tool.

- Save approved implementation plans and operational evidence under `worklog/`.
  Keep new plans, runtime paths, session records, logs, and credentials out of Git.
- Test cody with isolated command doubles and the credential-free loopback fixture.
  Tests must not update installed software or call external models.
- Interactive cody sessions use `--no-daemon`; preserve native argument handling,
  environment inheritance, platform support, and update consent.
- Do not stop existing Codex sessions or servers as part of wrapper maintenance.
- Commit implementation first, then current STATUS and a new DEVLOG entry in a
  separate commit. Push only when authorized. Preserve unrelated worktree changes.
