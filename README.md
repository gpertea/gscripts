# gscripts
This is a mishmash of utility scripts I have been using in my Genomics/Bioinformatics work.  Many of them may not make much sense outside of my work environment, some of them were project specific and I probably already forgot why I wrote them.. 

## Codex sessions by display

On Linux, `cody` recognizes numeric displays such as `:1`, `:2`, and `:2.0`.
It preserves the terminal's existing `TERM`, `COLORTERM`, and color-control
variables, while exporting the matching `LIN_BROWSER_USE_AUTH_DISPLAY` for
browser tools. It also supplies missing Xauthority and runtime bus defaults.

With a local numeric display, `cody` starts or reuses one app-server per live
X/VNC session. `cody --attach :2` attaches to an existing display server;
`--server-status`, `--server-endpoint`, and `--server-stop` manage it.

With `DISPLAY` unset or empty, `cody` starts the regular Codex CLI after its
usual update check. This path uses normal default app-server behavior and
requires no display manager, VNC, Xauthority, or desktop runtime setup. On
JHPCE compute and transfer nodes, it still loads the node module and puts
`~/.local/bin` first on PATH. Display-server reminders are omitted.

JHPCE's NFS home can prevent Codex from removing temporary helpers while
their lock files are open. For the installed workaround,
`~/.codex/tmp/arg0` links to `/tmp/cody-codex-arg0-UID` (the numeric user ID).
`cody` recreates this private directory on each node when that link is present.
`CODEX_HOME`, settings, authentication, and session history stay in home.

### Codex updates

Interactive `cody` launches check for a newer stable release before opening
the UI and ask before installing it. The check uses npm on Unix and
PowerShell on Windows, with a five-second request timeout. An unavailable
check still permits launching the installed CLI. When cody handles the
check, it disables Codex's duplicate startup update prompt for that launch.

On Windows/MSYS2, `cody` prefers the standalone CLI at
`%LOCALAPPDATA%\Programs\OpenAI\Codex\bin\codex.exe` and updates it with
OpenAI's Windows installer. The npm executable on PATH is only a fallback;
its version can differ from the CLI that `cody` launches. A successful check
is silent when the selected CLI is current. `--help` and `--version` skip
the automatic check, as do launches without terminal stdin and stdout.

Windows CLI launches enable unified execution for direct MSYS2 Bash selection,
disable login-shell execution, select inline mode for terminal scrollback,
and bind F10 to return from async questions. Interactive mintty launches set
`FORCE_COLOR=3` when neither `NO_COLOR` nor `FORCE_COLOR` is explicitly set.
The native CLI is executed directly, including in mintty, without winpty.
The Codex keymap includes `ctrl-enter` as a newline alias because Windows
ConPTY represents LF input as Ctrl+Enter in console records. The configuration
and input-path evidence are documented in `win-config/codex-cli/README.md`.
These overrides are passed to the CLI, rather than written to desktop settings.
The CLI inherits the MSYS2 tool PATH. Only the Windows installer receives a
Windows-first PATH for native `tar`; native PowerShell calls clear inherited
`PSModulePath` for their own process.

```bash
cody                       ## check, ask to update if newer, then launch
cody --noupdate             ## skip both wrapper and Codex startup checks
cody --update               ## update explicitly, report versions, then exit
cody --noupdate --docs      ## combine with other leading wrapper options
```

After a successful update, cody reports the installed CLI version, running
app-server version/status, daemon package version, and running Linux Codex
process versions. Executables already running before an installation can
remain older; their versions are read through `/proc`, not inferred from the
installed package.

If an app-server is running, cody offers to update its package to the installed
CLI and restart it. The default answer is no. Accepting stops the server,
pins the current CLI package with `app-server daemon update --from-cli --yes`,
and starts the updated server; active server work is interrupted. The version
report is then repeated. If package replacement fails, cody leaves the server
stopped and returns the error instead of starting an uncertain version.

An absent server is not started. Interactive CLI sessions are never killed.
Noninteractive launches do not run the wrapper's automatic check; explicit
noninteractive `--update` reports versions without restarting a server.
Remote clients can subsequently start a server again after `killall codex`;
stopping processes does not disconnect those clients.

## Named VNC sessions

This repository owns the generic VNC session layer. Any Linux host can manage
named desktops, view its own desktops, and be viewed from Linux, macOS, or
Windows. `lin-browser-use` consumes this layer but does not own or install it.

`vnc-get-display`, `vnc-session`, `vnc-self`, `vnc-session-view`, `vnc-lan`,
`vnc-win-lan`, `vnc-host`, `flameshot-vnc`, `vnc-xstartup`,
`vnc-secrets-forward`, `vnc-keyring-test`, and `jwm-vnc-session.xml` are
canonical here. Install the Linux entry points as links; do not maintain copies
in `~/bin`.

```bash
~/gscripts/install-vnc-tools install --dry-run
~/gscripts/install-vnc-tools install
~/gscripts/install-vnc-tools check
```

`install --dry-run` lists every link or copy it would change and changes
nothing. Installation moves an existing regular managed file into
`${VNC_BRU_BACKUP_DIR:-~/.local/state/vnc-bru-setup/backups/<run>}/<absolute path>`
before replacing it, prints each `>>> BACKUP:` path, and appends it to the
`MANIFEST` file there. `lin-config/browser-use/install-vnc-bru.sh` sets
`VNC_BRU_BACKUP_DIR` so one run keeps all backups together. Replaced symlinks
are reported with their previous target. It does not alter host-specific
`vnc-HOST` commands.
### Desktop startup, JWM, and D-Bus

Every managed desktop starts the same way. `vnc-session` sets only the
per-display parts, then runs the user's desktop startup:

1. A private D-Bus for the display (`dbus-run-session` with a generated
   `$XDG_RUNTIME_DIR/vnc-session/dbus-N/bus.conf`). Single-instance
   applications such as Flameshot therefore stay on their own display.
2. `DISPLAY`, `XAUTHORITY`, `RXVT_SOCKET` (one urxvtd per display), and
   `VNC_SESSION_JWM_CONFIG`; the session environment file for `vnc-session run`.
3. `~/.vnc/xstartup` when executable, otherwise the minimal `vnc-xstartup`.
   The startup loads `~/.profile` (so menu entries get the login `PATH`),
   merges `~/.Xresources`, starts `urxvtd` and `vncconfig`, and runs
   `jwm -f "$VNC_SESSION_JWM_CONFIG"`. It must not set up D-Bus or `DISPLAY`.

`~/.jwmrc-vnc` configures every VNC desktop; `~/.jwmrc` belongs only to a real
local `:0` desktop. Without `~/.jwmrc-vnc` the minimal `jwm-vnc-session.xml` is
used. The old `~/.config/vnc-session/jwm.xml` copy is no longer read;
`install-vnc-tools` reports it so it can be removed once desktops started
before this change are restarted.

The user's single gnome-keyring lives on the per-user bus
(`$XDG_RUNTIME_DIR/bus`), which a private bus cannot reach. Each private bus
therefore activates `vnc-secrets-forward` for `org.freedesktop.secrets` and
`org.gnome.keyring`. It relays Secret Service calls, replies, and signals to
the user bus, so Chromium, VS Code, Python `keyring`, and `secret-tool` work
on every desktop. Before relaying an unlock prompt it points the user bus
activation `DISPLAY` at its own desktop, so the dialog appears where it was
requested. It needs `python3-gi`; without it keyring calls fail at once
instead of waiting for the D-Bus timeout.

Check a desktop's keyring access from any terminal:

```bash
vnc-keyring-test                 # the desktop this terminal runs in
vnc-keyring-test :1              # or a running desktop by display or name
vnc-keyring-test --quick test    # no test item; nothing written to the keyring
vnc-keyring-test --unlock        # lock, then unlock through the dialog
```

It checks that a managed desktop has its own bus, that the keyring answers
within 2 s, and that a forwarder serves the display; by default it also stores,
reads, and deletes a test item while the keyring is unlocked. `--unlock`
verifies that the unlock dialog appears on the tested display. A desktop
started before this design fails with "keyring unreachable"; restart it.

Start desktops with `vnc-session start`; a desktop started from a terminal
in another desktop does not inherit that desktop's bus or terminal socket.

`flameshot-vnc` supports any numeric `DISPLAY`. Flameshot uses both D-Bus and a
Qt local socket for single-instance routing, neither of which is keyed by
`DISPLAY`. The wrapper gives every display a separate `TMPDIR`. It reuses the
private D-Bus already supplied by a managed `vnc-session`; for a legacy desktop
sharing the `:0` user bus, it creates a display-specific private bus. The
historical `:1` socket paths remain compatible with an already-running daemon.

Managed sessions also start an independent urxvtd for every display. JWM and
`vnc-session run` inherit a display-specific `RXVT_SOCKET` below
`$XDG_RUNTIME_DIR/vnc-session`, so stopping one display cannot terminate or
misroute terminals in another display. The runtime directory is cleared on
reboot, avoiding stale sockets from prior boots.

Each Linux VNC host requires Bash, TigerVNC, JWM, rxvt-unicode,
`dbus-run-session`, `flock`, `ss`, `pgrep`, `setsid`, `xrdb`, `xsetroot`, and a
private `~/.vnc/passwd`.
On Debian or Ubuntu install `bash`, `tigervnc-standalone-server`,
`tigervnc-tools`, `tigervnc-viewer`, `jwm`, `dbus-daemon`, `util-linux`,
`iproute2`, `procps`, `x11-xserver-utils`, and `rxvt-unicode`. Then verify
without starting a desktop:

```bash
vnc-session doctor
```

Create and inspect independent sessions on any Linux host:

```bash
vnc-get-display --session browser-a
vnc-get-display --list
vnc-get-display --status browser-a
```

The no-option `vnc-get-display [GEOMETRY]` interface remains compatible with
old `vnc-HOST` scripts: it always ensures direct-LAN display `:1` and prints the
scalar display number `1`. It never substitutes a named session on `:2` or
higher. Existing copied wrappers therefore continue working after the remote
host installs these tools.

To replace copied wrappers with one maintained implementation, link
`vnc-HOST` to `vnc-host`. Linux and macOS default to an SSH tunnel and display
`:1`; Windows retains the legacy direct-LAN default.

```bash
mv ~/bin/vnc-srv16 ~/bin/vnc-srv16.legacy
ln -s ~/gscripts/vnc-host ~/bin/vnc-srv16
vnc-srv16
vnc-srv16 2
vnc-srv16 --direct 2
```

Open local sessions independently. With no argument, `vnc-self` always targets
exact display `:1`; it never allocates or selects a named session:

```bash
vnc-self
vnc-self browser-a
vnc-self browser-b
vnc-self :2
vnc-self --list
```

An exact display such as `:2` attaches when running, restarts its retained
managed session when stopped, or creates `display-2` on exactly `:2`. Every
managed desktop uses JWM, a private D-Bus, and a loopback-only VNC listener.

Run GUI commands in a named desktop's recorded environment. `bru` derives its
profile and loopback CDP port from that desktop's `DISPLAY`; it also works when
launched normally on `:0`:

```bash
DISPLAY=:0 bru
vnc-session run browser-a bru
```

From another Linux host, use an explicit SSH target. Named sessions are
loopback-only and the viewer creates an SSH tunnel:

```bash
vnc-session-view linwks34
vnc-session-view linwks34:2
vnc-session-view linwks34:browser-a
```

With no suffix, the target is display `:1`. An exact numeric suffix ensures
that display exists on the target host. For example,
`vnc-session-view linwks34:2` creates retained session `display-2` when `:2` is
free, restarts it when stopped, or attaches when already running. A text suffix
selects a retained session name.

Tunnel mode is the default. SSH starts or discovers the desktop and carries the
VNC stream from its loopback listener:

```bash
vnc-session-view linwks34:2
vnc-session-view --tunnel linwks34:2
```

Direct LAN mode still uses SSH for session initiation, but exposes the managed
VNC listener on IPv4 `0.0.0.0` and connects the viewer directly:

```bash
vnc-session-view --direct linwks34:2
```

For a new or stopped session, `--direct` records LAN mode without destroying
processes. For a running loopback-only session, it refuses and leaves the
desktop unchanged because changing the bound listener requires replacing
Xtigervnc. Only the following explicit destructive request permits that:

```bash
vnc-session-view --direct --force-restart linwks34:2
```

Before restarting, the command prints warnings that identify the host, session,
display, listener change, and loss of every GUI process in that desktop. LAN
mode then persists. Later tunneled viewers reuse it without restarting or
disabling direct access.

Returning a running session to loopback-only mode has the same guard:

```bash
vnc-session ensure --tunnel :2                 # refuses while running
vnc-session ensure --tunnel --force-restart :2 # explicit destructive change
```

Stopping a running desktop terminates Xtigervnc and every GUI process inside
it. A normal stop therefore refuses and prints the exact destructive command:

```bash
vnc-session stop :2                  # refuses while running
vnc-session stop --force :2          # explicit destructive stop
vnc-self --stop --force :2           # same operation on this host
vnc-session-view --stop --force srv16:2 # same operation on another host
```

`vnc-HOST --stop --force 2` is the equivalent host-alias form. Stopping keeps
the managed display assignment, so a later open recreates the same desktop on
the same display. An explicit display also permits guarded stopping of an older
unmanaged desktop such as `:1`. A stopped desktop may be checked or stopped
again without `--force`.

Direct mode uses VNC password authentication but does not encrypt screen,
keyboard, pointer, or clipboard traffic. Use it only on trusted routed LANs.

On Linux, `vnc-session-view` allocates an available loopback port for each
new tunnel. Different hosts using the same remote display can be viewed
concurrently. Existing tunnels reuse the port from the verified SSH process
arguments, including tunnels created by older versions of the helper.
The viewer receives `localhost::PORT`, which specifies a TCP port explicitly.

`VNC_SESSION_LOCAL_DISPLAY` remains an optional fixed local display: `11`
requests port `5911`. A conflicting override or occupied fixed port produces
an error without replacing an active tunnel. Automatic allocation retries up
to five bind races; authentication and network failures are not retried.
Linux tunnel management requires `perl` (core `IO::Socket::INET`), `flock`,
`timeout`, and `/proc`. Locks, PID records, and socket cleanup are scoped to
the SSH endpoint and display. New tunnels use SSH server keepalives.

On macOS, `vnc-lan linwks34 --session browser-a` creates an SSH tunnel and
opens TurboVNC. Legacy direct-LAN usage remains
`vnc-lan linwks34 [vnc-host] [geometry]`.
The macOS tunnel port policy is unchanged by the Linux update.

On Windows/MSYS2, `vnc-win-lan` now accepts the same host/session target form
as `vnc-session-view` while retaining the established
`C:\util\vnc\tigervncviewer.exe` invocation:

```bash
vnc-win-lan gglin:2                  # direct LAN data connection
vnc-win-lan --direct gglin:browser-a
vnc-win-lan --tunnel srv16:2         # per-display SSH tunnel
vnc-win-lan --status srv16:2
vnc-win-lan --stop --force srv16:2
```

The legacy `vnc-win-lan HOST [GEOMETRY]` and
`vnc-win-lan HOST --session NAME [GEOMETRY]` forms remain accepted. Direct
mode resolves the VNC endpoint from `ssh -G HOST`; `VNC_DIRECT_HOST` and
`VNC_DIRECT_PORT_BASE` can override forwarded LAN endpoints. Tunnel control
sockets are keyed by the SSH endpoint and resolved remote display. Windows
allocates an available loopback port for each new tunnel, so multiple hosts
using `:1`, additional displays, and named sessions can coexist. A retained
tunnel reuses the port recorded in its live SSH process arguments. The former
`16500 + display` default is no longer used.

`VNC_LOCAL_PORT` remains an optional fixed-port request. An occupied explicit
port fails with a diagnostic; automatic allocation retries up to five times
if another process claims a selected port before SSH binds it. Authentication
and network failures are reported without allocation retries. Tunnel setup
uses `flock` to serialize launches for the same SSH endpoint and display.

Before opening the Windows viewer, `vnc-win-lan` normalizes the server's live
desktop name to `SHORT_HOST:DISPLAY` with `vncconfig`. This produces stable
window and taskbar titles such as `gglin:1`, `gglin:2`, and `srv16:3`, even
when the server originally advertised a session name or a fully-qualified
hostname. Changing the title does not restart the VNC desktop.

The title helper probes the running display for compatibility, trying the
default `vncconfig`/`tigervncconfig` first, then preserved Ubuntu utilities
(`/usr/bin/tigervncconfig.ubuntu-*` and `/usr/bin/vncconfig.ubuntu-*`). This
supports older desktops that remain running after a TigerVNC upgrade changes
the configuration extension from `VNC-EXTENSION` to `TIGERVNC`. Installed
package versions cannot identify the version of an already-running server.
An already-correct title is left alone; changes are verified by reading the
title back because some utility versions return success even on failure.
If no compatible utility is available, the launcher reports the failed probes
and stops before opening a viewer with an unverified singleton title.

The stable title also acts as a Windows-side singleton key. If a matching
TigerVNC window already exists, `vnc-win-lan` restores and focuses that viewer
instead of opening a second connection to the same remote display. It warns
when multiple matching windows already exist but never terminates a viewer
automatically. This avoids the indefinite multi-client pointer lock present in
TigerVNC server 1.13.1 when one viewer retains a pressed-button state.

Windows SSH tunnel control checks are bounded, and tunnel master PIDs are
recorded below the per-user control directory. If a retained master becomes
unresponsive after sleep or a network change, `vnc-win-lan` replaces only that
local forwarding process and its control socket. New masters use SSH server
keepalives so dead connections are normally removed automatically. Replacing a
tunnel does not stop the remote VNC desktop or any process running inside it.
Cleanup removes only the exact control socket expanded by SSH for that
endpoint and display. It does not remove another host's socket for the same
display number. The port allocator is `vnc-win-port.ps1`; keep it alongside
`vnc-win-lan` when updating the Windows client.

Regression checks:

```bash
bash tests/vnc-win-title.sh
bash tests/vnc-win-tunnel.sh
bash tests/vnc-session-view-tunnel.sh
```

`vnc-host` automatically dispatches to `vnc-win-lan` under MSYS2/Cygwin and to
`vnc-session-view` elsewhere. Its platform default follows that viewer:
direct LAN on Windows and an SSH tunnel elsewhere. Small host wrappers can set
`VNC_HOST_ALIAS` or override the default with `VNC_HOST_DEFAULT_TRANSPORT`.
Explicit `--direct` or `--tunnel` arguments always take precedence.

### Diagnose and recover VNC input

`vnc-fix HOST` is read-only by default and reports the active display, viewer
connections, TigerVNC input device and master state, JWM state, and only recent
relevant server-log errors. Use `--full` when the complete input and log detail
is needed.

```bash
vnc-fix gglin
vnc-fix --release-buttons gglin
vnc-fix --replace-master --force gglin
vnc-fix --guard-master gglin
```

The replacement-master action is an explicit fallback when synthetic button
releases cannot clear a wedged pointer. The guard keeps newly created X clients
on that replacement master for the lifetime of the current VNC server. Avoid
`--danger-disable-pointer`: disabling TigerVNC's input device can crash the
Ubuntu 24.04 TigerVNC 1.13.1 server and terminate X-attached processes.

## Update a VNC network

On every Linux host that may serve or view VNC:

```bash
git -C ~/gscripts pull --ff-only
~/gscripts/install-vnc-tools install
~/gscripts/install-vnc-tools check
vnc-session list
```

The installer does not stop or restart any VNC desktop. Existing sessions keep
their current processes and listener bindings. On hosts that use `bru`, update
the separate `lin-browser-use` skill after updating `codex-kit`.

On macOS, update the `gscripts` checkout and continue using `vnc-lan`. On
Windows/MSYS2, update the checkout and continue using `vnc-win-lan`. Existing
SSH aliases, proxy jumps, VNC password files, and host routing remain external
configuration.

## Claude Code sessions

`claudy` launches Claude Code with `--dangerously-skip-permissions` on every
host. It always runs the native per-user install at `~/.local/bin/claude`,
installing it with `https://claude.ai/install.sh` when missing. On JHPCE
compute and transfer nodes it loads the node module (for npm/npx MCP servers)
and then puts `~/.local/bin` first on PATH, because that module ships an old
claude that cannot be upgraded in place. There is no display server; with a
local numeric `DISPLAY` it only exports `LIN_BROWSER_USE_AUTH_DISPLAY`.

Launches compare the installed version with the `latest` release channel
(five-second timeout) and run `claude update` without asking when it is
newer. A failed or unavailable check still launches the installed CLI.
`--help` and `--version` skip the check. Windows/MSYS2 is not supported.

```bash
claudy                     ## install if missing, update if newer, then launch
claudy --noupdate          ## skip the wrapper update check
claudy --update            ## update explicitly, report versions, then exit
claudy -c                  ## other arguments pass through to claude
```
