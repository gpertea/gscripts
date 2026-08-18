# gscripts
This is a mishmash of utility scripts I have been using in my Genomics/Bioinformatics work.  Many of them may not make much sense outside of my work environment, some of them were project specific and I probably already forgot why I wrote them.. 

## Named VNC sessions

This repository owns the generic VNC session layer. Any Linux host can manage
named desktops, view its own desktops, and be viewed from Linux, macOS, or
Windows. `lin-browser-use` consumes this layer but does not own or install it.

`vnc-get-display`, `vnc-session`, `vnc-self`, `vnc-session-view`, `vnc-lan`,
`vnc-win-lan`, `vnc-host`, `flameshot-vnc`, and `jwm-vnc-session.xml` are
canonical here. Install the Linux entry points as links; do not maintain copies
in `~/bin`.

```bash
~/gscripts/install-vnc-tools install
~/gscripts/install-vnc-tools check
```

Installation backs up an existing regular managed file before replacing it
with a canonical link. It does not alter host-specific `vnc-HOST` commands.
The JWM configuration is a host-specific regular copy, not a repository link.
Installation refreshes `~/.config/vnc-session/jwm.xml` from `~/.jwmrc-vnc` when
present, otherwise from `~/.jwmrc`; the minimal repository configuration is
used only when neither file exists. `install-vnc-tools check` detects a stale
copy.

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
`vnc-HOST` to `vnc-host`. Its defaults match the legacy behavior: direct LAN
and display `:1`.

```bash
mv ~/bin/vnc-srv16 ~/bin/vnc-srv16.legacy
ln -s ~/gscripts/vnc-host ~/bin/vnc-srv16
vnc-srv16
vnc-srv16 2
vnc-srv16 --tunnel 2
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

On macOS, `vnc-lan linwks34 --session browser-a` creates the same tunnel and
opens TurboVNC. Legacy direct-LAN usage remains
`vnc-lan linwks34 [vnc-host] [geometry]`.

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
sockets and local ports are keyed by the resolved remote display, so `:2`,
`:3`, and named sessions can coexist.

Before opening the Windows viewer, `vnc-win-lan` normalizes the server's live
desktop name to `SHORT_HOST:DISPLAY` with `vncconfig`. This produces stable
window and taskbar titles such as `gglin:1`, `gglin:2`, and `srv16:3`, even
when the server originally advertised a session name or a fully-qualified
hostname. Changing the title does not restart the VNC desktop.

Windows SSH tunnel control checks are bounded, and tunnel master PIDs are
recorded below the per-user control directory. If a retained master becomes
unresponsive after sleep or a network change, `vnc-win-lan` replaces only that
local forwarding process and its control socket. New masters use SSH server
keepalives so dead connections are normally removed automatically. Replacing a
tunnel does not stop the remote VNC desktop or any process running inside it.

`vnc-host` automatically dispatches to `vnc-win-lan` under MSYS2/Cygwin and to
`vnc-session-view` on Linux. Small Windows host wrappers can set
`VNC_HOST_ALIAS`; off-LAN wrappers can additionally set
`VNC_HOST_DEFAULT_TRANSPORT=--tunnel`. Explicit `--direct` or `--tunnel`
arguments always override the wrapper default.

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
