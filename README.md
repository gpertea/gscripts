# gscripts
This is a mishmash of utility scripts I have been using in my Genomics/Bioinformatics work.  Many of them may not make much sense outside of my work environment, some of them were project specific and I probably already forgot why I wrote them.. 

## Named VNC sessions

`vnc-get-display`, `vnc-session`, `vnc-session-view`, `vnc-lan`, `vnc-win-lan`, and `jwm-vnc-session.xml` are the canonical VNC session tools. Link the scripts into `~/bin`; do not maintain separate copies there.

```bash
ln -s ~/gscripts/vnc-get-display ~/bin/vnc-get-display
ln -s ~/gscripts/vnc-session ~/bin/vnc-session
ln -s ~/gscripts/vnc-session-view ~/bin/vnc-session-view
mkdir -p ~/.config/vnc-session
ln -s ~/gscripts/jwm-vnc-session.xml ~/.config/vnc-session/jwm.xml
```

The Linux server requires TigerVNC, JWM, `dbus-run-session`, `flock`, `ss`, `pgrep`, `setsid`, `xrdb`, and `xsetroot`, plus a private `~/.vnc/passwd`. On Debian or Ubuntu these are supplied by `tigervnc-standalone-server`, `tigervnc-tools`, `jwm`, `dbus-daemon`, `util-linux`, `iproute2`, `procps`, and `x11-xserver-utils`.

Create and inspect independent sessions on the server:

```bash
vnc-get-display --session browser-a
vnc-get-display --list
vnc-get-display --status browser-a
vnc-get-display --stop browser-a
```

On macOS, `vnc-lan srv16 --session browser-a` creates an SSH tunnel to the loopback-only named session and opens TurboVNC. Legacy direct-LAN usage remains `vnc-lan srv16 [vnc-host] [geometry]`.

On Windows/MSYS2, `vnc-win-lan srv16 --session browser-a` provides the same named-session behavior while retaining the established `C:\util\vnc\tigervncviewer.exe` invocation. Legacy direct-LAN usage is `vnc-win-lan srv16 [geometry]`.
