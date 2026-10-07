# Approved cody display integration

Wire cody to one app-server per live local X/VNC session. The first local
invocation may start it; later clients attach to it. Provide attach-only remote
access, status, endpoint lookup, and explicit stop. Preserve update handling
and reject incompatible server configuration rather than silently changing
an existing display server. Test the actual cody wrapper in managed VNC/bru
sessions with tmux, plus concurrency and failure recovery in isolated sessions.
