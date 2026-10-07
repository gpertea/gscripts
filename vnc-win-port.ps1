$ErrorActionPreference = 'Stop'

# Ask Windows for an available loopback port. SSH must still confirm the bind
# because another process can claim it after this listener is released.
$listener = [Net.Sockets.TcpListener]::new([Net.IPAddress]::Loopback, 0)
try {
    $listener.Server.ExclusiveAddressUse = $true
    $listener.Start()
    $listener.LocalEndpoint.Port
}
finally {
    $listener.Stop()
}
