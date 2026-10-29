#!/usr/bin/env python3
"""Installed CLI/server regression with a local Responses fixture; no account used."""

import asyncio
import fcntl
import importlib.machinery
import importlib.util
import json
import os
import pty
import struct
import tempfile
import termios
import shlex
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
loader = importlib.machinery.SourceFileLoader("proxy", str(REPO / "cody-resume"))
m = importlib.util.module_from_spec(
    importlib.util.spec_from_loader(loader.name, loader)
)
loader.exec_module(m)


async def main():
    with tempfile.TemporaryDirectory(
        prefix="cody-live-", dir=Path.home() / ".cache"
    ) as tmp:
        root = Path(tmp)
        (root / "home").mkdir()
        (root / "workspace").mkdir()
        marker = root / "outside-workspace"
        run_tool = False

        async def model(reader, writer):
            nonlocal run_tool
            try:
                head = await reader.readuntil(b"\r\n\r\n")
                headers = dict(
                    line.split(b":", 1)
                    for line in head.split(b"\r\n")[1:]
                    if b":" in line
                )
                body = await reader.readexactly(
                    int(headers.get(b"content-length", b"0"))
                )
                json.loads(body)
                if run_tool:
                    run_tool = False
                    item = dict(
                        type="function_call",
                        id="fc_test",
                        call_id="call_test",
                        name="exec_command",
                        arguments=json.dumps(
                            dict(
                                cmd="git init -q && printf fixture > probe.txt && git add probe.txt && printf passed > "
                                + shlex.quote(str(marker)),
                                yield_time_ms=1000,
                            )
                        ),
                    )
                else:
                    item = dict(
                        type="message",
                        id="msg_test",
                        role="assistant",
                        content=[dict(type="output_text", text="done")],
                    )
                events = [
                    dict(
                        type="response.created",
                        response=dict(id="resp_test", status="in_progress", output=[]),
                    ),
                    dict(type="response.output_item.done", output_index=0, item=item),
                    dict(
                        type="response.completed",
                        response=dict(
                            id="resp_test",
                            status="completed",
                            output=[item],
                            usage=dict(input_tokens=1, output_tokens=1, total_tokens=2),
                        ),
                    ),
                ]
                payload = "".join(
                    "data: " + json.dumps(e) + "\n\n" for e in events
                ).encode()
                writer.write(
                    b"HTTP/1.1 200 OK\r\nContent-Type: text/event-stream\r\nContent-Length: "
                    + str(len(payload)).encode()
                    + b"\r\nConnection: close\r\n\r\n"
                    + payload
                )
                await writer.drain()
            finally:
                writer.close()

        http = await asyncio.start_server(model, "127.0.0.1", 0)
        port = http.sockets[0].getsockname()[1]
        (root / "home/config.toml").write_text(f"""model = "gpt-5.1-codex"
model_provider = "fixture"
check_for_update_on_startup = false
[model_providers.fixture]
name = "fixture"
base_url = "http://127.0.0.1:{port}/v1"
wire_api = "responses"
requires_openai_auth = false
[features]
unified_exec = true
code_mode = false
[projects."{root/'workspace'}"]
trust_level = "trusted"
""")
        env = dict(os.environ, CODEX_HOME=str(root / "home"))
        log = open(root / "server.log", "wb")
        child = await asyncio.create_subprocess_exec(
            "codex",
            "app-server",
            "--listen",
            "unix://" + tmp + "/server",
            env=env,
            stdout=log,
            stderr=log,
        )
        proxy = m.ResumeProxy("unix://" + tmp + "/server")
        front = await asyncio.start_unix_server(proxy.connect, tmp + "/proxy")
        streams = []
        try:
            for _ in range(100):
                if (root / "server").exists():
                    break
                await asyncio.sleep(0.05)

            async def connect(path):
                r, w = await asyncio.open_unix_connection(path)
                streams.append(w)
                w.write(
                    b"GET / HTTP/1.1\r\nHost: localhost\r\nUpgrade: websocket\r\nConnection: Upgrade\r\nSec-WebSocket-Key: dGhlIHNhbXBsZSBub25jZQ==\r\nSec-WebSocket-Version: 13\r\n\r\n"
                )
                await w.drain()
                await m.header(r)
                await rpc(
                    r,
                    w,
                    "initialize",
                    dict(
                        clientInfo=dict(name="cody-test", version="1"),
                        capabilities=dict(experimentalApi=True),
                    ),
                )
                return r, w

            seq = 0

            async def rpc(r, w, method, params):
                nonlocal seq
                seq += 1
                ident = seq
                w.write(m.encode(dict(id=ident, method=method, params=params), True))
                await w.drain()
                while True:
                    _, body, _ = await asyncio.wait_for(m.frame(r), 30)
                    msg = json.loads(body)
                    if msg.get("id") == ident:
                        if "error" in msg:
                            raise RuntimeError(msg["error"])
                        return msg["result"]

            r, w = await connect(tmp + "/server")
            start = await rpc(
                r,
                w,
                "thread/start",
                dict(
                    cwd=str(root / "workspace"),
                    permissions=":read-only",
                    approvalPolicy="never",
                ),
            )
            thread = start["thread"]["id"]
            print("saved restricted thread", start["sandbox"], flush=True)
            await rpc(
                r,
                w,
                "turn/start",
                dict(
                    threadId=thread,
                    input=[
                        dict(
                            type="text",
                            text="Seed the fixture thread.",
                            text_elements=[],
                        )
                    ],
                ),
            )
            while True:
                _, body, _ = await asyncio.wait_for(m.frame(r), 30)
                if json.loads(body).get("method") == "turn/completed":
                    break
            # a new thread has no persisted rollout; its start response verifies policy
            r, w = await connect(tmp + "/proxy")
            fresh = await rpc(
                r,
                w,
                "thread/start",
                dict(cwd=str(root / "workspace"), sandbox="read-only"),
            )
            assert fresh["sandbox"]["type"] == "dangerFullAccess"
            print("PASS: new thread before any saved rollout")
            # verify the warm-thread path before restarting the isolated server
            r, w = await connect(tmp + "/proxy")
            reply = await rpc(r, w, "thread/resume", dict(threadId=thread))
            assert reply["sandbox"]["type"] == "dangerFullAccess"
            run_tool = True
            await rpc(
                r,
                w,
                "turn/start",
                dict(
                    threadId=thread,
                    input=[
                        dict(
                            type="text",
                            text="Run the fixture command.",
                            text_elements=[],
                        )
                    ],
                ),
            )
            while True:
                _, body, _ = await asyncio.wait_for(m.frame(r), 30)
                if json.loads(body).get("method") == "turn/completed":
                    break
            assert marker.read_text() == "passed"
            marker.unlink()
            await rpc(
                r,
                w,
                "thread/settings/update",
                dict(threadId=thread, permissions=":read-only"),
            )
            r, w = await connect(tmp + "/proxy")
            reply = await rpc(r, w, "thread/resume", dict(threadId=thread))
            assert reply["sandbox"]["type"] == "readOnly"
            fork = await rpc(r, w, "thread/fork", dict(threadId=thread))
            assert fork["sandbox"]["type"] == "dangerFullAccess"
            side = await rpc(r, w, "thread/fork", dict(threadId=thread, ephemeral=True,
                                                     excludeTurns=True))
            assert side["sandbox"]["type"] == "dangerFullAccess"
            source = await rpc(r, w, "thread/resume", dict(threadId=thread))
            assert source["sandbox"]["type"] == "readOnly"
            print(
                "PASS: warm resume/tool write, explicit restriction across reconnect, fork isolation"
            )
            for stream in streams:
                stream.close()
                await stream.wait_closed()
            child.terminate()
            await child.wait()
            (root / "server").unlink(missing_ok=True)
            child = await asyncio.create_subprocess_exec(
                "codex",
                "app-server",
                "--listen",
                "unix://" + tmp + "/server",
                env=env,
                stdout=log,
                stderr=log,
            )
            for _ in range(100):
                if (root / "server").exists():
                    break
                await asyncio.sleep(0.05)
            # a fresh helper/TUI must restore the persisted restricted thread
            run_tool = True
            master, slave = pty.openpty()
            fcntl.ioctl(slave, termios.TIOCSWINSZ, struct.pack("HHHH", 40, 120, 0, 0))
            os.set_blocking(master, False)
            tui = await asyncio.create_subprocess_exec(
                "python3",
                str(REPO / "cody-resume"),
                "codex",
                "unix://" + tmp + "/server",
                "--no-alt-screen",
                "resume",
                thread,
                "Run the fixture command.",
                env=dict(env, TERM="xterm-256color"),
                cwd=root / "workspace",
                stdin=slave,
                stdout=slave,
                stderr=slave,
            )
            os.close(slave)
            screen = bytearray()
            side_tested = False
            try:
                for _ in range(300):
                    try:
                        part = os.read(master, 65536)
                        screen.extend(part)
                        if b"\x1b[6n" in part:
                            os.write(master, b"\x1b[1;1R")
                    except BlockingIOError:
                        pass
                    if marker.exists():
                        if side_tested:
                            break
                        marker.unlink()
                        # wait for the first turn to finish before entering the slash command
                        await asyncio.sleep(1)
                        os.write(master, b'/side')
                        await asyncio.sleep(.4)
                        os.write(master, b'\r')
                        await asyncio.sleep(2)
                        run_tool = True
                        os.write(master, b'Run the side fixture command.')
                        await asyncio.sleep(.4)
                        os.write(master, b'\r')
                        side_tested = True
                    if tui.returncode is not None:
                        break
                    await asyncio.sleep(0.1)
            finally:
                if tui.returncode is None:
                    tui.terminate()
                await tui.wait()
                os.close(master)
            if not marker.exists():
                print("fixture TUI:", screen.decode(errors="replace")[-4000:])
            assert side_tested
            assert b'Failed to start side conversation' not in screen
            assert marker.read_text() == "passed"
            print("PASS: actual /side TUI command and side tool outside-workspace write")
            assert (root / "workspace/.git/index").is_file()
            print(
                "PASS: cold saved-thread TUI resume, next tool git add and outside-workspace write"
            )
        finally:
            for w in streams:
                w.close()
            for t in list(proxy.tasks):
                t.cancel()
            await asyncio.gather(*proxy.tasks, return_exceptions=True)
            front.close()
            http.close()
            await front.wait_closed()
            await http.wait_closed()
            child.terminate()
            await child.wait()
            log.close()
            # only isolated fixture logs; no user credentials or histories are loaded
            if child.returncode not in (0, -15):
                print((root / "server.log").read_text()[-2500:])


if __name__ == "__main__":
    asyncio.run(asyncio.wait_for(main(), 90))
