#!/usr/bin/env python3
"""Embedded installed-CLI regression using only a local, credential-free fixture."""

import asyncio
import fcntl
import json
import os
from pathlib import Path
import pty
import shlex
import shutil
import struct
import tempfile
import termios


class Terminal:
    def __init__(self, child, master):
        self.child, self.master = child, master
        self.screen = bytearray()
        self.reader = asyncio.create_task(self.read())

    async def read(self):
        while self.child.returncode is None:
            try:
                part = os.read(self.master, 65536)
                if part:
                    self.screen.extend(part)
                    if b'\x1b[6n' in part:
                        os.write(self.master, b'\x1b[1;1R')
            except BlockingIOError:
                pass
            except OSError:
                break
            await asyncio.sleep(.02)

    async def send(self, text):
        os.write(self.master, text.encode())
        await asyncio.sleep(.15)
        os.write(self.master, b'\r')

    async def close(self):
        if self.child.returncode is None:
            self.child.terminate()
        await self.child.wait()
        await self.reader
        os.close(self.master)


async def until(check, label, seconds=20):
    deadline = asyncio.get_running_loop().time() + seconds
    while not check():
        if asyncio.get_running_loop().time() > deadline:
            raise AssertionError('deadline: ' + label)
        await asyncio.sleep(.05)


async def main():
    codex = shutil.which('codex')
    assert codex, 'installed codex is required'
    version = await asyncio.create_subprocess_exec(codex, '--version', stdout=asyncio.subprocess.PIPE)
    stdout, _ = await version.communicate()
    print(stdout.decode().strip(), flush=True)
    with tempfile.TemporaryDirectory(prefix='cody-embedded-') as temporary:
        root = Path(temporary)
        home, store, workspace = [root / name for name in ('home', 'codex', 'workspace')]
        for path in (home, store, workspace):
            path.mkdir()
        env = {'PATH': os.environ['PATH'], 'HOME': str(home), 'CODEX_HOME': str(store),
               'TERM': 'xterm-256color', 'LANG': 'C.UTF-8', 'SHELL': '/bin/bash'}
        requests = []
        tool = None
        sequence = 0

        async def model(reader, writer):
            nonlocal tool, sequence
            try:
                header = await reader.readuntil(b'\r\n\r\n')
                headers = dict(line.lower().split(b':', 1) for line in header.split(b'\r\n')[1:] if b':' in line)
                request = json.loads(await reader.readexactly(int(headers.get(b'content-length', b'0'))))
                requests.append(request)
                sequence += 1
                if tool:
                    command, tool = tool, None
                    item = dict(type='function_call', id=f'fc_{sequence}', call_id=f'call_{sequence}',
                                name='exec_command', arguments=json.dumps({'cmd': command, 'yield_time_ms': 1000}))
                else:
                    item = dict(type='message', id=f'msg_{sequence}', role='assistant',
                                content=[dict(type='output_text', text='fixture done')])
                events = [dict(type='response.created', response=dict(id=f'resp_{sequence}', status='in_progress', output=[])),
                          dict(type='response.output_item.done', output_index=0, item=item),
                          dict(type='response.completed', response=dict(id=f'resp_{sequence}', status='completed', output=[item],
                               usage=dict(input_tokens=1, output_tokens=1, total_tokens=2)))]
                data = ''.join('data: ' + json.dumps(event) + '\n\n' for event in events).encode()
                writer.write(b'HTTP/1.1 200 OK\r\nContent-Type: text/event-stream\r\nContent-Length: ' +
                             str(len(data)).encode() + b'\r\nConnection: close\r\n\r\n' + data)
                await writer.drain()
            finally:
                writer.close()

        server = await asyncio.start_server(model, '127.0.0.1', 0)
        port = server.sockets[0].getsockname()[1]
        (store / 'config.toml').write_text(f'''model = "gpt-5.1-codex"
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
[projects.{json.dumps(str(workspace))}]
trust_level = "trusted"
''')
        terminals = []

        def records():
            found = {}
            for path in (store / 'sessions').rglob('*.jsonl'):
                rows = []
                for line in path.read_text().splitlines():
                    try:
                        rows.append(json.loads(line))
                    except ValueError:
                        pass
                meta = next((r['payload'] for r in rows if r.get('type') == 'session_meta'), None)
                if meta:
                    found[meta['id']] = (path, rows)
            return found

        def contexts(thread):
            return [r['payload'] for r in records()[thread][1] if r.get('type') == 'turn_context']

        def has_text(thread, text):
            return text in json.dumps(records()[thread][1])

        async def launch(*arguments, restricted=False):
            master, slave = pty.openpty()
            fcntl.ioctl(slave, termios.TIOCSWINSZ, struct.pack('HHHH', 42, 140, 0, 0))
            os.set_blocking(master, False)
            if restricted:
                argv = [codex, '--no-daemon', '--no-alt-screen', '--sandbox', 'read-only',
                        '--ask-for-approval', 'never']
            else:
                # exercise cody's defaults; no permission or daemon override is supplied.
                argv = [str(Path(__file__).resolve().parents[1] / 'cody'), '--noupdate', '--no-alt-screen']
            child = await asyncio.create_subprocess_exec(*argv, *arguments, env=env, cwd=workspace,
                                                        stdin=slave, stdout=slave, stderr=slave)
            os.close(slave)
            terminal = Terminal(child, master)
            terminals.append(terminal)
            await asyncio.sleep(1)
            return terminal

        async def seed(name):
            existing = set(records())
            terminal = await launch('fictional seed ' + name, restricted=True)
            await until(lambda: len(set(records()) - existing) == 1, 'new seed ' + name)
            thread = next(iter(set(records()) - existing))
            await until(lambda: has_text(thread, 'fixture done'), 'seed answer ' + name)
            await terminal.close()
            policy = contexts(thread)[-1]
            assert policy['sandbox_policy']['type'] == 'read-only'
            assert policy['approval_policy'] == 'never'
            print('PASS: saved read-only seed ' + name, flush=True)
            return thread

        async def write_probe(terminal, thread, label):
            nonlocal tool
            marker = root / ('outside-' + label)
            text = 'fictional permission probe ' + label
            tool = 'git init -q && printf fixture > probe.txt && git add probe.txt && printf passed > ' + shlex.quote(str(marker))
            await terminal.send(text)
            await until(lambda: marker.exists(), 'unrestricted tool ' + label)
            await until(lambda: has_text(thread, text), 'saved target ' + label)
            await asyncio.sleep(.3)
            policy = contexts(thread)[-1]
            assert policy['sandbox_policy']['type'] == 'danger-full-access'
            assert policy['approval_policy'] == 'never'
            assert policy['permission_profile']['type'] == 'disabled'
            print('PASS: full-access target and actual tool writes ' + label, flush=True)
            assert marker.read_text() == 'passed'
            assert (workspace / '.git/index').is_file()

        try:
            first, second = await seed('alpha'), await seed('beta')
            third, fourth = await seed('gamma'), await seed('delta')
            existing = set(records())
            terminal = await launch('fictional fresh thread')
            await until(lambda: len(set(records()) - existing) == 1, 'fresh full-access thread')
            fresh = next(iter(set(records()) - existing))
            await until(lambda: has_text(fresh, 'fixture done'), 'fresh answer')
            await write_probe(terminal, fresh, 'new')
            await terminal.send('/resume')
            await asyncio.sleep(1)
            assert b'Resume a previous session' in terminal.screen
            # use the saved seed's distinctive text in the actual in-session picker.
            os.write(terminal.master, b'fictional seed alpha')
            await asyncio.sleep(.5)
            os.write(terminal.master, b'\r')
            await asyncio.sleep(1)
            await write_probe(terminal, first, 'in-session-resume')
            assert not has_text(fresh, 'fictional permission probe in-session-resume')
            await terminal.send('/resume')
            await asyncio.sleep(1)
            os.write(terminal.master, b'fictional seed beta')
            await asyncio.sleep(.5)
            os.write(terminal.master, b'\r')
            await asyncio.sleep(1)
            await write_probe(terminal, second, 'second-in-session-resume')
            assert not has_text(first, 'fictional permission probe second-in-session-resume')
            await terminal.close()
            terminal = await launch('resume', third, 'fictional explicit resume')
            await until(lambda: has_text(third, 'fictional explicit resume'), 'explicit resume target')
            await write_probe(terminal, third, 'resume')
            await terminal.close()
            existing = set(records())
            terminal = await launch('fork', fourth, 'fictional persisted fork')
            await until(lambda: len(set(records()) - existing) == 1, 'persisted fork identity')
            forked = next(iter(set(records()) - existing))
            await until(lambda: has_text(forked, 'fictional persisted fork'), 'fork prompt')
            await write_probe(terminal, forked, 'fork')
            await terminal.close()
            assert contexts(fourth)[-1]['sandbox_policy']['type'] == 'read-only'
            assert len(set(records())) == len(existing) + 1
            assert not list(store.rglob('auth.json'))
            print(json.dumps({'new': True, 'resume': True, 'fork': True, 'in_session_resume': True, 'in_session_switches': 2,
                              'local_fixture_requests': len(requests), 'external_model_requests': 0}), flush=True)
        finally:
            for terminal in terminals:
                if terminal.child.returncode is None:
                    print('failure TUI', terminal.screen.decode(errors='replace')[-4000:], flush=True)
                    await terminal.close()
            server.close()
            await server.wait_closed()


if __name__ == '__main__':
    asyncio.run(asyncio.wait_for(main(), 100))
