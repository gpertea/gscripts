#!/usr/bin/env python3
"""Exercise the installed CLI and shared-home swapping with fictional credentials."""

import base64
import fcntl
import json
import os
from pathlib import Path
import shutil
import shlex
import sqlite3
import subprocess
import tempfile
import threading
import time
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

CODY = Path(__file__).resolve().parents[1] / 'cody'
MANAGER = CODY.with_name('cody-account')


def jwt(claims):
    def encode(value):
        return base64.urlsafe_b64encode(json.dumps(value).encode()).rstrip(b'=').decode()
    return encode({'alg': 'none'}) + '.' + encode(claims) + '.fixture'


identity = jwt({'email': 'fixture@example.invalid', 'https://api.openai.com/auth':
                {'chatgpt_account_id': 'fixture-account', 'chatgpt_plan_type': 'plus'}})


class Handler(BaseHTTPRequestHandler):
    authorizations = []
    hold = False
    tool_probe = None
    entered = threading.Event()
    release = threading.Event()

    def log_message(self, *args):
        pass

    def do_POST(self):
        self.rfile.read(int(self.headers.get('content-length', '0')))
        if self.path == '/token':
            data = json.dumps({'id_token': identity, 'refresh_token': 'fictional-rotated',
                               'access_token': jwt({'exp': int(time.time()) + 7200})}).encode()
            self.send_response(200)
            self.send_header('content-type', 'application/json')
            self.send_header('content-length', str(len(data)))
            self.end_headers()
            self.wfile.write(data)
            return
        self.authorizations.append(self.headers.get('authorization', ''))
        if self.hold:
            self.entered.set()
            if not self.release.wait(timeout=20):
                raise RuntimeError('fixture response was not released')
        if Handler.tool_probe:
            command, Handler.tool_probe = Handler.tool_probe, None
            item = {'type': 'function_call', 'id': 'fc_fixture', 'call_id': 'call_fixture',
                    'name': 'exec_command', 'arguments': json.dumps({'cmd': command, 'yield_time_ms': 1000})}
        else:
            item = {'type': 'message', 'id': 'msg_fixture', 'role': 'assistant',
                    'content': [{'type': 'output_text', 'text': 'fixture ok'}]}
        events = [{'type': 'response.created', 'response': {'id': 'resp_fixture', 'status': 'in_progress', 'output': []}},
                  {'type': 'response.output_item.done', 'output_index': 0, 'item': item},
                  {'type': 'response.completed', 'response': {'id': 'resp_fixture', 'status': 'completed', 'output': [item],
                   'usage': {'input_tokens': 1, 'output_tokens': 1, 'total_tokens': 2}}}]
        data = ''.join('data: ' + json.dumps(event) + '\n\n' for event in events).encode()
        self.send_response(200)
        self.send_header('content-type', 'text/event-stream')
        self.send_header('content-length', str(len(data)))
        self.end_headers()
        self.wfile.write(data)


assert shutil.which('codex'), 'the installed CLI is required'
server = ThreadingHTTPServer(('127.0.0.1', 0), Handler)
threading.Thread(target=server.serve_forever, daemon=True).start()
with tempfile.TemporaryDirectory() as temporary:
    root = Path(temporary)
    home, work, state = root / 'home', root / 'workspace', root / 'home/.codex'
    for path in (home, work, state):
        path.mkdir()
    env = {'PATH': os.environ['PATH'], 'HOME': str(home), 'HOSTNAME': 'fixture', 'LANG': 'C.UTF-8'}
    subprocess.run([str(MANAGER), '--home', str(state), '--init'], env=env, check=True)
    (state / 'config.toml').write_text(f'''model = "gpt-5.1-codex"
model_provider = "fixture"
check_for_update_on_startup = false
[model_providers.fixture]
name = "fixture"
base_url = "http://127.0.0.1:{server.server_port}/v1"
wire_api = "responses"
requires_openai_auth = true
''')
    def run(*args, stdin=None):
        return subprocess.run([str(CODY), '--noupdate', *args], input=stdin, env=env,
                              cwd=work, capture_output=True, text=True, timeout=45)

    result = run('login', '--with-api-key', stdin='sk-fictional-work\n')
    assert result.returncode == 0, result.stderr
    result = run('exec', '--skip-git-repo-check', '--json', 'fictional first')
    assert result.returncode == 0, result.stderr
    thread = next(json.loads(line)['thread_id'] for line in result.stdout.splitlines()
                  if json.loads(line).get('type') == 'thread.started')
    result = run('-g', 'login', '--with-api-key', stdin='sk-fictional-g\n')
    assert result.returncode == 0, result.stderr
    result = run('-g', 'exec', 'resume', '--skip-git-repo-check', '--json', thread, 'fictional second')
    assert result.returncode == 0 and thread in result.stdout, result.stderr
    assert Handler.authorizations[0].endswith('sk-fictional-work')
    assert Handler.authorizations[-1].endswith('sk-fictional-g')
    assert (state / 'cody-accounts/default.json').is_file()
    assert not (home / '.codex-g').exists()

    # An active session must prevent a different-account swap.
    Handler.hold = True
    process = subprocess.Popen([str(CODY), '--noupdate', '-g', 'exec', '--skip-git-repo-check',
                                'fictional active'], env=env, cwd=work,
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    try:
        assert Handler.entered.wait(timeout=20), 'fixture session did not start'
        lock = os.open(state / 'cody-accounts/switch.lock', os.O_RDWR)
        try:
            try:
                fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError:
                pass
            else:
                raise AssertionError('active Codex session did not retain the account lock')
        finally:
            os.close(lock)
        before = (state / 'auth.json').read_bytes()
        blocked = run('-b', 'login', '--with-api-key', stdin='sk-fictional-b\n')
        assert blocked.returncode != 0 and 'close it before' in blocked.stderr
        assert (state / 'auth.json').read_bytes() == before
    finally:
        Handler.release.set()
        _, error = process.communicate(timeout=30)
        assert process.returncode == 0, error
    Handler.hold = False
    result = run('exec', 'resume', '--skip-git-repo-check', '--json', thread, 'fictional return')
    assert result.returncode == 0 and thread in result.stdout, result.stderr
    assert Handler.authorizations[-1].endswith('sk-fictional-work')
    assert (state / 'cody-accounts/g.json').is_file()
    assert not (state / 'cody-accounts/default.json').exists()
    # A refreshed ChatGPT credential follows its account through later swaps.
    primary_before = (state / 'auth.json').read_bytes()
    alternate = state / 'cody-accounts/g.json'
    alternate.write_text(json.dumps({'OPENAI_API_KEY': None, 'auth_mode': 'chatgpt',
        'last_refresh': '2020-01-01T00:00:00Z', 'tokens': {'id_token': identity,
        'access_token': jwt({'exp': 1}), 'refresh_token': 'fictional-original',
        'account_id': 'fixture-account'}}))
    refresh_env = {**env, 'CODEX_REFRESH_TOKEN_URL_OVERRIDE': f'http://127.0.0.1:{server.server_port}/token'}
    process = subprocess.Popen([str(CODY), '--noupdate', '-g', 'exec', '--skip-git-repo-check',
                                'fictional refresh'], env=refresh_env, cwd=work,
                               stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    try:
        deadline = time.monotonic() + 15
        refreshed = False
        while time.monotonic() < deadline:
            try:
                value = json.loads((state / 'auth.json').read_text())
                refreshed = value.get('tokens', {}).get('refresh_token') == 'fictional-rotated'
            except (OSError, ValueError):
                pass
            if refreshed or process.poll() is not None:
                break
            time.sleep(.05)
        assert refreshed, 'fictional credential did not refresh'
    finally:
        if process.poll() is None:
            process.terminate()
        process.wait(timeout=10)
    result = run('login', 'status')
    assert result.returncode == 0, result.stderr
    assert (state / 'auth.json').read_bytes() == primary_before
    assert json.loads(alternate.read_text())['tokens']['refresh_token'] == 'fictional-rotated'

    # A tool command must not inherit the session lock and keep it after exit.
    marker = root / 'lock-inherited'
    probe = ("import os,pathlib; links=[]\n"
             "for p in pathlib.Path('/proc/self/fd').iterdir():\n"
             " try: links.append(os.readlink(p))\n"
             " except OSError: pass\n"
             "pathlib.Path(" + repr(str(marker)) + ").write_text(str(any(x.endswith('/switch.lock') for x in links)))")
    Handler.tool_probe = 'python3 -c ' + shlex.quote(probe)
    result = run('exec', '--skip-git-repo-check', 'fictional tool probe')
    assert result.returncode == 0 and marker.read_text() == 'False', result.stderr

    for database in state.glob('*.sqlite'):
        with sqlite3.connect(database) as connection:
            assert connection.execute('PRAGMA integrity_check').fetchone()[0] == 'ok'
    print('PASS: actual CLI shared-home resume, account swaps, concurrent-session block, token refresh, lock lifetime, and SQLite integrity')
server.shutdown()
