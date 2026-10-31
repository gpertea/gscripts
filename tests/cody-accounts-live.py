#!/usr/bin/env python3
"""Native account-home isolation and shared-config checks without login."""
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import json
import threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

repo = Path(__file__).resolve().parents[1]
assert shutil.which('codex'), 'installed Codex CLI is required'
with tempfile.TemporaryDirectory(prefix='cody-accounts-', dir=Path.home() / '.cache') as temporary:
    root = Path(temporary)
    primary, personal, business = (root / name for name in ('.codex', '.codex-g', '.codex-b'))
    for directory in (primary, personal, business):
        directory.mkdir(mode=0o700)
    class Fixture(BaseHTTPRequestHandler):
        def log_message(self, *args):
            pass

        def do_POST(self):
            self.rfile.read(int(self.headers['Content-Length']))
            response = {'id': 'fixture-response', 'status': 'completed', 'output': [
                {'type': 'message', 'id': 'fixture-message', 'role': 'assistant',
                 'content': [{'type': 'output_text', 'text': 'fictional answer'}]}]}
            events = [{'type': 'response.created', 'response': {'id': 'fixture-response', 'status': 'in_progress', 'output': []}},
                      {'type': 'response.output_item.done', 'output_index': 0, 'item': response['output'][0]},
                      {'type': 'response.completed', 'response': response}]
            data = ''.join('data: ' + json.dumps(item) + '\n\n' for item in events).encode()
            self.send_response(200)
            self.send_header('Content-Type', 'text/event-stream')
            self.send_header('Content-Length', str(len(data)))
            self.end_headers()
            self.wfile.write(data)

    server = ThreadingHTTPServer(('127.0.0.1', 0), Fixture)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    (primary / 'config.toml').write_text(
        'model = "fictional"\nmodel_provider = "fixture"\n'
        '[model_providers.fixture]\nname = "fixture"\n'
        f'base_url = "http://127.0.0.1:{server.server_port}/v1"\n'
        'wire_api = "responses"\nrequires_openai_auth = false\n')
    (primary / 'auth.json').write_text('FICTIONAL_PRIMARY_AUTH_SENTINEL')
    (primary / 'sessions').mkdir()
    (primary / 'sessions' / 'sentinel').write_text('fictional')
    for directory in (personal, business):
        (directory / 'config.toml').symlink_to(primary / 'config.toml')
    env = {'PATH': os.environ['PATH'], 'HOME': temporary, 'HOSTNAME': 'fixture',
           'CODEX_HOME': str(primary), 'CODEX_SQLITE_HOME': str(primary)}
    for entry in ('gcody', 'bcody'):
        result = subprocess.run([str(repo / entry), '--noupdate', 'login', 'status'],
                                env=env, capture_output=True, text=True, timeout=20)
        assert result.returncode == 1 and 'not logged in' in (result.stdout + result.stderr).lower()
    result = subprocess.run([str(repo / 'gcody'), '--noupdate', 'mcp', 'add',
                             'fictional-shared-test', '--', '/bin/true'],
                            env=env, capture_output=True, text=True, timeout=20)
    assert result.returncode == 0, result.stderr[-500:]
    assert all((directory / 'config.toml').is_symlink() for directory in (personal, business))
    assert 'fictional-shared-test' in (primary / 'config.toml').read_text()
    assert 'fictional-shared-test' in (business / 'config.toml').read_text()
    assert (primary / 'auth.json').read_text() == 'FICTIONAL_PRIMARY_AUTH_SENTINEL'
    assert not (personal / 'auth.json').exists() and not (business / 'auth.json').exists()
    assert (primary / 'sessions' / 'sentinel').read_text() == 'fictional'
    try:
        result = subprocess.run([str(repo / 'gcody'), '--noupdate', 'exec',
                                 '--skip-git-repo-check', 'fictional prompt'],
                                env=env, cwd=root, capture_output=True, text=True, timeout=30)
        assert result.returncode == 0, result.stderr[-700:]
        assert list((personal / 'sessions').rglob('*.jsonl'))
        assert sorted(p.name for p in (primary / 'sessions').iterdir()) == ['sentinel']
        assert not list(primary.glob('state_*.sqlite'))
        assert not list(business.rglob('*.jsonl'))
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=3)
print('PASS: native CLI isolates login and session history, and writes through shared config links')
