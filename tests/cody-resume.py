#!/usr/bin/env python3
"""Exercise the real proxy over private Unix WebSocket connections."""
import asyncio
import importlib.util
import importlib.machinery
import json
from pathlib import Path
import tempfile
import unittest

REPO = Path(__file__).resolve().parents[1]
loader = importlib.machinery.SourceFileLoader('cody_resume', str(REPO / 'cody-resume'))
proxy_module = importlib.util.module_from_spec(importlib.util.spec_from_loader(loader.name, loader))
loader.exec_module(proxy_module)
FULL = dict(approvalPolicy='never', approvalsReviewer='user', sandbox={'type': 'dangerFullAccess'},
            activePermissionProfile={'id': ':danger-full-access'})
RESTRICTED = dict(approvalPolicy='on-request', approvalsReviewer='user', sandbox={'type': 'readOnly'},
                  activePermissionProfile={'id': ':read-only'})


class ProxyTests(unittest.IsolatedAsyncioTestCase):
    async def asyncSetUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.states = {}
        self.calls = []
        self.mode = ''
        self.streams = []
        self.handlers = set()
        self.server = await asyncio.start_unix_server(self.serve, self.temp.name + '/server')
        self.proxy = proxy_module.ResumeProxy('unix://' + self.temp.name + '/server')
        self.front = await asyncio.start_unix_server(self.proxy.connect, self.temp.name + '/proxy')
        self.reader, self.writer = await self.connect()

    async def connect(self):
        reader, writer = await asyncio.open_unix_connection(self.temp.name + '/proxy')
        self.streams.append(writer)
        writer.write(b'GET / HTTP/1.1\r\n\r\n')
        await writer.drain()
        await proxy_module.header(reader)
        return reader, writer

    async def serve(self, reader, writer):
        task = asyncio.current_task()
        self.handlers.add(task)
        try:
            await proxy_module.header(reader)
            writer.write(b'HTTP/1.1 101 Switching Protocols\r\n\r\n')
            parts = []
            while True:
                first, body, raw = await proxy_module.frame(reader)
                if first & 15 >= 8:
                    writer.write(raw)
                    continue
                parts.append(body)
                if not first & 128:
                    continue
                request = json.loads(b''.join(parts))
                parts = []
                self.calls.append(request)
                method, params = request['method'], request.get('params', {})
                thread = params.get('threadId', 'new')
                if method == 'thread/fork':
                    thread += '-fork'
                if method == 'thread/resume' and thread == 'ephemeral-fork':
                    writer.write(proxy_module.encode(dict(id=request['id'], error=dict(message='no rollout found'))))
                    continue
                self.states.setdefault(thread, dict(RESTRICTED))
                if self.mode != 'wrong' and method in ('thread/start', 'thread/fork') and params.get('permissions') == ':danger-full-access':
                    self.assertNotIn('sandbox', params)
                    self.states[thread] = dict(FULL)
                if method == 'thread/settings/update':
                    if self.mode == 'close':
                        writer.write(b'\x88\x00')
                        await writer.drain()
                        return
                    if self.mode == 'disconnect':
                        return
                    if self.mode == 'timeout':
                        await asyncio.Future()
                    if self.mode == 'reject':
                        writer.write(proxy_module.encode(dict(id=request['id'], error=dict(message='denied'))))
                        continue
                    if self.mode != 'wrong':
                        self.states[thread] = dict(FULL if params.get('permissions') == ':danger-full-access' else RESTRICTED)
                    result = {}
                elif method == 'turn/start':
                    if params.get('permissions') == ':read-only' or params.get('sandboxPolicy', {}).get('type') == 'readOnly':
                        self.states[thread] = dict(RESTRICTED)
                    result = dict(effective=self.states[thread])
                else:
                    result = dict(thread=dict(id=thread, turns=['history']), cwd='/original', model='original', **self.states[thread])
                writer.write(proxy_module.encode(dict(id=request['id'], result=result)))
                await writer.drain()
        except (asyncio.IncompleteReadError, ConnectionError):
            pass
        finally:
            writer.close()
            await writer.wait_closed()
            self.handlers.discard(task)

    async def rpc(self, method='thread/resume', thread='one', connection=None, **extra):
        reader, writer = connection or (self.reader, self.writer)
        writer.write(proxy_module.encode(dict(id=7, method=method, params=dict(threadId=thread, **extra)), True))
        await writer.drain()
        _, body, _ = await asyncio.wait_for(proxy_module.frame(reader), 2)
        reply = json.loads(body)
        self.assertEqual(reply['id'], 7)
        return reply

    async def test_repeated_resume(self):
        for thread in ('one', 'two', 'one'):
            reply = await self.rpc(thread=thread)
            self.assertEqual(reply['result']['sandbox'], FULL['sandbox'])
            self.assertEqual(reply['result']['thread']['turns'], ['history'])
            self.assertEqual(reply['result']['cwd'], '/original')
            self.assertEqual(reply['result']['model'], 'original')

    async def test_reconnect_new_thread(self):
        await self.rpc()
        reply = await self.rpc(thread='two', connection=await self.connect())
        self.assertEqual(reply['result']['sandbox'], FULL['sandbox'])

    async def test_new_then_in_session_resume_and_fork(self):
        for method in ('thread/start', 'thread/resume', 'thread/fork'):
            reply = await self.rpc(method, sandbox='read-only')
            self.assertEqual(reply['result']['sandbox'], FULL['sandbox'])

    async def test_ephemeral_fork_without_rollout(self):
        reply = await self.rpc('thread/fork', thread='ephemeral', ephemeral=True,
                               excludeTurns=True, sandbox='read-only')
        self.assertEqual(reply['result']['sandbox'], FULL['sandbox'])
        self.assertEqual([call['method'] for call in self.calls], ['thread/fork'])
        await self.rpc('thread/settings/update', thread='ephemeral-fork', permissions=':read-only')
        reply = await self.rpc('turn/start', thread='ephemeral-fork')
        self.assertEqual(reply['result']['effective'], RESTRICTED)

    async def test_ephemeral_fork_wrong_permissions_fails(self):
        self.mode = 'wrong'
        reply = await self.rpc('thread/fork', thread='ephemeral', ephemeral=True)
        self.assertIn('cannot restore Full Access', reply['error']['message'])

    async def test_readback_error_preserved(self):
        reply = await self.rpc('thread/fork', thread='ephemeral')
        self.assertIn('no rollout found', reply['error']['message'])

    async def test_prediction_fork_does_not_override_inherited_config(self):
        self.states['prediction-fork'] = dict(FULL)
        reply = await self.rpc('thread/fork', thread='prediction', ephemeral=True,
                               experimentalPredictionMode=True)
        self.assertEqual(reply['result']['sandbox'], FULL['sandbox'])
        self.assertNotIn('permissions', self.calls[0]['params'])
        self.assertEqual(len(self.calls), 1)

    async def test_explicit_choice_survives_reconnect_and_turn(self):
        await self.rpc()
        await self.rpc('thread/settings/update', permissions=':read-only')
        connection = await self.connect()
        reply = await self.rpc(connection=connection)
        self.assertEqual(reply['result']['sandbox'], RESTRICTED['sandbox'])
        reply = await self.rpc('turn/start', connection=connection)
        self.assertEqual(reply['result']['effective']['sandbox'], RESTRICTED['sandbox'])

    async def test_next_turn_and_turn_permission_choice(self):
        await self.rpc()
        reply = await self.rpc('turn/start')
        self.assertEqual(reply['result']['effective'], FULL)
        reply = await self.rpc('turn/start', permissions=':read-only')
        self.assertEqual(reply['result']['effective'], RESTRICTED)

    async def test_failed_restoration(self):
        self.proxy.timeout = .05
        for mode in ('reject', 'wrong', 'disconnect', 'close', 'timeout'):
            self.mode = mode
            reply = await self.rpc(thread=mode, connection=await self.connect())
            self.assertIn('cannot restore Full Access', reply['error']['message'])

    async def test_new_thread_inconsistent_settings_fail(self):
        self.mode = 'wrong'
        reply = await self.rpc('thread/start')
        self.assertIn('cannot restore Full Access', reply['error']['message'])

    async def test_concurrent_resumes_keep_rpc_ids_separate(self):
        for ident in (11, 12):
            self.writer.write(proxy_module.encode(dict(id=ident, method='thread/resume',
                              params=dict(threadId=str(ident))), True))
        await self.writer.drain()
        ids = set()
        for _ in range(2):
            _, body, _ = await asyncio.wait_for(proxy_module.frame(self.reader), 2)
            reply = json.loads(body)
            ids.add(reply['id'])
            self.assertEqual(reply['result']['sandbox'], FULL['sandbox'])
        self.assertEqual(ids, {11, 12})

    async def test_failed_update_can_be_retried(self):
        self.mode = 'reject'
        self.assertIn('error', await self.rpc())
        self.mode = ''
        reply = await self.rpc()
        self.assertEqual(reply['result']['sandbox'], FULL['sandbox'])

    async def test_fragmented_resume_with_ping(self):
        request = json.dumps(dict(id=8, method='thread/resume', params=dict(threadId='fragment'))).encode()
        # fragment a masked message around an interleaved control frame
        def fragment(first, body):
            mask = b'abcd'
            return bytes((first, 128 | len(body))) + mask + bytes(v ^ mask[i % 4] for i, v in enumerate(body))
        self.writer.write(fragment(1, request[:30]) + b'\x89\x00' + fragment(128, request[30:]))
        await self.writer.drain()
        first, _, _ = await asyncio.wait_for(proxy_module.frame(self.reader), 2)
        self.assertEqual(first, 137)
        _, body, _ = await asyncio.wait_for(proxy_module.frame(self.reader), 2)
        reply = json.loads(body)
        self.assertEqual(reply['id'], 8)
        self.assertEqual(reply['result']['sandbox'], FULL['sandbox'])

    async def test_two_connections_do_not_share_thread_choices(self):
        other = await self.connect()
        replies = await asyncio.gather(self.rpc(thread='one'), self.rpc(thread='two', connection=other))
        self.assertTrue(all(r['result']['sandbox'] == FULL['sandbox'] for r in replies))
        await self.rpc('thread/settings/update', thread='one', permissions=':read-only')
        reply = await self.rpc(thread='two', connection=other)
        self.assertEqual(reply['result']['sandbox'], FULL['sandbox'])
        self.assertEqual(self.states['one']['sandbox'], RESTRICTED['sandbox'])

    async def asyncTearDown(self):
        for writer in self.streams:
            writer.close()
            await writer.wait_closed()
        for server in (self.front, self.server):
            server.close()
        tasks = list(self.proxy.tasks | self.handlers)
        for task in tasks:
            task.cancel()
        await asyncio.gather(*tasks, return_exceptions=True)
        for server in (self.front, self.server):
            await server.wait_closed()


if __name__ == '__main__':
    unittest.main()
