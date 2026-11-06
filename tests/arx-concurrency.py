#!/usr/bin/env python3
"""Check generated arx commands without submitting jobs or running tasks."""
import os
from pathlib import Path
import shlex
import subprocess
import tempfile
import unittest


ARX = Path(__file__).resolve().parents[1] / 'arx'


class ConcurrencyTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='arx-concurrency-')
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.commands = self.root / 'commands.txt'
        self.commands.write_text('true\n' * 53)
        self.data = self.root / 'samples.txt'
        self.data.write_text(''.join('sample-%03d\n' % i for i in range(1, 54)))
        self.script = self.root / 'task.sh'
        self.script.write_text('#!/bin/sh\nset -eu\nexec linix "$1" "$SLURM_ARRAY_TASK_ID"\n')
        self.env = dict(os.environ)
        self.env['PATH'] = str(ARX.parent) + os.pathsep + self.env['PATH']

    def preview(self, options, command_file=True, manifest=False):
        args = ['-C', str(self.commands)] if command_file else [str(self.script)]
        if manifest:
            args.append(str(self.data))
        result = subprocess.run(['perl', str(ARX), 'sub', '-D', *options, *args],
                                cwd=self.root, env=self.env, capture_output=True,
                                text=True, check=True)
        if '-P' in options:
            lines = [s for s in result.stderr.splitlines() if ' --delay .01 ' in s]
        else:
            lines = [s for s in result.stderr.splitlines()
                     if s.startswith('sbatch ') and '--dependency=' not in s]
        self.assertTrue(lines, result.stderr)
        return shlex.split(lines[0])

    def test_slurm_command_file_defaults(self):
        for options in ([], ['-a1-'], ['-a', '1-'], ['--array=1-'],
                        ['-a', '1-53'], ['-j', '0']):
            with self.subTest(options=options):
                cmd = self.preview(options)
                self.assertEqual(cmd[cmd.index('-a') + 1], '1-53')

    def test_explicit_slurm_limits(self):
        cases = [(['-j', '7'], '1-53%7'),
                 (['-a', '1-%9'], '1-53%9'),
                 (['--array=1-53%9'], '1-53%9'),
                 (['-a', '1-53%9', '-j', '4'], '1-53%4'),
                 (['-a', '1-53%4', '-j', '9'], '1-53%4'),
                 (['-a', '2-40:2'], '2-40:2'),
                 (['-a', '2,4,8', '-j', '2'], '2,4,8%2'),
                 (['-j', '53'], '1-53')]
        for options, expected in cases:
            with self.subTest(options=options):
                cmd = self.preview(options)
                self.assertEqual(cmd[cmd.index('-a') + 1], expected)

    def test_explicit_array_without_manifest(self):
        for options, expected in [(['-a', '1-53'], '1-53'),
                                  (['-a', '1-53', '-j', '7'], '1-53%7'),
                                  (['-a', '1-53%7'], '1-53%7')]:
            with self.subTest(options=options):
                cmd = self.preview(options, command_file=False)
                self.assertEqual(cmd[cmd.index('-a') + 1], expected)

    def test_script_without_array(self):
        for manifest in (False, True):
            with self.subTest(manifest=manifest):
                cmd = self.preview([], command_file=False, manifest=manifest)
                self.assertNotIn('-a', cmd)
                self.assertIn(str(self.script), cmd)

    def test_script_with_manifest(self):
        for options, expected in [(['-a1-'], '1-53'),
                                  (['--array=1-'], '1-53'),
                                  (['-a', '1-53'], '1-53'),
                                  (['-a', '2-40:2'], '2-40:2'),
                                  (['-a', '2,4,8'], '2,4,8'),
                                  (['-a1-', '-j', '0'], '1-53'),
                                  (['-a', '1-', '-j', '7'], '1-53%7'),
                                  (['-a', '1-%9'], '1-53%9'),
                                  (['-a', '1-%9', '-j', '4'], '1-53%4')]:
            with self.subTest(options=options):
                cmd = self.preview(options, command_file=False, manifest=True)
                self.assertEqual(cmd[cmd.index('-a') + 1], expected)
                self.assertEqual(cmd[-2:], [str(self.script), str(self.data)])

    def test_script_selects_array_row(self):
        cmd = self.preview(['-a1-'], command_file=False, manifest=True)
        # emulate the task ID supplied by Slurm and run the generated script arguments.
        env = dict(self.env, SLURM_ARRAY_TASK_ID='17')
        result = subprocess.run(['sh', *cmd[-2:]], cwd=self.root, env=env,
                                capture_output=True, text=True, check=True)
        self.assertEqual(result.stdout.strip(), 'sample-017')

    def test_local_parallel_default_and_explicit_limit(self):
        for options, expected in [(['-P'], '6'), (['-P', '-a1-'], '6'),
                                  (['-P', '-j', '3'], '3')]:
            with self.subTest(options=options):
                cmd = self.preview(options)
                self.assertEqual(cmd[cmd.index('-j') + 1], expected)


if __name__ == '__main__':
    unittest.main()
