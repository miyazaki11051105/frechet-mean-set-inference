#!/usr/bin/env python3
"""Read-only checks for RELEASE_LOCK.json; never computes or repairs results.

python3 verify_release.py
python3 verify_release.py --recomputed /path/to/completed/run

Keep the release commit/archive or an independently retained lock SHA256 as the
trust anchor. Updating both files and their lock cannot be detected by a lock
stored beside those files alone. A hash check detects changes; it is not backup.
"""
from pathlib import Path, PurePosixPath
import argparse
import csv
import hashlib
import json
import os
import re
import subprocess
import sys


def fail(message):
    raise ValueError(message)


def relative_name(name):
    if not isinstance(name, str) or not name or '\\' in name:
        fail('Invalid relative path in lock')
    path = PurePosixPath(name)
    if path.is_absolute() or '..' in path.parts or ':' in name or str(path) != name:
        fail('Unsafe relative path in lock: ' + name)
    return name


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), ensure_ascii=True).encode('ascii')


def digest(path):
    before = path.stat()
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    after = path.stat()
    if (before.st_size, before.st_mtime_ns, before.st_ino) != (after.st_size, after.st_mtime_ns, after.st_ino):
        fail('File changed while being checked: ' + str(path))
    return h.hexdigest()


def check_files(root, files, label):
    failures = []
    for name, expected in sorted(files.items()):
        relative_name(name)
        if not re.fullmatch('[a-f0-9]{64}', expected):
            fail('Invalid SHA256 in ' + label)
        path = root / name
        try:
            if not path.resolve().is_relative_to(root.resolve()):
                fail('Path escapes checked directory: ' + name)
            if digest(path) != expected:
                failures.append({'path': name, 'reason': 'SHA256 mismatch'})
        except (FileNotFoundError, IsADirectoryError):
            failures.append({'path': name, 'reason': 'missing or not a regular file'})
    if failures:
        print(json.dumps({'status': 'FAIL', 'scope': label, 'failures': failures}), flush=True)
        fail(label + ' verification failed')


def check_source(root, identity):
    for name in identity['release_files']:
        path = PurePosixPath(relative_name(name))
        document_names = {'main', 'body', 'supplement', 'supplement-body', 'supplementary_material'}
        if (path.suffix.lower() in {'.pdf', '.bib', '.cls', '.bst'} or
                path.stem.lower() in document_names or
                any(re.search(r'(^|[-_])(manuscript|cover[-_ ]?letter|arxiv)([-_. ]|$)', part.lower())
                    for part in path.parts)):
            fail('Publication documents are outside this release: ' + name)
    check_files(root, identity['release_files'], 'release source')
    expected = set(identity['release_files']) | {'RELEASE_LOCK.json'}
    # Ignore local outputs, private audit records, Git metadata and installed libraries only at root.
    excluded = {'.git', '.cache', '.R-library', '.private-r-library', 'runs', 'build', 'data', 'evidence', '__pycache__'}
    found = set()
    for directory, folders, files in os.walk(root):
        folders[:] = [name for name in folders if name != '__pycache__']
        if Path(directory) == root:
            folders[:] = [name for name in folders if name not in excluded]
        for name in files:
            relative = (Path(directory) / name).relative_to(root).as_posix()
            found.add(relative)
    if found != expected:
        fail('Release file set differs: ' + json.dumps({'extra': sorted(found - expected), 'missing': sorted(expected - found)}))
    return len(expected)


def csv_fingerprint(path, keys, exclude, records=None):
    """Hash exact exported strings; optionally retain the same verified rows."""
    before = path.stat()
    rows, seen = [], set()
    with path.open(encoding='utf-8', newline='') as stream:
        reader = csv.reader(stream)
        headers = next(reader, [])
        if not headers or len(headers) != len(set(headers)) or not set(keys + exclude).issubset(headers):
            fail('CSV headers differ from the declared contract: ' + str(path))
        columns = [c for c in headers if c not in exclude]
        positions = [headers.index(c) for c in columns]
        key_positions = [headers.index(c) for c in keys]
        for row in reader:
            if len(row) != len(headers):
                fail('CSV row width differs: ' + str(path))
            key = hashlib.sha256(canonical([row[i] for i in key_positions])).digest()
            if key in seen:
                fail('Duplicate scientific row key: ' + str(path))
            seen.add(key)
            rows.append(hashlib.sha256(canonical([row[i] for i in positions])).hexdigest())
            if records is not None:
                records.append(dict(zip(headers, row)))
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        fail('CSV changed during verification: ' + str(path))
    value = {'columns': columns, 'row_sha256_sorted': sorted(rows)}
    return {'rows': len(rows), 'columns': columns, 'sha256': hashlib.sha256(canonical(value)).hexdigest()}


def check_reference(root, run=None):
    command = [sys.executable, '-B', str(root / 'export_results.py'), 'check']
    if run is not None:
        command.extend(['--run', str(run)])
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    if result.returncode != 0:
        fail('Reference output verification failed: ' + (result.stderr or result.stdout).strip())
    return json.loads(result.stdout)


def check_recomputed(run, identity):
    reference = identity['reference_exports']
    roots = {'run': run, 'ranking': Path(str(run) + '-ranking-candidates')}
    verified = []
    for name, spec in reference['csv'].items():
        path = roots[spec['location']] / relative_name(spec['path'])
        observed = csv_fingerprint(path, spec['keys'], spec['excluded_execution_columns'])
        if observed != spec['fingerprint']:
            fail('Scientific export differs from the fixed production reference: ' + name)
        verified.append({'component': name, 'rows': observed['rows']})
    check_files(run, reference['deterministic_tex_sha256'], 'reference TeX')
    return {'status': 'PASS', 'csv': verified, 'tex_files': len(reference['deterministic_tex_sha256']),
            'assurance': 'Exact scientific exports against the production references; the R run separately enforces numerical completion.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--recomputed', type=Path, help='Compare completed RUN exports to the fixed 1000-replication production reference')
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    lock_path = root / 'RELEASE_LOCK.json'
    lock_digest = digest(lock_path)
    lock = json.loads(lock_path.read_text(encoding='utf-8'))
    if set(lock) != {'schema', 'identity', 'identity_sha256'} or lock['schema'] != 'frechet-reproducibility-release-v1':
        fail('Unrecognized release lock schema')
    if set(lock['identity']) != {'release', 'license', 'release_files', 'reference_exports'}:
        fail('Unexpected release identity fields')
    if hashlib.sha256(canonical(lock['identity'])).hexdigest() != lock['identity_sha256']:
        fail('Release lock identity is corrupt')
    count = check_source(root, lock['identity'])
    results = check_reference(root, args.recomputed.resolve() if args.recomputed else None)
    recomputed = check_recomputed(args.recomputed.resolve(), lock['identity']) if args.recomputed else None
    check_source(root, lock['identity'])
    if digest(lock_path) != lock_digest:
        fail('Release lock changed during verification')
    print(json.dumps({'status': 'PASS', 'release_files': count, 'release_identity_sha256': lock['identity_sha256'],
                      'lock_sha256': lock_digest, 'reference_results': results, 'reference_exports': recomputed,
                      'assurance': 'Integrity relative to the retained release anchor; no simulation, repairs, writes, or universal cross-platform reproducibility claim.'}), flush=True)


if __name__ == '__main__':
    try:
        main()
    except (ValueError, KeyError, TypeError, OSError) as error:
        print(json.dumps({'status': 'FAIL', 'error': str(error)}), file=sys.stderr)
        raise SystemExit(1)
