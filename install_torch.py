#!/usr/bin/env python3
"""Install PyTorch 2.5.1 in the active conda environment (Linux x86_64)."""
import argparse
import importlib.metadata
import os
from pathlib import Path
import platform
import re
import subprocess
import sys

VERSION = '2.5.1'
VARIANTS = ('auto', 'cpu', 'cu118', 'cu121', 'cu124')


def select_variant(smi_output):
    # The reported CUDA version is the driver's supported maximum, not nvcc.
    match = re.search(r'CUDA Version:\s*(\d+)\.(\d+)', smi_output)
    if not match:
        raise RuntimeError('Cannot read driver CUDA support; specify cpu/cu118/cu121/cu124.')
    supported = tuple(map(int, match.groups()))
    for minimum, variant in (((12, 4), 'cu124'), ((12, 1), 'cu121'), ((11, 8), 'cu118')):
        if supported >= minimum:
            return variant
    raise RuntimeError('Driver CUDA support is below 11.8; update the driver or specify cpu.')


def detect_variant():
    try:
        devices = subprocess.check_output(
            ['nvidia-smi', '--query-gpu=name', '--format=csv,noheader'],
            text=True, stderr=subprocess.STDOUT, timeout=15)
        if not devices.strip() or 'No devices' in devices:
            raise RuntimeError('No GPU detected.')
        output = subprocess.check_output(
            ['nvidia-smi'], text=True, stderr=subprocess.STDOUT, timeout=15)
    except (OSError, subprocess.SubprocessError) as exc:
        raise RuntimeError('GPU detection failed. Run auto on a GPU compute node, '
                           'or explicitly specify cpu/cu118/cu121/cu124.') from exc
    return select_variant(output)


def verify(variant, require_gpu):
    # Use a fresh interpreter after pip installs torch.
    code = '''
import sys
import torch
variant, require_gpu = sys.argv[1], sys.argv[2] == 'True'
expected = {'cpu': None, 'cu118': '11.8', 'cu121': '12.1', 'cu124': '12.4'}[variant]
assert torch.__version__.split('+')[0] == '2.5.1', torch.__version__
assert torch.version.cuda == expected, (torch.version.cuda, expected)
print('PyTorch:', torch.__version__, 'CUDA runtime:', torch.version.cuda)
if variant == 'cpu':
    x = torch.ones((2, 2))
    assert (x @ x).sum().item() == 8
elif torch.cuda.is_available():
    for i in range(torch.cuda.device_count()):
        x = torch.ones((2, 2), device=f'cuda:{i}')
        assert (x @ x).sum().item() == 8
        torch.cuda.synchronize(i)
        print('GPU check passed:', torch.cuda.get_device_name(i))
elif require_gpu:
    raise RuntimeError('CUDA is not available. Check GPU allocation, visibility and driver compatibility.')
else:
    print('GPU check deferred: no usable GPU here. Run install_torch.py with the same variant on the compute node to verify.')
'''
    subprocess.run([sys.executable, '-c', code, variant, str(require_gpu)], check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('variant', choices=VARIANTS)
    parser.add_argument('--dry-run', action='store_true', help='Show the command without installing.')
    args = parser.parse_args()
    if platform.system() != 'Linux' or platform.machine() not in ('x86_64', 'AMD64'):
        parser.error('This installer targets Linux x86_64.')
    prefix = os.environ.get('CONDA_PREFIX')
    if not prefix or Path(prefix).resolve() != Path(sys.prefix).resolve():
        parser.error('Activate the prediction conda environment and use its python first.')
    if sys.version_info[:2] != (3, 10):
        parser.error('The supplied environment requires Python 3.10.')
    variant = detect_variant() if args.variant == 'auto' else args.variant
    target = f'{VERSION}+{variant}'
    command = [sys.executable, '-m', 'pip', 'install', f'torch=={target}',
               '--index-url', f'https://download.pytorch.org/whl/{variant}']
    print('Selected:', target, flush=True)
    if args.dry_run:
        print(' '.join(command))
        return
    try:
        installed = importlib.metadata.version('torch')
    except importlib.metadata.PackageNotFoundError:
        installed = None
    if installed is not None and installed != target:
        raise RuntimeError(f'This environment already contains torch {installed}. '
                           'Create a fresh environment for a different PyTorch build.')
    if installed is None:
        subprocess.run(command, check=True)
    verify(variant, require_gpu=(args.variant == 'auto'))


if __name__ == '__main__':
    try:
        main()
    except (RuntimeError, subprocess.SubprocessError) as exc:
        print(f'Error: {exc}', file=sys.stderr)
        sys.exit(1)
