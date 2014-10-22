import argparse
import textwrap

import yaml

from hexrd.indexing.orientations import find_orientations


def main_find_orientations():
    parser = argparse.ArgumentParser(
        description='Process diffraction data to find grain orientations',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent('''
            example:
            $ find_orientations -v configuration.yml
            ''')
        )
    parser.add_argument(
        'yml', type=str,
        help='YAML configuration file'
        )
    parser.add_argument(
        '-q', '--quiet', action='store_true',
        help="don't report progress in terminal"
        )
    parser.add_argument(
        '-f', '--force', action='store_true',
        help='overwrites existing analysis'
        )
    parser.add_argument(
        '--hkls', metavar='HKLs', type=str, default=None,
        help=textwrap.dedent(
            '''list hkl entries in the materials file to use for fitting
            if None, defaults to list specified in the yml file'''
            )
        )
    args = parser.parse_args()
    if args.hkls is not None:
        args.hkls = [int(i) for i in args.hkls.split(',') if i]

    with open(args.yml) as f:
        cfg = [cfg for cfg in yaml.load_all(f)][0]
    find_orientations(
        cfg, verbose=not args.quiet, hkls=args.hkls, force=args.force
        )
