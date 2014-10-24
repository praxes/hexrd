import argparse
import textwrap

from hexrd.coreutil import iter_cfg_sections
from hexrd.fitgrains import fit_grains


def main_fit_grains():
    parser = argparse.ArgumentParser(
        description='Extracts G vectors, grain position and strain',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent('''
            example:
            $ fit_grains -v configuration.yml
            ''')
        )
    parser.add_argument(
        'yml', type=str,
        help='YAML configuration file'
        )
    parser.add_argument(
        '-q', '--quiet', action='store_true',
        help='report progress in terminal'
        )
    parser.add_argument(
        '-f', '--force', action='store_true',
        help='force overwrite of existing data'
        )
    args = parser.parse_args()

    for cfg in iter_cfg_sections(args.yml):
        fit_grains(cfg, verbose=not args.quiet, force=args.force)
