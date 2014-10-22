import argparse
import textwrap

from hexrd.coreutil import iter_cfg_sections
from hexrd.peakextraction.extractgvecs import extract_g_vectors


def main_extract_g_vectors():
    parser = argparse.ArgumentParser(
        description='Extracts measured G vectors',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent('''
            example:
            $ extract_gvecs -v configuration.yml
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
        extract_g_vectors(cfg, verbose=not args.quiet, force=args.force)
