import argparse
import textwrap

from hexrd.indexing.find_orientations import find_orientations


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
        '-v', '--verbose', action='store_true',
        help='report progress in terminal'
        )
    parser.add_argument(
        '-f', '--force', action='store_true',
        help='force overwrite of existing data'
        )
    args = parser.parse_args()
    extract_g_vectors(args.yml, verbose=args.verbose, force=args.force)
