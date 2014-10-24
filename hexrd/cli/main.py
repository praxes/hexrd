from __future__ import print_function, absolute_import

import argparse
import sys

from . import hexrd_argparse
from . import main_grains
from . import main_gui
from . import main_help
from . import main_indexing

def main():
    if len(sys.argv) == 1:
        sys.argv.append('-h')

    import logging
    import hexrd

    p = hexrd_argparse.ArgumentParser(
        description='High energy diffraction data analysis'
    )
    p.add_argument(
        '-V', '--version',
        action = 'version',
        version = 'hexrd %s' % hexrd.__version__,
    )
    p.add_argument(
        "--debug",
        action = "store_true",
        help = argparse.SUPPRESS,
    )
    sub_parsers = p.add_subparsers(
        metavar = 'command',
        dest = 'cmd',
    )

    main_help.configure_parser(sub_parsers)
    main_gui.configure_parser(sub_parsers)
    main_indexing.configure_parser(sub_parsers)
    main_grains.configure_parser(sub_parsers)

    try:
        import argcomplete
        argcomplete.autocomplete(p)
    except ImportError:
        pass

    args = p.parse_args()

    if args.debug:
        logging.disable(logging.NOTSET)
        logging.basicConfig(level=logging.DEBUG)

    args_func(args, p)


def args_func(args, p):
    use_json = getattr(args, 'json', False)
    try:
        args.func(args, p)
    except RuntimeError as e:
        common.error_and_exit(str(e), json=use_json)
    except Exception as e:
        if e.__class__.__name__ not in ('ScannerError', 'ParserError'):
            message = """\
An unexpected error has occurred, please consider sending the
following traceback to the hexrd GitHub issue tracker at:

    https://github.com/praxes/hexrd/issues

Include the output of the command 'hexrd info' in your report.

"""
            print(message)
        raise  # as if we did not catch it


if __name__ == '__main__':
    main()
