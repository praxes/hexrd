from __future__ import print_function, absolute_import

import argparse
import logging
import sys

from . import grains
from . import gui
from . import help
from . import indexing

def main():
    if len(sys.argv) == 1:
        sys.argv.append('-h')

    import logging
    import hexrd

    p = argparse.ArgumentParser(
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

    help.configure_parser(sub_parsers)
    gui.configure_parser(sub_parsers)
    indexing.configure_parser(sub_parsers)
    grains.configure_parser(sub_parsers)

    try:
        import argcomplete
        argcomplete.autocomplete(p)
    except ImportError:
        pass

    args = p.parse_args()

    if args.debug:
        logging.disable(logging.NOTSET)
        logging.basicConfig(level=logging.DEBUG)

    args.func(args, p)
