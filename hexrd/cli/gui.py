from __future__ import print_function, division, absolute_import

import sys


help = "Launches the hexrd graphical user interface"


def configure_parser(sub_parsers):
    p = sub_parsers.add_parser('gui', description = help, help = help)
    #common.add_parser_json(p)
    p.set_defaults(func=execute)


def execute(args, parser):
    from hexrd.wx import mainapp

    # TODO: this should be improved to not draw directly on sys.argv
    mainapp.execute(*sys.argv[2:])
