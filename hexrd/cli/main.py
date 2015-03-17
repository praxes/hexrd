from __future__ import print_function, absolute_import

import argparse
import logging
import sys
import multiprocessing
import warnings

# These can't be relative imports on Windows because of the hack
# in main() for multiprocessing.freeze_support()
from hexrd.cli import cacheframes
from hexrd.cli import documentation
from hexrd.cli import findorientations
from hexrd.cli import fitgrains
from hexrd.cli import gui
from hexrd.cli import help
from hexrd.cli import test

try:
    from numbapro import nvtx
except ImportError:
    pass

try:
    import importlib
except ImportError:
    pass


def profile_instrument_function(fn_desc):
    """fn_desc contains:
    'fn' is the full path to the function
    """

    # we must split 'fn' into the module and the function itself. Note
    # this is not trivial as we may find several levels of objects inside
    # the containing module. Try sequentially...
    full_name = fn_desc['fn']
    color = fn_desc.get('color', 'black')
    color = getattr(nvtx.colors, color, nvtx.colors.black)
    parts = full_name.split('.')

    # last item will always be the function name
    fn_name = parts[-1]

    # number of parts in the path
    path_parts = len(parts) - 1

    # consume as many as possible with import (ignore last part that is the function name)
    pos = 0
    for i in range(1, path_parts+1):
        try:
            m = importlib.import_module('.'.join(parts[0:i]))
            pos = i
        except ImportError as e:
            break

    # at this point, i points at the starting of the dotted path to the function
    # to instrument... follow the parts till we get to the actual function
    try:
        o = m
        for i in range(pos, path_parts):
            o = getattr(o, parts[i])

        # instrument...
        original = getattr(o, fn_name)
        override = nvtx.profiled(full_name, color=color)(original)
        setattr(o, fn_name, override)
    except AttributeError:
        raise
        warnings.warn('Could not instrument "{0}"'.format(full_name))



def profile_parse_file(filename):
    try:
        import yaml
        with open(filename, 'r') as f:
            cfg = yaml.load(f)

        if 'profile' not in cfg:
            warnings.warn('profile file "{0}" missing a profile section'.format(filename))
            return

        profile_cfg = cfg['profile']
        if 'instrument' in profile_cfg:
            # instrument all
            [profile_instrument_function(fn_desc) for fn_desc in profile_cfg['instrument']]


    except Exception as e:
        raise
        warnings.warn('Failed to include profile file: {0}'.format(filename))
        warnings.warn(str(e))


def profile_instrument_all(args):
    """
    args contains a list of the yaml files containing instrumentation
    information
    """
    [profile_parse_file(filename) for filename in args]


def profile_dump_results(args):
    print(" STATS ".center(72, '='))
    fmt = "{2:>14}, {1:>8}, {0:<40}"
    print(fmt.format("FUNCTION", "CALLS", "TIME"))
    fmt = "{2:>14F}, {1:>8}, {0:<40}"
    sorted_by_time = sorted(nvtx.getstats().iteritems(), key=lambda tup: tup[1][1])
    for key, val in sorted_by_time:
        print(fmt.format(key, *val))


def main():
    if sys.platform.startswith('win'):
        # Hack for multiprocessing.freeze_support() to work from a
        # setuptools-generated entry point.
        if __name__ != "__main__":
            sys.modules["__main__"] = sys.modules[__name__]
        multiprocessing.freeze_support()

    if len(sys.argv) == 1:
        sys.argv.append('-h')

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
        help = 'verbose reporting',
    )
    p.add_argument(
        "--inst-profile",
        action="append",
        help='use the following files as source for functions to instrument',
    )
    sub_parsers = p.add_subparsers(
        metavar = 'command',
        dest = 'cmd',
    )

    help.configure_parser(sub_parsers)
    documentation.configure_parser(sub_parsers)
    gui.configure_parser(sub_parsers)
    findorientations.configure_parser(sub_parsers)
    fitgrains.configure_parser(sub_parsers)
    cacheframes.configure_parser(sub_parsers)
    test.configure_parser(sub_parsers)

    try:
        import argcomplete
        argcomplete.autocomplete(p)
    except ImportError:
        pass

    args = p.parse_args()

    log_level = logging.DEBUG if args.debug else logging.INFO
    logger = logging.getLogger('hexrd')
    ch = logging.StreamHandler()
    ch.setLevel(log_level)

    if hasattr(args, 'inst_profile'):
        profile_instrument_all(args.inst_profile)

    args.func(args, p)

    if hasattr(args, 'inst_profile'):
        profile_dump_results(args.inst_profile)
