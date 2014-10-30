from __future__ import print_function, division, absolute_import


descr = 'Process diffraction data to find grain orientations'
example = """
examples:
    hexrd find-orientations configuration.yml
"""


def configure_parser(sub_parsers):
    p = sub_parsers.add_parser(
        'find-orientations',
        description = descr,
        help = descr
        )
    p.add_argument(
        'yml', type=str,
        help='YAML configuration file'
        )
    p.add_argument(
        '-q', '--quiet', action='store_true',
        help="don't report progress in terminal"
        )
    p.add_argument(
        '-f', '--force', action='store_true',
        help='overwrites existing analysis'
        )
    p.add_argument(
        '--hkls', metavar='HKLs', type=str, default=None,
        help="""\
list hkl entries in the materials file to use for fitting
if None, defaults to list specified in the yml file"""
        )
    p.set_defaults(func=execute)


def execute(args, parser):
    import os
    import logging

    import yaml

    from hexrd import config
    from hexrd.findorientations import find_orientations


    if args.hkls is not None:
        args.hkls = [int(i) for i in args.hkls.split(',') if i]

    cfg = config.open(args.yml)[0]

    # now we know where to save the log file
    log_level = logging.DEBUG if args.debug else logging.INFO
    logger = logging.getLogger('hexrd')
    logger.setLevel(log_level)

    fn = os.path.join(
        cfg.working_dir,
        cfg.analysis_name,
        'find-orientations.log'
        )
    fh = logging.FileHandler(fn)
    fh.setLevel(log_level)
    fh.setFormatter(logging.Formatter())

    ch = logging.StreamHandler()
    ch.setLevel(log_level)
    ch.setFormatter(logging.Formatter())
    logger.addHandler(ch)
    logger.info("logging to %s", fn)
    logger.addHandler(fh)

    find_orientations(
        cfg, hkls=args.hkls, force=args.force
        )
