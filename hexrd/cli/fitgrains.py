from __future__ import print_function, division, absolute_import


descr = 'Extracts G vectors, grain position and strain'
example = """
examples:
    hexrd grains configuration.yml
"""


def configure_parser(sub_parsers):
    p = sub_parsers.add_parser('fit-grains', description = descr, help = descr)
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
    p.set_defaults(func=execute)


def execute(args, parser):
    import logging
    import os
    import sys

    import yaml

    from hexrd import config
    from hexrd.fitgrains import fit_grains


    log_level = logging.DEBUG if args.debug else logging.INFO
    if args.quiet:
        log_level = logging.ERROR
    logger = logging.getLogger('hexrd')
    logger.setLevel(log_level)
    ch = logging.StreamHandler()
    ch.setLevel(log_level)
    ch.setFormatter(
        logging.Formatter('%(asctime)s - %(message)s', '%y-%m-%d %H:%M:%S')
        )
    logger.addHandler(ch)

    for cfg in config.open(args.yml):
        if os.path.exists(cfg.analysis_dir) and not args.force:
            logger.error(
                'Analysis "%s" at %s already exists.'
                ' Change yml file or specify "force"',
                cfg.analysis_name, cfg.analysis_dir
                )
            sys.exit()
        # now we know where to save the log file
        if not os.path.exists(cfg.analysis_dir):
            os.makedirs(cfg.analysis_dir)
        logfile = os.path.join(
            cfg.working_dir,
            cfg.analysis_name,
            'fit-grains.log'
            )
        fh = logging.FileHandler(logfile, mode='w')
        fh.setLevel(log_level)
        fh.setFormatter(
            logging.Formatter(
                '%(asctime)s - %(name)s - %(message)s',
                '%m-%d %H:%M:%S'
                )
            )
        logger.info("logging to %s", logfile)
        logger.addHandler(fh)

        fit_grains(cfg, force=args.force)

        logger.removeHandler(fh)
        fh.close()
