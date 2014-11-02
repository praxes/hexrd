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


    # load the configuration settings
    cfgs = config.open(args.yml)

    # if find-orientations has not already been run, do so:
    quats_f = os.path.join(cfgs[0].working_dir, 'accepted_orientations.dat')
    if not os.path.exists(quats_f):
        from . import findorientations
        findorientations.execute(args, parser)

    # configure logging to the console:
    log_level = logging.DEBUG if args.debug else logging.INFO
    if args.quiet:
        log_level = logging.ERROR
    logger = logging.getLogger('hexrd')
    logger.setLevel(log_level)
    ch = logging.StreamHandler()
    ch.setLevel(log_level)
    cf = logging.Formatter('%(asctime)s - %(message)s', '%y-%m-%d %H:%M:%S')
    ch.setFormatter(cf)
    logger.addHandler(ch)
    logger.info('=== begin fit-grains ===')

    for cfg in config.open(args.yml):
        # prepare the analysis directory
        if os.path.exists(cfg.analysis_dir) and not args.force:
            logger.error(
                'Analysis "%s" at %s already exists.'
                ' Change yml file or specify "force"',
                cfg.analysis_name, cfg.analysis_dir
                )
            sys.exit()
        if not os.path.exists(cfg.analysis_dir):
            os.makedirs(cfg.analysis_dir)

        logger.info('*** begin analysis "%s" ***', cfg.analysis_name)

        # configure logging to file for this particular analysis
        logfile = os.path.join(
            cfg.working_dir,
            cfg.analysis_name,
            'fit-grains.log'
            )
        fh = logging.FileHandler(logfile, mode='w')
        fh.setLevel(log_level)
        ff = logging.Formatter(
                '%(asctime)s - %(name)s - %(message)s',
                '%m-%d %H:%M:%S'
                )
        fh.setFormatter(ff)
        logger.info("logging to %s", logfile)
        logger.addHandler(fh)

        # process the data
        fit_grains(cfg, force=args.force, show_progress=True)

        # stop logging for this particular analysis
        fh.flush()
        fh.close()
        logger.removeHandler(fh)

        logger.info('*** end analysis "%s" ***', cfg.analysis_name)

    logger.info('=== end fit-grains ===')
    # stop logging to the console
    ch.flush()
    ch.close()
    logger.removeHandler(ch)
