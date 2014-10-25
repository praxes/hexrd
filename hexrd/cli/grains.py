from __future__ import print_function, division, absolute_import


descr = 'Extracts G vectors, grain position and strain'
example = """
examples:
    hexrd grains configuration.yml
"""


def configure_parser(sub_parsers):
    p = sub_parsers.add_parser('grains', description = descr, help = descr)
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
    import yaml

    from hexrd.coreutil import iter_cfg_sections
    from hexrd.fitgrains import fit_grains

    for cfg in iter_cfg_sections(args.yml):
        fit_grains(cfg, verbose=not args.quiet, force=args.force)
