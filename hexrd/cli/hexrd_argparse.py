# hexrd cli draws heavily from ContinuumIO's conda effort
#
# (c) 2012-2013 Continuum Analytics, Inc. / http://continuum.io
# All Rights Reserved
#
# conda is distributed under the terms of the BSD 3-clause license.
# Consult LICENSE.txt or http://opensource.org/licenses/BSD-3-Clause.

from __future__ import print_function, division, absolute_import

import sys
import argparse

from difflib import get_close_matches

from .find_commands import find_commands


class ArgumentParser(argparse.ArgumentParser):
    def _get_action_from_name(self, name):
        """Given a name, get the Action instance registered with this parser.
        If only it were made available in the ArgumentError object. It is
        passed as it's first arg...
        """
        container = self._actions
        if name is None:
            return None
        for action in container:
            if '/'.join(action.option_strings) == name:
                return action
            elif action.metavar == name:
                return action
            elif action.dest == name:
                return action

    def error(self, message):
        import re
        import subprocess

        exc = sys.exc_info()[1]
        if exc:
            if hasattr(exc, 'argument_name'):
                argument = self._get_action_from_name(exc.argument_name)
            else:
                argument = None
        super(ArgumentParser, self).error(message)


    def print_help(self):
        super(ArgumentParser, self).print_help()

        if sys.argv[1:] in ([], ['help'], ['-h'], ['--help']):
            from .find_commands import help
            help()
