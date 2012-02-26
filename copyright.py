"""License information

This module provides two variables:

COPYRIGHT_FILE - the name of the file with copyright information
COPYRIGHT_TEXT - the text of that file

"""
import os

COPYRIGHT_FILE = 'COPYRIGHT'
with open(os.path.join(os.path.dirname(__file__), COPYRIGHT_FILE), 'r') as f:
     COPYRIGHT_TEXT = f.read()
