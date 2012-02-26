"""License information

This just reads the COPYRIGHT file in the same directory
as this file.  
"""
import os

COPYRIGHT_FILE = 'COPYRIGHT'
with open(os.path.join(os.path.dirname(__file__), COPYRIGHT_FILE), 'r') as f:
     COPYRIGHT_TEXT = f.read()
