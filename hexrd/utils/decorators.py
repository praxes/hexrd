# encoding: utf-8
"""Decorators that don't go anywhere else.

This module contains decorators that don't really go with another module
in :mod:`hexrd.utils`. Beore putting something here please see if it should
go into another topical module in :mod:`hexrd.utils`.
"""

def undoc(func):
    """Mark a function or class as undocumented.
    This is found by inspecting the AST, so for now it must be used directly
    as @undoc, not as e.g. @decorators.undoc
    """
    return func
