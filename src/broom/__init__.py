"""
Broom
=====

Provides a set of tools to interface with the DUST aerodynamic solver.

Available subpackages
---------------------
charm
    Functions to parse CHARM input files.
mesh
    Utilities to create mesh objects.
polartable
    Polar table class and parsers.
post
    DUST post-processing tools.
reference
    Reference frames and transformations.
settings
    DUST solver settings.
solver
    DUST solver case creation and interface.
utils
    Miscellaneous utilities.

"""


from importlib.metadata import version
__version__ = version("broom")

__all__ = [
    'charm',
    'mesh',
    'polartable',
    'post',
    'reference',
    'settings',
    'solver',
    'utils'
    ]
