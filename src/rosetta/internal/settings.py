import os
from .SLHA import CaseInsensitiveDict as CIdict
from .errors import ReadSettingsError
# from ..bases import __all__ as implemented_bases
from . import session
__doc__ ='''
    Reads config.txt and sets some variables relevant to running Rosetta.
'''
################################################################################
config = CIdict()
################################################################################
force = False
verbose = False
silent = False
################################################################################
import importlib.resources as resources

# Find the package directory
try:
    # Get the config.txt file from the package
    with resources.open_text("rosetta", "config.txt") as cfg:
        lines = [x.strip() for x in cfg.readlines()]
except IOError:
    err = 'error reading rosetta/config.txt'
    raise ReadSettingsError(err)
    
for i,l in enumerate(lines):
    if l and not l.startswith('#'):
        try:
            field, value = tuple(l.split())
            config[field] = value
        except Exception:
            session.verbose('Rosetta ignored line {} of config.txt'.format(i+1))
################################################################################
    