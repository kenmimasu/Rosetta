import os
from SLHA import CaseInsensitiveDict as CIdict
from errors import ReadSettingsError
from ..bases import __all__ as implemented_bases
import session
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
dirname = os.path.dirname(__file__)

# package base directory
rosetta_root = os.path.abspath(dirname+'/../')

# read config.txt
try:
    with open('{}/../config.txt'.format(dirname)) as cfg:
        lines = [x.strip() for x in cfg.readlines()]
except IOError:
    err = 'error reading {}/../config.txt'.format(dirname)
    raise ReadSettingsError(err)
    
for i,l in enumerate(lines):
    if l and not l.startswith('#'):
        try:
            field, value = tuple(l.split())
            config[field] = value
        except Exception:
            session.verbose('Rosetta ignored line {} of config.txt'.format(i+1))
################################################################################
    