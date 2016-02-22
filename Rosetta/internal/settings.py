import os
from SLHA import CaseInsensitiveDict as CIdict
from errors import ReadSettingsError
from ..bases import __all__ as implemented_bases
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
try:
    with open('{}/../config.txt'.format(dirname)) as cfg:
        lines = [x.strip() for x in cfg.readlines()]
except IOError:
    raise ReadSettingsError()
    
for i,l in enumerate(lines):
    if l and not l.startswith('#'):
        try:
            field, value = tuple(l.split())
        except Exception:
            print ('Rosetta ignored line {} of config.txt'.format(i+1))
            field, value = None, None
    else:
        field, value = None, None

    config[field] = value
    
################################################################################
    