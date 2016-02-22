import os
# print __name__
from .. import config
from errors import eHDECAYImportError, eHDECAYInterfaceError
# from ....bases import SILHBasis as SB

# eHDECAY executable
try:
    eHDECAY_dir = config['eHDECAY_dir']
except KeyError:
    err = ('Could not find option "eHDECAY_dir" in Rosetta/config.txt')
    raise eHDecayImportError(err)
    
executable = '{}/run'.format(eHDECAY_dir)

if not os.path.exists(executable):
    err = ('Could not find eHDECAY executable at {}'.format(executable))
    raise eHDecayImportError(err)
