import os
# from .. import config
from ...internal.settings import config

from errors import eHDECAYImportError, eHDECAYInterfaceError

# eHDECAY executable
try:
    eHDECAY_dir = config['eHDECAY_dir']
except KeyError:
    err = ('Could not find option "eHDECAY_dir" in Rosetta/config.txt')
    raise eHDECAYImportError(err)
    
executable = '{}/run'.format(eHDECAY_dir)

if not os.path.exists(executable):
    err = (('Could not find eHDECAY executable at {}: check option '
            '"eHDECAY_dir" in Rosetta/config.txt').format(executable))
    raise eHDECAYImportError(err)

from interface import eHDECAYInterface
