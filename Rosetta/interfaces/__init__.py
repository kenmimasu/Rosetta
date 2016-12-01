from importlib import import_module
from collections import OrderedDict
import os

from ..internal import session
from ..internal.errors import RosettaImportError

from errors import LoadInterfaceError

_all_interfaces = OrderedDict()

# Assumes all subdirectories are interface implementations
subdirs =  next(os.walk(os.path.dirname(__file__)))[1]

for intr in subdirs:
    try:
        intr_mod = import_module('.'+intr, 'Rosetta.interfaces')

        intr_class = '{0}Interface'.format(intr)
        
        _all_interfaces[intr_class] = getattr(intr_mod, intr_class)
                                            
    except ImportError as ee:
        msg = "Warning: Rosetta couldn't load {} interface.\n    Error:{}.\n".format(intr, ee)
        session.log(msg)
        
    except AttributeError as ee:
        msg = "Warning: Rosetta couldn't load {0}Interface Class from {0} module.\n    Error:{1}.\n".format(intr, ee)
        session.log(msg)

if 'TranslateInterface' not in _all_interfaces:
    msg = "Rosetta couldn't find the translate interface!"
    raise LoadInterfaceError(msg)

