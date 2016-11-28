from importlib import import_module
from collections import OrderedDict
import os

# from ..internal.settings import config
from ..internal import session

_all_interfaces = OrderedDict()

# Assumes all subdirectories are interface implementations
subdirs =  next(os.walk(os.path.dirname(__file__)))[1]

for intr in subdirs:
    try:
        intr_mod = import_module('.'+intr, 'Rosetta.interfaces')

        intr_class = '{0}Interface'.format(intr)
        
        _all_interfaces[intr_class] = getattr(intr_mod, intr_class)
                                            
    except ImportError as ee:
        msg = "Rosetta couldn't load {} interface. Error:{}.".format(intr, ee)
        session.log(msg)

