import sys
# from .. import config
from ...internal.settings import config

from .errors import LilithImportError, LilithInterfaceError
################################################################################
# check for NumPy >= 1.6.1 and Scipy >= 0.9.0
try:
    import numpy
    # if StrictVersion(numpy.__version__) < StrictVersion('1.6.1'):
    #     raise ImportError
except ImportError:
    err = ('NumPy version 1.6.1 or more recent '
           'must be installed to use Lilith interface')   
    raise LilithImportError(err)
except Exception as err:
    raise LilithImportError(err)

try:
    import scipy
    # if StrictVersion(scipy.__version__) < StrictVersion('0.9.0'):
    #     raise ImportError
except ImportError:
    err = ('SciPy version 0.9.0 or more recent '
           'must be installed to use Lilith interface')   
    raise LilithImportError(err)
except Exception as err:
    raise LilithImportError(err)
################################################################################
from .interface import LilithInterface
