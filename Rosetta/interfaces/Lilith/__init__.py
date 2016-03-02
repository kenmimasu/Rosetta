from distutils.version import StrictVersion
import sys
from .. import config 
from errors import LilithImportError, LilithInterfaceError
################################################################################
# check for NumPy >= 1.6.1 and Scipy >= 0.9.0
try:
    import numpy
    if StrictVersion(numpy.__version__) < StrictVersion('1.6.1'):
        raise ImportError
except ImportError:
    err = ('NumPy version 1.6.1 or more recent '
           'must be installed to use Lilith interface')   
    raise LilithImportError(err)

try:
    import scipy
    if StrictVersion(scipy.__version__) < StrictVersion('0.9.0'):
        raise ImportError
except ImportError:
    err = ('SciPy version 0.9.0 or more recent '
           'must be installed to use Lilith interface')   
    raise LilithImportError(err)

# Lilith path
try:
    import lilith
except ImportError:
    try:
        Lilith_dir = config['Lilith_dir']
        sys.path.append(Lilith_dir)
        import lilith
    except KeyError:
        err = ('Could not find option "Lilith_dir" in Rosetta/config.txt')
        raise LilithImportError(err)
    except ImportError:
        err = ('check Lilith_dir option in '
               'Rosetta/config.txt')
        raise LilithImportError(err)
################################################################################