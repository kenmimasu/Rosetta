
from errors import EWPOImportError, EWPOInterfaceError
# check for NumPy >= 1.6.1 and Scipy >= 0.9.0
try:
    import numpy
except ImportError:
    err = ('NumPy must be installed to use EWPO interface')   
    raise EWPOImportError(err)

try:
    import scipy
except ImportError:
    err = ('SciPy must be installed to use EWPO interface')   
    raise EWPOImportError(err)
    
from interface import EWPOInterface
