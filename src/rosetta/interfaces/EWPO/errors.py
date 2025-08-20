from ...internal.errors import RosettaImportError
from ..errors import RosettaInterfaceError

class EWPOInterfaceError(RosettaInterfaceError):
    '''
    Exception raised when a problem occurs when trying to import EWPO interface 
    package
    '''
    interface='EWPO'
    pass

class EWPOImportError(RosettaImportError):
    '''
    Exception raised when a problem occurs when importing EWPO interface
    '''
    interface='EWPO'
    pass