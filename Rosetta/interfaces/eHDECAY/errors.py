from ...internal.errors import RosettaImportError
from ..errors import RosettaInterfaceError

class eHDECAYImportError(RosettaImportError):
    '''
    Exception raised when a problem occurs when trying to import the eHDECAY 
    interface package
    '''
    interface='eHDECAY'
    pass
    
class eHDECAYInterfaceError(RosettaInterfaceError):
    '''Exception raised when a problem occurs within the eHDECAY interface'''
    interface='eHDECAY'
    pass