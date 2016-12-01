from ...internal.errors import RosettaImportError
from ..errors import RosettaInterfaceError

class LilithImportError(RosettaImportError):
    '''
    Exception raised when a problem occurs when trying to import the Lilith 
    interface package
    '''
    interface='Lilith'
    pass
    
class LilithInterfaceError(RosettaInterfaceError):
    '''Exception raised when a problem occurs within the Lilith interface'''
    interface='Lilith'
    pass