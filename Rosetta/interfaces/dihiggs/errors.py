from ...internal.errors import RosettaImportError
from ..errors import RosettaInterfaceError

class DiHiggsImportError(RosettaImportError):
    '''
    Exception raised when a problem occurs trying to import the 
    dihiggs interface package
    '''
    interface='dihiggs'
    pass
    
class DiHiggsInterfaceError(RosettaInterfaceError):
    '''
    Exception raised when a problem occurs within the dihiggs 
    interface
    '''
    interface='dihiggs'
    pass

class SqrtsError(DiHiggsInterfaceError):
    '''
    Exception for invalid value of sqrt(s) in TeV not in (7, 8, 13, 14, 100)
    '''
    interface='dihiggs'
    pass