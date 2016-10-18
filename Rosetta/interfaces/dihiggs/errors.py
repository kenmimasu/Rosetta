from ...internal.errors import RosettaImportError, RosettaInterfaceError

class dihiggsImportError(RosettaImportError):
    '''
    Exception raised when a problem occurs trying to import the 
    dihiggs interface package
    '''
    interface='dihiggs'
    pass
    
class dihiggsInterfaceError(RosettaInterfaceError):
    '''
    Exception raised when a problem occurs within the dihiggs 
    interface
    '''
    interface='dihiggs'
    pass

class SqrtsError(dihiggsInterfaceError):
    '''
    Exception for invalid value of sqrt(s) in TeV not in (7, 8, 13)
    '''
    interface='dihiggs'
    pass