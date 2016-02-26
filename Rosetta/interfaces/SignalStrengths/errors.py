from ...internal.errors import RosettaImportError, RosettaInterfaceError

class SignalStrengthsImportError(RosettaImportError):
    '''
    Exception raised when a problem occurs when trying to import the 
    SignalStrengths interface package
    '''
    interface='SignalStrengths'
    pass
    
class SignalStrengthsInterfaceError(RosettaInterfaceError):
    '''
    Exception raised when a problem occurs within the SignalStrengths 
    interface
    '''
    interface='SignalStrengths'
    pass

class SqrtsError(SignalStrengthsInterfaceError):
    '''
    Exception for invalid value of sqrt(s) in TeV not in (7, 8, 13)
    '''
    interface='SignalStrengths'
    pass