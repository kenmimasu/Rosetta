from ...internal.errors import RosettaImportError
from ..errors import RosettaInterfaceError,  RosettaInterfaceWarning

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


class SignalStrengthsNegativeWidthError(SignalStrengthsInterfaceError):
    '''Exception raised when total Higgs width is negative'''
    pass

class SignalStrengthsBrWarning(RosettaInterfaceWarning):
    interface='SignalStrengths'
    nmax_before_suppress=None
    

class SignalStrengthsBrGtOneWarning(SignalStrengthsBrWarning):
    '''Exception raised when calculated Higgs BRs sum to more than one'''
    pass


class SignalStrengthsBrNegativeWarning(SignalStrengthsBrWarning):
    '''Exception raised when negative Higgs BRs are encountered'''
    nmax_before_suppress = None 
    pass
