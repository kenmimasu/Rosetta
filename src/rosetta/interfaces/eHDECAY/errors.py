from ...internal.errors import RosettaImportError
from ..errors import RosettaInterfaceError,  RosettaInterfaceWarning

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

class eHDECAYNegativeWidthError(eHDECAYInterfaceError):
    '''Exception raised when total Higgs width is negative'''
    pass

class eHDECAYBrWarning(RosettaInterfaceWarning):
    interface='eHDECAY'
    nmax_before_suppress=None
    

class eHDECAYBrGtOneWarning(eHDECAYBrWarning):
    '''Exception raised when calculated Higgs BRs sum to more than one'''
    pass


class eHDECAYBrNegativeWarning(eHDECAYBrWarning):
    '''Exception raised when negative Higgs BRs are encountered'''
    nmax_before_suppress = None 
    pass
