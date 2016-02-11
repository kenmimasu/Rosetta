################################################################################
class RosettaError(Exception):
    '''Rosetta base exception class. Never raised'''
    pass
    
class TranslationPathError(RosettaError):
    '''Exception raised when a suitable translation path is not found'''
    pass
    
class eHDECAYInterfaceError(RosettaError):
    '''Exception raised when a problem occurs within the eHDECAY interface'''
    pass

class LilithInterfaceError(RosettaError):
    '''Exception raised when a problem occurs within the Lililth interface'''
    pass
################################################################################
# if __name__=='__main__':
#     pass
        
