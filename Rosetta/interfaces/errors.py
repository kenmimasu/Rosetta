from ..internal.errors import RosettaError

class RosettaInterfaceError(RosettaError):
    '''Base error class for Rosetta interfaces.'''
    interface=''
    def __init__( self, msg ):
        super(RosettaInterfaceError, self).__init__(
             'Error in Rosetta {} interface: {}'.format(self.interface, msg))
    pass

class ReadParamCardError(RosettaError):
    '''Exception raised inside RosettaInterface.read_param_card()'''
    pass

class LoadInterfaceError(RosettaError):
    '''Base error class for Rosetta interfaces.'''
    # def __init__( self, msg ):
    #     super(LoadInterfaceError, self).__init__(
    #          'Error loading Rosetta interface: {}'.format(self.interface, msg))
    pass