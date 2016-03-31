from ..errors import RosettaError

class FlavorMatrixError(RosettaError):
    '''Error in calling flavor_matrix()'''
    pass

class ParamCardReadError(RosettaError):
    '''Error in reading parameter card'''
    pass

class BasisNameError(RosettaError):
    '''Error reading basis name in parameter card'''
    pass