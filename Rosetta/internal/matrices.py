import SLHA
from itertools import product
from collections import OrderedDict
################################################################################
# Special 2D Matrices
class TwoDMatrix(SLHA.NamedMatrix):
    
    def __keymap__(self, key):
        if key in self._names:
            return self.__parse__(key)
        else:
            return key
            
    def __valuemap__(self, key, value):
        return value

    def __init__(self, *args):
        if not args:
            matrix = SLHA.NamedMatrix()
        elif len(args)>1:
            raise TypeError('TwoDMatrix constructor takes at most one argument.')
        else:
            matrix=args[0]
        self.name = matrix.name
        self.fmt = matrix.fmt
        self.keytype = matrix.keytype
        self.cast = matrix.cast
        self._data = OrderedDict(matrix._data.items())
        self.preamble = matrix.preamble
        self._names = {k:v for k,v in matrix._names.items()}
        self._numbers = SLHA.CaseInsensitiveDict({k:v for k,v in 
                                                  matrix._numbers.items()})
        for k,v in self.iteritems():
            self[k] = v
        # assert len(self.dimension())==2, 'TwoDMatrix must have dimension 2'

    def __setitem__(self, key, value):
        super(TwoDMatrix, self).__setitem__(self.__keymap__(key), 
                                            self.__valuemap__(key, value))
                                                 
    def __getitem__(self, key):
        value = super(TwoDMatrix, self).__getitem__(self.__keymap__(key))
        return self.__valuemap__(key, value)

    def __contains__(self, key):
        try:
            return super(TwoDMatrix, self).__contains__(self.__keymap__(key))
        except (ValueError, KeyError):
            return False

    def __delitem__(self, key):
        return super(TwoDMatrix, self).__delitem__(self.__keymap__(key))

class CTwoDMatrix(TwoDMatrix, SLHA.CNamedMatrix):
    
    def __keymap__(self, key):
        return self.__parse__(key)
    
    def __setreal__(self, real):
        return TwoDMatrix(real)
    
    def __setimag__(self, imag):
        return TwoDMatrix(imag)

    def __init__(self, *args):
        if not args:
            matrix = SLHA.CNamedMatrix()
        elif len(args)>1:
            raise TypeError('CTwoDMatrix constructor takes at most one argument.')
        else:
            matrix=args[0]
        self._re = self.__setreal__(matrix._re)
        self._im = self.__setimag__(matrix._im)
        super(CTwoDMatrix, self).__init__(matrix)

class SymmetricMatrix(TwoDMatrix):

    def __keymap__(self, key):
        key = super(SymmetricMatrix, self).__keymap__(key)
        i, j = key
        return (i, j) if i <= j else (j, i)
        
    def __valuemap__(self, key, value):
        return value

class AntisymmetricMatrix(TwoDMatrix):

    def __keymap__(self, key):
        key = super(AntisymmetricMatrix, self).__keymap__(key)
        i, j = key
        return (i, j) if i <= j else (j, i)
        
    def __valuemap__(self, key, value):
        i, j = key
        return 0. if i==j else -value if i > j else value

class CSymmetricMatrix(SymmetricMatrix, CTwoDMatrix):
    
    def __setreal__(self, real):
        return SymmetricMatrix(real)
    
    def __setimag__(self, imag):
        return SymmetricMatrix(imag)
        
class CAntisymmetricMatrix(CTwoDMatrix, AntisymmetricMatrix):
    
    def __setreal__(self, real):
        return AntisymmetricMatrix(real)
    
    def __setimag__(self, imag):
        return AntisymmetricMatrix(imag)
        
class HermitianMatrix(CTwoDMatrix):
    
    def __setreal__(self, real):
        return SymmetricMatrix(real)
    
    def __setimag__(self, imag):
        return AntisymmetricMatrix(imag)
        
    def __keymap__(self, key):
        key = self.__parse__(key)
        i, j = key
        return (i, j) if i <= j else (j, i)
        
    def __valuemap__(self, key, value):
        key = self.__parse__(key)
        i, j = key
        if i==j:
            return complex(value.real, 0.)
        elif i > j:
            return value.conjugate()
        else:
            return value

        
################################################################################
# Matrix operations
def matrix_mult(A, B, assign=None):
    '''
    Matrix operation A.B:
    Perform matrix multiplication of A and B. If the assign option is given, the 
    result will be assigned to its elements, otherwise, an object of identical 
    type to A will be returned.
    Objects should be SLHA.Matrix instances or be indexable as [i,j,..] and 
    posess the dimension() method which returns a tuple of array dimensions. 
    Assumes indexing from 1.
    '''
    
    dimA, dimB = A.dimension(), B.dimension()
    
    # matrix shape checks
    if dimA[-1] != dimB[0]:
        err = ("Size of last dimension of array A doesn't "
               "match the first of array B ")
        raise IndexError(err)
            
    if assign is None:
        dimC = dimA[:-1] + dimB[1:]
        if isinstance(A, SLHA.CBlock) or isinstance(B, SLHA.CBlock):
            C = CTwoDMatrix()
        else:
            C = TwoDMatrix()
    else:
        C = assign
        dimC = C.dimension()
        if dimC != dimA[:-1] + dimB[1:]:
            err = ("A.B -> C: Shape of array C doesn't match "
                   "the dimension of the result of A.B")
            raise IndexError(err)
                                           
    dim = dimA[-1]
    
    for k in product(*[range(1,x+1) for x in dimC]):
        C[k] = sum( [ A[k[:-1]+(x,)]*B[(x,)+k[1:]] 
                      for x in range(1,dim+1) ] )
        
    return C

def matrix_add(A, B, assign=None):
    '''
    Matrix operation A + B:
    Perform element-wise matrix addition of A and B. If a third argument is  
    given, the result will be assigned to its elements, otherwise, an object of 
    identical type to A will be returned.
    Objects should be SLHA.Matrix instances or be indexable as [i,j,..] and 
    posess the dimension() method which returns a tuple of array dimensions. 
    Assumes indexing from 1.
    '''
    
    dimA, dimB = A.dimension(), B.dimension()
    
    # matrix shape checks
    if dimA != dimB:
        err = ("Dimension of array A doesn't "
               "match that of array B")
        raise IndexError(err)
            
    if assign is None:
        dimC = dimA[:-1] + dimB[1:]
        if isinstance(A, SLHA.CBlock) or isinstance(B, SLHA.CBlock):
            C = CTwoDMatrix()
        else:
            C = TwoDMatrix()
    else:
        C = assign
        dimC = C.dimension()
        if dimC != dimA:
            err = ("A + B -> C: Shape of array C doesn't match "
                   "the dimension of the result of A + B")
            raise IndexError(err)
                                           
    dim = dimA[-1]
    
    for k in product(*[range(1,x+1) for x in dimC]):
        C[k] = A[k] + B[k]
        
    return C

def matrix_sub(A, B, assign=None):
    '''
    Matrix operation A - B:
    Perform element-wise matrix subtraction of A and B. If a third argument is  
    given, the result will be assigned to its elements, otherwise, an object of 
    identical type to A will be returned.
    Objects should be SLHA.Matrix instances or be indexable as [i,j,..] and 
    posess the dimension() method which returns a tuple of array dimensions. 
    Assumes indexing from 1.
    '''
    
    dimA, dimB = A.dimension(), B.dimension()
    
    # matrix shape checks
    if dimA != dimB:
        err = ("Dimension of array A doesn't "
               "match that of array B")
        raise IndexError(err)
            
    if assign is None:
        if isinstance(A, SLHA.CBlock) or isinstance(B, SLHA.CBlock):
            C = CTwoDMatrix()
        else:
            C = TwoDMatrix()
    else:
        C = assign
        dimC = C.dimension()
        if dimC != dimA:
            err = ("A - B -> C: Shape of array C doesn't match "
                   "the dimension of the result of A + B")
            raise IndexError(err)
                                           
    dim = dimA[-1]
    
    for k in product(*[range(1,x+1) for x in dimC]):
        C[k] = A[k] - B[k]
        
    return C

def matrix_eq(A, B):
    '''
    Matrix operation B -> A:
    Assigns the values of A element-wise to B.
    Objects should be SLHA.Matrix instances or be indexable as [i,j,..] and 
    posess the dimension() method which returns a tuple of array dimensions. 
    Assumes indexing from 1.
    '''
    
    dimA, dimB = A.dimension(), B.dimension()
    
    # matrix shape checks
    if dimA != dimB:
        err = ("Dimension of array A doesn't "
               "match that of array B")
        raise IndexError(err)
                                               
    for k in product(*[range(1,x+1) for x in dimA]):
        B[k] = A[k]

if __name__=='__main__':
    pass