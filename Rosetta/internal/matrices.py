import SLHA
from itertools import product
from collections import OrderedDict
import re

__doc__ = '''
Module implementing special cases of 2D matrices inheriting from the 
SLHA.Matrix structures in SLHA.py e.g. Hermitian, Symmetric etc. Some basic 
operations are also defined here for matrix multiplication, addition 
subtraction and element-wise assignment.
'''
################################################################################
# Special 2D Matrices

class TwoDMatrix(SLHA.NamedMatrix):
    mask = tuple()
    def __strkeymap__(self, key):
        if not re.match(r'.*\dx\d$', key):
            err = ('{}.__strkeymap__: Key {} '.format(self.__class__, key)
                  +'does not have format NAMEixj.')
            raise KeyError(err)
        intkeys = self.__parse__(key)
        i, j = self.__keymap__(intkeys)
        try:
            # return self._names[(j,i)]
            return self._names[(i,j)]
        except KeyError:
            # return key[:-3]+'{}x{}'.format(j,i)
            return key[:-3]+'{}x{}'.format(i,j)
        
    def __keymap__(self, key):
        '''
        Returns the modified key according to the properties of the matrix i.e. 
        if it is symmetric etc.
        '''
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
                                                  
        for k,v in self.items()[::-1]:
            self[k] = v

    def __setitem__(self, key, value):
        intkey = self.get_number(key) if type(key) is str else key
        if intkey in self.mask:
            if type(key) is str:
                newkey = self.__strkeymap__(key) 
            else:
                newkey = self.__keymap__(intkey)
        else:
            newkey = key
        super(TwoDMatrix, self).__setitem__(newkey,
                                            self.__valuemap__(intkey, value))
                                                 
    def __getitem__(self, key):
        intkey = self.get_number(key) if type(key) is str else key
        if intkey in self.mask:
            if type(key) is str:
                newkey = self.__strkeymap__(key) 
            else:
                newkey = self.__keymap__(key)
        else:
            newkey = key
        value = super(TwoDMatrix, self).__getitem__(newkey)
        return self.__valuemap__(intkey, value)

    def __contains__(self, key):
        try:
            intkey = self.get_number(key) if type(key) is str else key
            if intkey in self.mask:
                if type(key) is str:
                    newkey = self.__strkeymap__(key) 
                else:
                    newkey = self.__keymap__(key)
            else:
                newkey = key
            return super(TwoDMatrix, self).__contains__(newkey)
        except (ValueError, KeyError) as e:
            return False

    def __delitem__(self, key):
        key =  self.get_number(key) if type(key) is str else key
        return super(TwoDMatrix, self).__delitem__(key)
    
    def T(self):
        new = self.__class__(self)
        # transposed data
        newdata = OrderedDict([((j,i),v) for (i,j),v in new._data.items()])
        new._data = newdata
        return new

class CTwoDMatrix(TwoDMatrix, SLHA.CNamedMatrix):
    mask = tuple()    
    def __strkeymap__(self, key):
        if not re.match(r'.*\dx\d$', key):
            err = ('{}.__strkeymap__: Key {} '.format(self.__class__, key)
                  +'does not have format NAMEixj.')
            raise KeyError(err)
        
        if key in self._re or key in self._im:
            part = self._part(key)
        else:
            part = self
        
        intkeys = part.__parse__(key)
        i, j = self.__keymap__(intkeys)
        try:
            # return part._names[(j,i)]
            return part._names[(i,j)]
        except KeyError:
            return key[:-3]+'{}x{}'.format(i,j)
        
    def __keymap__(self, key):
        return key
            
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
        
    def __setitem__(self, key, value):
        return super(CTwoDMatrix, self).__setitem__(key, value)
        
    def __getitem__(self, key):
        if type(key) is str and (key in self._re or key in self._im):
            return self._part(key).__getitem__(key)
        else:
            return super(CTwoDMatrix, self).__getitem__(key)
    
    def __delitem__(self, key):
        return super(CTwoDMatrix, self).__delitem__(key)
    
    def __contains__(self, key):
        return super(CTwoDMatrix, self).__contains__(key)
        
    def T(self):
        new = self.__class__(self)
        # transposed data
        newdata = OrderedDict([((j,i),v) for (i,j),v in new._data.items()])
        re_data = OrderedDict([((j,i),v) for (i,j),v in new._re._data.items()])
        im_data = OrderedDict([((j,i),v) for (i,j),v in new._im._data.items()])
        new._data = newdata
        new._re._data = re_data
        new._im._data = im_data
        return new
    
    def dag(self):
        new = self.__class__(self)
        # transposed data
        newdata = OrderedDict([((j,i),v.conjugate()) 
                               for (i,j),v in new._data.items()])
        re_data = OrderedDict([((j,i),v) for (i,j),v in new._re._data.items()])
        im_data = OrderedDict([((j,i),-v) for (i,j),v in new._im._data.items()])
        new._data = newdata
        new._re._data = re_data
        new._im._data = im_data
        return new
        
class SymmetricMatrix(TwoDMatrix):
    mask = ((2,1),(3,1),(3,2))
    def __keymap__(self, key):
        key = super(SymmetricMatrix, self).__keymap__(key)
        i, j = key
        return (i, j) if i <= j else (j, i)
        
    def __valuemap__(self, key, value):
        return value

class AntisymmetricMatrix(TwoDMatrix):
    mask = ((2,1),(3,1),(3,2))
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
    mask = ((2,1),(3,1),(3,2))
    def __setreal__(self, real):
        return SymmetricMatrix(real)
    
    def __setimag__(self, imag):
        return AntisymmetricMatrix(imag)
        
    def __keymap__(self, key):
        i, j = key
        return (i, j) if i <= j else (j, i)
        
    def __valuemap__(self, key, value):
        i, j = key
        if i==j:
            return complex(value.real, 0.)
        elif i > j:
            return value.conjugate()
        else:
            return value
    def dag(self):
        return self

        
################################################################################
# Matrix operations
def matrix_mult(A, B, assign=None):
    '''
    Matrix operation A.B:
    Perform matrix multiplication of A and B. If the assign option is given, the 
    result will be assigned to its elements, otherwise, an TwoDMatrix or 
    CTwoDMatrix will be returned depending on the type of Matrices A & B.
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
        # dimC = dimA[:-1] + dimB[1:]
        dimC = dimA
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
        dimC = dimA
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
    
