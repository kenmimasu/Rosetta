################################################################################
import StringIO
import re
import sys
import os
import math
from collections import namedtuple, OrderedDict, MutableMapping
from itertools import product, combinations
from itertools import combinations_with_replacement as combinations2
################################################################################
from .. import SLHA
from ..SLHA import CaseInsensitiveDict
import checkers as check
from ..matrices import (TwoDMatrix, CTwoDMatrix, HermitianMatrix, 
                      SymmetricMatrix, AntisymmetricMatrix)
from ..constants import (PID, default_inputs, default_masses, input_names, 
                       default_ckm, particle_names, input_to_PID, 
                       PID_to_input, GammaZ, GammaW, Gammat)
from ..errors import TranslationError, TranslationPathError
from errors import FlavorMatrixError, ParamCardReadError
from .. import session
from io import read_param_card, write_param_card
from translator import translate
################################################################################
__doc__ = '''
Base class for Rosetta bases as well as some utility functions for defining the 
names of flavor block elements and sorting blocks names for writing SLHA output.
'''
# Base Basis class
class Basis(MutableMapping):
    '''
    Base class from which to derive other Higgs EFT basis classes. 
    
    Designed to be instantiated with an SLHA format parameter card using
    read_param_card(). The contents of the card are stored as an SLHA.Card 
    object and the special blocks "mass" and "sminputs" are stored as 
    SLHA.NamedBlock instances in self.mass and self.input respectively.
    
    self.name   - Unique basis identifier. This will be compared to the 0th 
                  element of the block basis in the SLHA parameter card to 
                  ensure that right basis class is being instantiated for a 
                  given input card.
    self.card   - SLHA.Card instance containing SLHA.NamedBlock and SLHA.Decay 
                  instances corresponding to those specified in the parameter 
                  card. Object names taken to be the first non-whitespace 
                  characters after a "#" character in a block or parameter 
                  definition are also stored.
    self.mass   - SLHA.NamedBlock instance for the "mass" block
    self.inputs - SLHA.NamedBlock instance for the "sminputs" block
    self.ckm    - Rosetta.matrices.CTwoDMatrix instance for the VCKM and 
                  IMVCKM blocks

    self.blocks, self.required_inputs and self.required_masses should be defined 
    in accordance with block structure of the SLHA parameter card, Blocks 
    "sminput" and "mass" respectively (the Z and Higgs masses are also stored 
    in self.inputs). Specifically, the block names specified as keys in 
    self.blocks should match those in the SLHA card as well as the names and 
    indices of the elements (taken to be the order in which the appear in the 
    list associated the the block in self.blocks). The blocks "mass" and 
    "sminputs" should minimally contain the entries specified in required_masses 
    and required_inputs respectively. A number of checks related to these 
    definitions are performed by check_param_data(), check_mass() and 
    check_sminputs() on the data read in from the parameter card. An example is 
    shown below where the required data members are defined outside the class 
    constructor so as to intrinsically belong to all instances of the class.
    
        >> from internal import Basis
        >> class MyBasis(Basis.Basis):
        >>    name = 'mybasis'
        >>    independent = {'A','B','C','1','2','3'}
        >>    required_inputs = {1,2,4}   # a_{EW}^{-1}, Gf and MZ required
        >>    required_masses = {23,25,6} # Z, Higgs and top masses required
        >>    blocks = {'letters':['A','B','C','D']
        >>              'numbers':['1','2','3','4']} # Expected block structure
    
    The list self.independent stored the basis parameters which should be read 
    in. Any other parameters declared in self.blocks are assumed to be set by 
    the user in the calculate_dependent() method. write_param_card() writes the 
    contents of the self.newcard into a new SLHA formatted file.
    
    Basis and any of its subclasses are designed to work similarly to a 
    dictionary in that parameter values can be referenced by name (duplicate 
    names in different blocks are not handled properly so try to avoid them). 
    A value can be referenced in various ways, see below example where the 
    parameter named 'D' is stored as entry 3 in the block 'letters' written in 
    'mycard.dat':
        
        >> instance = MyBasis(param_card='mycard.dat')
        >> instance['A'] = 0.5 # set value according to name in SLHA card
        >> instance.card['A'] = 0.5 # equivalent to the above
        >> instance.card.blocks['letters']['A'] = 0.5 # from block by name 
        >> instance.card.blocks['letters'][3] = 0.5 # from block by entry 
        
    The user can define calculate_dependent() to set any dependent parameters 
    (those listed in in self.blocks but not in self.independent) and any number 
    of additional functions, usually to translate the coefficients into a given 
    basis which sets the SLHA.Card object, self.newcard. Rosetta differentiates 
    from general utility functions and translation functions using the 
    "translation" decorator. Any function designed to take you to another 
    existing basis implementation should be decorated with the "translation" 
    decorator with an argument corresponding to the name of the target basis. 
    This name should match the unique name of an existing basis implementation 
    also contained in the Rosetta root directory. Below would be example of a 
    translation function from our example to the Warsaw Basis.
    
        >> @Basis.translation('warsaw') # identifies as translator to warsaw 
        >> def mytranslation(self, instance):
        >>     instance['cWW'] = 10.
        >>     instance['cpHl11'] = self['myparam']
        >>     return instance
    
    Any python module saved in the bases/ directory of the Rosetta package is 
    assumed to be a basis implementation. Furthermore, the class name of the 
    basis implementation should be the same as the file name of the python 
    module. In our example above, the class `MyBasis` would have to be saved in 
    a file named MyBasis.py. This is to help Rosetta automatically identify the 
    possible translation paths between bases. Any such class can then be used 
    by the command line script "translate". For example, Rosetta should be able 
    to figure out all possible multi-step translations.
    '''
    
    blocks = dict()
    independent = []
    # dependent = []
    flavored = dict()
    required_inputs, required_masses = set(), set()
    translate = translate
    
    def __init__(self, param_card=None, flavor = 'general', dependent=True):
        '''
        Keyword arguments:
            param_card   - SLHA formatted card conforming to the definitions in 
                           self.blocks, self.required_masses and 
                           self.required_inputs.
            flavor      - flavor structure of matrices: 'diagonal, 'universal' 
                           , 'MFV' or 'general'.
            dependent    - when param_card is None, whether or not to include 
                           dependent parameters in the SLHA.Card attribute of 
                           the basis instance.
        '''
        self.translations = CaseInsensitiveDict()
        
        self.flavor = flavor

        self.param_card = param_card

        if not hasattr(self, 'dependent'): self.dependent=[]

        self.set_dependents()

        self.set_fblocks(self.flavor)

        # read param card (sets self.inputs, self.mass, self.name, self.card)
        if param_card is not None: 
            if not os.path.exists(self.param_card):
                err = '{} does not exist!'.format(self.param_card)
                raise ParamCardReadError(err)
            
            read_param_card(self) 

            # various input checks
            check.sminputs(self, self.required_inputs) 
            check.masses(self, self.required_masses) 
            check.param_data(self)
            check.flavored_data(self)
            
            # generalises potentially reduced flavor structure
            self.set_flavor(self.flavor, 'general')
            
            # add dependent coefficients/blocks to self.card
            self.init_dependent()
            
            # generate internal OrderedDict() object for __len__, __iter__, 
            # items() and iteritems() method
            self._gen_thedict()
            
            # user defined function
            self.calculate_dependent()

        # if param_card option not given, instantiate with class name 
        # and all coeffs set to 0 (used for creating an empty basis 
        # instance for use in translate() method)
        else: 
            self.card = self.default_card(dependent=dependent) 
            self.inputs, self.mass = None, None           
            self._gen_thedict()
            
    # overloaded container (dict) methods for indexing etc.
    
    def __getitem__(self, key):
        try:
            return self.card.__getitem__(key)
        except KeyError:
            if hasattr(self,'_thedict'):
                return  self._thedict[key]
            else:
                raise KeyError

    def __setitem__(self, key, value):
        if hasattr(self,'_thedict'):
            self._thedict[key]=value
        return self.card.__setitem__(key, value)
    
    def __contains__(self, key):
        return self.card.__contains__(key)
        
    def __delitem__(self, key):
        if hasattr(self,'_thedict'):
            del self._thedict[key]
        return self.card.__delitem__(key)
    
    def __len__(self):
        return len(self._thedict)
        
    def __iter__(self):
        return iter(self._thedict)
    
    def _gen_thedict(self):
        # thedict = SLHA.CaseInsensitiveOrderedDict()
        thedict = CaseInsensitiveDict()
        for name, blk in self.card.blocks.iteritems():
            if name in self.blocks:
                for k, v in blk.iteritems():
                    thedict[blk.get_name(k)] = v
                    if isinstance(blk, SLHA.CBlock):
                        thedict[blk._re.get_name(k)] = v.real
                        thedict[blk._im.get_name(k)] = v.imag
                        
        for name, blk in self.card.matrices.iteritems():
            if name in self.fblocks:
                for k, v in blk.iteritems():
                    cname = blk.get_name(k)
                    if cname: thedict[cname] = v
                    if isinstance(blk, SLHA.CMatrix):
                        re_name = blk._re.get_name(k)
                        im_name = blk._im.get_name(k)
                        if re_name: thedict[re_name] = v.real
                        if im_name: thedict[im_name] = v.real
                        
        self._thedict = thedict
        
    
    def set_dependents(self):
        '''
        Populate self.independent and self.dependent lists according to 
        basis class definition.
        '''
        self.all_coeffs = [c for v in self.blocks.values() for c in v]
        # remove overlaps
        self.independent = [c for c in self.independent 
                            if c not in self.dependent]

        self.dependent.extend([c for c in self.all_coeffs if (c not in
                               self.independent and c not in self.dependent)])

        # check for block names in independent
        for k,v in self.blocks.iteritems():
            if k in self.independent and k not in self.dependent:
                for fld in v:
                    if fld not in self.independent:
                        self.independent.append(fld)
                    try:
                        self.dependent.remove(fld)
                    except ValueError:
                        pass
                        
        for name, opt in self.flavored.iteritems():
            coeffs = flavor_coeffs(name, **opt)

            if name not in (self.independent + self.dependent):
                dependents, independents = [], []
                for c in coeffs:
                    if c in self.dependent:
                        self.dependent.remove(c)
                        if c in self.independent: 
                            independents.append(c)
                        else:
                            dependents.append(c)
                    elif c in self.independent:
                        independents.append(c)
                    else:
                        dependents.append(c)
                
                if not independents:
                    self.dependent.append(name)
                elif not dependents:
                    self.independent.append(name)
                else:
                    self.dependent.extend(dependents)
                    self.independent.extend(independents)
        
    def set_fblocks(self, option='general'):
        self.fblocks = dict()
        for name, opt in self.flavored.iteritems():
            opt['flavor'] = option
            coeffs = flavor_coeffs(name, **opt)
            self.fblocks[name] = coeffs

    def default_card(self, dependent=True):
        '''
        Create a new default SLHA.Card instance according to the self.blocks 
        and self.flavored structure of the basis class specified in the 
        implementaion. The dependent option allows one to switch on or off the 
        inclusion of dependent parameters in the default card. By default they 
        are included.
        '''
        thecard = SLHA.Card(name=self.name)
        thecard.add_entry('basis', 1, self.name, name = 'translated basis')

        preamble = ('\n###################################\n'
            + '## INFORMATION FOR {} BASIS\n'.format(self.name.upper())
            + '###################################\n')
            
        thecard.blocks['basis'].preamble=preamble

        # default behaviour: create one 'newcoup' block, ignoring flavored
        if not self.blocks: 
            omit = ([] if not self.flavored else 
                    [c for (k,v) in self.fblocks.items() for c in v+[k]])
            
            all_coeffs = [c for c in self.independent if c not in omit]
            self.blocks = {'newcoup':all_coeffs}
            
        # otherwise follow self.blocks structure
        for blk, flds in self.blocks.iteritems():                
            for i,fld in enumerate(flds):                    
                if dependent or (blk not in self.dependent 
                                 and fld not in self.dependent):
                    thecard.add_entry(blk, i+1, 0., name = fld)

        # deal with flavored
        for blk, flds in self.fblocks.iteritems():
            for fld in flds:
                index = (int(fld[-3]), int(fld[-1])) # XBcoeff(I)x(J)               
                if dependent or (blk not in self.dependent 
                                 and fld not in self.dependent):
                    if self.flavored[blk]['domain']=='complex':              
                        thecard.add_entry(blk, index, 0., name = 'R'+fld[1:])
                        thecard.add_entry('IM'+blk, index, 0., name = 'I'+fld[1:])
                    else:
                        thecard.add_entry(blk, index, 0., name = fld)
        
        vckm = default_ckm
        thecard.add_block(vckm)
        thecard.ckm = vckm
        thecard.set_complex()
        self.fix_matrices(card = thecard)
        
        for ptcl, wid in zip((23, 24, 6),(GammaZ, GammaW, Gammat)):
            thecard.add_decay(SLHA.Decay(ptcl, wid))
        
        return thecard
    
    def set_flavor(self, _from, to):
        if _from == to: return

        self.set_fblocks(to) # reset fblocks according to flavor option
        # newcard = self.default_card(dependent=False)
        newcard = self.default_card(dependent=True)
        
        if (_from, to) in (('general', 'universal'), ('general', 'diagonal')):
            blks_to_del = []
            for bname, blk in self.card.matrices.iteritems():
                # only consider declared flavor matrices
                if bname not in self.flavored: continue
                
                # only consider independent blocks
                if not newcard.has_matrix(bname):
                    blks_to_del.append(bname)
                    continue
                
                # delete elements not present in default card
                to_del = []
                no_del = []
                for k in blk.keys():
                    cname = blk.get_name(k)
                    if k not in newcard.matrices[bname]:
                        if abs(blk[k]) < 1e-6:
                            to_del.append(k)
                        else:                            
                            no_del.append(k)
                
                # only keep 2,2 and 3,3 elements as anomalous if they 
                # sufficiently different from the 1,1 element in the universal 
                # case.
                if no_del and (_from, to) == ('general', 'universal'): 
                    for diag in ((2,2), (3,3)):
                        if diag in no_del:
                            val = abs(blk[diag])
                            if val - abs(blk[1,1]) < 1e-4*val:
                                no_del.remove(diag)
                                to_del.append(diag)
                
                for k in to_del:
                    del blk[k]

                if no_del:
                    no_del_names = [blk.get_name(k) for k in no_del]
                    no_del_values = [blk[k] for k in no_del]
                    session.verbose(
                    '    Warning in {}.set_flavour():\n'.format(self.__class__)+
                    '    Reduction in flavour structure ' +
                    'from "{}" to "{}" '.format(_from, to) +
                    'encountered some unexpected non-zero elements ' +
                    ' which were not deleted.\n    Not deleted: ' +
                    '{}\n'.format(
                    ', '.join(['{}={}'.format(x,y) for x,y in 
                               zip(no_del_names,no_del_values)])
                    )
                    )
                
            for blk in blks_to_del:
                del self.card.matrices[blk]
                        
        elif to=='general':
            for bname, blk in newcard.matrices.iteritems():
                # only consider declared flavor matrices
                if bname not in self.flavored: continue
                # Add blocks absent in self.card but present in default card
                if not self.card.has_matrix(bname):
                    self.card.add_block(blk)
                    continue
                
                oneone = self.card.matrices[bname].get((1,1), 0.)
                for k, v in blk.iteritems():
                    # only add value if element doesn't already exist in block
                    if k in self.card.matrices[bname]: continue

                    diag = (k[0] == k[1])
                    if diag: v = oneone

                    cname = blk.get_name(k)
                    self.card.add_entry(bname, k, v, name = cname)


    def init_dependent(self):
        '''
        Adds entries defined as dependent to the corresponding block of 
        self.card so that they can be assigned values in calculate_dependent().
        '''
        for bname, fields in self.blocks.iteritems():
            theblock = self.card.blocks.get(bname,[])
            to_add = [f for f in fields if f in self.dependent 
                                        and f not in theblock]
            for entry in to_add:
                self.card.add_entry(bname, fields.index(entry)+1, 0., name=entry)
        
        for bname, fields in self.fblocks.iteritems():
            theblock = self.card.matrices.get(bname,[])
            to_add = [f for f in fields if 
                      (f in self.dependent or bname in self.dependent)
                      and f not in theblock]
            
            for entry in to_add:
                index = (int(entry[-3]), int(entry[-1]))
                if self.flavored[bname]['domain']=='complex':      
                    theblock = self.card.matrices.get(bname, None)
                    if not (isinstance(theblock, SLHA.CMatrix) 
                         or isinstance(theblock, SLHA.CNamedMatrix)):
                        self.card.add_entry(bname, index, 0., 
                                            name='R'+entry[1:])
                        self.card.add_entry('IM'+bname, index, 0.,
                                            name='I'+entry[1:])
                    else:
                        self.card.add_entry(bname, index, complex(0.,0.), 
                                            name=entry)
                else:
                    value = 0.
        
        self.card.set_complex()
        self.fix_matrices()
        

    def fix_matrices(self, card=None):
        if card is None:
            card = self.card
        for name, matrix in card.matrices.iteritems():
            if name.lower() == 'vckm':
                card.matrices['vckm'] = CTwoDMatrix(matrix)
            elif name not in self.flavored:
                continue
            else:
                opt = self.flavored[name]
                kind, domain = opt['kind'], opt['domain']
                if (kind, domain) == ('hermitian', 'complex'):
                    MatrixType = HermitianMatrix
                elif (kind, domain) == ('symmetric', 'complex'):
                    MatrixType = CSymmetricMatrix
                elif ((kind, domain) == ('hermitian', 'real') or
                      (kind, domain) == ('symmetric', 'real')):
                    MatrixType = SymmetricMatrix
                elif (kind, domain) == ('general', 'real'):
                    MatrixType = TwoDMatrix                    
                elif (kind, domain) == ('general', 'complex'):
                    MatrixType = CTwoDMatrix
                else: 
                    continue
                card.matrices[name] = MatrixType(matrix)
    
    def reduce_hermitian_matrices(self):
        '''
        Deletes imaginary parts of the diagonal elements of HermitianMatrix 
        instances belonging to the basis instance.
        '''
        for matrix in self.card.matrices.values():
            if isinstance(matrix, HermitianMatrix):
                for i,j in matrix.keys():
                    if i==j: del matrix._im._data[i,j]
    
    def delete_dependent(self):
        '''
        Deletes all named coefficients present in self.dependent from their 
        respective containers.
        '''
        for container in (self.card.blocks, self.card.matrices):
            for name, blk in container.iteritems():
                if name in self.dependent:
                    del container[name]
                else:    
                    for n in blk._numbers:
                        if n in self.dependent:
                            del blk[n]
                    if len(blk)==0:
                        del container[name]
        
    def calculate_dependent(self):
        '''
        Default behaviour of calculate_dependent(). Called if a subclass of 
        Basis has not implemented the function.
        '''
        session.verbose('Nothing done for {}.calculate_'\
              'dependent()\n'.format(self.__class__.__name__))
              
    def modify_inputs(self):
        '''
        Default behaviour of modify_inputs(). Called if a subclass of 
        Basis has not implemented the function.
        '''
        session.verbose('Nothing done for {}.modify_'\
              'inputs()\n'.format(self.__class__.__name__))
        
################################################################################
    
def flavor_coeffs(name, kind='hermitian', domain='real', flavor='general', 
                         cname = None):
    '''
    Function to create flavor components of a coefficient according to its 
    properties. Takes a parameter name as an argument and returns a tuple of 
    lists corresponding to the real and imaginary parts of the matrix elements. 
    The naming convention for real coefficients is to suffix the coefficient 
    name with the matrix indices separates by "x". 
    '''
    index = (1, 2, 3)
    if cname is None: cname = name
    
    if flavor.lower() == 'diagonal':
        include = lambda i,j: i == j 
    elif flavor.lower() in ('universal','mfv'):
        include = lambda i,j: i == 1 and j == 1
    else:
        include = lambda i,j: True
        
    if (kind, domain) == ('hermitian', 'complex'):
        real = ['{0}{1}x{1}'.format(cname,i) for i in index if include(i,i)]
        cplx = ['{}{}x{}'.format(cname,i,j) for i,j in 
                combinations(index,2) if include(i,j)]

    elif (kind, domain) == ('symmetric', 'complex'):
        real = []
        cplx = ['{}{}x{}'.format(cname,i,j) for i,j in 
                combinations2(index,2) if include(i,j)]
        
    elif ((kind, domain) == ('hermitian', 'real') or
          (kind, domain) == ('symmetric', 'real')):
        real = ['{}{}x{}'.format(cname,i,j) for i,j in 
                combinations2(index,2) if include(i,j)]
        cplx = []
        
    elif (kind, domain) == ('general', 'real'):
        real = ['{}{}x{}'.format(cname,i,j) for i,j in 
                product(index,index) if include(i,j)]
        cplx = []
        
    elif (kind, domain) == ('general', 'complex'):
        real = []
        cplx = ['{}{}x{}'.format(cname,i,j) for i,j in 
                product(index,index) if include(i,j)]
        
    else:
        err = ('flavor_matrix function got and unrecognised combination of '
               '"kind" and "domain" keyword arguments')
        raise FlavorMatrixError(err)
        # return [name]
    
    if (not cplx and domain!='complex'):
        return real
    else:
        return ['C'+c for c in real+cplx]

    