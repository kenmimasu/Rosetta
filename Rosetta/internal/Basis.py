################################################################################
import StringIO
import re
import sys
import os
import math
import datetime
from collections import namedtuple, OrderedDict, MutableMapping
from itertools import product, combinations
from itertools import combinations_with_replacement as combinations2
################################################################################
import SLHA
from query import query_yes_no as Y_or_N
from decorators import translation
from matrices import (TwoDMatrix, CTwoDMatrix, HermitianMatrix, 
                      SymmetricMatrix, AntisymmetricMatrix)
from constants import (PID, default_inputs, default_masses, input_names, VCKM, 
                      IMVCKM, VCKMele, particle_names, input_to_PID, 
                      PID_to_input)
from __init__ import RosettaError
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
    "sminput" and "mass" respectively (the Z ahd Higgs masses are also stored 
    in self.inputs). Specifically, the block names specified as keys in 
    self.blocks should match those in the SLHA card as well as the names and 
    indices of the elements (taken to be the order in which the appear in the 
    list associated the the block in self.blocks). The blocks "mass" and 
    "sminputs" should minimally contain the entries specified in required_masses 
    and required_inputs respectively. A number of checks related to these 
    definitions are performed by check_param_data(), check_mass() and 
    check_sminputs() on the data read in from the parameter card. An example is 
    shown below where the required data members are defined outside the class 
    construtor so as to instrinsically belong to all instances of the class.
    
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
    basis which sets the SLHA.Card object, self.newcard. Rosetta differentiaties 
    from general utility functions and translation functions using the 
    "translation" decorator. Any function designed to take you to another 
    existing basis implemenation should be decorated with the "translation" 
    decorator with an argument corresponding to the name of the target basis. 
    This name should match the unique name of an existing basis implementation 
    also contained in the Rosetta root directory. Below would be example of a 
    translation function from our example to the Warsaw Basis.
    
        >> @Basis.translation('warsaw') # identifies as translator to warsaw 
        >> def mytranslation(self, instance):
        >>     instance['cWW'] = 10.
        >>     instance['cpHl11'] = self['myparam']
        >>     return instance
    
    Any python module saved in the root directory of the Rosetta package is 
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
    required_inputs, required_masses = set(),set()
    
    def __init__(self, param_card=None, output_basis='bsmc', flavor = 'general',
                 ehdecay=False, silent=False, translate=True,
                 dependent=True, modify_inputs=True):
        '''
        Keyword arguments:
            param_card   - SLHA formatted card conforming to the definitions in 
                           self.blocks, self.required_masses and 
                           self.required_inputs.
            output_basis - target basis to which to translate coefficients.
            ehdecay      - whether to run the eHDECAY interface to calculate the 
                           Higgs width and BRs.
            silent       - suppress all warnings and take default answers for 
                           all questions.
            translate    - Whether to call the translate method when reading in 
                           from a parameter card.
            flavor      - flavor structure of matrices: 'diagonal, 'universal' 
                           , 'MFV' or 'general'.
            dependent    - when param_card is None, whether or not to include 
                           dependent parameters in the SLHA.Card attribute of 
                           the basis instance.
        '''
        
        self.flavor = flavor

        self.param_card = param_card
        self.output_basis = output_basis
        self.newname = 'Basis'
        if not hasattr(self, 'dependent'): self.dependent=[]

        self.set_dependents()

        self.set_fblocks(self.flavor)

        # read param card (sets self.inputs, self.mass, self.name, self.card)
        if param_card is not None: 
            assert os.path.exists(self.param_card), \
                   '{} does not exist!'.format(self.param_card)
            
            self.read_param_card() 

            # various input checks
            self.check_sminputs(self.required_inputs) 
            self.check_masses(self.required_masses, silent=silent) 
            self.check_param_data(silent=silent)
            self.check_flavored_data(silent=silent)
            # generalises potentially reduced flavor structure
            
            self.set_flavor(self.flavor, 'general')
            # add dependent coefficients/blocks to self.card
            self.init_dependent()
            # generate internal OrderedDict() object for __len__, __iter__, 
            # items() and iteritems() method
            self._gen_thedict()
            
            # user defined function
            self.calculate_dependent()

            if translate:
                # translate to new basis (User defined) 
                # return an instance of a class derived from Basis
                self.newbasis = self.translate() 
            else:
                # do nothing
                self.newbasis = self
            
            if modify_inputs:
                self.newbasis.modify_inputs()
                self.newbasis.check_modified_inputs()
            
            # delete imaginary parts of diagonal elements in hermitian matrices
            if self.output_basis != 'bsmc':
                self.newbasis.reduce_hermitian_matrices()
                
            # set new SLHA card
            self.newcard = self.newbasis.card 
            self.newname = self.newbasis.name
            
            preamble = ('###################################\n'
                      + '## DECAY INFORMATION\n'
                      + '###################################')
            for decay in self.card.decays.values():
                decay.preamble = preamble
                break
                        
            try: 
                if ehdecay: 
                    self.run_eHDECAY()
            except TranslationError as e:
                print e
                print 'Translation to SILH Basis required, skipping eHDECAY.'
            
        # if param_card option not given, instantiate with class name 
        # and all coeffs set to 0 (used for creating an empty basis 
        # instance for use in translate() method)
        else: 
            self.card = self.default_card(dependent=dependent)            
            self._gen_thedict()
            
    # overloaded container (dict) methods for indexing etc.
    
    def __getitem__(self, key):
        return self.card.__getitem__(key)

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
        thedict = SLHA.CaseInsensitiveOrderedDict()
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
        
        vckm = self.default_ckm()
        thecard.add_block(vckm)
        thecard.ckm = vckm
        thecard.set_complex()
        self.fix_matrices(card = thecard)
        return thecard
    
    def set_flavor(self, _from, to):
        if _from == to: return

        self.set_fblocks(to) # reset fblocks according to flavor option
        newcard = self.default_card(dependent=False)
        
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
                    print  'Warning in {}.set_flavour():'.format(self.__class__)
                    print ('    Reduction in flavour structure ' +
                           'from "{}" to "{}" '.format(_from, to) +
                           'encountered some unexpected non-zero elements ' +
                           '(> 1e-6) which were not deleted.')
                    print 'Not deleted: {}\n'.format(', '.join(no_del_names))
                
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

    def read_param_card(self):
        '''
        Call SLHA.read() and set up the Basis instance accordingly.
        '''
        self.card = SLHA.read(self.param_card, set_cplx=False)
        
        try:
            card_name = self.card.blocks.get('basis',[''])[1]
            if self.name.lower()!=card_name.lower():
                err = ('Rosetta was expecting to read an instance of ' 
                     + '{}, named "{}", '.format(self.__class__.name, self.name)
                     + 'but read the name '
                     + '"{}" in block "basis" of {}.'.format(card_name,
                                                             self.param_card))
                raise RosettaError(err)
        except KeyError:
            raise RosettaError('Formatting error for block basis. '
                               'Check input card, {}.'.format(self.param_card))
        
        if not self.blocks:
            # if self.blocks not defined, default behaviour is to automatically
            # create self.blocks structure for independent parameters
            for par in self.independent:
                block = self.card._parent_block(par)
                name = block.name
                if name not in self.blocks:
                    self.blocks[name] = [par]
                else:
                    self.blocks[name].append(par)
        
        to_add = []
        for bname, blk in self.card.matrices.iteritems():
            is_cplx = self.flavored.get(bname,{}).get('domain','')=='complex'
            if is_cplx:
                if bname.lower().startswith('im'):
                    other_part = bname[2:]
                else: 
                    other_part ='IM' + bname
                if other_part not in self.card.matrices:
                    to_add.append((bname,other_part))

        for part, other_part in to_add:
            blk = self.card.matrices[part]
            for k, v in blk.iteritems():
                if other_part.lower().startswith('im'):
                    imname = 'R' + blk.get_name(k)[1:]
                else:
                    imname = 'I' + blk.get_name(k)[1:]
                self.card.add_entry(other_part, k, 0., name=imname)
                
        self.inputs = self.card.blocks.get('sminputs', None)
        self.mass = self.card.blocks.get('mass', None)
        
        self.card.set_complex()
        self.fix_matrices()
            
    def check_sminputs(self, required_inputs, message='Rosetta'):
        '''
        Check consistency of sminputs w.r.t required_inputs. Any inconsistencies
        found will raise a warning and a question of whether the user wishes 
        to continue assuming default values for offending inputs. Higgs and Z 
        masses will also be stored as SM inputs.
        '''
        # Check all required inputs have been set
        input_list = ', '.join( ['{} ({})'.format(x,input_names[x]) 
                                for x in required_inputs] )
        if required_inputs:
            if self.inputs is None:
                repr_default = ['{}={: .5e}'.format(input_names[k], 
                                 default_inputs.get(k,0.)) for 
                                 k in required_inputs]
                
                print 'Warning: Block "sminputs" not found. '\
                      'Assume default values for required inputs?'
                print 'Required inputs: {}'.format(input_list)
                carry_on = Y_or_N('Continue with default values '\
                                 'for unspecified inputs? '\
                                 '({})'.format(', '.join(repr_default)))
                if carry_on:
                    theblock = SLHA.NamedBlock(name='sminputs')
                    for m in required_inputs: 
                        theblock.new_entry(m, default_inputs[m], 
                                           name=input_names[m])
                    self.card.add_block(theblock)
                    self.inputs = theblock
                else:
                    print 'Exit'
                    sys.exit()
            else:
                if self.mass:
                    for k,v in [(i,j) for i,j in self.mass.iteritems() 
                                if i in (23,25)]:
                        i = PID_to_input[k]
                        if i in self.inputs:
                            v2 = float(self.inputs[i])
                            if v!=v2:
                                print ('Warning: {} '.format(input_names[i])
                                + 'specified in block sminput[{}] '.format(i)
                                + '({:.5e}) not consistent with '.format(v2)
                                + 'value specified in block mass '
                                + '[{}] ({:.5e}).\n'.format(k,float(v))
                                + 'Rosetta will keep value from sminputs.')
                                self.mass[k]=v2
                        else:
                            mstring = 'M{}'.format(particle_names[k])
                            self.inputs.new_entry(i, v,name=mstring)

                required = set(required_inputs)
                inputs = set(self.inputs.keys())
                missing_inputs = required.difference(inputs)
                missing_values = {k:default_inputs.get(k,0.) 
                                  for k in missing_inputs}
                repr_default = ['{}={: .5e}'.format(input_names[k], 
                                 default_inputs.get(k,0.)) for 
                                 k in missing_inputs]
                if missing_inputs: # Deal with unassigned SM inputs
                    print 'Warning: Not all required SM inputs are '\
                          'defined for {}.'.format(message)
                    print 'Required inputs: {}'.format(input_list)
                    missing_list = ', '.join([str(x) for x in missing_inputs])
                    print 'Missing inputs: {}'.format(missing_list)
                    carry_on = Y_or_N(
                        'Continue with default values '\
                        'for unspecified inputs? '\
                        '({})'.format(', '.join(repr_default))
                                     )
                    if carry_on: # sets unspecified inputs
                        for m in missing_inputs: 
                            self.inputs.new_entry(m, default_inputs[m], 
                                                  name=input_names[m])
                    else:
                        print 'Exit'
                        sys.exit()
            # ensure presence of CKM matrix
            vckm = self.card.matrices.get('vckm', None)
            if vckm is None:
                print 'Block "VCKM" not found, will use default values.\n'
                vckm = self.default_ckm()
                self.card.add_block(vckm)
            self.ckm = vckm
                

    def default_ckm(self):
        preamble = ('\n###################################\n'
                + '## CKM INFORMATION\n'
                + '###################################\n')
                
        ckm = SLHA.NamedMatrix(name='VCKM', preamble=preamble)
        for i,j in product((1,2,3),(1,2,3)):
            cname = 'RV{}{}x{}'.format(VCKMele[(i,j)], i, j)
            ckm.new_entry((i,j), VCKM[i][j], name=cname)
            
        imckm = SLHA.NamedMatrix(name='IMVCKM')
        for i,j in product((1,2,3),(1,2,3)):
            cname = 'IV{}{}x{}'.format(VCKMele[(i,j)], i, j)
            imckm.new_entry((i,j), IMVCKM[i][j], name=cname) 
            
        # vckm = SLHA.CNamedMatrix(ckm, imckm)
        vckm = CTwoDMatrix(SLHA.CNamedMatrix(ckm, imckm))
        return vckm

        
    def check_masses(self, required_masses, silent=False, message='Rosetta'):
        '''
        Check consistency of particle masses w.r.t required_masses. Any 
        inconsistencies found will raise a warning and a question of whether 
        the user wishes to continue assuming default values for masses. Higgs 
        and Z masses are also stored as SM inputs.
        '''
        PID_list = ', '.join([str(x) for x in required_masses])
        if required_masses:
            if self.mass is None:

                repr_default = ['{}={:.5e}'.format(particle_names[k], 
                                 default_masses.get(k,0.)) for 
                                 k in required_masses]
                mass_list = ', '.join( ['{} (M{})'.format(x,particle_names[x]) 
                                        for x in required_masses] )
                if not silent:
                    print 'Warning: Block "mass" not found. '\
                          'Assume default values for required masses?'
                    print 'Required PIDs: {}'.format(mass_list)

                    carry_on = Y_or_N(
                        ' values '\
                        'for unspecified masses? '\
                        '({})'.format(', '.join(repr_default))
                                     )
                else:
                    carry_on=True
                    
                if carry_on:
                    theblock = SLHA.NamedBlock(name='mass')
                    for m in required_masses: 
                        theblock.new_entry(m, default_masses[m], 
                                           name='M%s' % particle_names[m])
                    self.card.add_block(theblock)
                    self.mass = theblock
                else:
                    print 'Exit'
                    sys.exit()
            else:
                if self.inputs:
                    for k,v in [(i,j) for i,j in self.inputs.iteritems() 
                                if i in (4,25)]:
                        i = input_to_PID[k]
                        if i in self.mass:
                            v2 = float(self.mass[i])
                            if v!=v2:
                                if not silent:
                                    print ('Warning: M{} '.format(particle_names[i])
                                + 'specified in block sminput[{}] '.format(k)
                                + '({:.5e}) not consistent with '.format(v)
                                + 'value specified in block mass '
                                + '[{}] ({:.5e}).\n'.format(i,float(v2))
                                + 'Rosetta will keep value from sminputs.')
                                self.mass[i]=v
                        else:
                            mstring = 'M{}'.format(particle_names[i])
                            self.mass.new_entry(i, v, name=mstring)
                            
                masses = set(self.mass.keys())
                required = set(required_masses)
                missing_masses = required.difference(masses)
                missing_values = {k:default_masses.get(k,0.) 
                                  for k in missing_masses}
                if missing_masses: # Deal with unassigned fermion masses
                    repr_default = ['M{}={: .5e} GeV'.format(particle_names[k],v) 
                                    for k,v in missing_values.items()]
                    mass_list = ', '.join(['{} (M{})'.format(x,particle_names[x]) 
                                            for x in missing_masses])
                    if not silent:
                        print 'Warning: Not all required masses are '\
                              'defined for {}.'.format(message)
                        print 'Required PIDs: {}'.format(PID_list)
                        print 'Missing PIDs: {}'.format(mass_list)
                        carry_on = Y_or_N(
                            'Continue assuming default values '\
                            'for unspecified masses? '\
                            '({})'.format(', '.join(repr_default))
                            )
                    else:
                        carry_on = True
                    if carry_on: 
                        for m in missing_masses: 
                            self.mass.new_entry(m, default_masses[m], 
                                                name='M%s' % particle_names[m])
                            
                    else:
                        print 'Exit'
                        sys.exit()

    def check_param_data(self, silent=False, 
                               do_unknown=True, 
                               do_consistency=True,
                               do_dependent=True,
                               do_independent=True):
        '''
        Cross check of parameter data read in the SLHA formatted card.
        Compares lists of coefficients declared in self.independent and 
        self.dependent to those read in from self.param_card.
           1)  Deals with unrecognised names by removing them from the block
           2)  Check consistency of name & numbering for the EFT basis blocks 
               declared in self.blocks. Renames coefficients according to their 
               position in each block's list.
           2)  Prints a warning if coefficients declared as dependent are 
               assigned values in param_card. The user is asked whether or not 
               they wish to continue knowing these values could be overwritten 
               by calculate_dependent().
           3)  Checks if all coefficients declared as independent are 
               assigned values. If not, the user is given the option to 
               continue with them set to 0.
        '''
        for bname, defined in self.blocks.iteritems():
            # collect block info
            inputblock = self.card.blocks.get(bname, None)
            if inputblock is None:
                theblock = SLHA.NamedBlock(name=bname)
                self.card.add_block(theblock)
                inputblock = theblock
            
            input_eles = set(inputblock.keys())
            
            defined_block = {i+1:v for i,v in enumerate(defined)}
            defined_eles = set(defined_block.keys())
            
            independent = {i:v for i,v in defined_block.iteritems() 
                           if v in self.independent}
            dependent = {i:v for i,v in defined_block.iteritems() 
                           if v in self.dependent}
                           
            # check for unrecognised coefficient numbers              
            unknown = input_eles.difference(defined_eles)
            if unknown and do_unknown:
                unknown_names = {i:inputblock.get_name(i,'none') 
                                 for i in unknown}
                if not silent:
                    print 'Warning: you have declared coefficients '\
                          'undefined in {}, block:{}.'.format(self.__class__,
                                                              bname)
                    print 'The following will be ignored - '\
                          '{}'.format(', '.join(['{}:"{}"'.format(k,v) for k,v 
                                                 in unknown_names.iteritems()]))
                                             
                for x in unknown: del inputblock[x]
                                             
            # check that name, number pairs match
            mismatched = []
            unnamed = []
            
            for index, name in independent.items():
                input_name = inputblock.get_name(index, None)
                if input_name is None: 
                    inputblock._names[index] = name
                elif input_name.lower()!=name.lower():
                    mismatched.append((index,name,input_name))
                    
            if mismatched and do_consistency:
                if not silent:
                    print 'Warning: Mismatch of coefficient names '\
                          'in {}, block "{}".'.format(self.__class__,
                                                              bname)
                for index, name, input_name in mismatched:
                    if not silent:
                        print ('    Coefficient ' +
                               '{}, named {} '.format(index, input_name) +
                               'will be renamed to {}'.format(name))
                    inputblock._names[index] = name
                    inputblock._numbers[name] = index
            
            # check if coefficients defined as dependent 
            # have been assigned values
            defined_dependent = set(dependent).intersection(input_eles)
            if defined_dependent and do_dependent: 
                if not silent:
                    print 'Warning: you have assigned values to some '\
                          'coefficients defined as dependent '\
                          'in {}, block "{}".'.format(self.__class__, bname)
                    print 'Coefficients: {}'.format(', '.join(
                                                    ['{}:"{}"'.format(k,v) 
                                                    for k,v in dependent.items() 
                                                    if k in defined_dependent]
                                                    ))
                    print 'These may be overwritten by an implementation of '\
                          '{}.calculate_dependent()'.format(self.__class__.__name__)
                    carry_on = Y_or_N('Continue?')
                else:
                    carry_on = True
                if carry_on:
                    pass
                else:
                    print 'Exit'
                    sys.exit()
            
            # check if all independent coefficients are assigned values
            missing = set(independent).difference(input_eles).difference(set(dependent))
            if missing and do_independent: # Deal with unassigned independent coefficients
                if not silent:
                    print 'Warning: some coefficients defined as independent '\
                          'in {}, block "{}", have not been assigned values'\
                          '.'.format(self.__class__,bname)
                    print 'Undefined: {}'.format(', '.join(
                                                    ['{}:"{}"'.format(k,v) 
                                                     for k,v in independent.items() 
                                                     if k in missing]))
                    carry_on = Y_or_N('Continue assuming unspecified '\
                                      'coefficients are Zero?')
                else:
                    carry_on = True
                if carry_on:
                    for m in missing: 
                        inputblock.new_entry(m, 0., independent[m])
                else:
                    print 'Exit'
                    sys.exit()
                    
            if len(inputblock)==0:
                del self.card.blocks[inputblock.name]
    
    
    def check_flavored_data(self, silent=False):
        '''
        Cross check of flavored part of parameter data read in the SLHA 
        formatted card. Compares lists of coefficients declared in 
        self.independent and self.dependent to those read in from self.param_card.
           1)  Deals with unrecognised names by removing them from the block
           2)  Check consistency of name & numbering for the EFT basis blocks 
               declared in self.blocks. Renames coefficients according to their 
               position in each block's list.
           2)  Prints a warning if coefficients declared as dependent are 
               assigned values in param_card. The user is asked whether or not 
               they wish to continue knowing these values could be overwritten 
               by calculate_dependent().
           3)  Checks if all coefficients declared as independent are 
               assigned values. If not, the user is given the option to 
               continue with them set to 0.
        '''

        for bname, defined in self.fblocks.iteritems():
            # collect block info
            inputblock = self.card.matrices.get(bname, None)
            if inputblock is None:
                theblock = SLHA.NamedMatrix(name=bname)
                self.card.add_block(theblock)
                inputblock = theblock
                
            input_eles = set(inputblock.keys())
            
            defined_block = {(int(v[-3]),int(v[-1])):v for 
                                i,v in enumerate(defined)}
            defined_eles = set(defined_block.keys())
            
            independent = {k:v for k,v in defined_block.iteritems() 
                           if v in self.independent or bname in self.independent}
            dependent = {k:v for k,v in defined_block.iteritems() 
                           if v in self.dependent}
            # check for unrecognised coefficient numbers              
            unknown = input_eles.difference(defined_eles)
            if unknown:
                unknown_names = {i:inputblock.get_name(i,'none') 
                                 for i in unknown}
                if not silent:
                    print 'Warning: you have declared coefficients '\
                          'undefined in {}, matrix:{}.'.format(self.__class__,
                                                              bname)
                    print 'The following will be ignored - '\
                          '{}'.format(', '.join(['{}:"{}"'.format(k,v) for k,v 
                                                 in unknown_names.iteritems()]))
                                             
                for x in unknown: del inputblock[x]
                                             
            # check that name, number pairs match
            mismatched = []
            unnamed = []
            
            for index, name in independent.items():
                input_name = inputblock.get_name(index, None)
                if input_name is None: 
                    inputblock._names[index] = name
                elif input_name.lower()!=name.lower():
                    mismatched.append((index,name,input_name))
                    
            if mismatched:
                if not silent:
                    print 'Warning: Mismatch of coefficient names '\
                          'in {}, matrix "{}".'.format(self.__class__,
                                                              bname)
                for index, name, input_name in mismatched:
                    if not silent:
                        print ('    Coefficient ' +
                               '{}, named {} '.format(index, input_name) +
                               'will be renamed to {}'.format(name))
                    inputblock._names[index] = name
                    inputblock._numbers[name] = index

            # check if coefficients defined as dependent 
            # have been assigned values
            defined_dependent = set(dependent).intersection(input_eles)
            if defined_dependent: 
                if not silent:
                    print 'Warning: you have assigned values to some '\
                          'coefficients defined as dependent '\
                          'in {}, matrix "{}".'.format(self.__class__, bname)
                    print 'Coefficients: {}'.format(', '.join(
                                                   ['{}:"{}"'.format(k,v) 
                                                    for k,v in dependent.items() 
                                                    if k in defined_dependent]
                                                    ))
                    print 'These may be overwritten by an implementation of '\
                          '{}.calculate_dependent()'.format(self.__class__.__name__)
                    carry_on = Y_or_N('Continue?')
                else:
                    carry_on = True
                if carry_on:
                    pass
                else:
                    print 'Exit'
                    sys.exit()

            # check if all independent coefficients are assigned values
            missing = set(independent).difference(input_eles).difference(set(dependent))
            if missing: # Deal with unassigned independent coefficients
                if not silent:
                    print 'Warning: some coefficients defined as independent '\
                          'in {}, matrix "{}", have not been assigned values'\
                          '.'.format(self.__class__,bname)
                    print 'Undefined: {}'.format(', '.join(
                                                ['{}:"{}"'.format(k,v) 
                                                 for k,v in independent.items() 
                                                 if k in missing]))
                    carry_on = Y_or_N('Continue assuming unspecified '\
                                      'coefficients are Zero?')
                else:
                    carry_on = True
                if carry_on:
                    for m in missing: 
                        inputblock.new_entry(m, 0., independent[m])
                else:
                    print 'Exit'
                    sys.exit()
            
            if len(inputblock)==0:
                del self.card.matrices[inputblock.name]

    def write_param_card(self, filename, overwrite=False, 
                                         include_dependent=True):
        '''Write contents of self.newcard to filename'''

        preamble = ('\n###################################\n'
                    + '## INFORMATION FOR {} BASIS\n'.format(self.newname.upper())
                    + '###################################\n')
        if 'basis' in self.newcard.blocks:
            self.newcard.blocks['basis'].preamble = preamble
            self.newcard.blocks['basis'][1] = self.newname

        dec_preamble = ('\n###################################\n'
                    + '## DECAY INFORMATION\n'
                    + '###################################\n')
        for decay in self.newcard.decays.values():
            decay.preamble = dec_preamble
            break

        mass_preamble = ('\n###################################\n'
                       + '## INFORMATION FOR MASS\n'
                       + '###################################\n')
        if 'mass' in self.newcard.blocks:
            self.newcard.blocks['mass'].preamble = mass_preamble

        sm_preamble = ('\n###################################\n'
                     + '## INFORMATION FOR SMINPUTS\n'
                     + '###################################\n')
        if 'sminputs' in self.newcard.blocks:
            self.newcard.blocks['sminputs'].preamble = sm_preamble

        ckm_preamble = ('\n###################################\n'
                      + '## CKM INFORMATION\n'
                      + '###################################\n')
        if 'vckm' in self.newcard.matrices:
            self.newcard.matrices['vckm'].preamble = ckm_preamble

        card_preamble = ('################################################'
                    +'######################\n'
                    +'############# COEFFICIENTS TRANSLATED BY ROSETTA'
                    +' MODULE  #############\n'
                    +'########### PARAM CARD GENERATED {}  ##########'\
                     '#\n'.format(datetime.datetime.now().ctime().upper())
                    +'################################################'
                    +'######################\n\n')

        if os.path.exists(filename) and not overwrite:
            print '{} already exists.'.format(filename)
            carry_on = Y_or_N('Overwrite?', default='no')
        else:
            carry_on=True
        
        if carry_on:
            special_blocks = ['loop','mass','sminputs','yukawa','vckm','basis']
            coefforder = sortblocks(self.newcard, ignore = special_blocks)
            self.newcard.write(filename, blockorder=special_blocks + coefforder,
                                         preamble=card_preamble)
            print 'Wrote "{}".'.format(filename)
            return True
        else:
            return False
    
    def write_template_card(self, filename, value=0.):
        from machinery import bases
        try: 
            val = float(value)
            rand = False
        except ValueError:
            if value.lower() == 'random':
                import random
                rand = True
            else:
                print ('In write_template_card: "value" keyword argument '
                       'must either be a number or the string, "random".')  
                sys.exit()          
        newinstance = bases[self.name](flavor=self.flavor, dependent=False)
        for k in newinstance.keys():            
                try:
                    if rand:
                        newinstance[k] = complex(random.uniform(-1.,1.),
                                                 random.uniform(-1.,1.))
                    else:
                        newinstance[k] = complex(val, val)
                except TypeError as e:
                    if rand:
                        newinstance[k] = random.uniform(-1.,1.)
                    else:
                        newinstance[k] = val
        
        newinstance.reduce_hermitian_matrices()
        
        SLHA_card = newinstance.card
        
        mass_preamble = ('\n###################################\n'
                       + '## INFORMATION FOR MASS\n'
                       + '###################################\n')

        massblock = SLHA.NamedBlock(name='mass', preamble=mass_preamble)
        for m in self.required_masses: 
            massblock.new_entry(m, default_masses[m], 
                               name='M%s' % particle_names[m])
        SLHA_card.add_block(massblock)
                    
        sm_preamble = ('\n###################################\n'
                     + '## INFORMATION FOR SMINPUTS\n'
                     + '###################################\n')
                    
        inputblock = SLHA.NamedBlock(name='sminputs', preamble=sm_preamble)
        for m in self.required_inputs: 
            inputblock.new_entry(m, default_inputs[m], 
                                 name='%s' % input_names[m])
        SLHA_card.add_block(inputblock)
        
        title = ' ROSETTA: {} BASIS INPUT CARD '.format(self.name.upper())
        time = ' GENERATED {} '.format(datetime.datetime.now().ctime().upper())
        nleft = int((80-len(title))/2.)
        nright = 80-len(title)-nleft
        preamble = ( '#'*80 + '\n'
                    +'#'*nleft + title + '#'*nright +'\n'
                    +'#'*22 + time + '#'*22 + '\n'
                    +'#'*80 +'\n')
        
        special_blocks = ['mass','sminputs','vckm','basis']
                
        theorder = sortblocks(SLHA_card, ignore = special_blocks)
        SLHA_card.write(filename, blockorder=special_blocks + theorder, 
                        preamble = preamble, postamble = ('\n'+'#'*80)*2)
        print ('wrote {}'.format(filename))
        
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
        
    def check_calculated_data(self):
        '''
        Compares self.dependent and self.card to see if all dependent 
        coefficients have been calculated. If not, asks whether the user wants 
        to continue assuming they are zero.
        '''
        missing_dependents = set(self.dependent).difference(self.par_dict.keys())
        if missing_dependents and self.dependent:
            print 'Warning: Set of dependent coefficients calculated '\
                  'by {0}.calculate_dependent() does not match those '\
                  'specified in {0}.'.format(self.__class__)
            print 'Undefined: {}'.format(', '.join(missing_dependents))
            carry_on = Y_or_N('Continue assuming coefficients are Zero?')
            if carry_on:
                for m in missing_dependents: self.par_dict[m]=0.
            else:
                print 'Exit'
                sys.exit()   
            print 'Calculated coefficients match those defined '\
                  'in {}.dependent.'.format(self.__class__.__name__) 
    
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
    
    def translate(self, target=None, verbose=True): 
        '''
        Translation function. Makes use of existing translations and attempts to 
        find a path between the basis instance and the target basis. If a path 
        exists, successive translations are performed to get from A to B.
        Keyword arguments:
            target - target basis, must be defined in implemented.bases.
            verbose - print some status messages.
        '''
        from machinery import get_path, relationships, bases
                
        # default target
        target = target if target is not None else self.output_basis
        
        if target == self.name: return self
        
        # find the chain of bases and translation functions 
        # to get from self.name to target
        chain = get_path(self.name.lower(), target, relationships)

        if not chain: # no translation possible
            inputbasis = self.__class__.__name__
            outputbasis = bases[target].__name__
            err = ('Sorry, Rosetta cannot translate from ' +
                  '{} to {}'.format(inputbasis, outputbasis))
            raise TranslationError(err)
        
        names = [self.name.lower()]+[x[0] for x in chain]
        if verbose:
            print 'Rosetta will be performing the translation:'
            print '    '+' -> '.join([bases[x].__name__ 
                                                    for x in names])
        
        # perform succesive translations, checking for  
        # required SM inputs/masses along the way

        current = self
        
        required_inputs = current.required_inputs
        required_masses = current.required_masses

        
        for target, translate_function in chain:
            # new target basis instance
            # instance = bases[target](dependent=False)
            instance = bases[target](dependent=True)
            # ensure required inputs are present
            message = 'translation ({} -> {})'.format(current.name, 
                                                      instance.name)
                                                      
            current.check_sminputs(required_inputs, message=message)
            current.check_masses(required_masses, message=message)
            
            new = translate_function(current, instance)
            
            # update new basis instance with non-EFT blocks, decays
            all_coeffs = (current.blocks.keys() + current.flavored.keys())
            new.get_other_blocks(current.card, all_coeffs)
            
            message = ' {}'.format(instance.__class__)
            # checks & calculate dependent
            print 'Checking for required SM inputs'
            new.check_sminputs(new.required_inputs, message=message)
            print 'Checking for required masses'
            new.check_masses(new.required_masses, message=message)
            print 'Checking EFT coefficients'
            # NOTE cannot check for presence of dependent coefficients as there 
            # is no difference between the coefficient existing and it having 
            # been assigned a value.
            new.check_param_data(do_dependent=False)
            print 'Calling {}.calculate_dependent()'.format(new.__class__)
            new.init_dependent()
            new.calculate_dependent()
            
            # prepare for next step
            required_inputs = set(instance.required_inputs)
            required_masses = set(instance.required_masses)
            current = new
        
        if current.name =='bsmc':
            current.expand_matrices()
        else:
            current.flavor = self.flavor
            # reduce flavor structure back to user set option
            current.set_flavor('general', self.flavor)            
            
        if verbose:
            print 'Translation successful.\n'
        
        
        return current

    def get_other_blocks(self, card, ignore):
        ignore = [x.lower() for x in ignore] 
        
        other_blocks, other_matrices = {}, {}
        for k, v in card.blocks.iteritems():
            if k.lower() != 'basis' and k.lower() not in ignore:
                other_blocks[k]=v
        for k, v in card.matrices.iteritems():
            if k.lower() != 'basis' and k.lower() not in ignore:
                other_matrices[k]=v

        for block in other_blocks:
            theblock = card.blocks[block]
            self.card.add_block(theblock)
        
        for matrix in other_matrices:
            theblock = card.matrices[matrix]
            self.card.add_block(theblock)
        
        for decay in card.decays.values():
            self.card.add_decay(decay, preamble = decay.preamble)
        
        if card.has_block('mass'):
            self.mass=self.card.blocks['mass']
        if card.has_block('sminputs'):
            self.inputs=self.card.blocks['sminputs']
            
        self.ckm = card.matrices['vckm']
    
    def expand_matrices(self):
        '''
        Special function to populate redundant elements of matrix blocks when 
        translating to the bsmc Lagrangian so that values for all 9 entries are 
        explicitly stored before writing out the parameter card. This is to 
        stay in accordance with the SLHA format.
        The function directly modifies the _data and  _names attributes of the 
        matrices since matrices with special properties i.e. Hermitian, 
        Symmetric etc. do not grant direct access to the redundant keys such as 
        the lower triangle of a Hermitian matrix.
        '''
        all_keys = [(1,1), (1,2), (1,3),
                    (2,1), (2,2), (2,3),
                    (3,1), (3,2), (3,3)]
                    
        for matrix in self.card.matrices.values():
            # list of missing elements in _data member of matrix instance
            missing_keys = [k for k in all_keys if k not in matrix._data]
            
            if missing_keys:
                # randomly select parameter name since they all should have 
                # the same structure: (R|I)NAMEixj
                elename = matrix._names.values()[0]
                cname = elename[1:-3] # matrix name
                pref = elename[0] 
                for k in missing_keys:
                    tail = cname + '{}x{}'.format(*k)
                    matrix._data[k] = matrix[k]
                    matrix._names[k] = pref + tail
                    matrix._numbers[pref+tail] = k
                    try:
                        matrix._re._data[k] = matrix._re[k]
                        matrix._re._names[k] = 'R' + tail
                        matrix._re._numbers['R'+tail] = k
                    except AttributeError:
                        pass
                    try:
                        matrix._im._data[k] = matrix._im[k]
                        matrix._im._names[k] = 'I' + tail
                        matrix._im._numbers['I'+tail] = k
                    except AttributeError:
                        pass
    
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
        
        
    def run_eHDECAY(self):
        '''
        Interface Rosetta with eHDECAY to calculate Higgs widths and branching 
        fractions. 
        '''
        import eHDECAY
        
        BRs = eHDECAY.run(self, electroweak=True)
        
        sum_BRs = sum([v for k,v in BRs.items() if k is not 'WTOT'])
        
        # sometimes eHDECAY gives a sum of BRs slightly greater than 1. 
        # for now a hacky global rescaling is implemented to deal with this.
        if sum_BRs > 1: 
            if sum_BRs - 1. > 1e-2: # complain if its too wrong
                raise RuntimeError('Sum of branching fractions > 1 by more than 1%')
            else:
                for channel, BR in BRs.iteritems():
                    if channel!='WTOT':
                        BRs[channel] = BR/sum_BRs
                
        totalwidth = BRs.pop('WTOT')
        
        if totalwidth < 0.:
            print ('eHDECAY: Negative total Higgs width. Check your EFT inputs.')
            return
        
        hdecays = {}
        
        # sometimes eHDECAY gives negative branching fractions. 
        for channel, BR in BRs.iteritems(): 
            # n1, n2 = particle_names[channel[0]], particle_names[channel[1]]
            # comment = 'H -> {}{}'.format(n1,n2)
            # if BR < 0.:
            #     print ('eHDECAY: Negative branching fraction encountered ' +
            #           'for {}. Rosetta will ignore it.'.format(comment))
            #     totalwidth -= BR # adjust total width
            #     continue
            # elif BR == 0.:
            #     continue
            if BR==0.:
                continue
            else:
                hdecays[channel] = BR

        # credit
        preamble = ('# Higgs widths and branching fractions '
                    'calculated by eHDECAY.\n# Please cite '
                    'arXiv:hep-ph/974448 & arXiv:1403.3381.')
        # new SLHA.Decay block
        decayblock = SLHA.Decay(25, totalwidth, data=hdecays, preamble=preamble)
        
        if 25 not in self.newcard.decays:
            self.newcard.add_decay(decayblock)
        else:
            self.newcard.decays[25] = decayblock
        print '#############################\n'
        
    def calculate_dependent(self):
        '''
        Default behavoiur of calculate_dependent(). Called if a subclass of 
        Basis has not implemented the function.
        '''
        print 'Nothing done for {}.calculate_'\
              'dependent()\n'.format(self.__class__.__name__)
              
    def modify_inputs(self):
        '''
        Default behavoiur of modify_inputs(). Called if a subclass of 
        Basis has not implemented the function.
        '''
        print 'Nothing done for {}.modify_'\
              'inputs()\n'.format(self.__class__.__name__)
    
    def check_modified_inputs(self, silent=False):
        for k,v in [(i,j) for i,j in self.inputs.iteritems() 
                    if i in (4,5,6,7,25)]:
            i = input_to_PID[k]
            if i in self.mass:
                v2 = float(self.mass[i])
                if v!=v2:
                    if not silent:
                        print ('Warning: M{} '.format(particle_names[i])
                    + 'in block sminputs[{}] '.format(k)
                    + '({:.5e}) not consistent with '.format(v)
                    + 'value specified in block mass'
                    + '[{}] ({:.5e}) '.format(i,float(v2))
                    + 'after modify_inputs().\n')
                    
                        keep_from_input = Y_or_N('Keep value from sminputs?')
                    else:
                        keep_from_input = True
                        
                    if keep_from_input:
                        if not silent:
                            print ('Modified M{} '.format(particle_names[i]) +
                                   'in block mass[{}]'.format(i))
                        self.mass[i]=v
                    else:
                        if not silent:
                            print ('Modified M{} '.format(particle_names[i]) +
                                   'in block sminputs[{}]'.format(k))
                        self.inputs[k]=v2
                        
            else:
                mstring = 'M{}'.format(particle_names[i])
                self.mass.new_entry(i, v,name=mstring)
        
        
################################################################################
class TranslationError(Exception):
    '''Fancy error name.'''
    pass
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
        print ('flavor_matrix function got and unrecognised combination of '
               '"kind" and "domain" keyword arguments')
        return [name]
    
    if (not cplx and domain!='complex'):
        return real
    else:
        return ['C'+c for c in real+cplx]

def sortblocks(card, ignore = []):
    normal_blocks = sorted([k for k in card.blocks.keys() 
                            if k.lower() not in ignore],
                            key=str.lower)
    flav_blocks = sorted([k for k in card.matrices.keys() 
                          if k.lower() not in ignore],
                          key=str.lower)
    flav_cplx = []
    for im in [fl for fl in flav_blocks if fl.lower().startswith('im')]:
        re = im[2:]
        flav_cplx.append(re)
        flav_cplx.append(im)
        flav_blocks.remove(re)
        flav_blocks.remove(im)

    return normal_blocks + flav_blocks + flav_cplx
    
################################################################################
if __name__=='__main__':
    pass