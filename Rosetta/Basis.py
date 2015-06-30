import StringIO
import re
import sys
import os
import math
import datetime
from collections import namedtuple,OrderedDict
from itertools import product, combinations
from itertools import combinations_with_replacement as combinations2
# from implemented import bases
import SLHA
from eHDECAY import eHDECAY
from query import query_yes_no as Y_or_N
from __init__ import (PID, default_inputs, default_masses, input_names, 
                      particle_names, input_to_PID, PID_to_input)
########################################################################
# Base Basis class
class Basis(object):
    '''
    Base class from which to derive other Higgs EFT basis classes. 
    
    Designed to be instantiated with an SLHA format parameter card using
    read_param_card(). The contents of the card are stored as an SLHA.Card 
    object and the special blocks "mass" and "sminputs" are stored as 
    SLHA.NamedBlock instances in self.mass and self.input respectively.
    
    self.name   - string storing the value of the 0th element of 
                     "block" basis.
    self.card   - SLHA.Card instance containing SLHA.NamedBlock and SLHA.Decay 
                 instances corresponding to those specified in the parameter card. 
                 Object names taken to be the first non-whitespace characters 
                 after a "#" character in a block or parameter definition are 
                 also stored.
    self.mass   - SLHA.NamedBlock instance for the "mass" block
    self.inputs - SLHA.NamedBlock instance for the "sminputs" block

    self.blocks, self.required_inputs and self.required_masses should be defined
     in accordance with block structure of the SLHA parameter card, Blocks 
    "sminput" and "mass" respectively (the Z ahd Higgs masses are also stored 
    in self.inputs). Specifically, the block names specified as keys in 
    self.blocks should match those in the SLHA card as well as the names and 
    indices of the elements (taken to be the order in which the appear in the 
    list associated the the block in self.blocks). The blocks "mass" and 
    "sminputs" should minimally contain the entries specified in required_masses
     and required_inputs respectively.A number of checks related to these 
    definitions are performed by check_param_data(), check_mass() and 
    check_inputs() on the data read in from the parameter card. An example is 
    shown below where the data members are defined outside the class construtor 
    so as to instrinsically belong to all instances of the class.
    
        >> independent = {'A','B','C','1','2','3'}
        >> dependent = {'D','4'}
        >> required_inputs = {1,2,4}   # a_{EW}^{-1}, Gf and MZ required
        >> required_masses = {24,25,6} # Z, Higgs and top masses required
        >> translate_to = {'mass'}     # translation to mass basis implemented
        >> blocks = {'letters':['A','B','C','D]
                  'number':['1','2','3','4']} # Expected block structure
    
    The lists self.independent and self.dependent classify the basis parameters 
    into those which should be read in and those which should be derived by 
    calling the calculate_dependent() method. The set of parameters read in 
    can then by translated into a different basis using the translate() method.
    
    Basis is designed to work similarly to a dictionary in that parameter values
    can be referenced by name (duplicate names are not handled properly so try 
    to avoid them). A value can be referenced in various ways, see below example
    where the parameter named 'myparam' is stored as entry 16 in the block 
    'myblock' written in 'mycard.dat':
        
        >> instance = MyBasis(param_card='mycard.dat')
        >> instance['myparam'] = 0.5 # set value according to name in SLHA card
        >> instance.card['myparam'] = 0.5 # equivalent to the above
        >> instance.card.blocks['myblock']['myparam'] = 0.5 # from block by name 
        >> instance.card.blocks['myblock'][16] = 0.5 # from block by entry 
        
    The user can define methods calculate_dependent() and translate() 
    to calculate any dependent parameters and then translate the 
    coefficients into a given basis (MassBasis by default) which sets
    the SHLA.Card object, self.newcard.
        
        from MassBasis import MassBasis
        def translate():
            instance = MassBasis() # Empty instance with all coeffs 0
            instance['myparam'] = 10.
            return instance
    
    write_param_card() writes the contents of the self.newcard into a new SLHA 
    formatted file. Any derived class can then be used by the command line 
    script "translate", provided Rosetta knows about it by including it in the 
    bases dict of implemented.py.
    '''
    
    independent, dependent=[], []
    required_inputs, required_masses = set(),set()
    translate_to = {'mass'}
    blocks = dict()
    
    def __init__(self, param_card=None, output_basis='mass', 
                 ehdecay=False, silent=False):
        '''
        Keyword arguments:
            param_card - SLHA formatted card conforming to the definitions in 
                         self.blocks, self.required_masses and 
                         self.required_inputs.
            output_basis - target basis to which to translate coefficients.
            ehdecay - whether to run the eHDECAY interface to calculate the 
                      Higgs width and BRs.
            silent  - suppress all checks and warnings
        '''
        
        # Check that target basis has a translation implemented
        if output_basis in self.translate_to:
            self.target_basis = output_basis
        else:
            error ='''
            {}.translate_to does not contain "{}". 
            Rosetta doesn't know of an implemented translation in {}
            '''.format(self.__class__.__name__, output_basis, 
                       self.__class__)
            raise ValueError(error)
            
        self.param_card = param_card
        self.newname = 'Basis'
        self.all_coeffs = self.independent + self.dependent
        # make blocks case insensitive
        self.blocks = {k.lower():v for k,v in self.blocks.iteritems()} 
        # read param card (sets self.inputs, self.mass, self.name, self.card)
        if param_card is not None: 
            assert os.path.exists(self.param_card), \
                   '{} does not exist!'.format(self.param_card)
            
            self.read_param_card() 
            # various input checks 
            self.check_sminputs() 
            self.check_masses(silent=silent) 
            self.check_param_data(silent=silent)
            # add dependent coefficients/blocks to self.card
            self.init_dependent()
            # (User defined function)
            self.calculate_dependent()
            # translate to new basis (User defined) 
            # return an instance of a class derived from Basis
            newbasis = self.translate() 
            # set new SLHA card
            self.newcard = newbasis.card 
            # fill remaining block & decay info from original
            other_blocks = {k:v for k,v in self.card.blocks.iteritems() 
                            if not (k in self.blocks or k.lower()=='basis')}
                            
            for block in other_blocks:
                theblock = self.card.blocks[block]
                self.newcard.add_block(theblock)
            
            preamble = ('###################################\n'
                        + '## DECAY INFORMATION\n'
                        + '###################################')
                        
            for i,decay in enumerate(self.card.decays.values()):
                if i !=0: preamble = ''
                self.newcard.add_decay(decay, preamble = preamble)
            
            # eHDECAY interface BROKEN!!
            try: 
                if ehdecay: self.BRs = eHDECAY(self) 
            except NotImplementedError:
                pass
            
        # if param_card option not given, instantiate with class name 
        # and all coeffs set to 0 (used for creating an empty basis 
        # instance for use in translate() method)
        else: 
            self.name=self.__class__.__name__
            self.card = SLHA.Card(name=self.name)
            self.card.add_entry('basis', 0, self.name, name = 'translated basis')
            
            preamble = ('###################################\n'
                        + '## INFORMATION FOR {}\n'.format(self.name.upper())
                        + '###################################')
            self.card.blocks['basis'].preamble=preamble
            
            for blk, flds in self.blocks.iteritems():
                theblock = SLHA.NamedBlock(name=blk)
                for i,fld in enumerate(flds):
                    theblock.new_entry(i,0., name=fld)
                self.card.add_block(theblock)

    # overloaded container (dict) methods for indexing etc.
    def __getitem__(self, key):
        return self.card.__getitem__(key)

    def __setitem__(self, key, value):
        return self.card.__setitem__(key, value)
    
    def __contains__(self,key):
        return self.card.__contains__(key)
        
    def __delitem__(self,key):
        return self.card.__delitem__(key)
        
    def read_param_card(self):
        '''
        Call SHLA.read() and set up the Basis instance accordingly.
        '''
        self.card = SLHA.read(self.param_card)
        self.inputs = self.card.blocks.get('sminputs',None)
        self.mass = self.card.blocks.get('mass',None)
        try:
            self.name = self.card.blocks.get('basis',[''])[0]
        except KeyError:
            print 'Formatting error for block basis. Check input card.'
        
    def check_sminputs(self):
        '''
        Check consistency of sminputs w.r.t required_inputs. Any inconsistencies
        found will raise a warning and a question of whether the user wishes 
        to continue assuming default values for offending inputs. Higgs and Z 
        masses will also be stored as SM inputs.
        '''
        # Check all required inputs have been set
        input_list = ', '.join( ['{} ({})'.format(x,input_names[x]) 
                                for x in self.required_inputs] )
        if self.required_inputs:
            if self.inputs is None:
                repr_default = ['{}={: .5e}'.format(input_names[k], 
                                 default_inputs.get(k,0.)) for 
                                 k in self.required_inputs]
                
                print 'Warning: Block "sminputs" not found. '\
                      'Assume default values for required inputs?'
                print 'Required inputs: {}'.format(input_list)
                carry_on = Y_or_N('Continue with default values '\
                                 'for unspecified inputs? '\
                                 '({})'.format(', '.join(repr_default)))
                if carry_on:
                    theblock = SLHA.NamedBlock(name='sminputs')
                    for m in self.required_inputs: 
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

                required = set(self.required_inputs)
                inputs = set(self.inputs.keys())
                missing_inputs = required.difference(inputs)
                missing_values = {k:default_inputs.get(k,0.) for k in missing_inputs}
                repr_default = ['{}={: .5e}'.format(input_names[k], 
                                 default_inputs.get(k,0.)) for 
                                 k in missing_inputs]
                if missing_inputs: # Deal with unassigned SM inputs
                    print 'Warning: Not all required SM inputs are '\
                          'defined.'
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
        print 'SM inputs are OK.'
        
    def check_masses(self, silent=False):
        '''
        Check consistency of particle masses w.r.t required_masses. Any 
        inconsistencies found will raise a warning and a question of whether 
        the user wishes to continue assuming default values for masses. Higgs 
        and Z masses are also stored as SM inputs.
        '''
        if silent: return
        mass_list = ', '.join( ['{} (M{})'.format(x,particle_names[x]) 
                                for x in self.required_masses] )
        PID_list = ', '.join([str(x) for x in self.required_masses])
        if self.required_masses:
            if self.mass is None:
                repr_default = ['{}={: .5e}'.format(particle_names[k], 
                                 default_inputs.get(k,0.)) for 
                                 k in self.required_masses]
                print 'Warning: Block "mass" not found. '\
                      'Assume default values for required masses?'
                print 'Required PIDs: {}'.format(mass_list)

                carry_on = Y_or_N(
                    'Continue with default values '\
                    'for unspecified masses? '\
                    '({})'.format(', '.join(repr_default))
                                 )
                if carry_on:
                    self.card.add_block('sminputs')
                    for m in self.required_inputs: 
                        self.inputs[m]=default_inputs[m]
                else:
                    print 'Exit'
                    sys.exit()
            else:
                masses = set(self.mass.keys())
                required = set(self.required_masses)
                missing_masses = required.difference(masses)
                missing_values = {k:default_masses.get(k,0.) 
                                  for k in missing_masses}
                if missing_masses: # Deal with unassigned fermion masses
                    repr_default = ['{}={: .5e}'.format(particle_names[k], 
                                     default_inputs.get(k,0.)) for 
                                     k in missing_masses]
                    print 'Warning: Not all required masses are '\
                          'defined.'
                    print 'Required PIDs: {}'.format(PID_list)
                    print 'Missing PIDs: {}'.format(mass_list)
                    carry_on = Y_or_N(
                        'Continue assuming default values '\
                        'for unspecified masses? '\
                        '({})'.format(', '.join(repr_default))
                        )
                    if carry_on: 
                        for m in missing_inputs: 
                            self.mass.new_entry(m, default_masses[m], 
                                                name='M%s' % particle_names[m])
                    else:
                        print 'Exit'
                        sys.exit()
        print 'Particle masses are OK.'

    def check_param_data(self, silent=False):
        '''
        Cross check of parameter data read in the SLHA formatted card.
        Compares lists of coefficients declared in self.independent and 
        self.dependent to those read in from self.param_card.
           1)  Deals with unrecognised names by removing them from 
               self.par_dict
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
        if silent: return        
        for bname, defined in self.blocks.iteritems():
            
            # collect block info
            inputblock = self.card.blocks.get(bname,SLHA.NamedBlock())
            input_eles = set(inputblock.keys())
            
            defined_block = {i:v for i,v in enumerate(defined)}
            defined_eles = set(defined_block.keys())
            
            independent = {i:v for i,v in defined_block.iteritems() 
                           if v in self.independent}
            dependent = {i:v for i,v in defined_block.iteritems() 
                           if v in self.dependent}
                           
            # check for unrecognised coefficient numbers              
            unknown = input_eles.difference(defined_eles)
            if unknown:
                unknown_names = {i:inputblock.get_name(i,'none') 
                                 for i in unknown}
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
                elif input_name!=name:
                    mismatched.append((index,name,input_name))
                    
            if mismatched:
                print 'Warning: Mismatch of coefficient names '\
                      'in {}, block "{}".'.format(self.__class__,
                                                          bname)
                for index, name, input_name in mismatched:
                    print ('    Coefficient {}, named {} '.format(index, input_name) +
                          'will be renamed to {}'.format(name))
                    inputblock._names[index] = name
                    inputblock._numbers[name] = index
            
            # check if coefficients defined as dependent 
            # have been assigned values
            defined_dependent = set(dependent).intersection(input_eles)
            if defined_dependent: 
                print 'Warning: you have assigned values to some '\
                      'coefficients defined as dependent '\
                      'in {}, block "{}".'.format(self.__class__, bname)
                print 'Coefficients: {}'.format(', '.join(
                                                  ['{}:"{}"'.format(k,v) 
                                                   for k,v in dependent.items() 
                                                   if k in defined_dependent]
                                                   ))
                print 'These may be overwritten by an implementation of '\
                      '{}.translate()'.format(self.__class__.__name__)
                carry_on = Y_or_N('Continue?')
                if carry_on:
                    pass
                else:
                    print 'Exit'
                    sys.exit()
            
            # check if all independent coefficients are assigned values
            missing = set(independent).difference(input_eles)
            if missing: # Deal with unassigned independent coefficients
                print 'Warning: some coefficients defined as independent '\
                      'in {}, block "{}", have not been assigned values'\
                      '.'.format(self.__class__,bname)
                print 'Undefined: {}'.format(', '.join(
                                                ['{}:"{}"'.format(k,v) 
                                                 for k,v in independent.items() 
                                                 if k in missing]))
                carry_on = Y_or_N('Continue assuming unspecified '\
                                  'coefficients are Zero?')
                if carry_on:
                    for m in missing: 
                        inputblock.new_entry(m, 0., independent[m])
                else:
                    print 'Exit'
                    sys.exit()
        
    def write_param_card(self,filename,overwrite=False):
        '''Write contents of self.newcard to filename'''
        preamble = ('\n###################################\n'
                    + '## INFORMATION FOR {}\n'.format(self.name.upper())
                    + '###################################\n')
        if 'basis' in self.card.blocks:
            self.card.blocks['basis'].preamble=preamble
            self.card.blocks['basis'][0]=self.newname
            
        preamble_2 = ('\n###################################\n'
                    + '## DECAY INFORMATION\n'
                    + '###################################\n')
        for decay in self.card.decays.values():
            decay.preamble=preamble_2
            break
        
        preamble = ( '################################################'
                    +'######################\n'
                    +'############# COEFFICIENTS TRANSLATED BY ROSETTA'
                    +' MODULE  #############\n'
                    +'########### PARAM CARD GENERATED {}  ##########'\
                     '#\n'.format(datetime.datetime.now().ctime().upper())
                    +'################################################'
                    +'######################\n\n')
        if os.path.exists(filename) and not overwrite:
            print '{} already exists.'.format(filename)
            carry_on = Y_or_N('Overwrite?')
        else:
            carry_on=True
        if carry_on:
            order = ['loop','mass','sminputs','yukawa','basis']
            self.newcard.write(filename, blockorder=order, preamble=preamble)
            print 'Wrote "{}".'.format(filename)
        else: 
            return False

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
                self.card.add_entry(bname, fields.index(entry), 0., name=entry)

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
    
    def calculate_dependent(self):
        print 'Nothing done for {}.calculate_'\
              'dependent()'.format(self.__class__.__name__)
    
    def translate(self): # default behaviour for translate()
        return self
        
    def eHDECAY_inputs(self):
        raise NotImplementedError
        
########################################################################

def flavour_matrix(name, kind='hermitian', domain='real', 
                   dimension=3):
    '''
    Function to declare flavour components of a coefficient according to its 
    properties. Takes a parameter name as an arguement and returns a list of 
    coefficients with the string suffixed by:
        'ij' indices for matrix coefficients (e.g. 12, 22, 13...)
        '_Im' and '_Re' for complex parameters
    '''
    index = range(1,dimension+1)
    
    if (kind, domain) == ('hermitian', 'complex'):
        real = ['{0}{0}'.format(i) for i in index]
        cplx = ['{}{}'.format(i,j) for i,j in combinations(index,2)]

    elif (kind, domain) == ('symmetric', 'complex'):
        real = []
        cplx = ['{}{}'.format(i,j) for i,j in combinations2(index,2)]
        
    elif ((kind, domain) == ('hermitian', 'real') or
          (kind, domain) == ('symmetric', 'real')):
        real = ['{}{}'.format(i,j) for i,j in combinations2(index,2)]
        cplx = []
        
    elif (kind, domain) == ('general', 'real'):
        real = ['{}{}'.format(i,j) for i,j in product(index,index)]
        cplx = []
        
    elif (kind, domain) == ('general', 'complex'):
        real = []
        cplx = ['{}{}'.format(i,j) for i,j in product(index,index)]
        
    else:
        print ('flavour_matrix function got and unrecognised combination of '
               '"kind" and "domain" keyword arguments')
        return [name]
        
    coefficients = []
    for indices in real:
        coefficients.append('{}{}'.format(name,indices))
    for indices in cplx:
        coefficients.append('{}{}_Re'.format(name,indices))
        coefficients.append('{}{}_Im'.format(name,indices))

    return coefficients
        
