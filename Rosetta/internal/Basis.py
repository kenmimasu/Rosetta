import StringIO
import re
import sys
import os
import math
import datetime
from collections import namedtuple,OrderedDict
from itertools import product, combinations
from itertools import combinations_with_replacement as combinations2
import SLHA
from query import query_yes_no as Y_or_N
from __init__ import (PID, default_inputs, default_masses, input_names, 
                      particle_names, input_to_PID, PID_to_input, RosettaError)
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
    blocks = dict()
    independent = []
    required_inputs, required_masses = set(),set()
    translate_to = {'mass'}
    
    def __init__(self, param_card=None, output_basis='mass', 
                 ehdecay=False, silent=False, translate=True):
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
            
        self.param_card = param_card
        self.output_basis = output_basis
        self.newname = 'Basis'
        
        self.all_coeffs = [c for v in self.blocks.values() for c in v ]
        self.dependent = [c for c in self.all_coeffs 
                          if c not in self.independent]
        
        # make blocks case insensitive
        self.blocks = {k.lower():v for k,v in self.blocks.iteritems()} 
        # read param card (sets self.inputs, self.mass, self.name, self.card)
        if param_card is not None: 
            
            assert os.path.exists(self.param_card), \
                   '{} does not exist!'.format(self.param_card)
            
            self.read_param_card() 
            # various input checks 
            self.check_sminputs(self.required_inputs) 
            self.check_masses(self.required_masses, silent=silent) 
            self.check_param_data(silent=silent)
            # add dependent coefficients/blocks to self.card
            self.init_dependent()
            # user defined function
            self.calculate_dependent()
            if translate:
                # translate to new basis (User defined) 
                # return an instance of a class derived from Basis
                newbasis = self.translate() 
                # set new SLHA card
                self.newcard = newbasis.card 
            
                preamble = ('###################################\n'
                          + '## DECAY INFORMATION\n'
                          + '###################################')
                for i,decay in enumerate(self.card.decays.values()):
                    decay.preamble=preamble
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
            self.card = SLHA.Card(name=self.name)
            self.card.add_entry('basis', 0, self.name, name = 'translated basis')
            
            preamble = ('###################################\n'
                        + '## INFORMATION FOR {}\n'.format(self.name.upper())
                        + '###################################')
            self.card.blocks['basis'].preamble=preamble
            
            # default behaviour: create one 'newcoup' block
            if not self.blocks: 
                theblock = SLHA.NamedBlock(name='newcoup')
                for i,fld in enumerate(self.all_coeffs):
                    theblock.new_entry(i,0., name=fld)
                self.card.add_block(theblock)
            # otherwise follow self.blocks structure
            else: 
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
        
        if not self.blocks:
            # create self.blocks structure for independent parameters
            for par in self.independent:
                block = self.card._parent_block(par)
                name = block.name
                if name not in self.blocks:
                    self.blocks[name] = [par]
                else:
                    self.blocks[name].append(par)
            # default behaviour, create new block for dependent parameters
            self.blocks['dependent']=self.dependent

        self.inputs = self.card.blocks.get('sminputs',None)
        self.mass = self.card.blocks.get('mass',None)
        try:
            self.name = self.card.blocks.get('basis',[''])[0].lower()
        except KeyError:
            raise RosettaError('Formatting error for block basis.'
                                'Check input card.')
        
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
                missing_values = {k:default_inputs.get(k,0.) for k in missing_inputs}
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
        # print 'SM inputs are OK.'
        
    def check_masses(self, required_masses, silent=False, message='Rosetta'):
        '''
        Check consistency of particle masses w.r.t required_masses. Any 
        inconsistencies found will raise a warning and a question of whether 
        the user wishes to continue assuming default values for masses. Higgs 
        and Z masses are also stored as SM inputs.
        '''
        if silent: return
        PID_list = ', '.join([str(x) for x in required_masses])
        if required_masses:
            if self.mass is None:
                repr_default = ['{}={: .5e}'.format(particle_names[k], 
                                 default_inputs.get(k,0.)) for 
                                 k in required_masses]
                mass_list = ', '.join( ['{} (M{})'.format(x,particle_names[x]) 
                                        for x in required_masses] )
                print 'Warning: Block "mass" not found. '\
                      'Assume default values for required masses?'
                print 'Required PIDs: {}'.format(mass_list)

                carry_on = Y_or_N(
                    'Continue with default values '\
                    'for unspecified masses? '\
                    '({})'.format(', '.join(repr_default))
                                 )
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
                                if i in (4,8)]:
                        i = input_to_PID[k]
                        if i in self.mass:
                            v2 = float(self.mass[i])
                            if v!=v2:
                                print ('Warning: M{} '.format(particle_names[i])
                                + 'specified in block sminput[{}] '.format(k)
                                + '({:.5e}) not consistent with '.format(v)
                                + 'value specified in block mass '
                                + '[{}] ({:.5e}).\n'.format(i,float(v2))
                                + 'Rosetta will keep value from sminputs.')
                                self.mass[k]=v
                        else:
                            mstring = 'M{}'.format(particle_names[i])
                            self.mass.new_entry(i, v,name=mstring)
                            
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
                    print 'Warning: Not all required masses are '\
                          'defined for {}.'.format(message)
                    print 'Required PIDs: {}'.format(PID_list)
                    print 'Missing PIDs: {}'.format(mass_list)
                    carry_on = Y_or_N(
                        'Continue assuming default values '\
                        'for unspecified masses? '\
                        '({})'.format(', '.join(repr_default))
                        )
                    if carry_on: 
                        for m in missing_masses: 
                            self.mass.new_entry(m, default_masses[m], 
                                                name='M%s' % particle_names[m])
                            
                    else:
                        print 'Exit'
                        sys.exit()
        # print 'Particle masses are OK.'

    def check_param_data(self, silent=False):
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
        if 'basis' in self.newcard.blocks:
            self.newcard.blocks['basis'].preamble=preamble
            self.newcard.blocks['basis'][0]=self.newname
            
        preamble_2 = ('\n###################################\n'
                    + '## DECAY INFORMATION\n'
                    + '###################################\n')
        for decay in self.newcard.decays.values():
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
            return True
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
        '''
        Default behavoiur of calculate_dependent(). Called if a subclass of 
        Basis has not implemeneted the function.
        '''
        print 'Nothing done for {}.calculate_'\
              'dependent()'.format(self.__class__.__name__)
    
    def translate(self, target=None, verbose=True): 
        '''
        Translation function. Makes use of known translations in implemented.py 
        and attempts to find a path between the basis instance and the target 
        basis. If a path exists, successive translations are performed to get 
        from A to B.
        Keyword arguments:
            target - target basis, must be defined in implemented.bases.
            verbose - print some status messages.
        '''
        
        from .. import implemented
        from relate import get_path
        
        # default target
        target = target if target is not None else self.output_basis
        
        if target == self.name: return self
        
        # find the chain of bases and translation functions 
        # to get from self.name to target
        chain = get_path(self.name.lower(), target, implemented.relationships)

        if not chain: # no translation possible
            inputbasis = self.__class__.__name__
            outputbasis = implemented.bases[target].__name__
            err = ('Sorry, Rosetta cannot translate from ' +
                  '{} to {}'.format(inputbasis, outputbasis))
            raise TranslationError(err)
        
        names = [self.name.lower()]+[x[0] for x in chain]
        if verbose:
            print 'Rosetta will be performing the translation:'
            print '\t'+' -> '.join([implemented.bases[x].__name__ 
                                                    for x in names])
        
        # perform succesive translations, checking for  
        # required SM inputs along the way
        current = self
        required_inputs = current.required_inputs
        required_masses = current.required_masses
        for target, function in chain:
            # new target basis instance
            instance = implemented.bases[target]()
            # ensure required inputs are present
            message = 'translation ({} -> {})'.format(current.name, 
                                                      instance.name)
            current.check_sminputs(required_inputs, message=message)
            current.check_masses(required_masses, message=message)

            # call translation function
            new = function(current, instance)
            # update new basis instance with non-EFT blocks, decays
            new.inputs, new.mass = current.inputs, current.mass
            
            other_blocks = {k:v for k,v in current.card.blocks.iteritems() 
                            if not (k in current.blocks or k.lower()=='basis')}
                            
            for block in other_blocks:
                theblock = current.card.blocks[block]
                new.card.add_block(theblock)
            
            for i,decay in enumerate(current.card.decays.values()):
                new.card.add_decay(decay, preamble = decay.preamble)
            
            # new.calculate_dependent()
            current.newname=new.name
            # prepare for next step
            required_inputs = set(instance.required_inputs)
            required_masses = set(instance.required_masses)
            current = new
        if verbose:
            print 'Translation successful.'

        return current
            
    def run_eHDECAY(self):
        '''
        Interface Rosetta with eHDECAY to calculate Higgs widths and branching 
        fractions. 
        '''
        import eHDECAY
        BRs = eHDECAY.run(self,electroweak=True)
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
        hdecays = {}
        
        # sometimes eHDECAY gives negative branching fractions. 
        # print warning and leave out that channel.
        for channel, BR in BRs.iteritems(): 
            n1, n2 = particle_names[channel[0]], particle_names[channel[1]]
            comment = 'H -> {}{}'.format(n1,n2)
            if BR < 0.: 
                print ('eHDECAY: Negative branching fraction encountered ' +
                      'for {}. Rosetta will ignore it.'.format(comment))
                continue
            elif BR == 0.:
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
        
########################################################################
class TranslationError(Exception):
    '''Fancy error name.'''
    pass
    
def flavour_matrix(name, kind='hermitian', domain='real', 
                   dimension=3):
    '''
    Function to declare flavour components of a coefficient according to its 
    properties. Takes a parameter name as an arguement and returns a list of 
    coefficients with the string suffixed by:
        'ij' indices for matrix coefficients (e.g. 12, 22, 13...)
        'Im' and 'Re' for complex parameters
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
        coefficients.append('{}{}Re'.format(name,indices))
        coefficients.append('{}{}Im'.format(name,indices))

    return coefficients
