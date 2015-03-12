import StringIO
import re, sys, os, math, datetime
from collections import namedtuple,OrderedDict
from query import query_yes_no as Y_or_N
from Rosetta import PID, default_inputs, default_masses

####################################################################################################
# Base Basis class
class Basis(object):
    '''
    Base class from which to derive other Higgs EFT basis classes. 
    
    Designed to be instantiated with an SLHA format parameter card using
    read_param_card(). The contents of the parameter card are stored in 
    self.param_card, while relevant information from blocks "basis", 
    "newcoup", "mass" and "sminput" are stored in self.name, self.par_dict, 
    self.mass and self.input respectively.
    
    self.name     - string storing the value of the 0th element of Block basis.
    self.par_dict - ``OrderedDict`` of (name,value) pairs with the name taken to be 
                    the first non-whitespace characters after a "#" character
                    in a parameter definition within Block newcoup.
    self.mass     - ``OrderedDict`` of (PID,value) pairs within Block mass
    self.input    - ``OrderedDict`` as for self.par_dict but within Block sminput.
    
    self.independent and self.dependent should be defined in accordance 
    with the contents of Block newcoup, as well as self.required_inputs and 
    self.required_masses with Blocks sminput and mass respectively 
    (the Z ahd Higgs masses are also stored in SM input). A number of checks
    related to these definitions are performed by check_param_data() on the 
    data read in from the parameter card.
    
    Basis also stores the coefficients in self.coeffs as a ``namedtuple`` 
    which stores the coefficients as data members and from which an ``OrderedDict``
    can also be recovered  with self.coeffs._asdict()
        
        type(self.par_dict)==type(self.coeffs._asdict()) # True
        
        self.par_dict['myvar'] = 10.
        self.coeffs._asdict()['myvar'] == 10. # True
        self.coeffs.myvar == 10. # True
    
    The user should define methods calculate_dependent() and translate() to 
    firstly calculate any dependent parameters and then translate the coefficients 
    into the mass basis and populate self.newpar, a template of which can be 
    obtained by creating an instance of HEFTRosetta.MassBasis.MassBasis() without 
    the optional param_card argument. 
        
        from MassBasis import MassBasis
        new_instance = MassBasis() # Empty MassBasis instance with all coeffs set to 0
        new_dict = new_instance.coeffs._asdict()
    
    set_newcard() uses the contents of self.newpar to write a modified parameter card 
    in the mass basis and write_param_card() saves it to file.
    
    Users should call translate() and set_newcard() in the __init__ of their derived class.
    
    Derived classes can be used by the command line script "translate"
    '''
    independent, dependent=[], []
    required_inputs, required_masses = set(),set()
    translate_to = {'mass'}
    def __init__(self, param_card=None, block_in='newcoup', block_out='newcoup', output_basis='mass', keep_old=True):
        # Check that target basis has a translation implemented
        if output_basis in self.translate_to:
            self.target_basis = output_basis
        else:
            raise ValueError('''
            {}.translate_to does not contain "{}". 
            Rosetta doesn't know of an implemented translation in {}
            '''.format(self.__class__.__name__,output_basis, self.__class__))
            
        self.param_card = param_card
        self.block_in, self.block_out, self.keep_old = block_in, block_out, keep_old # other options
        self.card, self.newcard = StringIO.StringIO(), StringIO.StringIO()
        self.input, self.mass, self.par_dict, self.newpar, self.newmass = OrderedDict(),OrderedDict(),OrderedDict(),OrderedDict(),OrderedDict()
        self.newname = 'Basis'
        self.all_coeffs = self.independent + self.dependent
        
        if param_card is not None: # read param card (sets self.par_dict, self.input, self.mass, self.name, self.card)
            assert os.path.exists(self.param_card), '{} does not exist!'.format(self.param_card)
            self.read_param_card()
            self.check_param_data()
            self.calculate_dependent()
            self.check_calculated_data()
            coeffclass = namedtuple(self.name, self.all_coeffs)
            self.coeffs = coeffclass(**self.par_dict)
        else: # if param_card option not given, instatiate with class name and all coeffs set to 0 (used for creating an empty MassBasis instance for use in translate() method)
            self.name=self.__class__.__name__
            coeffclass = namedtuple(self.name, self.all_coeffs)
            self.coeffs = coeffclass(**{k:0. for k in self.all_coeffs})
            
        self.translate()
        self.set_newcard()
        
    def read_param_card(self):
        '''
        Reads SLHA style param card and stores values of parameters.
        Looks for Blocks "basis", "mass", "newcoup" and "sminput" and 
        updates self.name, self.mass (by PDG number), self.par_dict (name, value) 
        and self.input (name, value) respectively. 
        The whole param_card are also stored in self.card for later use.
        '''
        def read_pattern(plines, patt, pdict, ikey=2, ival=1, convert_key=str, convert_val=float):
            '''Takes a list of lines and looks for regexp pattern "patt" assuming the pattern captures 2 groups.
               "key" and "val" refer to the index of the resulting group (1 or 2).
               "convert_val" and "convert_key" specify a functions with which to treat the resulting string match e.g. convert to float.
               If pdict is None, function returns convert_val( match.group(val) ) (used to assign one variable).
               Otherwise, pdict is appended with a key value pair ( convert_key(match.group(key)), convert_val(match.group(val)) ).
               Also prints each line to self.card'''
            for pline in plines[:-1]:
                comment_or_block = (pline and not pline.startswith('#') and not 'block' in pline.lower())
                print >> self.card, '    ' + pline.strip('\n') if comment_or_block else pline.strip('\n')
                try:
                    match = re.match(patt,pline)
                    if pdict is not None:
                        key, val = convert_key(match.group(ikey)), convert_val(match.group(ival))
                        if pdict.has_key(key): # check if variable has already been assigned
                            label = 'PDG mass' if type(key)==int else 'variable'
                            print 'Warning: {} "{}" assigned more than once, kept value {}'.format( label,key,val )
                        pdict[key]= val
                    else:
                        return convert_val(match.group(ival))
                except AttributeError:
                    pass
            return None
                
        with open(self.param_card,'r') as card:
            relevant_blocks =['mass','basis',self.block_in,'sminputs']
            lines = iter(card)
            for line in lines:
                line, block = self.line_is_block(line)
                while block in relevant_blocks:
                    if block=='mass': # Get masses
                        print >> self.card, line.strip('\n')
                        param_lines = self.read_until(lines,'Block','DECAY') # lines in block mass
                        read_pattern(param_lines, r'\s*(\d+)\s+(\S+)\s+.*', self.mass, ikey=1, ival=2, convert_key=int) # by PDG number
                        relevant_blocks.remove('mass')
                        line, block = self.line_is_block(param_lines[-1]) # set last line read to current line
                    
                    elif block=='basis': # Get basis name
                        print >> self.card, line.strip('\n')
                        param_lines = self.read_until(lines,'Block','DECAY') # lines in block basis
                        self.name = read_pattern(param_lines,r'\s*0\s+(\S+).*', None, convert_val=str).strip()
                        relevant_blocks.remove('basis')
                        line, block = self.line_is_block(param_lines[-1]) # set last line read to current line
                    
                    elif block==self.block_in: # Get basis coefficients
                        print >> self.card, line.strip('\n')
                        param_lines = self.read_until(lines,'Block', 'DECAY') # lines in block newcoup
                        read_pattern(param_lines, r'\s*\d+\s+(\S+)\s+#+\s+(\S+)', self.par_dict)
                        relevant_blocks.remove(self.block_in)
                        line, block = self.line_is_block(param_lines[-1]) # set last line read to current line

                    elif block=='sminputs': # Get SM inputs
                        print >> self.card, line.strip('\n')
                        param_lines = self.read_until(lines,'Block', 'DECAY') # lines in block sminputs
                        read_pattern(param_lines, r'\s*\d+\s+(\S+)\s+#+\s+(\S+)', self.input)
                        relevant_blocks.remove('sminputs')
                        line, block = self.line_is_block(param_lines[-1]) # set last line read to current line
                    
                if not block or (block not in relevant_blocks): print >> self.card, line.strip('\n')
            if relevant_blocks:
                if relevant_blocks==['basis']: # only optional basis block unread
                    print 'Warning: block basis not found.'
                    carry_on = Y_or_N('Continue assuming default value "{}"?'.format(self.__class__.__name__))
                    if carry_on:
                        self.name=self.__class__.__name__
                    else:
                        print 'Exit'
                        sys.exit()
                else:
                    raise IOError("Required block(s) ({}) not found in {}".format(', '.join(relevant_blocks), self.param_card))
            
        if self.mass.has_key(23): self.input['MZ']=self.mass[23] # MZ counts as SM input
        if self.mass.has_key(25): self.input['MH']=self.mass[25] # MH counts as SM input
        
    def check_param_data(self):
        '''
        Compares lists of coefficients declared in self.independent and self.dependent to those read in from self.param_card.
           1) Deals with unrecognised names by removing them from self.par_dict
           2) Prints a warning if coefficients declared as dependent are assigned values in param_card.
              If so, the user is asked whether or not they wish to continue.
           3) Checks if all coefficients declared as independent are assigned values. 
              If not, the user is given the option to continue with them set to 0.
           4) Ensures all fermion masses declared in self.required_masses are defined.
              If not, the user is given the option to continue with them set to 0.
           5) Ensures all sm imputs declared in self.required_inputs are defined.
              If not, the user is given the option to continue with them set to default 
              values defined in HEFTrosetta.default_inputs.
        '''

        unknown_coeffs = set(self.par_dict.keys()).difference( self.all_coeffs )
        if unknown_coeffs: # Check for unrecognised coeff names
            print 'Warning: you have declared coefficients undefined in {}.'.format(self.__class__)
            print 'The following will be ignored: {}'.format(','.join(unknown_coeffs))
            for c in unknown_coeffs: del self.par_dict[c]
            
        defined_dependent = set(self.dependent).intersection( set(self.par_dict.keys()) )
        if defined_dependent: # Check if coefficients defined as dependent have been assigned values
            print 'Warning: you have assigned values to some coefficients defined as dependent in {}.'.format(self.__class__)
            print 'Coefficients: {}'.format(','.join(defined_dependent))
            print 'These may be overwritten by an implementation of {}.translate()'.format(self.__class__.__name__)
            carry_on = Y_or_N('Continue?')
            if carry_on:
                pass
            else:
                print 'Exit'
                sys.exit()
                
        missing_coeffs = set(self.independent).difference( set(self.par_dict.keys()) )
        if missing_coeffs: # Deal with unassigned independent coefficients
            print 'Warning: Set of independent coefficients read from {} does not match those specified in {}.'.format(self.param_card,self.__class__)
            print 'Undefined: {}'.format(', '.join(missing_coeffs))
            carry_on = Y_or_N('Continue assuming unspecified coefficients are Zero?')
            if carry_on:
                for m in missing_coeffs: self.par_dict[m]=0.
            else:
                print 'Exit'
                sys.exit()
                
        missing_masses = set(self.required_masses).difference(self.mass.keys())
        repr_default_masses = ['{}={:.5e}'.format(k,default_masses[k]) for k in missing_masses]
        if missing_masses: # Deal with unassigned fermion masses
            print 'Warning: Not all required fermion masses are defined in {}.'.format(self.param_card)
            print 'Required PIDs: {}'.format(', '.join([str(x) for x in self.required_masses]))
            print 'Missing PIDs: {}'.format(', '.join([str(x) for x in missing_masses]))
            carry_on = Y_or_N('Continue assuming default values for unspecified masses? ({})'.format(', '.join(repr_default_masses)))
            if carry_on:
                for m in missing_masses: self.mass[m]=0.
            else:
                print 'Exit'
                sys.exit()
                
        missing_inputs = set(self.required_inputs).difference(self.input.keys())
        repr_default_inputs = ['{}={:.5e}'.format(k,default_inputs[k]) for k in missing_inputs]
        if missing_inputs: # Deal with unassigned SM inputs
            print 'Warning: Not all requires SM inputs are defined in {}.'.format(self.param_card)
            print 'Required inputs: {}'.format(', '.join([str(x) for x in self.required_inputs]))
            print 'Missing inputs: {}'.format(', '.join([str(x) for x in missing_inputs]))
            carry_on = Y_or_N('Continue with default values for unspecified inputs? ({})'.format(', '.join(repr_default_inputs)))
            if carry_on:
                for m in missing_inputs: self.input[m]=default_inputs[m]
            else:
                print 'Exit'
                sys.exit()
        print 'Param card data are OK.'
    
    def check_calculated_data(self):
        missing_dependents = set(self.dependent).difference(self.par_dict.keys())
        if missing_dependents and self.dependent:
            print 'Warning: Set of dependent coefficients calulcated by {0}.calculate_dependent() does not match those specified in {0}.'.format(self.__class__)
            print 'Undefined: {}'.format(', '.join(missing_dependents))
            carry_on = Y_or_N('Continue assuming coefficients are Zero?')
            if carry_on:
                for m in missing_dependents: self.par_dict[m]=0.
            else:
                print 'Exit'
                sys.exit()   
            print 'Calculated coefficients match those defined in {}.dependent.'.format(self.__class__.__name__)         
    
    def set_newcard(self):
        '''
        Generates a new param card in self.newcard from self.card, adding the new set of coefficients in self.newpar after "Block newcoup".
        If self.keep_old is True, the original names and values of newcoup variables are included, commented out.
        '''
        print >> self.newcard, '######################################################################'
        print >> self.newcard, '############# COEFFICIENTS TRANSLATED BY ROSETTA MODULE  #############'
        print >> self.newcard, '########### PARAM_CARD GENERATED {}  ###########'.format(datetime.datetime.now().ctime().upper())
        print >> self.newcard, '######################################################################'
        blocks_to_modify =('basis',self.block_in, 'mass')
        lines = iter(self.card.getvalue().splitlines()) # Lines of old param card
        for line in lines:
            line, block = self.line_is_block(line)
            if block==self.block_in:
                print >> self.newcard, line.strip('\n').replace(self.block_in,self.block_out) # Add line to self.newcard, modify block name
            else:
                print >> self.newcard, line.strip('\n') # Add line to self.newcard
            while block in blocks_to_modify:
                
                if block=='basis': # When Block basis is reached
                    print >> self.newcard, '    0 {} # translated basis'.format(self.newname)
                    param_lines = self.read_until(lines,'Block') 
                    for pline in param_lines[:-1]:
                        if self.keep_old:
                            comment_out = (pline and not pline.startswith('#'))
                            print >> self.newcard, '#    ' + pline.strip('\n') if comment_out else pline.strip('\n')
                        else:
                            if not pline.strip() or pline.startswith('#'): print >> self.newcard, pline.strip('\n')
                    line, block = self.line_is_block(param_lines[-1]) # set last line read to current line
                    
                elif block==self.block_in: # When Block self.block_in is reached,
                    for i,(par,val) in enumerate(self.newpar.items()):  # write out new couplings
                        print >> self.newcard, '    {} {:.5e} # {}'.format(i,val,par)
                    param_lines = self.read_until(lines,'Block') 
                    if self.keep_old:
                        print >> self.newcard, ''
                        print >> self.newcard, '###################################'
                        print >> self.newcard, '## COEFFICIENTS IN {} BASIS'.format(self.name.upper())
                        print >> self.newcard, '###################################'
                        for pline in param_lines[:-1]:
                            comment_out = (pline and not pline.startswith('#'))
                            print >> self.newcard, '#    ' + pline.strip('\n') if comment_out else pline.strip('\n')
                    else:
                        for pline in param_lines[:-1]:
                            if not pline.strip() or pline.startswith('#'): print >> self.newcard, pline.strip('\n')
                    line, block = self.line_is_block(param_lines[-1]) # set last line read to current line
                
                elif block=='mass': # When block mass is reached
                    param_lines = self.read_until(lines,'Block') 
                    for pline in param_lines[:-1]:
                        match = re.match(r'\s*(\d+)\s+(\S+)\s+.*',pline.strip()) # read old mass
                        try:
                            ID, mass = int(match.group(1)), match.group(2)
                            if self.newmass.has_key(ID):  # if a new value exists
                                print >> self.newcard, pline.replace(mass, '{:.5e}'.format(self.newmass[ID])).strip('\n')
                                if self.keep_old: print >> self.newcard, '# '+pline.strip('\n')+' # old value'
                            else:
                                print >> self.newcard, pline.strip('\n')
                        except AttributeError:
                            print >> self.newcard, pline.strip('\n')
                    line, block = self.line_is_block(param_lines[-1]) # set last line read to current line
                     
                if block==self.block_in:
                    print >> self.newcard, line.strip('\n').replace(self.block_in,self.block_out) # Add line to self.newcard, modify block name
                else:
                    print >> self.newcard, line.strip('\n') # Add line to self.newcard     
        return None
        
    def write_param_card(self,filename,overwrite=False):
        '''Write contents of self.newcard to filename'''
        contents = self.newcard.getvalue()
        if not contents: 
            print '{}.newcard is empty, nothing done.'.format(self)
            return None
        if os.path.exists(filename) and not overwrite:
            print '{} already exists.'.format(filename)
            carry_on = Y_or_N('Overwrite?')
        else:
            carry_on=True
        if carry_on:
            with open(filename,'w') as param_card:
                param_card.write(contents)
            return True
        else: 
            return False

    @staticmethod
    def read_until(lines, here, *args):
        '''Loops through an iterator of strings by calling next() until 
           it reaches a line starting with a particular string.
           Case insensitive.
           Args:
               lines - iterator of strings
               here - string (plus any further argumnts). Reading will end if the line matches any of these.
           Return:
               lines_read - list of lines read
               line - last line that was read (containing string "here")
        '''
        end_strings = [here.lower()]+[a.lower() for a in args]
        lines_read = []
        line = ''
        while not any([line.lower().startswith(x) for x in end_strings]): 
            line = lines.next().strip()
            lines_read.append(line.strip())
        return lines_read
        
    @staticmethod
    def line_is_block(line):
        '''
        Looks for "BLOCK XXXX" pattern in line.
        returns (line, XXXX) if line matches.
        returns (line, False) if line doesn't match
        '''
        try:
            return line, re.match(r'block\s+(\S+).*',line.strip().lower()).group(1)
        except AttributeError:
            return line, False
            
    def calculate_dependent(self):
        print 'Nothing done for {}.calculate_dependent()'.format(self.__class__.__name__)
    
    def translate(self): # default behaviour for translate()
        self.keep_old=False
        self.newpar = self.coeffs._asdict()
####################################################################################################