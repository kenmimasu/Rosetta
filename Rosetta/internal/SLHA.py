from collections import OrderedDict
import sys
import re


class Block(OrderedDict):
    '''    
    Container class for SLHA block. A subclass of `collections.OrderedDict`with 
    a restriction on integer keys and a modified __repr__(). Block can't be 
    initialised with any positional arguments but rather with the 'data' 
    keyword argument. It can optionally be named using the 'name' keyword 
    argument. __repr__ and __str__ functions are also defined to output the 
    SLHA formatted block.
    '''

    def __checkkey__(self,key):
        '''Forces key to be of integer type as in SLHA convention.'''
        if type(key) is not int: 
            raise TypeError( self.__class__.__name__ +
                             ' only accepts integer keys.' ) 
        else:
            return key
    
    def __init__(self, name=None, data=None, 
                 decimal=5, dtype=lambda x: x, preamble=''):
        '''    Intialisation keyword arguments:
            name - A name for the block that will appear 
                   in the __str__ and __repr__ methods.
            data - A dict type object supporting iteritems() 
                   to initialize Block values.
            decimal - Number of decimal points with which 
                      to write out parameter values.
            dtype - a function to cast parameter values if 
                    read in as strings i.e. int, float.
            preamble - Some text to print before the __str__ output.
        '''
        super(Block, self).__init__()
        self.name = name
        self.cast = dtype # value casting function
        self.fmt = ':+.{}e'.format(decimal) # format string
        self.preamble = preamble
        if data is not None:
            for k,v in data.iteritems():
                self[k]=v
                
    def __setitem__(self,key,value):
        return super(Block,self).__setitem__(self.__checkkey__(key),
                                             self.cast(value))

    def __repr__(self):
        return ( '<SHLA Block: "{}"; {} entries.>'.format(self.name,len(self)) )

    def __str__(self):
        content = []
        for k,v in self.iteritems():
            try:
                val = float(v)
                fmt = self.fmt
            except ValueError:
                fmt = ''
            
            line = ('    {{: <4}} {{{}}}\n'.format(fmt)).format
            content.append(line(k,val)) 
            
        string = (self.preamble+'\n'
                  +'BLOCK {}\n'.format(self.name)
                  + ''.join(content) + '\n')
        return string
    
    def dict(self):
        '''Return SHLA Block object as a regular python dict.'''
        return {k:v for k,v in self.iteritems()}

class NamedBlock(Block):
    '''
    Class derived from 'Block' with the added functionality of assigning a name 
    to each key via the 'names' keyword argument. The values can then also be 
    referenced by name via a lookup in self._numbers. 
    '''
    
    def __parse__(self, key):
        '''    
        Allows for elements to additionally be referred to by name. 
        Looks up self._numbers for an assigned name and returns the 
        associated key. Rasies KeyError if name is not booked. 
        Used in overloading __setitem__, __getitem__, __contains__ 
        and __delitem__.
        '''
        if type(key) is str:
            try:
                return self._numbers[key.lower()]
            except KeyError:
                return None
                # err = ('Name "{}" has not been assigned to a '.format(key) +
               #         'key in SLHA.NamedBlock {}._names'.format(self.name))
               #  raise KeyError(err)
        else:
            return key
    
    def __init__(self, name=None, data=None, 
                 comment='', decimal=5, dtype=lambda x:x, preamble=''):
        '''    
        Same as the SLHA.Block constructor but additionally checks if the data 
        keyword argument has a "_names" attribute (i.e. if it is an existing 
        instance of NamedBlock) and stores it. The "comment" keyword argument 
        prints a comment to the right of the block declaration in __str__().
        '''        
        if data is not None and hasattr(data, '_names'):
            self._names = data._names
            try:
                self._numbers = {v:k.lower() for k,v 
                                      in self._names.iteritems()}
            except AttributeError:
                raise ValueError('"names" keyword argument must be a'
                                 'dictionary or support iteritems() method.')
        else:
            self._names, self._numbers = {}, {}

        self.comment=comment

        return super(NamedBlock, self).__init__(name=name, data=data, 
                                                decimal=decimal, dtype=dtype,
                                                preamble = preamble)
                                                
    def __setitem__(self, key, value):
        super(NamedBlock,self).__setitem__(self.__parse__(key),
                                             self.cast(value))

    def __getitem__(self, key):
        return super(NamedBlock,self).__getitem__(self.__parse__(key))
    
    def __delitem__(self, key):
        # Additional cleanup needed for _names and _numbers lookup dicts.
        if type(key) is str:
            key = key.lower()
            try:
                del self._names[self._numbers[key]]
                del self._numbers[key]
            except KeyError:
                pass
        else:
            try:
                del self._numbers[self._names[key]]
                del self._names[key]
            except KeyError:
                pass
                
        return super(NamedBlock,self).__delitem__(self.__parse__(key))
                     
    def __contains__(self, key):
        return super(NamedBlock,self).__contains__(self.__parse__(key))
        
    def __repr__(self):
        return ('<SHLA NamedBlock: "{}"; {} entries, {} named.>'.format(
                                        self.name, len(self), len(self._names)))
    def __str__(self):
        content = []
        for k,v in self.iteritems():
            try:
                strval = ('{{{}}}'.format(self.fmt)).format(float(v))
            except ValueError:
                strval = v

            if k in self._names:
                content.append('    {: <4} {} # {}'.format(k, strval,
                                                           self._names[k]))
            else:
                content.append('    {: <4} {}'.format(k,strval))

        string = (self.preamble + '\n'
                  + 'BLOCK {} # {}\n'.format(self.name, self.comment)
                  + '\n'.join(content) + '\n')
        return string

    def get_name(self, key, default=''):
        return self._names.get(key, default)
        
    def get_number(self, name, default=-1):
        return self._numbers.get(name.lower(),default)
    
    def new_entry(self, key, value, name=None):
        '''
        Add a new named block entry. Ensures that the key hasn't already 
        been set and the name given hasn't already been booked.
        '''
           
        if key in self:
            err = ("Key '{}' already belongs to NamedBlock ".format(key.lower())
                 +  "'{}', mapped to name '{}'.".format(self.name,
                                                        self.get_name(key)))
            raise KeyError(err)
            
        if name in self._numbers:
            err = ("Name '{}' already booked in ".format(name) +
                   "NamedBlock {}, mapped to key '{}'.".format(
                                   self.name, self._numbers[name]))
            raise KeyError(err)
        else:
            if name is not None:
                self._names[key] = name.lower()
                self._numbers[name.lower()] = key
            self[key]=value
            
    def namedict(self):
        '''Return python dict of name:key pairs.'''
        return {self._names[k]:v for k,v in self.iteritems() 
                                 if k in self._names }

class Matrix(Block):
    '''
    Caontainer class for multi-counter SLHA Blocks. Inherits from Block and is 
    intended to take further instances of Matrix or Block. items() and 
    iteritems() methods are implemented to recursively iterate through the 
    structure until non-Block values are encountered and return/yield a tuple 
    of keys referring to the indices along with the values e.g. 
    ((i,j,k,..), val).'''

    def __checkval__(self,val):
        '''Forces value to be an instance of SLHA.Block (or any of its subclasses).'''
        if not isinstance(val,Block): 
            raise TypeError( self.__class__.__name__ +
                             ' only accepts SLHA.Block type as values.' ) 
        else:
            return val
                
    def __setitem__(self,keys,value):
        if type(keys) is int:
            return super(Matrix,self).__setitem__(self.__checkkey__(keys),
                                                  self.__checkval__(value))
        else:
            key, subkeys = keys[0], keys[1:]
            tmp = self.__getitem__(key)
            if len(subkeys)==1: subkeys = subkeys[0]
            tmp[subkeys] = value

    def __getitem__(self,keys):
        getitem = super(Matrix,self).__getitem__
        if type(keys) is int:
            return getitem(keys)
        else:
            tmp=self
            for i,k in enumerate(keys):
                try:
                    tmp = tmp[k]
                except KeyError as e:
                    err = ', '.join([str(k) for k in keys]) + ('At index '+
                        '{} of {}'.format(i+1, self))
                    raise KeyError(err)
            return tmp
    
    def __len__(self):
        length = 0
        for val in self.values():
            length += len(val)
        return length
    
    def __repr__(self):
        return '<SHLA Matrix: "{}"; {} entries.>'.format(self.name, len(self.items()))        

    def __str__(self):
        content = []
        for k,v in self.iteritems():
            try:
                v = float(v)
                fmt = self.fmt
            except ValueError:
                fmt=''
            ind = '{: <4}'*len(k)
            line = ('    {} {{{}}}\n'.format(ind,fmt)).format
            args = list(k)+[v]
            content.append(line(*args))
        string = (self.preamble+'\n'
                  +'BLOCK {}\n'.format(self.name)
                  + ''.join(content) + '\n')
        return string
        
    def items(self):
        itms = []
        for k,v in zip(self.keys(),self.values()):
            if isinstance(v, Block):
                for sk, sv in v.items():
                    if type(sk) is int: 
                        sk = [sk]
                    newkeys = tuple([k] + list(sk))
                    itms.append((newkeys,sv))
            else:
                itms.append((k,v))
        return itms
    
    def iteritems(self):
        for x in self.items():
            yield x
    
    def dimension(self):
        for i,(k,v) in enumerate(self.iteritems()):
            if i!=0: subdim_prev = subdim
            
            if isinstance(v,Matrix):
                subdim = v.dimension()
            else:
                subdim = [len(v)]

            if i!=0 and subdim_prev != subdim:
                err = ('Inconsistent dimension of SLHA Matrix '
                      '"{}" encountered'.format(self))
                raise SLHAError(err)
        
        return [len(self)] + [d for d in subdim]
    
    def array(self):
        array = []
        for k,v in self.iteritems():
            if isinstance(v,Matrix):
                array.append(v.array())
            else:
                array.append([v for v in v.values()])
        return array
        
        
    def dict(self):
        thedict = {}
        for k,v in self.iteritems():
            thedict[k] = v.dict()
        return thedict

class NamedMatrix(NamedBlock, Matrix):
    '''
    Class derived from 'Block' with the added functionality of assigning a name 
    to each key via the 'names' keyword argument. The values can then also be 
    referenced by name via a lookup in self._numbers. 
    '''

    def __setitem__(self, key, value):
        return Matrix.__setitem__(self, self.__parse__(key), self.cast(value))

    def __getitem__(self, key):
        return Matrix.__getitem__(self, self.__parse__(key))
        
    def __len__(self):
        return Matrix.__len__(self)
   
    def __repr__(self):
        return ('<SHLA NamedMatrix: "{}"; {} entries, {} named.>'.format(
                                       self.name, len(self), len(self._names)))
   
    def __str__(self):
        content = []
        for k,v in self.iteritems():
            try:
                v = float(v)
                fmt = self.fmt
            except ValueError:
                fmt=''
            vname = self._names.get(k,None)
            comment = '' if vname is None else '# {}'.format(vname)
            
            ind = '{: <4}'*len(k)
            line = ('    {} {{{}}} {{}}\n'.format(ind,fmt)).format
            args = list(k)+[v, comment]
            content.append(line(*args))
        string = (self.preamble+'\n'
                  +'BLOCK {}\n'.format(self.name)
                  + ''.join(content) + '\n')
        return string

    def items(self):
        return Matrix.items(self)
    
    def iteritems(self):
        return Matrix.iteritems(self)

    def get_name(self, key, default=''):
        return self._names.get(key, default)
        
    def get_number(self, name, default=-1):
        return self._numbers.get(name.lower(),default)
    
    def new_entry(self, key, value, name=None):
        '''
        Add a new named block entry. Ensures that the key hasn't already 
        been set and the name given hasn't already been booked.
        '''
           
        if key in self:
            err = ("Key '{}' already belongs to NamedBlock ".format(key.lower())
                 +  "'{}', mapped to name '{}'.".format(self.name,
                                                        self.get_name(key)))
            raise KeyError(err)
            
        if name in self._numbers:
            err = ("Name '{}' already booked in ".format(name) +
                   "NamedBlock {}, mapped to key '{}'.".format(
                                   self.name, self._numbers[name]))
            raise KeyError(err)
        else:
            if name is not None:
                self._names[key] = name.lower()
                self._numbers[name.lower()] = key
            try:
                self[key] = value
            except KeyError:
                subkey = key[:-2]
                penult = key[-2]
                submatrix = self
                for i,sk in enumerate(subkey):
                    if sk not in submatrix:
                        submatrix[sk] = NamedMatrix()
                    submatrix = submatrix[sk]
                submatrix[penult] = NamedBlock()
                self[key] = value
                
            
    def namedict(self):
        '''Return python dict of name:key pairs.'''
        return {self._names[k]:v for k,v in self.iteritems() 
                                 if k in self._names }


class Decay(OrderedDict):
    '''
    Container class for SLHA Decay blocks. A subclass of `collections.OrderedDict` 
    with a restriction on tuple keys and float values less than 1. A modified 
    __repr__() function is implemented for easy writing to file. Decay is 
    initialised with a PID argument to specify the particle to which it refers 
    as well as its total width. The sum of branching ratios is kept and a 
    ValueError will be raised if the total exceeds 1. It can optionally be 
    named using the 'name' keyword argument. Finally a __str__() function is 
    also defined to output the SLHA formatted block.
    '''
    
    def __checkkey__(self, key):
        '''    Forces key to be a tuple and casts the elements to integers.'''
        if type(key) is not tuple: 
            raise TypeError( self.__class__.__name__ +
                             ' only accepts tuple keys: (PID1, PID2,...).' ) 
        else:
            return tuple(map(int,key))
    
    def __checkval__(self, val):
        '''    Forces values (i.e. branching ratios) to be a float less than 1.'''
        try:
            fval = float(val)
        except ValueError:
            raise TypeError( self.__class__.__name__ +
                         ' only accepts floats or values castable via float()' )
        if fval > 1.:
            raise TypeError("SLHA Decay object for PID = {}. ".format(self.PID)+
                            "Branching ratio > 1 encountered : {}".format(fval))
        return fval
    
    def __init__(self, PID, total, 
                 data=None, comment='', decimal=5, preamble= ''):
        '''    Positional arguments:
            PID - integer to denote particle whose decay is being described
            total - total width
            
            Keyword arguments:
            data - A dict type object supporting iteritems() 
                   to initialize Decay values.
            comment - prints a comment to the right of the block 
                      declaration in __str__().
            decimal - Number of decimal points with which 
                      to write out width and BRs.
            preamble - Some text to print before the __str__ output.'''
            
        super(Decay, self).__init__()
        try:
            self.PID=int(PID)
        except ValueError:
            err = ("SLHA Decay object for PID = {}. ".format(PID) + 
                   "'PID' argument must be an integer or be castable via int()")
            raise TypeError(err)
            
        try:
            self.total=float(total)
        except ValueError:
            err = ("SLHA Decay object for PID = {}. ".format(PID) + 
                  "'total' argument must be a float or be castable via float()")
            raise TypeError(err)

        self._BRtot = 0.
        self._fmt = ':+.{}e'.format(decimal)
        self._comment = comment
        self._decimal = decimal
        self._comments = {}
        self.preamble = preamble
        
        if data is not None:
            for k,v in data.iteritems():
                self[k]=v
                
    def __setitem__(self, key, value):
        super(Decay,self).__setitem__(self.__checkkey__(key),
                                      self.__checkval__(value))
        self._BRtot+=self.__checkval__(value)
        if self._BRtot > 1.:
            print '!!! ERROR !!!'
            print self
            raise ValueError("SLHA Decay object for PID = {}. ".format(self.PID)
                            + "Sum of branching ratios > 1!")            
    
    def __delitem__(self, key):
        super(Decay,self).__delitem__(self.__checkkey__( key ),
                                      self.__checkval__(value))
        self._BRtot-=self.__checkval__(self[key])
        
    def __repr__(self):
        return ( '<SHLA Decay: PID={}; {} entries.>'.format(self.PID, 
                                                            len(self)) )
    
    def __str__(self):
        above = '#        PDG        Width\n'
        title = ('DECAY    {{:<10}} {{{}}} # {{}}\n'.format(self._fmt)).format
        below = '#    BR{}    NDA       ID1       ID2...\n'.format(' '*self._decimal)
        if len(self)==0: below=''
        content = []
        for k,v in self.iteritems():
            nparts = len(k)
            idfmt = nparts*'{: <10}'
            line = ('{{{}}}    {{: <10}}{}'.format(self._fmt,idfmt)).format

            cmnt = ('# {}'.format(self._comments[k]) 
                     if k in self._comments else '')
            
            content.append(line(v, nparts, *k) + cmnt)

        string = (self.preamble + '\n' + above 
                 + title(self.PID, self.total, self._comment) 
                 + below + '\n'.join(content) )
        return string
    
    def new_channel(self, PIDs, BR, comment=None):
        '''
        Add a new decay channel. 
        PIDs - tuple of integers.
        BR - the branching ratio into that channel.
        comment - optional comment to be written to the 
                  right of that channel in __str__().
        '''
        self[PIDs]=BR
        if comment is not None: self._comments[PIDs]=comment

class Card(object):
    '''
    SLHA card object: a container for storing multiple SLHA.Block, 
    SLHA.NamedBlock and SLHA.Decay objects. Blocks and decays are stored in 
    OrderedDicts self.blocks and self.decays respectively. 
    Index syntax mycard[key] can be used to get or set elements of blocks or 
    decays possessed by the Card instance. If passed a string key, the card 
    instance will look up the first block in self.blocks possessing a parameter 
    with the given name while an int key will return the SLHA.Decay object for 
    that PID.
    '''
    def _parent_block(self, key):
        '''
        Finds the parent block of a given string key. Raises Key error 
        if the key is not booked in any of the card's blocks.
        '''
        if key in self.matrices: 
            return self.matrices[key]
        for block in self.matrices.values():
            if key in block:
                return block
        for block in self.blocks.values():
            if key in block:
                return block
        err = 'Key/Matrix "{}" not found in {}'.format(key, self)
        raise KeyError(err)
    
    def _container(self, key):
        '''
        Returns the container containing the indexed key. Used in overloading 
        __setitem__, __getitem__, __contains__ and __delitem__.
        '''
        if type(key) is int:
            return self.decays
        elif type(key) is str:
            return self._parent_block(key)
        else:
            err = ('SLHA Card has integer keys for '
                   'DECAY and string keys for BLOCK.')
            raise ValueError(err)
        
    def __init__(self, blocks=None, decays=None, name=None):
        '''    Keyword arguments:
        blocks - a dictionary of items with which to initialise 
                 individual blocks in the card instance i.e. 
                 name:block pairs with block an SLHA.NamedBlock.
        decays - a dictionary of SLHA.Decay objects.
        matrices - a dictionary of SLHA.NamedMatrix objects.
        name - card name'''
        
        self.blocks = OrderedDict()
        self.decays = OrderedDict()
        self.matrices = OrderedDict()
        self.name = name if name is not None else ''
        if blocks is not None:
            for bname, block in blocks.iteritems():
                self.blocks[bname] = NamedBlock(name=bname, data=block)  
        if decays is not None:
            for PID, decay in decays.iteritems():
                self.decays[PID]= Decay(PID,decay.total, data=decay)
            
                
    def __repr__(self):
        return ('<SHLA Card "{}": {} blocks, {} decays.>'.format(self.name,
                                        len(self.blocks), len(self.decays)))
    def __contains__(self, key):
        return key in self._container(key)
    
    def __getitem__(self, key):
        return self._container(key)[key]
    
    def __delitem__(self, key):
        del self._container(key)[key]
    
    def __setitem__(self, key, value):
        return self._container(key).__setitem__(key,value)
        
    def add_block(self, block, preamble=''):
        '''Append an SLHA.Block or NamedBlock to self.blocks.'''
        block.preamble = preamble
        if isinstance(block,Matrix):
            self.matrices[block.name] = block
        else:
            self.blocks[block.name] = block
        
    def add_decay(self, decay, preamble=''):
        '''Append an SLHA.Decay to self.decays.'''
        decay.preamble = preamble
        self.decays[decay.PID] = decay
    
    def add_entry(self, blockname, key, value, name=None):
        '''
        Add a new field in a given block. If the Card instance already  
        has such a block, a new field is appended. If the Card doesn't, a new 
        block is created.
        '''
        if self.has_block(blockname):
            self.blocks[blockname].new_entry(key, value, name=name)
        elif self.has_matrix(blockname):
            self.matrices[blockname].new_entry(key, value, name=name)
        else:
            if type(key) is tuple:
                container = Matrix
                namedcontainer = NamedMatrix
            else:
                container = Block
                namedcontainer = NamedBlock
                
            if name is not None:
                theblock = namedcontainer(name=blockname)
                theblock.new_entry(key, value, name=name)
            else:
                theblock = container(name=blockname)
                theblock.new_entry(key, value)
                
            self.add_block(theblock)
                
    def add_channel(self, PID, channel, partial, comment=None):
        '''
        Add a new channel in a given decay. The Decay for the given PID 
        must already exist in the Card instance.
        '''
        assert(self.has_decay(PID)),('Tried adding a decay channel '
                                     'for non existent decay PID')
        self.decays[PID].new_channel(channel, partial, comment=comment)
    
    def new_channel(self, PIDs, BR, comment=None):
        '''Append a new channel to a given Decay block.'''
        self[PIDs]=BR
        if comment is not None: self._comments[PIDs]=comment
        
    def has_block(self, name):
        return name in self.blocks 
    
    def has_matrix(self, name):
        return name in self.matrices
    
    def has_decay(self, PID):
        return PID in self.decays
    
    def write(self, filename, blockorder = [], preamble='',postamble=''):
        '''
        Write contents of Card in SLHA formatted style to "filename". 
        Makes use of the __str__ methods written for the SLHA elements. 
        Keyword arguments:
            blockorder - Specify an ordering for printing the block names. 
                         Names given in blockorder will be printed first in the 
                         order specified and others will be printed as ordered 
                         in self.blocks
            preamble - Some text to write before the Block and Decay information
            postamble - Some text to write after the Block and Decay information
        '''
        
        with open(filename,'w') as out:
            out.write(preamble)
            
            allblocks = self.blocks.keys() + self.matrices.keys()
            other_blocks = [b for b in allblocks if b not in blockorder]
            
            for block in blockorder:
                if self.has_block(block):
                    out.write(str(self.blocks[block]))
                elif self.has_matrix(block):
                    out.write(str(self.matrices[block]))
                    
            for block in other_blocks:
                if self.has_block(block):
                    out.write(str(self.blocks[block]))
                elif self.has_matrix(block):
                    out.write(str(self.matrices[block]))
                                
            for decay in self.decays.values():
                out.write(str(decay))
            out.write(postamble)

class SLHAReadError(Exception): # Custom error name for reader
    pass

def read(card):
    '''
    SLHA formatted card reader. Blocks and Decay structures are read 
    into an SLHA.Card object which is returned. Comments are specified by the # 
    character. Comments to the right of the structures as wel as the individual 
    fields are stored.
    '''
    
    def get_comment(line):
        '''
        Returns characters following a "#" in line.
        returns empty string if match is not found.
        '''
        match = re.match(r'.*#\s*(.*)\s*', line)
        if match: return match.group(1).strip()
        else: return ''

    thecard = Card()
    pcard = open(card,'r')
    lines = iter(pcard)
    counter = 0
    
    try:
        while True:
            counter+=1 # keep track of line number
            stop=False
            try: ll=(last_line.strip()).lower()
            except NameError: ll = (lines.next()).strip().lower()
        
            # ll = (line.strip()).lower()
            if not ll: continue

            first_chars = re.match(r'\s*(\S+).*',ll).group(1)
            
            if first_chars=='block':
                try:
                    block_details = re.match(r'\s*block\s+([^\n]+)',ll,
                                             re.IGNORECASE).group(1)
                except AttributeError:
                    err = ('Invalid block format encountered ' + 
                          'in line {} of {}.'.format(counter, card) )
                    raise SLHAReadError(err)
                    
                bname = block_details.split('#')[0].strip()
                
                comment = get_comment(ll)
                
                block_data, last_line, stop = read_until(lines,'block','decay')
                
                elements = []
                
                for datum in block_data:
                    counter +=1                    
                    is_comment = re.match(r'\s*#.*',datum)
                
                    if (not datum.strip() or is_comment): continue
                    
                    info = re.match(r'\s*((?:\d+\s+)+)\s*(\S+).*',
                                    datum)
                    if not info:
                        print datum
                        print ('Ignored datum in block '+
                                '{},'.format(theblock.name) +
                                ' (line {} of {})'.format(counter, card))
                        continue
                    
                    key, value = info.group(1), info.group(2)
                    
                    try:
                        key = int(key)
                    except ValueError:
                        key = tuple( map(int, key.split()) )

                    dname = get_comment(datum)
                    
                    try:
                        entry = float(value)
                    except ValueError:
                        entry = value
                    finally:
                        elements.append((key,value,dname))
                
                
                if type(elements[0][0]) is tuple:
                    theblock = NamedMatrix(name=bname.lower(), comment=comment)
                else:
                    theblock = NamedBlock(name=bname.lower(), comment=comment)
            
                for ele in elements:
                    theblock.new_entry(*ele)
                
                thecard.add_block(theblock)

            elif first_chars=='decay':
                try:
                    decay_details = re.match(r'\s*decay\s+(.+)',ll,
                                             re.IGNORECASE).group(1)
                except AttributeError:
                    err = ('Invalid decay format encountered' + 
                          'in line {} of {}.'.format(counter, card) )
                    raise SHLAReadError(err)
                    
                info = re.match(r'\s*(\d+)\s+(\S+)\s+.*', decay_details)
                PID, total = info.group(1), info.group(2)
                
                comment = get_comment(decay_details)
                    
                thedecay = Decay(PID=int(PID), total=float(total), 
                                 comment=comment)
                
                decay_data, last_line, stop = read_until(lines,'block','decay')
                
                for datum in decay_data:
                    counter +=1                    
                    is_comment = re.match(r'\s*#.*',datum)
                
                    if ((not datum.strip()) or (is_comment)): continue

                    info = re.match(r'\s*(\S+)\s+(\d+)\s+(.+)',datum)
                    if not info:
                        print datum
                        print ('Ignored above datum in decay '+
                                '{},'.format(thedecay.PID) +
                                ' (line {} of {})'.format(counter, card))
                        continue
                        
                    BR, nout = info.group(1), info.group(2)
                    PIDinfo = info.group(3).split('#')[0]
                    
                    PIDs = tuple( map(int, PIDinfo.split()) )
                    
                    if len(PIDs)!=int(nout):
                        print ("Number of external particles in column 2 doesn't "
                               "match number of subsequent columns:")
                        print datum
                        print ('Ignored above datum in decay '+
                                '{},'.format(thedecay.PID) +
                                ' (line {} of {})'.format(counter, card))
                        continue

                    comment = get_comment(datum)
                        
                    thedecay.new_channel(PIDs, float(BR), comment=comment)
                    
                thecard.add_decay(thedecay)
                
            if stop: raise StopIteration
            
    except StopIteration:
        # print 'Finished reading "{}".'.format(card)
        pcard.close()
        # pass
                
    
    return thecard


def read_until(lines, here, *args):
    '''
    Loops through an iterator of strings by calling next() until 
    it reaches a line starting with a particular string.
    Case insensitive.
    Args:
        lines - iterator of strings
        here - string (plus any further argumnts). 
               Reading will end if the line matches any of these.
    Return:
        lines_read - list of lines read
        line - last line that was read (containing string "here")
    '''
    end_strings = [here.lower()]+[a.lower() for a in args]
    lines_read = []
    line = ''
    stopiter = False
    while not any([line.strip().lower().startswith(x) for x in end_strings]):
        try:
            line = lines.next()
            lines_read.append(line.strip('\n'))
        except StopIteration:
            stopiter=True
            break
    try:
        if stopiter:
            return lines_read, '', stopiter
        else:
            return lines_read[:-1],lines_read[-1], stopiter
    except IndexError:
        return [],'',stopiter


if __name__=='__main__':
    # from itertools import product
    # mm = NamedMatrix(name='wehey')
    # for i,j,k in product((1,2,3),(1,2,3),(1,2,3)):
    #     mm.new_entry((i,j,k),i+j+k, name = '{}{}{}'.format(i,j,k))
    #     # print '{} {} {}'.format(i,j,k),m.dict()[i][j][k]
    #
    # print mm
    # a = Block( name='a', data={1:'1',2:'seven'})
    # b = Block( name='b', data={1:'1',2:'2'}, dtype=float)
    # c = Block( name='c', data={1:'1',2:'2'}, dtype=float)
    # d = Block( name='d', data={1:'1',2:'2'}, dtype=float)
    # e = Block( name='e', data={1:'1',2:'2'}, dtype=float)
    # f = Block( name='f', data={1:'1',2:'2'}, dtype=float)
    # g = Block( name='g', data={1:'1',2:'2'}, dtype=float)
    # h = Block( name='h', data={1:'1',2:'2'}, dtype=float)
    # i = Block( name='i', data={1:'1',2:'2'}, dtype=float)
    #
    # aa = Matrix( name='aa', data={1:a,2:b,3:c})
    # dd = Matrix( name='dd', data={1:d,2:e,3:f})
    # gg = Matrix( name='gg', data={1:g,2:h,3:i})
    #
    # m = Matrix( name='m', data={1:aa,2:dd,3:gg})
    #
    # print m
    # for k,v in m.iteritems():
        # print k,v
    # for i,j,k in product((1,2,3),(1,2,3),(1,2,3)):
    #     print '{} {} {}'.format(i,j,k),m.dict()[i][j][k]
    card = read('testmatrix.dat')
    
    print card.matrices['wehey']
    print card.blocks['wehey2']