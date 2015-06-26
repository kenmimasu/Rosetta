from collections import OrderedDict
import sys
import re

# import warnings
#
class Block(OrderedDict):
    '''    Container class for SLHA block. A subclass of `collections.OrderedDict`
    with a restriction on integer keys and a modified __repr__(). 
    Block can't be initialised with any positional arguments but rather 
    with the 'data' keyword argument. 
    It can optionally be named using the 'name' keyword argument. 
    Finally a _print() function is also defined to output the SLHA 
    formatted block.'''
        
    def __check__(self,key):
        if type(key) is not int: 
            raise TypeError( self.__class__.__name__ +
                             ' only accepts integer keys.' ) 
        else:
            return key
        
    def __init__(self, name=None, data=None, decimal=5, dtype=lambda x: x):
        super(Block, self).__init__()
        self.name = name
        self.cast = dtype
        self.fmt= ':+.{}e'.format(decimal)
        if data is not None:
            for k,v in data.iteritems():
                self[k]=dtype(v)
                
    def __setitem__(self,key,value):
        super(Block,self).__setitem__(self.__check__(key),self.cast(value))
    
    def __repr__(self):
        return ( '<SHLA Block: "{}"; {} entries.>'.format(self.name,len(self)) )
    
    def __str__(self):
        line = ('    {{: <4}} {{{}}}'.format(self.fmt)).format
        content = '\n'.join([line(k,v) 
                             for k,v in self.iteritems()])
        string = ('BLOCK {}\n'.format(self.name)
                   + content+ '\n\n')
        return string
    
    def dict(self):
        return {k:v for k,v in self.iteritems()}

class NamedBlock(Block):
    '''    Class derived from 'Block' with the added functionality of assigning a 
    name to each key via the 'names' keyword argument. The values can then 
    also be referenced by name via a lookup in self._numbers. '''
    
    def __parse__(self, key):
        if type(key) is str:
            try:
                return self._numbers[key]
            except KeyError:
                err = ('Name "{}" has not been assigned '.format(key) +
                       'to a key in {}.names'.format(self.name))
                raise KeyError(err)
        else:
            return key
    
    def __init__(self, name=None, data=None, 
                 comment='',decimal=5, dtype=lambda x:x):
        super(NamedBlock, self).__init__(name=name, data=data, 
                                         decimal=decimal, dtype=dtype)
        if data is not None:
            self._names = data._names
            try:
                self._numbers = {v:k for k,v 
                                      in self._names.iteritems()}
            except AttributeError:
                raise ValueError('"names" keyword argument must be a'
                                 'dictionary or support iteritems() method.')
        else:
            self._names, self._numbers = {}, {}
        self.cast=dtype
        self.comment=comment

    def __setitem__(self,key, value):
        return super(Block,self).__setitem__(self.__parse__(key),
                                             self.cast(value))

    def __getitem__(self, key):
        return super(Block,self).__getitem__(self.__parse__(key))
    
    def __delitem__(self, key):
        if type(key) is int:
            try:
                del self._numbers[self._names[key]]
                del self._names[key]
            except KeyError:
                pass
        elif type(key) is str:
            key = key
            try:
                del self._names[self._numbers[key]]
                del self._numbers[key]
            except KeyError:
                pass
        return super(Block,self).__delitem__(self.__parse__(key))
                     
    def __contains__(self, key):
        if type(key) is str:
            return key in self._numbers 
        elif type(key) is int:
            return super(Block,self).__contains__(key)
        
    def __repr__(self):
        return ('<SHLA NamedBlock: "{}"; {} entries, {} named.>'.format(
                                        self.name,len(self), len(self._names)))
    def __str__(self):
        content = []
        for k,v in self.iteritems():
            try:
                strval = ('{{{}}}'.format(self.fmt)).format(float(v))
            except ValueError:
                strval = v
            if k in self._names:
                content.append('    {: <4} {} # {}'.format(k,strval,
                                                           self._names[k]))
            else:
                content.append('    {: <4} {}'.format(k,strval))
        string = ('BLOCK {} # {}\n'.format(self.name,self.comment)
                   + '\n'.join(content) + '\n\n')       
        return string

    def get_name(self, key, default=''):
        return self._names.get(key,default)
        
    def get_number(self, name):
        return self._numbers[name]
    
    def new_entry(self,key,value,name=None):
        if key in self:
            err = ("Key '{}' already belongs to NamedBlock ".format(key) +
                   "'{}', mapped to name '{}'.".format(self.name,
                                                    self._names[key]))
            raise KeyError(err)
        elif name in self._numbers:
            err = ("Name '{}' already booked in ".format(name) +
                   "NamedBlock {}, mapped to key '{}'.".format(
                                   self.name,self._numbers[name]))
            raise KeyError(err)
        else:
            if name is not None:
                self._names[key] = name
                self._numbers[name] = key
            self[key]=value
            
    def namedict(self):
        return {self._names.get(k,''):v for k,v in self.iteritems()}
        

class Decay(OrderedDict):
    '''Container class for SLHA Decay. A subclass of `collections.OrderedDict`
    with a restriction on tuple keys and float values less than 1. A modified 
    __repr__() funciton is implemented for easy writing to file. 
    Decay is initialised with a PID argument to specify the particle to which 
    it refers as well as its total width. 
    It can optionally be named using the 'name' keyword argument. 
    Finally a _print() function is also defined to output the SLHA 
    formatted block.'''
    
    def __check__(self,key):
        if type(key) is not tuple: 
            raise TypeError( self.__class__.__name__ +
                             ' only accepts tuple keys: (PID1, PID2,...).' ) 
        else:
            return tuple(map(int,key))
    
    def __checkval__(self,val):
        try:
            fval = float(val)
        except ValueError:
            raise TypeError( self.__class__.__name__ +
                         ' only accepts floats or values castable via float()' )
        if fval > 1.:
            raise TypeError("SLHA Decay object for PID = {}. ".format(self.PID)+
                            "Branching ratio > 1 encountered : {}".format(fval))
        return fval
    
    def __init__(self, PID, total, data=None, comment="", decimal=5):
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

        if data is not None:
            for k,v in data.iteritems():
                self[k]=v
        
        self._BRtot = 0.
        
        self._fmt= ':+.{}e'.format(decimal)
        self._comment=comment
        self._decimal=decimal
        self._comments={}
                
    def __setitem__(self,key,value):
        super(Decay,self).__setitem__(self.__check__( key ),self.__checkval__(value))
        self._BRtot+=self.__checkval__(value)
        if self._BRtot > 1.:
            print '!!! ERROR !!!'
            print self
            raise TypeError("SLHA Decay object for PID = {}. ".format(self.PID)
                            + "Sum of branching ratios > 1!")            
    
    def __delitem__(self, key):
        self._BRtot-=self.__checkval__(self[key])
        super(Decay,self).__delitem__(self.__check__( key ),self.__checkval__(value))
        
    def __repr__(self):
        return ( '<SHLA Decay: PID={}; {} entries.>'.format(self.PID,len(self)) )
    
    def __str__(self):
        above = '#\tPDG\t\tWidth\n'
        title = ('DECAY\t{{:<10}} {{{}}} # {{}}\n'.format(self._fmt)).format
        below = '#    BR{}\tNDA\tID1       ID2...\n'.format(' '*self._decimal)
        if len(self)==0: below=''
        content = []
        for k,v in self.iteritems():
            nparts = len(k)
            idfmt = nparts*'{: <10}'
            line = ('{{{}}}\t{{}}\t{}'.format(self._fmt,idfmt)).format

            cmnt = ('# {}'.format(self._comments[k]) 
                     if k in self._comments else '')
            
            content.append(line(v, nparts, *k) + cmnt)

        string = (above + title(self.PID, self.total, self._comment) 
                + below + '\n'.join(content) + '\n')
        return string
    
    def new_channel(self,PIDs,partial,comment=None):
        self[PIDs]=partial
        if comment is not None: self._comments[PIDs]=comment
        
class Card(object):
    
    def _parent_block(self, key):
        for block in self.blocks.values():
            if key in block:
                return block
        return OrderedDict()
    
    def _container(self, key):
        if type(key) is int:
            return self.decays
        elif type(key) is str:
            return self._parent_block(key)
        else:
            err = ('SLHA Card has integer keys for '
                   'DECAY and string keys for BLOCK.')
            raise ValueError(err)
        
    def __init__(self, blocks=None, decays=None):
        self.blocks=OrderedDict()
        self.decays=OrderedDict()
        if blocks is not None:
            for bname, block in blocks.iteritems():
                self.blocks[bname] = NamedBlock(name=bname, data=block)  
        if decays is not None:
            for PID, decay in decays.iteritems():
                self.decays[PID]= Decay(PID,decay.total, data=decay)
                
    
    def __repr__(self):
        return ('<SHLA Card: {} blocks, {} decays.>'.format(
                                        len(self.blocks), len(self.decays)))
    def __contains__(self,key):
        return key in self._container(key)
    
    def __getitem__(self, key):
        return self._container(key)[key]
    
    def __delitem__(self,key):
        del self._container(key)[key]
    
    def __setitem__(self,key,value):
        self._container(key)[key]=value
        
            
    def add_block(self, block):
        self.blocks[block.name] = block
        
    def add_decay(self, decay):
        self.decays[decay.PID] = decay
    
    def add_entry(self,block,key,value,name=None):
        if self.has_block(block):
            self.blocks[block].new_entry(key,value,name=name)
        else:
            if name is not None:
                theblock = NamedBlock(name=block)
                theblock.new_entry(key,value,name=name)
            else:
                theblock = Block(name=block)
                theblock.new_entry(key,value)
                
            self.add_block(theblock)
                
    def add_channel(self,PID,channel,partial,comment=None):
        assert(self.has_decay(PID)),('Tried adding a decay channel '
                                     'for non existent decay PID')
        self.decays[PID].new_channel(channel,partial,comment=comment)
    
    def new_channel(self,PIDs,partial,comment=None):
        self[PIDs]=partial
        if comment is not None: self._comments[PIDs]=comment
        
    def has_block(self, name):
        return name in self.blocks
    
    def has_decay(self, PID):
        return PID in self.decay
    
    def write(self,filename):
        with open(filename,'w') as out:
            for block in self.blocks.values():
                out.write(str(block))
            for decay in self.decays.values():
                out.write(str(decay))
            
    

class SLHAReadError(Exception):
    pass

def read(card):
    thecard = Card()
    pcard = open(card,'r')
    lines = iter(pcard)
    counter = 0
    try:
        while True:
            counter+=1
        
            try: line=(last_line.strip()).lower()
            except NameError: line = lines.next()
        
            ll = (line.strip()).lower()
            if not ll: continue

            first_chars = re.match(r'\s*(\S+).*',ll).group(1)
            if first_chars=='block':
                try:
                    block_details = re.match(r'\s*\S+\s+(.+)',ll).group(1)
                except AttributeError:
                    err = ('Invalid block format encountered' + 
                          'in line {} of {}.'.format(counter, card) )
                    raise SHLAReadError(err)
                try:
                    bname, comment = map(str.strip, block_details.split('#'))
                except ValueError:
                    bname, comment = block_details.strip(), ''
            
                theblock = NamedBlock(name=bname.lower(), comment=comment)
            
                block_data, last_line = read_until(lines,'block','decay')

                for datum in block_data:
                    counter +=1                    
                    is_comment = re.match(r'\s*#.*',datum)
                
                    if (not datum.strip() or is_comment): continue
                
                    try:
                        info = re.match(r'\s*(\S+)\s+(\S+)\s*#\s*(\S+).*',
                                        datum)
                        index, value, dname = tuple([info.group(x) 
                                                    for x in xrange(1,4)])
                    except AttributeError:
                        info = re.match(r'\s*(\d+)\s+(\S+).*',
                                        datum)
                        if not info:
                            print datum
                            print ('Ignored datum in block '+
                                    '{},'.format(theblock.name) +
                                    ' (line {} of {})'.format(counter, card))
                            continue
                        index, value = info.group(1), info.group(2)  
                        dname=None
                    try:
                        entry = float(value)
                    except ValueError:
                        entry = value
                    finally:
                        theblock.new_entry(int(index),entry,name=dname)
                        
                
                thecard.add_block(theblock)
                
            elif first_chars=='decay':
                try:
                    decay_details = re.match(r'\s*\S+\s+(.+)',ll).group(1)
                except AttributeError:
                    err = ('Invalid decay format encountered' + 
                          'in line {} of {}.'.format(counter, card) )
                    raise SHLAReadError(err)
                try:
                    info = re.match(r'\s*(\d+)\s+(\S+)\s*#\s*(\S+).*',
                                    decay_details)
                    PID, total, comment = tuple([info.group(x) 
                                                    for x in xrange(1,4)])
                except AttributeError:
                    info = re.match(r'\s*(\d+)\s+(\S+).*',decay_details)
                    PID, total = info.group(1), info.group(2)
                    comment=''

                thedecay = Decay(PID=int(PID), total=float(total), 
                                 comment=comment)
                
                decay_data, last_line = read_until(lines,'block','decay')

                for datum in decay_data:
                    counter +=1                    
                    is_comment = re.match(r'\s*#.*',datum)
                
                    if (not datum.strip() or is_comment): continue
                
                    try:
                        info = re.match( (r'\s*(\S+)\s+(\d+)\s+([\d+-]+)'+
                                          r'\s+([\d+-]+)\s*#\s*(.*)'), datum)
                        partial, nout = info.group(1), info.group(2)  
                        ID1, ID2 = info.group(3), info.group(4) 
                        comment = info.group(5)
                    except AttributeError:
                        info = re.match( (r'\s*(\S+)\s+(\d+)\s+([\d+-]+)'+
                                          r'\s+([\d+-]+).*'), datum)
                        partial, nout = info.group(1), info.group(2)  
                        ID1, ID2 = info.group(3), info.group(4) 
                        comment=''
                        if not info:
                            print datum
                            print ('Ignored above datum in decay '+
                                    '{},'.format(thedecay.PID) +
                                    ' (line {} of {})'.format(counter, card))
                            continue
                        index, value = info.group(1), info.group(2)  
                        dname=None
                        
                    thedecay.new_channel( (int(ID1),int(ID2)), float(partial) )
                    
                thecard.add_decay(thedecay)
                
    except StopIteration:
        print 'Finished reading.'
        print 
                
    pcard.close()
    
    return thecard
    
def read_until(lines, here, *args):
    '''Loops through an iterator of strings by calling next() until 
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
    while not any([line.strip().lower().startswith(x) 
                   for x in end_strings]): 
        line = lines.next()
        lines_read.append(line.strip('\n'))
    return lines_read[:-1],lines_read[-1]

if __name__=='__main__':
    mycard = read('param_card_test.dat')
    print mycard
    sys.exit()
    _input = OrderedDict([(0,0.),(1,1.),(2,2.)])
    _entries = {0:'zero', 1:'one', 2:'two'}
    myblock = NamedBlock(data=_input, name='test', names=_entries, decimal=3)
    print myblock[1], myblock['one']
    del myblock[2]
    myblock.new_entry(3,3.,'three')
    print myblock['three']
    myblock['three']=3.3
    # myblock.new_entry(700,-7,'seven')
    # myblock.new_entry(700,-7,'argh')
    # print 'seven' in myblock
    mydecay = Decay(7,100., comment='bums')
    mydecay[ (1,2) ]=0.78
    mydecay[ (1,2,40006) ]= 0.3
    # del mydecay[(1,2)]
    print mydecay
    
# class Decay(object):
#     def __init__(self, name, PID, total, data)