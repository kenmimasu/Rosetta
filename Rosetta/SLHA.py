from collections import OrderedDict
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
        
    def __init__(self, name=None, data=None, decimal=5):
        super(Block, self).__init__()
        self.name = name
        self.fmt= ':+.{}e'.format(decimal)
        if data is not None:
            for k,v in data.iteritems():
                self[k]=v
                
    def __setitem__(self,key,value):
        super(Block,self).__setitem__(self.__check__(key),value)
    
    def __repr__(self):
        return ( '<SHLA Block: "{}"; {} entries.>'.format(self.name,len(self)) )
    
    def __str__(self):
        line = ('    {{: <4}} {{{}}}'.format(self.fmt)).format
        content = '\n'.join([line(k,v) 
                             for k,v in self.iteritems()])
        string = ('BLOCK {}\n'.format(self.name)
                   + content)
        return string
        

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
                       names=None, decimal=5):
        super(NamedBlock, self).__init__(name=name, 
                                         data=data, decimal=decimal)
        if names is not None:
            self._names = names
            try:
                self._numbers = data={v:k for k,v 
                                      in self._names.iteritems()}
            except AttributeError:
                raise ValueError('"names" keyword argument must be a'
                                 'dictionary or support iteritems() method.')
        else:
            self._names, self._numbers = {}, {}

    def __setitem__(self,key, value):
        return super(Block,self).__setitem__(self.__parse__(key),value)

    def __getitem__(self, key):
        return super(Block,self).__getitem__(self.__parse__(key))
    
    def __delitem__(self, key):
        if type(key) is int:
            try:
                del self._numbers[self._names[key]]
                del self._names[key]
            except KeyError:
                pass
        if type(key) is str:
            try:
                del self._names[self._numbers[key]]
                del self._numbers[key]
            except KeyError:
                pass
        return super(Block,self).__delitem__(self.__parse__(key))
                     
    def __contains__(self, key):
        return super(Block,self).__contains__(self.__parse__(key))
        
    def __repr__(self):
        return ('<SHLA NamedBlock: "{}"; {} entries, {} named.>'.format(
                                        self.name,len(self), len(self._names)))
    def name(key):
        return self._names[key]
        
    def number(name):
        return self._numbers[name]
    
    def new_entry(self,key,value,name=None):
        if key in self:
            err = ("Key '{}' already belongs to NamedBlock ".format(key) +
                   "'{}', mapped to name '{}'.".format(self.name,
                                                    self._names[key]) )
            raise KeyError(err)
        elif name in self._numbers:
            err = ("Name '{}' already booked in ".format(name) +
                   "NamedBlock {}, mapped to key '{}'.".format(
                                   self.name,self._numbers[name]) )
            raise KeyError(err)
        else:
            if name is not None:
                self._names[key] = name
                self._numbers[name] = key
            self[key]=value

    def __str__(self):
        line = ('    {{: <4}} {{{}}} # {{}}'.format(self.fmt)).format
        content = '\n'.join([line(k, v, self._names[k]) 
                             for k,v in self.iteritems()])
        string = ('BLOCK {}\n'.format(self.name)
                   + content)
        return string        

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
                             ' only accepts tuple keys: (PID1, PID2).' ) 
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
        
        self.fmt= ':+.{}e'.format(decimal)
        self.comment=comment
        self.decimal=decimal
                
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
        below = '#    BR{}\tNDA\tID1       ID2\n'.format(' '*self.decimal)
        line = ('{{{}}}\t{{}}\t{{: <10}}{{: <10}}'.format(self.fmt)).format
        content = '\n'.join([line(v,len(k),*k) for k,v in self.iteritems()])
        title = ('DECAY\t{{:<10}} {{{}}} # {{}}\n'.format(self.fmt)).format
        string = (above + title(self.PID, self.total, self.comment) 
                + below + content)
        return string
    
if __name__=='__main__':
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