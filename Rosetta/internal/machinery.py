# Imports all modules in the Rosetta base directory
from ..bases import __all__ as basisnames
from ..bases import *
from errors import TranslationPathError, BasesError, RelationshipsError
import inspect
################################################################################
__doc__ = '''
Rosetta's machinery for generating the translation paths available for the set 
of Basis implementations contained in its root directory. relate() builds
the possible translation paths from one basis to all of the others by looking 
for translation functions belonging to the class tagged by the translation 
decorator. The relationships dictionary contains all of this information and 
is used in the translate() function of the class Basis.
'''

def relate(mydict):
    '''
    Recursively build the relations tree for translation functions between 
    bases. 
    '''
    def traverse(base, state={}, chain=[]):
        
        subset = [(x,f) for x,f in mydict[base].iteritems() 
                  if (x!=k and x not in state)]
                  
        for ele in subset: 
            state[ele[0]] = chain

        for ele in subset:
            subchain = chain[:]
            subchain.append(ele)
            return traverse(ele[0], chain=subchain, state=state)
        
        return state
    
    relations = {}
               
    for k in mydict.keys():
        relations[k] = traverse(k, state={}, chain=[])
        
    for k, v in mydict.iteritems():
        for kk, vv in v.iteritems():
            relations[k][kk] = [(kk,vv)]
            
    return relations

def get_path(start, end, rels):
    '''traverse relations tree and return the path needed'''

    avail = rels[start]
    try:
        path = avail[end][:]
    except KeyError:
        err = 'No translation path available to "{}"'.format(end)
        raise TranslationPathError(err)
        
        # return []

    step = path[-1][0]

    while step!=end:
        try:
            intermediate = rels[step][end]
            path += intermediate
            step = intermediate[-1][0]
        except KeyError:
            err = ('\tFailed to find translation step '+
                   'between "{}" & "{}"'.format(step,end))
            raise TranslationPathError(err)
            # return []
    return path
    
################################################################################
# Constructs the relationship dictionary from all .py files in Rosetta/bases dir

# Collect basis classes in base directory of Rosetta according to their name
modules = {b:v for b,v in globals().iteritems() if b in basisnames}
bases = {}
for bname, module in modules.iteritems():
    try:
        bclass = module.__dict__[bname]
        if not inspect.isclass(bclass):
            raise KeyError
        bases[bclass.name] = bclass
    except KeyError:
        print ('Warning: Rosetta did not find a class named ' +
               '{0} in {0}.py. File ignored.'.format(bname) )

if not bases: raise BasesError('No valid basis implementations found.')
# Build dictionary of all basis classes and their implemented translation 
# functions
translations = {}
for basis, instance in bases.iteritems():
    functions = [i for i in instance.__dict__.values() 
                 if hasattr(i,'_target')]
    tmap = {f._target:f for f in functions}
    translations[basis] = tmap

# Feed dictionary to relate function to build all possible paths between bases
relationships = relate(translations)
if not relationships: 
    raise RelationshipsError('No valid translation relationships between '
                             'basis implementations found.')

