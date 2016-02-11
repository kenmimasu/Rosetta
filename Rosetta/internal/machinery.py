# Imports all modules in the Rosetta base directory
from .. import __all__ as basisnames
from .. import *
# from .. import HiggsBasis
from errors import TranslationPathError
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
        err = 'No translation path available for "{}"'.format(end)
        raise TranslationPathError(err)
        
        # return []

    step = path[-1][0]

    while step!=end:
        try:
            intermediate = rels[step][end]
            path += intermediate
            step = intermediate[-1][0]
        except KeyError:
            err = ('\tProblem finding translation step '+
                   'between "{}" & "{}"'.format(step,end))
            raise TranslationPathError(err)
            # return []
    return path
    
################################################################################
# Constructs the relationship dictionary from all .py files in Rosetta/ dir

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

# Build dictionary of all basis classes and their implemented translation 
# functions
translations = {}
for basis, instance in bases.iteritems():
    functions = [i for i in instance.__dict__.values() 
                 if hasattr(i,'_target')]
    tmap = {f._target:f for f in functions}
    translations[basis] = tmap

# Feed dictionary to relate funciton to build all possible paths between bases
relationships = relate(translations)

# for k,v in  relationships.items():
#     print k
#     for kk, vv in v.items():
#         print '\t',kk,vv
# import sys
# sys.exit()
################################################################################
if __name__=='__main__':
    pass
    # mydict = {
    #          'A':{'B':'A.Bfunc','C':'A.Cfunc'},
    #          'B':{'C':'B.Cfunc','D':'B.Dfunc'},
    #          'C':{'A':'C.Afunc','D':'C.Dfunc','B':'CBfunc'},
    #          'D':{'A':'D.Afunc'},
    #          'E':{'D':'E.Dfunc'},
    #          'F':{'C':'F.Cfunc','E':'F.Efunc'},
    #          'G':{'E':'G.Efunc'}
    #          }
    #
    # relations =  relate(mydict)
    # for k,v in relations.items():
    #     print k,v
    # print '\n\n'
    # print get_path('A','C',relations)