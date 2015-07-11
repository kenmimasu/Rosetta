# Imports all modules in the Rosetta bas directory
from Rosetta import __all__ as basisnames
from Rosetta import *
import inspect
################################################################################
def relate(mydict):
    '''build the relations tree for implemented.py'''
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
        
    for k,v in mydict.iteritems():
        for kk,vv in v.iteritems():
            relations[k][kk] = [(kk,vv)]
            
    return relations

def get_path(start,end,rels):
    '''traverse relations tree and return the path needed'''
    avail = rels[start]
    try:
        path = avail[end]
    except KeyError:
        return []
    step = path[-1][0]
    while step!=end:
        intermediate = rels[step][end]
        path += intermediate
        step = intermediate[-1][0]
    return path
################################################################################
# Constructs the relationship dictionary from all .py files in
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

translations = {}

for bas,inst in bases.iteritems():
    functions = [i for i in inst.__dict__.values() if hasattr(i,'_target')]
    tmap = {f._target:f for f in functions}
    translations[bas] = tmap

relationships = relate(translations)
################################################################################
if __name__=='__main__':
    mydict = {
             'A':{'B':'A.Bfunc','C':'A.Cfunc'},
             'B':{'C':'B.Cfunc','D':'B.Dfunc'},
             'C':{'A':'C.Afunc','D':'C.Dfunc','B':'CBfunc'},
             'D':{'A':'D.Afunc'},
             'E':{'D':'E.Dfunc'},
             'F':{'C':'F.Cfunc','E':'F.Efunc'},
             'G':{'E':'G.Efunc'}
             }

    relations =  relate(mydict)
    for k,v in relations.items():
        print k,v
    print '\n\n'
    print get_path('A','C',relations)