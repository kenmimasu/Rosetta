# Imports all modules in the Rosetta base directory
from operator import itemgetter
import inspect

from ..bases import __all__ as basisnames
from ..bases import *
from errors import TranslationPathError, BasesError, RelationshipsError
################################################################################
__doc__ = '''
Rosetta's machinery for generating the translation paths available for the set 
of Basis implementations contained in its root directory. relate() builds
the possible translation paths from one basis to all of the others by looking 
for translation functions belonging to the class tagged by the translation 
decorator. The relationships dictionary contains all of this information and 
is used in the translate() function of the class Basis.
'''
################################################################################
def djik(graph, source):
    '''
    Variant of Djikstra's algorithm to find the shortest path on a graph of 
    all nodes from a specified source node. Here all paths have length one 
    in order to find the path with the least number of steps. 
    Returns a dictionary of pointers to the next node on the shortest 
    path to the source node.
    https://en.wikipedia.org/wiki/Dijkstra%27s_algorithm
    '''
    nodes = []
    dist, prev = {}, {}
            
    for v in graph.keys():
        nodes.append(v)  # Collect all nodes on graph
        dist[v] = float('inf') # Set all distances to infinity
        prev[v] = None # Previous node in shortest path from source
    
    # distance to self is zero
    dist[source] = 0 
    
    while nodes:
        # sort nodes by distance to source
        avail = [(k,v) for k,v in dist.items() if k in nodes]
        
        # find closest node to source and remove it from nodes
        u = sorted(avail, key=itemgetter(1))[0][0]
        nodes.remove(u)
        
        # iterate over neighbours of u
        for v in graph[u].keys():
            nsteps = dist[u] + 1
            if nsteps < dist[v]:
                dist[v] = nsteps
                prev[v] = u

    return prev

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
# functions, this graph will be fed into djik()
translations = {}
for basis, instance in bases.iteritems():
    functions = [i for i in instance.__dict__.values() 
                 if hasattr(i,'_target')]
    tmap = {f._target:f for f in functions if f._target in bases}
    translations[basis] = tmap

relationships = {k:djik(translations, k) for k in translations}

def get_path(start, end, rels=relationships, graph=translations):
    djik = rels[start]
    
    u = end
    path = []

    while djik[u] is not None:
        path = [u] + path
        u = djik[u]
    path = [u] + path
    if not path:
        err = 'No translation path available from {} to "{}"'.format(start, end)
        raise TranslationPathError(err)
    return [(y,graph[x][y]) for x,y in zip(path,path[1:])]
    