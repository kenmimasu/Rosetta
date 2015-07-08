################################################################################
def relate(mydict):
    '''build the relations tree for implemented.py'''
    def traverse(base, state={}, chain=[]):
        
        subset = [(x,f) for x,f in mydict[base].iteritems() 
                  if (x!=k and x not in state)]
                  
        for ele in subset: 
            # state[ele[0]] = [(base,ele[1])]
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
    
    