from .. import session
import checkers as check
from ..errors import TranslationPathError

def translate(basis, target=None, cache=True, avoid_cache = False, 
                      minimal=False): 
    '''
    Translation function. Makes use of existing translations and attempts to 
    find a path between the basis instance and the target basis. If a path 
    exists, successive translations are performed to get from A to B. The 
    resulting basis instance of every translation performed is stored as a 
    value in the basis.translations dictionary with the basis name as a key. If 
    an instance of the target basis is cached, the function returns this 
    instance and prints a message.
    
    Keyword arguments:
        target      - Target basis, must be implemented in Rosetta.bases.
        cache       - Store the result of the translation in the 'translated'
                      attribute of the basis instance to avoid repeating the same
                      translations.
        avoid_cache - Do not use the cached version of translation.
        minimal     - Forego set checks and calculate_dependent() for the target 
                      basis instance. Useful if you want things to be a bit 
                      faster and don't need to propagate the SM input parameters 
                      etc. to the new basis instance
    '''
    from ..machinery import get_path, bases
    # default target
    target = target if target is not None else basis.output_basis
    
    if target == basis.name: return basis
    if target in basis.translations and not avoid_cache:
        msg = 'translation ({} -> {}) already performed.\n'.format
        session.once(msg(basis.name, target))
        return basis.translations[target]
    
    # find the chain of bases and translation functions 
    # to get from basis.name to target
    chain = get_path(basis.name.lower(), target)            

    if not chain: # no translation possible
        inputbasis = basis.__class__.__name__
        outputbasis = bases[target].__name__
        err = ('Sorry, Rosetta cannot translate from ' +
              '{} to {}'.format(inputbasis, outputbasis))
        raise TranslationPathError(err)
    
    names = [basis.name.lower()]+[x[0] for x in chain]
    session.log('Rosetta will be performing the translation:\n'+
                    '    ' + ' -> '.join([bases[x].__name__ 
                                          for x in names]))
    
    # perform succesive translations, checking for  
    # required SM inputs/masses along the way
    current = basis
    
    # required_inputs = current.required_inputs
    # required_masses = current.required_masses

    # ensure required inputs are present
    check.sminputs(current, current.required_inputs)
    check.masses(current, current.required_masses)
    
    for trgt, translate_function in chain:
        # Check if this particular translation path has been saved previously
        if trgt in basis.translations and not avoid_cache:
            new = basis.translations[trgt]
            msg = 'translation ({} -> {}) already performed.\n'.format
            session.once(msg(current.name, trgt))
            required_inputs = set(new.required_inputs)
            required_masses = set(new.required_masses)
            current = new
            
            continue

        instance = bases[trgt](dependent=True, flavor='general')
        
        # update new basis instance with non-EFT blocks, decays
        all_coeffs = (current.blocks.keys() + current.flavored.keys())
        get_other_blocks(instance, current.card, all_coeffs)
        
        message = 'translation ({} -> {})'.format(current.name, 
                                                  instance.name)
        session.verbose(message)
        
        new = translate_function(current, instance)
        
        if (minimal and new.name == target):
            current = new
            break
        
        message = '{}'.format(instance.__class__)
        # checks & calculates dependent parameters
        session.verbose('    Checking required SM inputs '
                        'for "{}"'.format(new.name))
        # new.check_sminputs(new.required_inputs, message=message)
        check.sminputs(new, new.required_inputs, message=message)
        session.verbose('    Checking for required masses '
                        'for "{}"'.format(new.name))
        check.masses(new, new.required_masses, message=message)
        session.verbose('    Checking EFT coefficients '
                        'for "{}"'.format(new.name))
        # NOTE switch off check for presence of dependent coefficients as 
        # there is no difference between the coefficient existing and it 
        # having been assigned a value.
        check.param_data(basis, do_dependent=False)
        
        session.verbose('    Calling ' +
                        '{}.calculate_dependent()\n'.format(new.__class__))
        new.calculate_dependent()

        if cache: basis.translations[trgt] = new
        # prepare for next step
        required_inputs = set(new.required_inputs)
        required_masses = set(new.required_masses)
        current = new
    
    if current.name =='bsmc' and not minimal:
        # expand matrix blocks in the case of bsmc output to respect SLHA 
        # convention
        expand_matrices(current)
    else:
        # reduce flavor structure back to user set option
        current.flavor = basis.flavor
        current.set_flavor('general', basis.flavor)            
        
    session.verbose('\nTranslation successful.\n')
        
    return current


def get_other_blocks(basis, card, ignore):
    ignore = [x.lower() for x in ignore] 
    
    other_blocks, other_matrices = {}, {}
    for k, v in card.blocks.iteritems():
        if k.lower() != 'basis' and k.lower() not in ignore:
            other_blocks[k]=v
    for k, v in card.matrices.iteritems():
        if k.lower() != 'basis' and k.lower() not in ignore:
            other_matrices[k]=v

    for block in other_blocks:
        theblock = card.blocks[block]
        basis.card.add_block(theblock)
    
    for matrix in other_matrices:
        theblock = card.matrices[matrix]
        basis.card.add_block(theblock)
    
    for decay in card.decays.values():
        basis.card.add_decay(decay, preamble = decay.preamble)
    
    if card.has_block('mass'):
        basis.mass=basis.card.blocks['mass']
    if card.has_block(basis.inputs_blockname):
        basis.inputs=basis.card.blocks[basis.inputs_blockname]
        
    basis.ckm = card.matrices['vckm']
    

def expand_matrices(basis):
    '''
    Special function to populate redundant elements of matrix blocks when 
    translating to the bsmc Lagrangian so that values for all 9 entries are 
    explicitly stored before writing out the parameter card. This is to 
    stay in accordance with the SLHA format.
    The function directly modifies the _data and  _names attributes of the 
    matrices since matrices with special properties i.e. Hermitian, 
    Symmetric etc. do not grant direct access to the redundant keys such as 
    the lower triangle of a Hermitian matrix.
    '''
    all_keys = [(1,1), (1,2), (1,3),
                (2,1), (2,2), (2,3),
                (3,1), (3,2), (3,3)]
                
    for matrix in basis.card.matrices.values():
        # list of missing elements in _data member of matrix instance
        missing_keys = [k for k in all_keys if k not in matrix._data]
        
        if missing_keys:
            # randomly select parameter name since they all should have 
            # the same structure: (R|I)NAMEixj
            elename = matrix._names.values()[0]
            cname = elename[1:-3] # matrix name
            pref = elename[0] 
            for k in missing_keys:
                tail = cname + '{}x{}'.format(*k)
                matrix._data[k] = matrix[k]
                matrix._names[k] = pref + tail
                matrix._numbers[pref+tail] = k
                try:
                    matrix._re._data[k] = matrix._re[k]
                    matrix._re._names[k] = 'R' + tail
                    matrix._re._numbers['R'+tail] = k
                except AttributeError:
                    pass
                try:
                    matrix._im._data[k] = matrix._im[k]
                    matrix._im._names[k] = 'I' + tail
                    matrix._im._numbers['I'+tail] = k
                except AttributeError:
                    pass