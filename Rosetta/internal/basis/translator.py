from .. import session
import checkers as check

def translate(basis, target=None, verbose=True): 
    '''
    Translation function. Makes use of existing translations and attempts to 
    find a path between the basis instance and the target basis. If a path 
    exists, successive translations are performed to get from A to B. The 
    resulting basis instance of every translation performed is stored as a 
    value in the basis.translations dictionary with the basis name as a key. 
    Keyword arguments:
        target - target basis, must be implemented in Rosetta.bases.
        verbose - print some status messages.
    '''
    from ..machinery import get_path, relationships, bases

    # default target
    target = target if target is not None else basis.output_basis
    
    if target == basis.name: return basis
    
    # find the chain of bases and translation functions 
    # to get from basis.name to target
    chain = get_path(basis.name.lower(), target, relationships)            

    if not chain: # no translation possible
        inputbasis = basis.__class__.__name__
        outputbasis = bases[target].__name__
        err = ('Sorry, Rosetta cannot translate from ' +
              '{} to {}'.format(inputbasis, outputbasis))
        raise TranslationPathError(err)
    
    names = [basis.name.lower()]+[x[0] for x in chain]
    session.log('Rosetta will be performing the translation:\n'+
                    '    ' + ' -> '.join([bases[x].__name__ 
                                          for x in names])+'\n')
    
    # perform succesive translations, checking for  
    # required SM inputs/masses along the way

    current = basis
    
    required_inputs = current.required_inputs
    required_masses = current.required_masses

    for target, translate_function in chain:
        try:
            new = basis.translations[target]
            msg = 'translation ({} -> {}) already performed.\n'.format
            session.once(msg(current.name, target))
            required_inputs = set(new.required_inputs)
            required_masses = set(new.required_masses)
            current = new
            
            continue
        except KeyError:
            pass
            
        instance = bases[target](dependent=True)
        # ensure required inputs are present
        message = 'translation ({} -> {})'.format(current.name, 
                                                  instance.name)
                                                  
        # current.check_sminputs(required_inputs, message=message)
        check.sminputs(current, required_inputs, message=message)
        # current.check_masses(required_masses, message=message)
        check.masses(current, required_masses, message=message)
        
        session.verbose(message)
        new = translate_function(current, instance)
        
        # update new basis instance with non-EFT blocks, decays
        all_coeffs = (current.blocks.keys() + current.flavored.keys())
        new.get_other_blocks(current.card, all_coeffs)
        
        message = '{}'.format(instance.__class__)
        # checks & calculate dependent
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
        session.verbose('    Calling {}.calculate_dependent()\n'.format(new.__class__))
        new.init_dependent()
        new.calculate_dependent()
        
        basis.translations[target] = new
        # prepare for next step
        required_inputs = set(new.required_inputs)
        required_masses = set(new.required_masses)
        current = new
    
    if current.name =='bsmc':
        current.expand_matrices()
    else:
        current.flavor = basis.flavor
        # reduce flavor structure back to user set option
        current.set_flavor('general', basis.flavor)            
        
    if verbose:
        session.log('\nTranslation successful.\n')
    
    return current
