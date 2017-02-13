import datetime
import os

from .. import SLHA
from .. import session
from errors import BasisNameError, ParamCardReadError

def read_param_card(basis, SLHAcard = None):
    '''
    Call SLHA.read() and set up the Basis instance accordingly.
    '''
    if SLHAcard is not None:
        basis.card = SLHAcard
    else:
        basis.card = SLHA.read(basis.param_card, set_cplx=False)
    
    try:
        card_name = basis.card.blocks.get('basis',[''])[1]
        if basis.name.lower()!=card_name.lower():
            err = ('Rosetta was expecting to read an instance of ' 
                 + '{}, named "{}", '.format(basis.__class__.name, basis.name)
                 + 'but read the name '
                 + '"{}" in block "basis" of {}.'.format(card_name,
                                                         basis.param_card))
            raise BasisNameError(err)
    except KeyError:
        raise  ReadParamCardError('Formatting error for block basis. '
                           'Check input card, {}.'.format(basis.param_card))
    
    if not basis.blocks:
        # if basis.blocks not defined, default behaviour is to automatically
        # create basis.blocks structure for independent parameters
        for par in basis.independent:
            block = basis.card._parent_block(par)
            name = block.name
            if name not in basis.blocks:
                basis.blocks[name] = [par]
            else:
                basis.blocks[name].append(par)
    
    to_add = []
    for bname, blk in basis.card.matrices.iteritems():
        is_cplx = basis.flavored.get(bname,{}).get('domain','')=='complex'
        if is_cplx:
            if bname.lower().startswith('im'):
                other_part = bname[2:]
            else: 
                other_part ='IM' + bname
            if other_part not in basis.card.matrices:
                to_add.append((bname,other_part))

    for part, other_part in to_add:
        blk = basis.card.matrices[part]
        for k, v in blk.iteritems():
            if other_part.lower().startswith('im'):
                imname = 'R' + blk.get_name(k)[1:]
            else:
                imname = 'I' + blk.get_name(k)[1:]
            basis.card.add_entry(other_part, k, 0., name=imname)
        
    basis.inputs = basis.card.blocks.get(basis.inputs_blockname, None)
    basis.mass = basis.card.blocks.get('mass', None)
    
    basis.card.set_complex()
    basis.fix_matrices()
    # Deletes redundant elements of hermitian matrices
    # basis.reduce_hermitian_matrices()

def write_param_card(card, filename, overwrite=False):
    '''Write contents of card to filename'''

    dec_preamble = ('\n###################################\n'
                + '## DECAY INFORMATION\n'
                + '###################################\n')
    for decay in card.decays.values():
        decay.preamble = dec_preamble
        break

    mass_preamble = ('\n###################################\n'
                   + '## INFORMATION FOR MASS\n'
                   + '###################################\n')
    if 'mass' in card.blocks:
        card.blocks['mass'].preamble = mass_preamble

    sm_preamble = ('\n###################################\n'
                 + '## INFORMATION FOR SMINPUTS\n'
                 + '###################################\n')
    if basis.inputs_blockname in card.blocks:
        card.blocks[basis.inputs_blockname].preamble = sm_preamble

    ckm_preamble = ('\n###################################\n'
                  + '## CKM INFORMATION\n'
                  + '###################################\n')
    if 'vckm' in card.matrices:
        card.matrices['vckm'].preamble = ckm_preamble

    card_preamble = ('#'*80 +'\n'
                +'############# COEFFICIENTS TRANSLATED BY ROSETTA'
                +' MODULE  #############\n'
                +'#### See Eur.Phys.J. C75 (2015) 12, 583 (arXiv:1508.05895) '\
                 'for more details ####\n'
                +'########### {} BASIS PARAM CARD GENERATED {}  ##########'\
                 '#\n'.format(card.blocks['basis'][1].upper(),
                              datetime.datetime.now().ctime().upper())
                +'#'*80 +'\n\n')
        
    if os.path.exists(filename) and not overwrite:
        session.log('\n{} already exists.'.format(filename))
        carry_on = session.query('Overwrite?', default='no')
    else:
        carry_on=True
    
    if carry_on:
        special_blocks = ['loop','mass',basis.inputs_blockname,'yukawa','vckm','basis']
        coefforder = SLHA.sortblocks(card, ignore = special_blocks)
        card.write(filename, blockorder=special_blocks + coefforder,
                                     preamble=card_preamble)
        session.log('Wrote "{}".'.format(filename))
        session.log('')
        
        return True
    else:
        return False

def write_template_card(basis, filename, value=0.):
    from ..machinery import bases
    from ..constants import default_masses, particle_names, default_inputs, input_names
    try: 
        val = float(value)
        rand = False
    except ValueError:
        if value.lower() == 'random':
            import random
            rand = True
        else:
            session.log('In write_template_card: "value" keyword argument '
                        'must either be a number or the string, "random".')  
            sys.exit()          
    # newinstance = bases[basis.name](flavor=basis.flavor, dependent=False)
    newinstance = basis.__class__(flavor=basis.flavor, dependent=False)
    for k in newinstance.keys():            
            try:
                if rand:
                    newinstance[k] = complex(random.uniform(-1.,1.),
                                             random.uniform(-1.,1.))
                else:
                    newinstance[k] = complex(val, val)
            except TypeError as e:
                if rand:
                    newinstance[k] = random.uniform(-1.,1.)
                else:
                    newinstance[k] = val
    
    newinstance.reduce_hermitian_matrices()
    
    SLHA_card = newinstance.card
    
    mass_preamble = ('\n###################################\n'
                   + '## INFORMATION FOR MASS\n'
                   + '###################################\n')

    massblock = SLHA.NamedBlock(name='mass', preamble=mass_preamble)
    for m in basis.required_masses: 
        massblock.new_entry(m, default_masses[m], 
                           name='M%s' % particle_names[m])
    SLHA_card.add_block(massblock)
                
    sm_preamble = ('\n###################################\n'
                 + '## INFORMATION FOR {}\n'
                 + '###################################\n'.format(
                     basis.inputs_blockname.upper()
                 ))
                
    inputblock = SLHA.NamedBlock(name=basis.inputs_blockname, 
                                 preamble=sm_preamble)
    for m in basis.required_inputs: 
        inputblock.new_entry(m, default_inputs[m], 
                             name='%s' % input_names[m])
    SLHA_card.add_block(inputblock)
    
    title = ' ROSETTA: {} BASIS INPUT CARD '.format(basis.name.upper())
    time = ' GENERATED {} '.format(datetime.datetime.now().ctime().upper())
    nleft = int((80-len(title))/2.)
    nright = 80-len(title)-nleft
    preamble = ( '#'*80 + '\n'
                +'#'*nleft + title + '#'*nright +'\n'
                +'#### See Eur.Phys.J. C75 (2015) 12, 583 (arXiv:1508.05895) '\
                 'for more details ####\n'
                +'#'*22 + time + '#'*22 + '\n'
                +'#'*80 +'\n')
    
    special_blocks = ['mass',basis.inputs_blockname,'vckm','basis']
            
    theorder = SLHA.sortblocks(SLHA_card, ignore = special_blocks)
    SLHA_card.write(filename, blockorder=special_blocks + theorder, 
                    preamble = preamble, postamble = ('\n'+'#'*80)*2)
    session.log('wrote {}'.format(filename))
    session.log('')
    