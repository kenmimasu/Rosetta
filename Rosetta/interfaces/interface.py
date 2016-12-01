from ..internal import SLHA
from .. import session, settings, implemented_bases
import argparse
import os

allowed_flav = ('general', 'diagonal', 'universal')

class RosettaInterface(object):

    description = ''
    
    def __init__(self, parsed_args):
        self(parsed_args)
    
    def __call__(self):
        pass
        
    def from_param_card(self, param_card, flavor='general'):
        # Read basis class from parameter card, must be defined in 
        # Rosetta.implemented_bases i.e. a valid basis class implementation 
        # exists in the Rosetta/bases directory.
        
        session.log('Reading "{}"'.format(param_card))
        
        thecard = SLHA.read(param_card, set_cplx = False)
        if 'basis' not in thecard.blocks:
            raise ReadParamCardError('\n    Parameter card ' +
                                     '{} '.format(param_card) +
                                     'does not have a block "BASIS".')
        else:
            try:
                basis = thecard.blocks['basis'][1]
            except KeyError:
                raise ReadParamCardError('\n    Formatting error for block basis. '
                                       'Check input card, {}.'.format(param_card))
                                       
        try:
            mybasis = implemented_bases[basis.lower()]
        except KeyError:
            err = ('\n    Element 1 of block "BASIS", "{}", '.format(basis) +
                   'not recognised. Rosetta accepts one of: ' +
                   '{}.'.format(', '.join(implemented_bases.keys())))
            raise ReadParamCardError(err)
        
        session.log('\nBasis class used to read in param card:\n'+
                    '    {}\n'.format(mybasis))
        
        # create instance of mybasis, automatically translating to target basis
        return  mybasis(param_card=param_card, flavor=flavor)
