from ..internal import SLHA
from ..internal.errors import RosettaError
from .. import session, settings, implemented_bases
import argparse
import os
# for translate
from ..internal.basis import checkers as check
from ..internal.basis import write_param_card
from ..interfaces.eHDECAY.eHDECAY import SLHAblock
from ..internal.errors import TranslationError
# for defaultcard
from ..internal.basis.io import write_template_card

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
        
        session.log('Reading "{}"\n'.format(param_card))
        
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
        
        session.log('Basis class used to read in param card:\n'+
                    '    {}\n'.format(mybasis))
        
        # create instance of mybasis, automatically translating to target basis
        return  mybasis(param_card=param_card, flavor=flavor)
        
class TranslateInterface(RosettaInterface):
    
    interface = 'translate'
    description = ("Read in an SLHA format parameter card in a particular "
                   "basis and write a new card in another implemented basis.")
    helpstr = ("Translate in SLHA input card from one basis to another")
    
    parser_args = {
        ('param_card',):{
            'metavar':'PARAMCARD', 'type':str, 
            'help':'Input parameter card.'
        },
        ('-o','--output'):{
            'metavar':'OUTPUT', 'type':str, 
            'help':'Output file name. Default: [PARAMCARD]_new'
        },
        ('-w','--overwrite'):{
            'action':'store_true',
            'help':'Overwrite any pre-existing output file.'
        },
        ('--target',):{
            'type':str, 'default':'bsmc', 'choices':implemented_bases.keys(),
            'metavar':'',
            'help':('Basis into which to translate. Allowed values are: '+
                    ', '.join(implemented_bases.keys()) + ' (default = bsmc)')
        },
        ('--flavor',):{
            'type':str, 'default':'general', 'choices':allowed_flav, 'metavar':'',
            'help':('Specify flavor structure. Allowed values are: '+
                    ', '.join(allowed_flav)+' (default = general).')
        },
        ('--dependent',):{
            'action':'store_true', 
            'help':'Also write out dependent parameters to output card'
        },
        ('--ehdecay',):{
            'action':'store_true', 
            'help':'Interface with eHDECAY for Higgs branching fractions.'
        }
    }

    def __call__(self, args):
        
        basis_instance = self.from_param_card(args.param_card)
        
        newbasis = basis_instance.translate(target=args.target)
        
        newbasis.modify_inputs()
        check.modified_inputs(newbasis)
        
        # delete imaginary parts of diagonal elements in hermitian matrices
        if args.target != 'bsmc':
            newbasis.reduce_hermitian_matrices()
        
        if not args.dependent and args.target.lower()!='bsmc':
            newbasis.delete_dependent()
        
        preamble = ('###################################\n'
                  + '## DECAY INFORMATION\n'
                  + '###################################')
        for decay in newbasis.card.decays.values():
            decay.preamble = preamble
            break
        
        # run eHDECAY
        try:
            if args.ehdecay:
                decayblock = SLHAblock(basis_instance)
                
                if 25 not in newbasis.card.decays:
                    newbasis.card.add_decay(decayblock)
                else:
                    newbasis.card.decays[25] = decayblock
                session.log('#############################\n')
                
        except TranslationError as e:
            print e
            print 'Translation to SILH Basis required, skipping eHDECAY.'
    
        # write param card
        if not args.output: # output file name
            inputfile = os.path.basename(args.param_card)
            if '.' in args.param_card:
                new_name = inputfile.split('.')
                new_name.insert(-1,'_new.')
                new_param_card = ''.join(new_name)
            else:
                new_param_card = args.param_card+'_new'
        else:
            new_param_card = args.output
            
        if write_param_card(newbasis.card, new_param_card, 
                            overwrite=args.overwrite): 
            session.log('#############################')
            session.exit(0)
        else:
            session.log('#############################')
            session.exit(0)


class DefaultCardInterface(RosettaInterface):
    
    interface = 'defaultcard'
    description = ("Generate a parameter card for an implemented basis")
    helpstr = ("Generate a parameter card for an implemented basis")
    
    parser_args = {
        ('basis',):{
            'metavar':'BASIS', 'type':str, 'choices':implemented_bases.keys(),
            'help':('Basis class for which to generate the parameter card. ' +
            'Allowed values are: ' + ', '.join(implemented_bases.keys()) + 
            ' (default = bsmc)')
        },
        ('-o','--output'):{
            'metavar':'OUTPUT', 'type':str, 
            'help':'Output file name. Default: [BASIS]_[FLAVOR]_default.dat'
        },
        ('-w','--overwrite'):{
            'action':'store_true',
            'help':'Overwrite any pre-existing output file'
        },
        ('--flavor',):{
            'type':str, 'default':'general', 'choices':allowed_flav, 
            'metavar':'FLAVOR',
            'help':('Specify flavor structure. Allowed values are: '+
                    ', '.join(allowed_flav)+' (default = general)')
        },
        ('--value',):{
            'type':float, 'metavar':'VALUE', 'default': 0.,
            'help':'Set value of all parameters to VALUE'
        }
        # ('--dependent',):{
        #     'action':'store_true',
        #     'help':'Also write out dependent parameters to output card'
        # },
        # ('--ehdecay',):{
        #     'action':'store_true',
        #     'help':'Interface with eHDECAY for Higgs branching fractions.'
        # }
    }

    def __call__(self, args):
        
        instance = implemented_bases[args.basis.lower()](flavor = args.flavor)
        
        if not args.output: # output file name
            outputfile = '{}_{}_default.dat'.format(args.basis, args.flavor)
        else:
            outputfile = args.output
        
        if os.path.exists(outputfile) and not args.overwrite:
            session.log('{} already exists.'.format(outputfile))
            carry_on = session.query('Overwrite?', default='no')
        else:
            carry_on=True
        
        if carry_on:
            write_template_card(instance, 
                                outputfile,
                                value = args.value)
                                
        session.log('#############################')
        session.exit(0)
        

class ReadParamCardError(RosettaError):
    '''Exception raised inside RosettaInterface.read_param_card()'''
    pass