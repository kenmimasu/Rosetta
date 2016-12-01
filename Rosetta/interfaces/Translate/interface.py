# from ..internal import SLHA
# from ..internal.errors import RosettaError
from ..interface import RosettaInterface, allowed_flav
from ... import session, settings, implemented_bases
import os
# for translate
from ...internal.basis import checkers as check
from ...internal.basis import write_param_card
from ...internal.errors import TranslationError

allowed_flav = ('general', 'diagonal', 'universal')

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
        # Initialise starting basis instance from card
        basis_instance = self.from_param_card(args.param_card, flavor=args.flavor)

        # Translate to target basis instance
        newbasis = basis_instance.translate(target=args.target)
        
        # Call modify_inputs() to apply any changes in e.g. mW & check
        # Only ever do this for the FINAL output basis instance to avoid 
        # multiple modifications of parameters
        newbasis.modify_inputs()
        check.modified_inputs(newbasis)
        
        # delete imaginary parts of diagonal elements in hermitian matrices
        if args.target != 'bsmc':
            newbasis.reduce_hermitian_matrices()
        
        # Remove occurrences of parameters defined as dependent for output 
        # unless --dependent option specified OR output basis is e.g. tied to a 
        # MC implementation that would like to see all parameters
        if not args.dependent and args.target.lower()!='bsmc':
            newbasis.delete_dependent()
        
        # Preamble for decay blocks
        preamble = ('###################################\n'
                  + '## DECAY INFORMATION\n'
                  + '###################################')
        for decay in newbasis.card.decays.values():
            decay.preamble = preamble
            break
        
        # run eHDECAY for Higgs branching fractions
        try:
            if args.ehdecay:
                from ...interfaces.eHDECAY.eHDECAY import create_SLHA_block

                decayblock = create_SLHA_block(basis_instance)
                # Create Higgs decay block if not present
                if 25 not in newbasis.card.decays:
                    newbasis.card.add_decay(decayblock)
                else:
                    newbasis.card.decays[25] = decayblock
                # session.drawline()
                
        except TranslationError as e:
            # Catch translation error in map to SILH
            session.log('')
            session.log('Translation to SILH Basis required, skipping eHDECAY.')
            session.log('TranslationError: ' + (e))
            session.log('')
            
        except ImportError as e:
            session.log('')
            session.log('eHDECAY interface not loaded, skipping eHDECAY.')
            session.log('ImportError: ' + str(e))
            session.log('')
        
    
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
            session.drawline()
            session.exit(0)
        else:
            session.drawline()
            session.exit(0)