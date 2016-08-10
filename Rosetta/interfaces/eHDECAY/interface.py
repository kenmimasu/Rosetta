from ..interface import RosettaInterface, allowed_flav
from ...internal.basis import checkers as check
from ...internal.errors import TranslationError
from eHDECAY import create_SLHA_block
from ...internal import session

#
class eHDECAYInterface(RosettaInterface):
    
    interface = 'ehdecay'
    description = ("Run the eHDECAY interface to obtain "
                   "the Higgs branching fractions")
    helpstr = "Standalone eHDECAY interface"
    
    parser_args = {
        ('param_card',):{
            'metavar':'PARAMCARD', 'type':str, 
            'help':'Input parameter card'
        },
        ('-o','--output'):{
            'metavar':'OUTPUT', 'type':str, 
            'help':'Write out a new SLHA card containing the decay block'
        },
        ('-w','--overwrite'):{
            'action':'store_true',
            'help':'Overwrite any pre-existing output file'
        },
        ('--flavor',):{
            'type':str, 'default':'general', 'choices':allowed_flav, 'metavar':'',
            'help':('Specify flavor structure. Allowed values are: '+
                    ', '.join(allowed_flav)+' (default = general)')
        },            
        ('--dependent',):{
            'action':'store_true', 
            'help':'Also write out dependent parameters to output card'
        },
        ('--EW',):{
            'action':'store_true', 
            'help':'switch on electroweak corrections in eHDECAY'
        }
    }
    def __call__(self, args):

        basis_instance = self.from_param_card(args.param_card)

        basis_instance.modify_inputs()
        check.modified_inputs(basis_instance)

        # run eHDECAY
        try:
            decayblock = create_SLHA_block(basis_instance, electroweak=args.EW)
        except TranslationError as e:
            print e
            print ('Translation to modified-SILH Basis '
                   'required for eHDECAY interface.')
            session.exit(1)

        if not args.output:
            session.stdout('###### eHDECAY results ######')
            session.stdout(str(decayblock))
            session.stdout('#############################')
            session.exit(0)
        else:

            if not args.dependent:
                basis_instance.delete_dependent()
                
            preamble = ('###################################\n'
                      + '## DECAY INFORMATION\n'
                      + '###################################')
            for decay in basis_instance.card.decays.values():
                decay.preamble = preamble
                break
                
            if 25 not in basis_instance.card.decays:
                basis_instance.card.add_decay(decayblock)
            else:
                basis_instance.card.decays[25] = decayblock
            session.log('#############################\n')

            if write_param_card(basis_instance.card, args.output,
                                overwrite=args.overwrite):
                session.log('#############################')
                session.exit(0)
            else:
                session.log('#############################')
                session.exit(0)