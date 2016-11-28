from itertools import product
from ...internal.basis import checkers as check
from ...internal.errors import TranslationError
from ...internal import session
from ..interface import RosettaInterface, allowed_flav
from Lilith import compute_likelihood
#
channels = {'bb':(5,-5),'mumu':(13,-13), 'tautau':(15,-15), 
            'gammagamma':(22,22), 'ZZ':(23,23), 'WW':(24,-24)}

class LilithInterface(RosettaInterface):
    
    interface = 'lilith'
    description = ('Run the Lilith interface to obtain the likelihood value '
                   'with respect to the latest Higgs signal strength data for '
                   'a particular point in EFT parameter space')
    helpstr = 'Standalone Lilith interface'
    
    parser_args = {
        ('param_card',):{
            'metavar':'PARAMCARD', 'type':str,
            'help':'Input parameter card.'
        },
        # ('--sqrts',):{
        #     'type':int, 'choices':(7,8,13), 'default':8,
        #     'help':'Specify pp collider centre-of mass energy in TeV'
        # },
        ('--flavor',):{
            'type':str, 'default':'general', 'choices':allowed_flav, 'metavar':'',
            'help':('Specify flavor structure. Allowed values are: '+
                    ', '.join(allowed_flav)+' (default = general)')
        },
        ('--squares',):{
            'action':'store_true',
            'help':('Retain quadratic order in EFT coefficients in the '
                    'SignalStrengths interface (NOT IMPLEMENTED)')
        }
    }
    def __call__(self, args):

        basis_instance = self.from_param_card(args.param_card)
        basis_instance.modify_inputs()
        check.modified_inputs(basis_instance)

        # if not args.dependent:
        #     basis_instance.delete_dependent()

        # run SignalStrengths
        
        likelihood = compute_likelihood(basis_instance, sqrts=8)
        
        session.drawline(text='Lilith results', ignore_silent=True)
        session.stdout('Likelihood = {:.3f}'.format(likelihood))         
        session.drawline()
        session.exit(0)





