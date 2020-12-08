from itertools import product
from ...internal.basis import checkers as check
from ...internal.errors import TranslationError
from ...internal import session
from ..interface import RosettaInterface, allowed_flav
from production import production
from decay import decay
#
channels = {'bb':(5,-5),'mumu':(13,-13), 'tautau':(15,-15), 
            'gammagamma':(22,22), 'ZZ':(23,23), 'WW':(24,-24)}

class SignalStrengthsInterface(RosettaInterface):
    
    interface = 'signalstrengths'
    description = ("Run the SignalStrengths interface to obtain "
                   "the mu's for all of the Higgs "
                   "production and decay channels.")
    helpstr = "Standalone SignalStrengths interface"
    
    parser_args = {
        ('param_card',):{
            'metavar':'PARAMCARD', 'type':str, 
            'help':'Input parameter card.'
        },
        ('--sqrts',):{
            'type':int, 'choices':(7,8,13), 'default':8, 'metavar':'',
            'help':('Specify pp collider centre-of mass energy in TeV. ' +
                   'Allowed values are: '+', '.join(('7','8','13'))+' (default = 8).')
        },
        ('--flavor',):{
            'type':str, 'default':'general', 'choices':allowed_flav, 'metavar':'',
            'help':('Specify flavor structure. Allowed values are: '+
                    ', '.join(allowed_flav)+' (default = general).')
        },
        # ('--squares',):{
        #     'action':'store_true',
        #     'help':'Retain quadratic order in EFT coefficients (NOT IMPLEMENTED)'
        # }
    }
    def __call__(self, args):

        basis_instance = self.from_param_card(args.param_card,flavor=args.flavor)
        basis_instance.modify_inputs()
        check.modified_inputs(basis_instance)

        # if not args.dependent:
        #     basis_instance.delete_dependent()

        # run SignalStrengths
        
        # ratios of decay partial widths and total width
        decays = decay(basis_instance, 
                       electroweak=True, SM_BRs=None, ratio=True)
        # ratios of production cross sections
        prods = production(basis_instance, sqrts=args.sqrts)
        
        session.drawline(text='SignalStrengths results', ignore_silent=True)
        session.stdout('  production  decay           mu ') 
        session.stdout('  ---------------------------------') 
        for kp, (kd, vd) in product(prods.keys(), channels.items()):
            mu = prods[kp]*decays[vd]/decays['WTOT']
            mustr = '  {:<8}    {:<10}      {:.3f}'.format(kp, kd, mu)
            session.stdout(mustr)
        
        session.drawline()
        session.exit(0)





