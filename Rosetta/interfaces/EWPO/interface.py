from itertools import product
from ...bases import HiggsBasis
from ...internal.basis import checkers as check
from ...internal.errors import TranslationError
from ...internal import session
from ..interface import RosettaInterface, allowed_flav
from chisq import chisq_and_pvalue
#

class EWPOInterface(RosettaInterface):
    
    interface = 'ewpo'
    description = ("Run the EWPO interface to obtain the compatibility of a "
                   "parameter point with a fit to electroweak precision "
                   "observables")
    helpstr = "EWPO interface"
    
    parser_args = {
        ('param_card',):{
            'metavar':'PARAMCARD', 'type':str, 
            'help':'Input parameter card.'
        },
        ('--flavor',):{
            'type':str, 'default':'general', 'choices':allowed_flav, 'metavar':'',
            'help':('Specify flavor structure. Allowed values are: '+
                    ', '.join(allowed_flav)+' (default = general).')
        }
    }
    def __call__(self, args):

        basis_instance = self.from_param_card(args.param_card,
                                              flavor=args.flavor)
        
        HB_instance = basis_instance.translate(target='higgs')
        # basis_instance.modify_inputs()
        # check.modified_inputs(basis_instance)
        if HB_instance.flavor == 'universal':
            DOF = 23
        else:
            DOF = 36

        inpt = create_input(HB_instance)
        
        csq, pv = chisq_and_pvalue(inpt, flavor=HB_instance.flavor)
        
        
        session.drawline(text='EWPO results', ignore_silent=True)
        session.stdout('delta Chi^2 ({} d.o.f):'.format(DOF))
        session.stdout('{:.2f}'.format(csq))
        session.stdout('Corresponding p-value:')
        session.stdout('{:.2e}'.format(pv))
        
        session.drawline()
        session.exit(0)

def create_input(H):
    '''
    Create input arrays for delta chi-squared function of EWPO interface from 
    HiggsBasis instance. Based on flavour structure of instance.
    '''
    if H.flavor=='universal':
        return [H['HBxdGLwl'][1,1].real, H['HBxdGLze'][1,1].real,
                H['HBxdGRze'][1,1].real, H['HBxdGLzu'][1,1].real,
                H['HBxdGRzu'][1,1].real, H['HBxdGLzd'][1,1].real,
                H['HBxdGRzd'][1,1].real,
                H['dG1z'], H['dKa'], H['Lz'],
                H['cll1111'], H['cle1111'], H['cee1111'],
                H['cll1221'], H['cll1122'], H['cle1122'],
                H['cle2211'], H['cee1122'], H['cll1331'],
                H['cll1133'], (H['cle1133']+H['cle3311']),
                H['cee1133'], H['cll2332']]
    else:
        return [H['HBxdGLwl'][1,1].real, H['HBxdGLwl'][2,2].real,
                H['HBxdGLwl'][3,3].real, H['HBxdGLze'][1,1].real,
                H['HBxdGLze'][2,2].real, H['HBxdGLze'][3,3].real,
                H['HBxdGRze'][1,1].real, H['HBxdGRze'][2,2].real,
                H['HBxdGRze'][3,3].real, H['HBxdGLzu'][1,1].real,
                H['HBxdGLzu'][2,2].real, H['HBxdGLzu'][3,3].real,
                H['HBxdGRzu'][1,1].real, H['HBxdGRzu'][2,2].real,
                H['HBxdGLzd'][1,1].real, H['HBxdGLzd'][2,2].real,
                H['HBxdGLzd'][3,3].real, H['HBxdGRzd'][1,1].real,
                H['HBxdGRzd'][2,2].real, H['HBxdGRzd'][3,3].real,
                H['dG1z'], H['dKa'], H['Lz'],
                H['cll1111'], H['cle1111'], H['cee1111'],
                H['cll1221'], H['cll1122'], H['cle1122'],
                H['cle2211'], H['cee1122'], H['cll1331'],
                H['cll1133'], (H['cle1133']+H['cle3311']) ,
                H['cee1133'], H['cll2332']]



