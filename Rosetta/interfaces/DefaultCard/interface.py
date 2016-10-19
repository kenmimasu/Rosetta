from ..interface import RosettaInterface, allowed_flav
from ...internal import SLHA
from ...internal.errors import RosettaError
from ... import session, settings, implemented_bases
import argparse
import os

# for defaultcard
from ...internal.basis.io import write_template_card

allowed_flav = ('general', 'diagonal', 'universal')

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
            'type':str, 'metavar':'VALUE', 'default': '0',
            'help':('Set value of all parameters to VALUE. The value "random" '
                    'will set random coefficients between -1. and 1.')
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
        
        if args.value.lower()!='random':
            val = float(args.value)
        else:
            val = args.value.lower()
        
        if carry_on:
            write_template_card(instance, 
                                outputfile,
                                value = val)
                                
        session.drawline()
        session.exit(0)
        