from math import sqrt
import os
import sys
from itertools import product
from .errors import LilithImportError, LilithInterfaceError
from ..SignalStrengths.production import production
from ..SignalStrengths.decay import decay
from ...internal import session
from ...internal.settings import config

################################################################################
reference = 'J. Bernon & B. Dumont, Eur. Phys. J. C75 (2015) 9, 440'
################################################################################

channels = {'bb':(5,-5),'mumu':(13,-13), 'tautau':(15,-15), 
            'gammagamma':(22,22), 'ZZ':(23,23), 'WW':(24,-24),}
            
def _get_lillith():
    """Load lillith module"""
    try:
        from . import lilith
    except ImportError:
        try:
            Lilith_dir = config['Lilith_dir']
            sys.path.append(Lilith_dir)
            from . import lilith
        except KeyError:
            err = ('Could not find option "Lilith_dir" in Rosetta/config.txt')
            raise LilithImportError(err)
        except ImportError:
            err = ('check Lilith_dir option in '
                  'Rosetta/config.txt')
            raise LilithImportError(err)
    return lilith

_cached_calc = None

def lilithcalc():
    global _cached_calc
    if _cached_calc is None:
        lilith = _get_lillith()
        _cached_calc = lilith.Lilith(verbose=False,timer=False)
    return _cached_calc
################################################################################

def compute_likelihood(basis, sqrts=8):
    session.cite('Lilith', reference)
    
    # ratios of decay partial widths and total width
    session.verbose('Calculating Higgs decay branching fractions.')
    decays = decay(basis, electroweak=True, SM_BRs=None, ratio=True)
    # ratios of production cross sections
    session.verbose('Calculating Higgs production cross sections')
    prods = production(basis, sqrts=sqrts)
    
    xml_input = generate_input(basis.mass[25], prods, decays)
    
    lilithcalc.computelikelihood(userinput=xml_input)
    
    return lilithcalc.l    
        
def generate_input(MH, prod, decay):
    mus = []
    
    for kp, (kd, vd) in product(list(prod.keys()), list(channels.items())):
        mu = prod[kp]*decay[vd]/decay['WTOT']
        mustr = '<mu prod="{}" decay="{}">{}</mu>'.format(kp, kd, mu)
        mus.append(mustr)
        
    return \
'''<?xml version="1.0"?>

<lilithinput>
  <signalstrengths part="h">
    <mass>{}</mass>
    {}
  </signalstrengths>
</lilithinput>
'''.format(MH,'\n    '.join(mus))
    
