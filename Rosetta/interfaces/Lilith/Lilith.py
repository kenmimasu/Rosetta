from math import sqrt
import os
import sys
from itertools import product
from . import lilith, LilithInterfaceError
from ..SignalStrengths.production import production
from ..SignalStrengths.decay import decay
from ...internal import session
################################################################################
reference = 'J. Bernon & B. Dumont, Eur. Phys. J. C75 (2015) 9, 440'
################################################################################

channels = {'bb':(5,-5),'mumu':(13,-13), 'tautau':(15,-15), 
            'gammagamma':(22,22), 'ZZ':(23,23), 'WW':(24,-24),}
            
lilithcalc = lilith.Lilith(verbose=False,timer=False)
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
    
    for kp, (kd, vd) in product(prod.keys(), channels.items()):
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
    
