import __init__ as Rosetta
from math import sqrt
import os
import sys
from distutils.version import StrictVersion
from tempfile import mkdtemp
import subprocess as sub
from collections import namedtuple
from ..internal.settings import settings 
################################################################################
# check for NumPy >= 1.6.1 and Scipy >= 0.9.0
try:
    import numpy
    if StrictVersion(numpy.__version__) < StrictVersion('1.6.1'):
        raise ImportError
except ImportError:
    err = ('NumPy version 1.6.1 or more recent '
           'must be installed to use Lilith interface')   
    raise RuntimeError('Lilith interface error: {}'.format(err))

try:
    import scipy
    if StrictVersion(scipy.__version__) < StrictVersion('0.9.0'):
        raise ImportError
except ImportError:
    err = ('SciPy version 0.9.0 or more recent '
           'must be installed to use Lilith interface')   
    raise RuntimeError('Lilith interface error: {}'.format(err))

# Lilith path
try:
    import lilith
except ImportError:
    try:
        Lilith_dir = settings['Lilith_dir']
        sys.path.append(Lilith_dir)
        import lilith
    except KeyError:
        err = ('Could not find option "Lilith_dir" in Rosetta/config.txt')
        raise RuntimeError('Lilith interface error: {}'.format(err)) 
    except ImportError:
        err = ('Failed to import Lilith')
        raise RuntimeError('Lilith interface error: {}'.format(err))
################################################################################

def generate_input(prod, decay, width):
    
    lilith_input = '''
<?xml version="1.0"?>

<lilithinput>
  <!-- signal strengths in theory space, like mu(gg -> H -> ZZ), as input -->
  <signalstrengths part="h">
    <mass>125</mass>
    <!-- optional:
    if not given, Higgs mass of 125 GeV is assumed
    valid Higgs masses are in the [123,128] GeV range
    -->

    <!--
    "VV" decay can be split into "WW" and "ZZ" in all that follows
    -->
    <mu prod="ggH" decay="gammagamma">1.0</mu>
    <mu prod="ggH" decay="WW">1.0</mu>
    <mu prod="ggH" decay="ZZ">1.0</mu>
    <mu prod="ggH" decay="bb">1.0</mu>
    <mu prod="ggH" decay="tautau">1.0</mu>
    <mu prod="ggH" decay="mumu">1.0</mu>

    <!--
    if necessary, possible to specify "VBF", "WH" and "ZH" production
    instead of a common "VVH"
    -->
    <mu prod="VBF" decay="gammagamma">1.0</mu>
    <mu prod="VBF" decay="ZZ">1.0</mu>
    <mu prod="VBF" decay="WW">1.0</mu>
    <mu prod="VBF" decay="bb">1.0</mu>
    <mu prod="VBF" decay="tautau">1.0</mu>
    <mu prod="VBF" decay="mumu">1.0</mu>
    <mu prod="WH" decay="gammagamma">1.0</mu>
    <mu prod="WH" decay="ZZ">1.0</mu>
    <mu prod="WH" decay="WW">1.0</mu>
    <mu prod="WH" decay="bb">1.0</mu>
    <mu prod="WH" decay="tautau">1.0</mu>
    <mu prod="WH" decay="mumu">1.0</mu>
    <mu prod="ZH" decay="gammagamma">1.0</mu>
    <mu prod="ZH" decay="ZZ">1.0</mu>
    <mu prod="ZH" decay="WW">1.0</mu>
    <mu prod="ZH" decay="bb">1.0</mu>
    <mu prod="ZH" decay="tautau">1.0</mu>
    <mu prod="ZH" decay="mumu">1.0</mu>
    
    <!--
    ttH is optional: if not provided, SM-like ttH is assumed
    -->
    <mu prod="ttH" decay="gammagamma">1.0</mu>
    <mu prod="ttH" decay="ZZ">1.0</mu>
    <mu prod="ttH" decay="WW">1.0</mu>
    <mu prod="ttH" decay="bb">1.0</mu>
    <mu prod="ttH" decay="tautau">1.0</mu>
    <mu prod="ttH" decay="mumu">1.0</mu>

    <!--
    the following is optional: if not given, no decay into invisible particles
    used for the ZH/VBF->ll+invisible constraint
    corresponds to (sigma(ZH)/sigma(ZH_SM))*BR(H->invisible)
    and similarly for VBF
    -->
    <redxsBR prod="ZH" decay="invisible">0.0</redxsBR>
    <redxsBR prod="VBF" decay="invisible">0.0</redxsBR>
  </signalstrengths>
</lilithinput>
'''
    
def compute_likelihood(prod, decay, width):
    xml_input = generate_input(prod, decay)
    lilithcalc.computelikelihood(userinput=xml_input)
    return lilithcalc.l    














