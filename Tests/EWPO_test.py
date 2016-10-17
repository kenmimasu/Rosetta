#!/usr/bin/env python
import tempfile
import os
import sys
import re
import random

sys.path.append('../')

from Rosetta import HiggsBasis as HB
from Rosetta import WarsawBasis as WB
from Rosetta import SILHBasis as SB
from Rosetta import BSMCharacterisation as MB
from Rosetta import TemplateBasis as TB
from Rosetta import HISZ as HZ
from Rosetta.internal import SLHA, session

from Rosetta.interfaces.EWPO.chisq import chisq, chisq_and_pvalue

# dGLwl1x1 dGLwl2x2 dGLwl3x3 
# dGLze1x1 dGLze2x2 dGLze3x3 dGRze1x1 dGRze2x2 dGRze3x3 
# dGLzu1x1 dGLzu2x2 dGLzu3x3 dGRzu1x1 dGRzu2x2 
# dGLzd1x1 dGLzd2x2 dGLzd3x3 dGRzd1x1 dGRzd2x2 dGRzd3x3 
# dG1z dKa Lz 
# cll1111 cle1111 cee1111 cll1221 cll1122 cle1122  cle2211 cee1122
# cll1331 cll1133 (cle1133+cle3311) cee1133 cll2332


def create_input(H):
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

# H = HB.HiggsBasis(flavor='universal', param_card='Cards/HiggsBasis_universal_1e-3.dat')
# H = HB.HiggsBasis(flavor='universal', param_card='Cards/HiggsBasis_universal.dat')
S = SB.SILHBasis(flavor='universal', param_card='Cards/SILHBasis_universal_1e-3.dat')
H = S.translate(target='higgs')
inp = create_input(H)
print chisq_and_pvalue(inp, flavor=H.flavor)


# instance = HZ.HISZ(flavor='universal', param_card = '../HISZ_universal_1e-3.dat', translate=False)
# instance = HZ.HISZ(flavor='universal', param_card = '../HISZ_universal_0.dat', translate=False)

# bsmc = basis.translate(target='bsmc')


# lik = Lilith.compute_likelihood(instance)

# session.log('Lilith Likelihood: '+str(lik))
# session.log('#############################')
# session.log('')