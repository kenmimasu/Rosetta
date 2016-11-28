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

from Rosetta.interfaces.Lilith import Lilith

# instance = HB.HiggsBasis(flavor='universal', param_card = '../HiggsBasis_universal_1e-3.dat', translate=False)
# instance = HZ.HISZ(flavor='universal', param_card = '../HISZ_universal_1e-3.dat', translate=False)
instance = HZ.HISZ(flavor='universal', param_card = 'Cards/HISZ_universal.dat')

lik = Lilith.compute_likelihood(instance)

session.log('Lilith Likelihood: '+str(lik))
session.log('#############################')
session.log('')