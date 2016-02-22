#!/usr/bin/env python
import tempfile
import os
import sys
import re
import random

sys.path.append('../')
# sys.path.append('Users/Ken/Work/Projects/Rosetta/dev/rosetta/')

from Rosetta import HiggsBasis as HB
from Rosetta import WarsawBasis as WB
from Rosetta import SILHBasis as SB
from Rosetta import BSMCharacterisation as MB
from Rosetta import TemplateBasis as TB
from Rosetta import HISZ as HZ
# from Rosetta import MufBasis as MUF
from Rosetta.internal import SLHA

from Rosetta.interfaces.Lilith import Lilith

instance = HB.HiggsBasis(flavor='universal', param_card = '../HiggsBasis_universal_1e-3.dat', translate=False)



print Lilith.compute_likelihood(instance)
