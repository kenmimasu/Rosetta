#!/usr/bin/env python
import tempfile
import os
import sys
import re
import random

import pytest

from rosetta import HiggsBasis as HB
from rosetta import WarsawBasis as WB
from rosetta import SILHBasis as SB
from rosetta import BSMCharacterisation as MB
from rosetta import TemplateBasis as TB
from rosetta import HISZ as HZ
from rosetta.internal import SLHA, session

from rosetta.interfaces.Lilith import Lilith

@pytest.mark.xfail(reason="Requires Lilith to be installed, though also HISZ_universal.dat seems to have Lambda=0")
def test_lilith(request):
    cards_dir = request.path.parent / 'Cards'
    # instance = HB.HiggsBasis(flavor='universal', param_card = '../HiggsBasis_universal_1e-3.dat', translate=False)
    # instance = HZ.HISZ(flavor='universal', param_card = '../HISZ_universal_1e-3.dat', translate=False)
    instance = HZ.HISZ(flavor='universal', param_card = cards_dir / 'HISZ_universal.dat')

    lik = Lilith.compute_likelihood(instance)

    session.log('Lilith Likelihood: '+str(lik))
    session.log('#############################')
    session.log('')