import sys

sys.path.append('../')

from Rosetta import HiggsBasis as HB
from Rosetta.interfaces.dihiggs import dihiggs, production_xs

instance = HB.HiggsBasis(flavor='universal', 
                         param_card = 'Cards/HiggsBasis_universal_1e-3.dat')

print dihiggs.get_xs_times_br(instance)

