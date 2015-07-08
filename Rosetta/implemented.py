import MassBasis as MB
import HiggsBasis as HB
import WarsawBasis as WB
import SILHBasis as SB
import TemplateBasis as TB
from internal.relate import relate
################################################################################
# dict of implemented bases
bases = {'warsaw':WB.WarsawBasis,'higgs':HB.HiggsBasis,
         'mass':MB.MassBasis,'template':TB.TemplateBasis,
         'silh':SB.SILHBasis}
################################################################################
# dict representing the translation relationships between bases via the function
# object that implements the translation in question
translations = {
    # translate function not implemented in MassBasis making it default to the
    # trivial translate function defined in the parent class, Basis.
    'mass':{'mass':MB.MassBasis.translate},
    # Also a trivial translation, here HiggsBasis.calculate_dependent() does
    # everything needed
    'higgs': {'mass':HB.HiggsBasis.translate},
    # to_mass implemented in WarsawBasis
    'warsaw': {'mass':WB.WarsawBasis.to_mass,
               'silh':WB.WarsawBasis.to_silh},
    # to_mass, to_warsaw implemented in SILHBasis
    'silh': {'mass':SB.SILHBasis.to_mass,
             'warsaw':SB.SILHBasis.to_warsaw},
    # toy translations
    'template':{'mass':TB.TemplateBasis.to_mass,
                'silh':TB.TemplateBasis.to_silh},
}
# translations = {
#          'A':{'B':'A.Bfunc','C':'A.Cfunc'},
#          'B':{'C':'B.Cfunc','D':'B.Dfunc'},
#          'C':{'A':'C.Afunc','D':'C.Dfunc','B':'CBfunc'},
#          'D':{'A':'D.Afunc'},
#          'E':{'D':'E.Dfunc'},
#          'F':{'C':'F.Cfunc','E':'F.Efunc'},
#          'G':{'E':'G.Efunc'}
#          }


relationships= relate(translations)
# for k,v in rels.items():
#     print k,v# print relations

