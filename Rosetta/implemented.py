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
    'higgs': {'mass':HB.HiggsBasis.translate,
              'silh':HB.HiggsBasis.to_silh,
              'warsaw':HB.HiggsBasis.to_warsaw},
    # to_mass implemented in WarsawBasis
    'warsaw': {'mass':WB.WarsawBasis.to_mass,
               'higgs':WB.WarsawBasis.to_mass,
               'silh':WB.WarsawBasis.to_silh},
    # to_mass, to_warsaw implemented in SILHBasis
    'silh': {'mass':SB.SILHBasis.to_mass,
             'higgs':SB.SILHBasis.to_mass,
             'warsaw':SB.SILHBasis.to_warsaw},
    # toy translations
    'template':{'mass':TB.TemplateBasis.to_mass,
                'silh':TB.TemplateBasis.to_silh},
}

relationships= relate(translations)


