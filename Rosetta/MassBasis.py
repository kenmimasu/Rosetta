from internal import Basis
import HiggsBasis as HB
################################################################################
class MassBasis(Basis.Basis):
    '''
    The main output basis for the Rosetta translation module. Based entirely on 
    the Higgs basis parametrisation without assuming the SU(2)xU(1) preserving 
    relations. It is a general collection of most possible Lorentz and charge 
    conserving structures at dimension 6 excluding those not described in the 
    current incarnation of the Higgs basis. Most importantly, a FeynRules/UFO 
    implementation of the interactions in this model is associated with Rosetta 
    and can be found at https://feynrules.irmp.ucl.ac.be/wiki/EFTmassbasis.
    '''
    
    name ='mass'
    # all coefficients stored in HiggsBasis.blocks
    independent = [c for v in HB.HiggsBasis.blocks.values() for c in v]
    # rename Higgs Basis blocks, replacing HB with MB
    blocks = {k.replace('HB','MB'):v for k,v in HB.HiggsBasis.blocks.iteritems()}
    # all other undefined behavoiur inherited from Basis.Basis by default
################################################################################
