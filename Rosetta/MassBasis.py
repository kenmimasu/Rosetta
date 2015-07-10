from internal import Basis
import HiggsBasis as HB

################################################################################
# Mass basis class
class MassBasis(Basis.Basis):
    name ='mass'
    # all coefficients stored in HiggsBasis.blocks
    independent = [c for v in HB.HiggsBasis.blocks.values() for c in v]
    # independent = HiggsBasis.independent + HiggsBasis.dependent
    dependent=[]
    blocks = {k.replace('HB','MB'):v for k,v in HB.HiggsBasis.blocks.iteritems()}
    # def __init__(self,*args,**kwargs):
    #     super(MassBasis, self).__init__(*args,**kwargs)
################################################################################
