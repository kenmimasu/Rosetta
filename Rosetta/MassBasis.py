from Basis import Basis
from HiggsBasis import HiggsBasis
################################################################################
# Mass basis class
class MassBasis(Basis):
    independent = HiggsBasis.independent + HiggsBasis.dependent
    dependent=[]
    blocks = {k.replace('HB','MB'):v for k,v in HiggsBasis.blocks.iteritems()}
    # def __init__(self,*args,**kwargs):
    #     super(MassBasis, self).__init__(*args,**kwargs)
################################################################################
