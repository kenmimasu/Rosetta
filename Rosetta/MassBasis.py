from Basis import Basis
from HiggsBasis import HiggsBasis
################################################################################
# Mass basis class
class MassBasis(Basis):
    independent = HiggsBasis.independent + HiggsBasis.dependent
    dependent=[]
    # def __init__(self,*args,**kwargs):
    #     super(MassBasis, self).__init__(*args,**kwargs)
################################################################################
