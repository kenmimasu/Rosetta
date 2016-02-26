from __future__ import division
from math import pi
from decay import Hgg
from ...internal.basis import checkers as check
from errors import SqrtsError
################################################################################
# required info
masses = {25,5,6} # H, b, t masses
################################################################################
def production_ratios(basis, sqrts):
    '''
    Calculate ratios of all Higgs production rates w.r.t the SM prediction and 
    returns a dictionary of rescaling factors.
    '''
    bsmc = basis.translate(target='bsmc')
    check.masses(bsmc, masses, message='Signal strength calculation')
    
    # Store relevant coefficients
    MH, MT, MB = bsmc.mass[25], bsmc.mass[6], bsmc.mass[5]
    dCz, Czbx, Czz = bsmc['dCz'], bsmc['Czbx'], bsmc['Czz']
    dCw, Cwbx, Cww = bsmc['dCw'], bsmc['Cwbx'], bsmc['Cww']
    Cgg, Cza, Caa, Cabx = bsmc['Cgg'], bsmc['Cza'], bsmc['Caa'], bsmc['Cabx']
    
    if basis.flavor == 'universal':
        dYu, dYd = bsmc['BCxdYu'][1,1], bsmc['BCxdYd'][1,1]
    else: 
        dYu, dYd = bsmc['BCxdYu'][3,3], bsmc['BCxdYd'][3,3]
    
    ratios = {
        'ggH':Hgg(MH, MT, MB, Cgg, dYu, dYd),
        'VBF':VBF(dCz, Czbx, Czz, Cza, Caa, Cabx, dCw, Cwbx, Cww, sqrts=sqrts),
        'WH':WH(dCw, Cwbx, Cww, sqrts = sqrts),
        'ZH':ZH(dCz, Czbx, Czz, Cza, Cabx, sqrts = sqrts),
        'ttH':ttH(dYu)
    }
    
    return ratios
    
    
################################################################################

def VBF(dCz, Czbx, Czz, Cza, Caa, Cabx, dCw, Cwbx, Cww, sqrts=8):
    '''
    Return the approximate ratio of vector boson fusion production cross 
    section of the Higgs to the SM prediction as a function of the EFT 
    parameters in the BSMC Lagrangian.
    '''
    if sqrts == 7:
        Awbx, Aww, Azbx, Azz, Aza, Aabx, Aaa, AdCw, AdCz = \
        -1.07716, -0.099027, -0.345889, -0.0453625, -0.0213649,\
        -0.100547, -0.00445525, 1.48299, 0.517011
    elif sqrts == 8:
        Awbx, Aww, Azbx, Azz, Aza, Aabx, Aaa, AdCw, AdCz = \
        -1.11064, -0.0963591, -0.354127, -0.0437853, -0.0213434,\
        -0.100675, -0.00243252, 1.48635, 0.513653
    elif sqrts == 13:
        Awbx, Aww, Azbx, Azz, Aza, Aabx, Aaa, AdCw, AdCz = \
        -1.23262, -0.0935994, -0.402271, -0.0426703, -0.017894,\
        -0.11287, -0.00481762, 1.4862, 0.513799
    else:
        raise SqrtsError(srt(sqrts))
    
    lin = (1. + AdCz*dCz + Azbx*Czbx + Azz*Czz + Aza*Cza + Aaa*Caa 
          + Aabx*Cabx + AdCw*dCw + Awbx*Cwbx + Aww*Cww)
    
    return lin
    
def WH(dCw, Cwbx, Cww, sqrts = 8):
    '''
    Return the approximate ratio of W boson associated production cross section 
    of the Higgs to the SM prediction as a function of the EFT parameters in 
    the BSMC Lagrangian.
    '''
    if sqrts == 7:
        Awbx, Aww = 6.38857, 1.48504
    elif sqrts == 8:
        Awbx, Aww = 6.51141, 1.48866
    elif sqrts == 13:
        Awbx, Aww = 6.95706, 1.50497
    else:
        raise SqrtsError(srt(sqrts))
    
    lin = 1. + 2.*dCw + Awbx*Cwbx + Aww*Cww
    
    return lin

def ZH(dCz, Czbx, Czz, Cza, Cabx, sqrts = 8):
    '''
    Return the approximate ratio of Z boson associated production cross section 
    of the Higgs to the SM prediction as a function of the EFT parameters in the 
    BSMC Lagrangian.
    '''
    if sqrts == 7:
        Azbx, Azz, Aza, Aabx = 5.29969, 1.79327, 0.220717, 0.798601
    elif sqrts == 8:
        Azbx, Azz, Aza, Aabx = 5.40375, 1.80297, 0.220526, 0.816691
    elif sqrts == 13:
        Azbx, Azz, Aza, Aabx = 5.72384, 1.82393, 0.222767, 0.86938
    else:
        raise SqrtsError(srt(sqrts))
    
    lin = 1. + 2.*dCz + Azbx*Czbx + Azz*Czz + Aza*Cza + Aabx*Cabx 
    
    return lin 

def ttH(dYu):
    '''
    Return the approximate ratio of top associated production cross section of 
    the Higgs to the SM prediction as a function of the EFT parameters in the 
    BSMC Lagrangian.
    '''
    lin = 1. + 2.*dYu
    
    return lin
    
