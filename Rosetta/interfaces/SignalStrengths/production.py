from __future__ import division
from math import pi
from decay import Hgg
from ...internal.basis import checkers as check
################################################################################
# required info
masses = {25,5,6} # H, b, t masses
################################################################################
def production_ratios(basis, sqrts):
    '''
    Calculate ratios of all Higgs production rates w.r.t the SM prediction and 
    returns a dictionary of rescaling factors.
    '''
    bsmc = basis.translate(target='bsmc',verbose=False)
    check.masses(bsmc, masses, message='Signal strength calculation')
    
    # Store relevant coefficients
    MH, MT, MB = bsmc.mass[25], bsmc.mass[6], bsmc.mass[5]
    dCz, Czbx, Czz = bsmc['dCz'], bsmc['Czbx'], bsmc['Czz']
    Cgg, Cza, Caa = bsmc['Cgg'], bsmc['Cza'], bsmc['Caa']
    
    if basis.flavor == 'universal':
        dYu, dYd = bsmc['BCxdYu'][1,1], bsmc['BCxdYd'][1,1]
    else: 
        dYu, dYd = bsmc['BCxdYu'][3,3], bsmc['BCxdYd'][3,3]
    
    ratios = {
        'ggH':Hgg(MH, MT, MB, Cgg, dYu, dYd),
        'VBF':VBF(dCz, Czbx, Czz, Cza, Caa),
        'WH':WH(dCz, Czbx, Czz, Cza, Caa, sqrts = sqrts),
        'ZH':ZH(dCz, Czbx, Czz, Cza, Caa, sqrts = sqrts),
        'ttH':ttH(dYu)
    }
    
    return ratios
    
    
################################################################################
def VBF(dCz, Czbx, Czz, Cza, Caa):
    '''
    Return the approximate ratio of vector boson fusion production cross 
    section of the Higgs to the SM prediction as a function of the EFT 
    parameters in the Higgs Basis.
    '''
    return 1. + 2.*dCz - 2.25*Czbx - 0.83*Czz + 0.30*Cza + 0.12*Caa
    
    
def WH(dCz, Czbx, Czz, Cza, Caa, sqrts = 8):
    '''
    Return the approximate ratio of W boson associated production cross section 
    of the Higgs to the SM prediction as a function of the EFT parameters in 
    the Higgs Basis.
    '''
    Azbx = {7:9.26, 8:9.43, 13:10.08}
    Azz = {7:4.35, 8:4.41, 13:4.63}
    Aza = {7:0.81, 8:0.81, 13:0.93}
    Aaa = {7:0.43, 8:0.44, 13:0.48}
    
    return (1. + 2.*dCz - Azbx[sqrts]*Czbx - Azz[sqrts]*Czz 
           + Aza[sqrts]*Cza + Aaa[sqrts]*Caa)
           
    
def ZH(dCz, Czbx, Czz, Cza, Caa, sqrts = 8):
    '''
    Return the approximate ratio of Z boson associated production cross section 
    of the Higgs to the SM prediction as a function of the EFT parameters in the 
    Higgs Basis.
    '''
    Azbx = {7:7.61, 8:7.77, 13:8.24}
    Azz = {7:3.31, 8:3.35, 13:3.47}
    Aza = {7:0.58, 8:0.60, 13:0.65}
    Aaa = {7:0.27, 8:0.28, 13:0.30}
    
    return (1. + 2.*dCz - Azbx[sqrts]*Czbx - Azz[sqrts]*Czz 
           + Aza[sqrts]*Cza + Aaa[sqrts]*Caa)    

def ttH(dYu):
    '''
    Return the approximate ratio of top associated production cross section of 
    the Higgs to the SM prediction as a function of the EFT parameters in the 
    Higgs Basis.
    '''
    return 1. + 2.*dYu
    
