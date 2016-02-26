from math import pi
import os

from ...internal.basis import checkers as check
from loopfunctions import Af

from . import use_eHDECAY
if use_eHDECAY:
    from . import SM_BR

################################################################################
# required info
masses = {25,4,5,6,15} # H, b, t masses
################################################################################
def decay(basis, electroweak=True, SM_BRs=None, ratio=False):
    '''
    Return the new Higgs BRs and Total widths based on the rescaling factors 
    computed in ratio(). SM BRs can be provided through the SM_BRs keyword 
    argument, otherwise they will be calculated using eHDECAY.
    Arguments:
        basis - Rosetta.internal.Basis instance
    Options:
        SM_BRs - Provide dict for SM Higgs branching fractions and total width 
        electroweak - if SM_BRs is None and SignalStrengths.use_eHDECAY is True 
                      the options is fed into the eHDECAY interface
    '''
    # get partial width rescaling factors
    rscl = partial_width_ratios(basis)

    # calculate SM Higgs Branching ratios
    if SM_BRs is None:
        if use_eHDECAY:
            BRs = SM_BR(basis=basis, electroweak=electroweak)
        else: 
            BRs = get_BR(basis.mass[25])
    else:
        BRs = SM_BRs
    
    # SM total width
    SMwid = BRs['WTOT']
    # Compute total width rescaling factor
    rscl['WTOT'] = 0.
    for k, v in rscl.iteritems():
        rscl['WTOT'] += v*BRs[k]
    
    if not ratio:
    # Compute new width, branching fractions & return
        BRs['WTOT'] *= rscl['WTOT']
        for k,fact in rscl.iteritems():
            PW = BRs[k]*SMwid*fact
            BRs[k] = PW/BRs['WTOT']
        return BRs
    # Return rescaling factors as they are
    else:
        return rscl
    
def partial_width_ratios(basis):
    '''
    Calculate ratios of all Higgs partial widths w.r.t the SM prediction and 
    returns a dictionary of rescaling factors.
    Arguments:
        basis - Rosetta.internal.basis.Basis instance
    '''
    bsmc = basis.translate(target='bsmc')
    check.masses(bsmc, masses, message='Signal strength calculation')

    # Store relevant coefficients
    MH, MT, MB = bsmc.mass[25], bsmc.mass[6], bsmc.mass[5]
    dCz, Czbx, Czz = bsmc['dCz'], bsmc['Czbx'], bsmc['Czz']
    dCw, Cwbx, Cww = bsmc['dCw'], bsmc['Cwbx'], bsmc['Cww']
    Cgg, Cza, Caa, Cabx = bsmc['Cgg'], bsmc['Cza'], bsmc['Caa'], bsmc['Cabx']
    
    if basis.flavor == 'universal':
        dYt, dYb = bsmc['BCxdYu'][1,1], bsmc['BCxdYd'][1,1]
        dYc, dYs = dYt, dYb
        dYta = bsmc['BCxdYe'][1,1]
        dYmu = dYta
    else: 
        dYt, dYb = bsmc['BCxdYu'][3,3], bsmc['BCxdYd'][3,3]
        dYc, dYs = bsmc['BCxdYu'][2,2], bsmc['BCxdYd'][2,2]
        dYta, dYmu = bsmc['BCxdYe'][3,3], bsmc['BCxdYe'][2,2]    
    
    ratios = { 
        (3,-3):Hff(dYs), 
        (4,-4):Hff(dYc), 
        (5,-5):Hff(dYb), 
        (6,-6):Hff(dYt),
        (15,-15):Hff(dYta), 
        (13,-13):Hff(dYmu), 
        (21,21):Hgg(MH, MT, MB, Cgg, dYt, dYb),
        (22,22):Haa(Caa), 
        (23,22):Hza(Cza), 
        (24,-24):H2l2v(dCw, Cwbx, Cww),
        (23,23):H4l(dCz, Czbx, Czz, Cza, Caa, Cabx)
    }

    return ratios
################################################################################
def Hgg(MH, MT, MB, Cgg, dYu, dYd):
    '''
    Return the approximate ratio of gluon-gluon partial width of the Higgs to 
    the SM prediction as a function of the EFT parameters in the Higgs/BSMC 
    Basis. Also used in SignalStrengths.production module to compute the 
    gluon-gluon fusion production cross section ratio w.r.t the SM.
    '''
    At, Ab = 3.*Af(MH**2/4./MT**2),  3.*Af(MH**2/4./MB**2)
    
    chat = Cgg + 1./(12.*pi**2)*( At*dYu + Ab*dYd )
    cSM =  1./(12.*pi**2)*( At + Ab )

    lin = 1. + 2.*(chat/cSM).real
    quad = abs( (1. + chat/cSM) )**2
    
    return lin

def Haa(Caa):
    '''
    Return the approximate ratio of gamma-gamma partial width of the Higgs to 
    the SM prediction as a function of the EFT parameters in the Higgs/BSMC 
    Basis. 
    '''
    cSM = -8.3e-2

    lin = 1. + 2.*Caa/cSM
    quad = (1. + Caa/cSM)**2
    
    return lin
    
def Hza(Cza):
    '''
    Return the approximate ratio of Z-gamma partial width of the Higgs to 
    the SM prediction as a function of the EFT parameters in the Higgs/BSMC 
    Basis. 
    '''
    cSM = -5.9e-2

    lin = 1. + 2.*Cza/cSM
    quad = (1. + Cza/cSM)**2
    
    return lin
    
    
def Hff(dYf):
    '''
    Return the approximate ratio of fermion-antifermion partial width of the 
    Higgs to the SM prediction as a function of the EFT parameters in the 
    Higgs/BSMC Basis. 
    '''
    lin = 1. + 2.*dYf
    quad = 1. + 2.*dYf * dYf**2
    
    return lin

def H2l2v(dCw, Cwbx, Cww):
    '''
    Return the approximate ratio of 2 lepton - 2 neutrino (WW) partial width of 
    the Higgs to the SM prediction as a function of the EFT parameters in the 
    Higgs/BSMC Basis. 
    '''
    # return 1 + 2.*dCz + 0.67*Czbx + 0.05*Czz - 0.17*Cza - 0.05*Caa
    lin = 1 + 2.*dCw + 0.46074*Cwbx - 0.153652*Cww
    
    return lin
    
def H2e2mu(dCz, Czbx, Czz, Cza, Caa, Cabx):
    '''
    Return the approximate ratio of 2 electron - 2 muon (ZZ) partial width of 
    the Higgs to the SM prediction as a function of the EFT parameters in the 
    Higgs/BSMC Basis. 
    '''
    # return 1 + 2.*dCz + 0.35*Czbx + 0.19*Czz - 0.09*Cza - 0.01*Caa
    lin = (1 + 2.*dCz + 0.40640*Czbx - 0.14865*Czz - 0.06883*Cza - 0.000366*Caa 
          - 0.02108*Cabx)

    return lin

def H4l(dCz, Czbx, Czz, Cza, Caa, Cabx):
    '''
    Return the approximate ratio of 4 lepton (ZZ) partial width of the 
    Higgs to the SM prediction as a function of the EFT parameters in the 
    Higgs/BSMC Basis. 
    '''
    
    lin = (1 + 2.*dCz + 0.398148*Czbx - 0.144664*Czz + 0.060691*Cza 
          - 0.0226033*Cabx - 0.012563*Caa)

    return lin
    # return 1 + 2.*dCz + 0.32*Czbx + 0.19*Czz - 0.08*Cza - 0.02*Caa
    
################################################################################
def get_datum(MH, file):
    '''
    Read line from a tabulated data file, "file", from the LHCHXSWG 
    corresponding to a Higgs mass of MH. Information is linearly interpolated 
    between the two closest available mass points.
    '''
    last_data = [1e10]
    for line in open(file, 'r'):
        try:
            data = [float(x) for x in line.split()]
            MHi = data[0]
        except ValueError:
            continue
                
        if MHi == MH:  
            break

        if MHi > MH:
            # linear interpolation
            dMH0 = abs(MHi-last_data[0])
            dMH= abs(MHi-MH)
            data = [x + (y-x)*dMH/dMH0 for x,y in zip(last_data, data)]
            break
            
        last_data= data
    
    BRs = [d for i,d in enumerate(data[1:]) if i%3==0 ]
    
    return BRs
            
    
def get_BR(MH):
    '''
    Get SM Higgs branching fractions and total width for a given Higgs mass, MH.
    Data taken from tabulated values provided by the LHCHXSWG at:
    
    https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageBR2014
    
    Information is linearly interpolated between the two closest available mass 
    points.
    '''
    basedir = os.path.dirname(__file__)
    
    BRff = get_datum(MH, '{}/BR/ff.dat'.format(basedir))
    BB, TATA, MUMU, CC, SS, TT = BRff
    
    BRVV = get_datum(MH, '{}/BR/VV.dat'.format(basedir))
    GG, AA, ZA, WW, ZZ, WTOT = BRVV

    return { 
        (5,-5):BB, (15,-15):TATA, (13,-13):MUMU, 
        (3,-3):SS, (4,-4):CC, (6,-6):TT , (21,21):GG,
        (22,22):AA, (23,22):ZA, (24,-24):WW, (23,23):ZZ, 
        'WTOT':WTOT
    }


