################################################################################
from math import sqrt,pi
import os
from tempfile import mkdtemp
import subprocess as sub
from collections import namedtuple

from ...internal.constants import (particle_names, default_masses, 
                                   default_ckm, default_inputs, 
                                   GammaZ, GammaW)

from ...internal.basis import checkers as check
from ...internal import session, SLHA

from . import executable
from .errors import (eHDECAYInterfaceError, eHDECAYImportError, 
                     eHDECAYNegativeWidthError, eHDECAYBrGtOneWarning,
                     eHDECAYBrNegativeWarning)
################################################################################
# required info
masses = {25,3,4,5,6,15,13,24,23} # H, c, b, t, tau, mu, Z, W masses
inputs = {1,2,3} # aEWM1, Gf, aS@MZ
################################################################################
__doc__='''
Interface with eHDECAY program (arXiv:1403.3381) to calculate new Higgs width 
and branching ratio to SM particles. Relies on the existence of a translation of 
ones basis to the SILH basis. The minimal set of SILH coefficients that should 
be translated to are:

    'sH','sT','sW','sB','sHW','sHB','sa','g',
    'RsHe3x3','RsHe2x2','RsHu3x3','RsHu2x2','RsHd3x3','RsHd2x2'
    
Required SM inputs are: 

    'MH','aSMZ','MC','MB','MT','MTAU','MMU','aEWM1','Gf','MZ', 'MW'. 
    
The widths of the W and Z are also looked for in the Decay blocks of the silh 
instance and are set to default PDG values if not found. 
The absolute path to the local eHDECAY directory containing the executable 
should be specified in config.txt as:

eHDECAY_dir     /FULL/PATH/TO/eHDECAY
'''

reference = ('R. Contino et al., Comput.Phys.Commun. 185 (2014) 3412\n'
             'A. Djouadi, J. Kalinowski, M. Spira et al., '
             'Comput.Phys.Commun. 108 (1998) 56 \n')
            
################################################################################
def create_SLHA_block(basis, electroweak=True):
    '''
    Interface Rosetta with eHDECAY to calculate Higgs widths and branching
    fractions and return an SLHA decay block for the Higgs.
    '''
    
    try:
        BRs = run(basis, electroweak=True)
        # BR2 = run(basis, interpolate=True)
    except eHDECAYInterfaceError:
        print e
        return

    sum_BRs = sum([v for k,v in BRs.items() if k is not 'WTOT'])

    # sometimes eHDECAY gives a sum of BRs slightly greater than 1.
    # for now a hacky global rescaling is implemented to deal with this.
    if sum_BRs > 1:
        if abs(sum_BRs - 1.) > 1e-2: # complain if its too wrong
            msg = ('Sum of branching fractions = '
                   '{} ; > 1 by more than 1%'.format(sum_BRs))
            session.warnings.warn(msg, eHDECAYBrGtOneWarning)
        for channel, BR in BRs.iteritems():
            if channel!='WTOT':
                BRs[channel] = BR/sum_BRs

    totalwidth = BRs.pop('WTOT')
    if totalwidth < 0.:
        # session.log('eHDECAY: Negative total Higgs width. Check your EFT inputs.')
        raise eHDECAYNegativeWidthError('eHDECAY: Negative total Higgs width. '
                                        'Check your EFT inputs.')
        return

    hdecays = {}

    # sometimes eHDECAY gives negative branching fractions.
    for channel, BR in BRs.iteritems():

        if BR < 0.:
            n1, n2 = particle_names[channel[0]], particle_names[channel[1]]
            comment = 'H -> {}{}'.format(n1,n2)
            msg = ('eHDECAY: Negative branching fraction '
                  'encountered for {}.'.format(comment))
            session.warnings.warn(msg, eHDECAYBrNegativeWarning)
        #     totalwidth -= BR # adjust total width
        #     continue
        # elif BR == 0.:
        #     continue
        if BR==0.:
            continue
        else:
            hdecays[channel] = BR

    # credit
    preamble = ('# Higgs widths and branching fractions '
                'calculated by eHDECAY.\n')
                
    return SLHA.Decay(25, totalwidth, data=hdecays, preamble=preamble)

def run(basis, electroweak=True, interpolate=False, SM_BRs=None):
    '''
    Translate basis instance to SILH and convert to eHDECAY input. 
    return the dictionary returned by execute().
    Keyword arguments:  
        electroweak - switch for electroweak corrections, IELW
        interpolate - use the numerical formulae given in the eHDECAY paper
        SM_BRs - if interpolate==True provide the SM Higgs branching fractions 
                 to rescale as a dict formatted as {(PID1, PID2):BR,...}.
    '''

    session.cite('eHDECAY', reference)

    # ensure required masses & inputs
    check.masses(basis, masses, message='eHDECAY interface')
    check.sminputs(basis, inputs, message='eHDECAY interface')
    
    # translate to silh instance
    thesilh = basis.translate(target='silh')
    
    # generalise flavour matrices
    thesilh.set_flavor(thesilh.flavor, 'general')
    
    # convert silh to eHDECAY inputs
    inpt = from_silh(thesilh, ew = electroweak) 
    
    if interpolate:
        return interpolated(inpt, electroweak=electroweak, SM_BRs=SM_BRs)
    else:
        return execute(inpt)

def execute(inpt):
    '''
    Run local installation of eHDECAY and return the resulting Higgs width 
    and branching fraction information. Takes a dictionary as input containing 
    the following keys:
        'MH','aSMZ','MC','MB','MT','MTAU','MMU','aEWM1','Gf','MZ', 'MW'.
    Keyword arguments:  
        electroweak - switch for electroweak corrections, IELW
    '''
    # if not os.path.exists(executable):
    #     err = ('Rosetta: could not find eHDECAY ' +
    #            'executable in {}'.format(eHDECAY_dir))
    #     raise eHDECAYInterfaceError(err)
        
    session.cite('eHDECAY', reference)
    
    tmpdir = mkdtemp(prefix='eHDECAY_',dir = os.getcwd())
    
    # write out eHDECAY input file
    with open('{}/ehdecay.in'.format(tmpdir),'w') as infile:
        infile.write( create_input(inpt) )
        
    process = sub.Popen(executable, stdout = sub.PIPE, 
                        stderr = sub.PIPE, cwd = tmpdir)
    out, err = process.communicate()
    
    if err: raise eHDECAYInterfaceError(err)
    
    session.verbose('eHDECAY output:\n{}'.format(out))
    session.drawline()
    # read BRs and total width
    result = read_output(tmpdir)
    
    # clean up temp directory
    sub.call(['cp', '{}/ehdecay.in'.format(tmpdir), '.'])
    sub.call(['rm', '-r', tmpdir])
    
    return result

def SM_BR(basis=None, inputs={}, electroweak=True):
    '''
    Calculate the Higgs width and BRs in the SM using the input parameters 
    specified in a basis instance or the defaults.
    '''
    session.cite('eHDECAY', reference)
    wid ={}
    if basis is not None:
        # ensure required masses & inputs
        check.masses(basis, masses, message='eHDECAY interface')
        mass = {k:nonzero_mass(basis,k) for k in basis.mass}
        
        check.sminputs(basis, inputs, message='eHDECAY interface')
        sminputs = basis.inputs
        
        # get W,Z widths if present
        if basis.card.has_decay(23):
            wid['GAMZ'] = basis.card.decays[23].total
        if basis.card.has_decay(24):
            wid['GAMW'] = basis.card.decays[24].total
    else:
        mass = default_masses
        sminputs = default_inputs
        
        idicts = {1:(sminputs, 'aEWM1'), 2:(sminputs, 'Gf'), 
                  3:(sminputs, 'aSMZ'), 23:(mass, 'MZ'), 
                  24:(mass, 'MW'), 25:(mass, 'MH')}
        for k, (idict, name) in idicts.iteritems():
            try:
                idict[k]=inputs[name]
            except KeyError:
                pass
        
        try:
            mass[24]=inputs['MW']
        except KeyError:
            s2w = (1.- sqrt(1. - (4.*pi*sminputs[3])/
                  (sqrt(2.)*sminputs[2]*masses[23]**2)))/2. # sin^2(theta_W)
            masses[24] = masses[23]*(1. - s2w)
            
    # create SM input for eHDECAY
    coeffs = ['CHbar','CTbar','CWbar','CBbar','CHWbar','CHBbar','Cgambar',
              'Cgbar','Ctaubar','Cmubar','Ctbar','Ccbar','Cbbar','Csbar']
    einpts = {c:0. for c in coeffs}
    
    # EW option
    einpts['IELW'] = 1 if electroweak else 0
    
    # SM inputs & masses
    einpts['MH'] = mass[25]
    einpts['aSMZ'] = sminputs[3]
    einpts['MC'] = mass[4]
    einpts['MB'] = mass[5]
    einpts['MT'] = mass[6]
    einpts['MTAU'] = mass[15]
    einpts['MMU'] = mass[13]
    einpts['aEWM1'] = sminputs[1]
    einpts['Gf'] = sminputs[2]
    einpts['MZ'] = mass[23]
    einpts['MW'] = mass[24]
    
    # W & Z widths
    einpts['GAMZ'] = wid.get('GAMZ', GammaZ)
    einpts['GAMW'] = wid.get('GAMW', GammaW)
    
    return execute(einpts)

def read_output(workdir):
    '''
    Read eHDECAY output files br.eff1 and br.eff2
    '''
    with open('{}/br.eff1'.format(workdir),'r') as br1, open('{}/br.eff2'.format(workdir),'r') as br2:
        br1_dat, br2_dat = br1.readlines()[3], br2.readlines()[3]
        
    MH, BB, TATA, MUMU, SS, CC, TT   = tuple(map(float, br1_dat.split()))
    __, GG, AA,   ZA,   WW, ZZ, WTOT = tuple(map(float, br2_dat.split()))
    
    BR = { (5,-5):BB, (15,-15):TATA, (13,-13):MUMU, 
           (3,-3):SS, (4,-4):CC, (6,-6):TT , (21,21):GG,
           (22,22):AA, (23,22):ZA, (24,-24):WW, (23,23):ZZ, 
           'WTOT':WTOT}
           
    return BR

def from_silh(silh_instance, ew=True):
    '''
    Set up eHDECAY inputs from SILH basis instance.
    '''
    si = silh_instance

    SILH = {}
    mass = {k:nonzero_mass(si,k) for k in masses}
    
    # EW option
    SILH['IELW'] = 1 if ew else 0
    
    # SM inputs & masses
    SILH['MH'] = mass[25]
    SILH['aSMZ'] = si.inputs[3]
    SILH['MC'] = mass[4]
    SILH['MB'] = mass[5]
    SILH['MT'] = mass[6]
    SILH['MTAU'] = mass[15]
    SILH['MMU'] = mass[13]
    SILH['aEWM1'] = si.inputs[1]
    SILH['Gf'] = si.inputs[2]
    SILH['MZ'] = mass[23]
    SILH['MW'] = mass[24]
    
    # W & Z widths
    if si.card.has_decay(23):
        SILH['GAMZ'] = si.card.decays[23].total
    else:
        SILH['GAMZ'] = GammaZ # default
    if si.card.has_decay(24):
        SILH['GAMW'] = si.card.decays[24].total
    else:
        SILH['GAMW'] = GammaW # default
        
    # gauge & Higgs coefficients
    SILH['CHbar'] = si['sH']
    SILH['CTbar'] = si['sT']
    SILH['CWbar'] = si['sW']
    SILH['CBbar'] = si['sB']
    SILH['CHWbar'] = si['sHW']
    SILH['CHBbar'] = si['sHB']
    SILH['Cgambar'] = si['sa']
    SILH['Cgbar'] = si['sg']
    
    # Flavor diagonal Yukawa operators for 2nd & 3rd generations
    SILH['Ctaubar'] = si['SBxHe'][3,3].real
    SILH['Cmubar'] = si['SBxHe'][2,2].real
    SILH['Ctbar'] = si['SBxHu'][3,3].real
    SILH['Ccbar'] = si['SBxHu'][2,2].real
    SILH['Cbbar'] = si['SBxHd'][3,3].real
    SILH['Csbar'] = si['SBxHd'][2,2].real
    
    return SILH

def nonzero_mass(basis, PID):
    '''
    Force non-zero values for all relevant masses. eHDECAY gives nan otherwise.
    '''
    themass = basis.mass[PID]
    if themass == 0.:
        name = particle_names[PID]
        default = default_masses[PID]
        session.log('eHDECAY requires nonzero mass for {}. '.format(name) +
                'Default value of {} GeV used.'.format(default))
        return default
    else:
        return themass


def create_input(inp):
    '''
    Write out input file for eHDECAY.
    '''

    return \
'''SLHAIN   = 0
SLHAOUT  = 0
COUPVAR  = 1
HIGGS    = 0
SM4      = 0
FERMPHOB = 0
MODEL    = 1
TGBET    = 1.D0
MABEG    = {MH}
MAEND    = 1000.D0
NMA      = 1
ALS(MZ)  = {aSMZ}
MSBAR(2) = 0.100D0
MC       = {MC}
MB       = {MB}
MT       = {MT}
MTAU     = {MTAU}
MMUON    = {MMU}
1/ALPHA  = {aEWM1}
GF       = {Gf}
GAMW     = {GAMW}
GAMZ     = {GAMZ}
MZ       = {MZ}
MW       = {MW}
VUS      = 0.2253D0
VCB      = 0.0410D0
VUB/VCB  = 0.0846D0
********************* 4TH GENERATION *************************************
  SCENARIO FOR ELW. CORRECTIONS TO H -> GG (EVERYTHING IN GEV):
  GG_ELW = 1: MTP = 500    MBP = 450    MNUP = 375    MEP = 450
  GG_ELW = 2: MBP = MNUP = MEP = 600    MTP = MBP+50*(1+LOG(M_H/115)/5)

GG_ELW   = 1
MTP      = 500.D0
MBP      = 450.D0
MNUP     = 375.D0
MEP      = 450.D0
**************************************************************************
SUSYSCALE= 1000.D0
MU       = 1000.D0
M2       = 1000.D0
MGLUINO  = 1000.D0
MSL1     = 1000.D0
MER1     = 1000.D0
MQL1     = 1000.D0
MUR1     = 1000.D0
MDR1     = 1000.D0
MSL      = 1000.D0
MER      = 1000.D0
MSQ      = 1000.D0
MUR      = 1000.D0
MDR      = 1000.D0
AL       = 1000.D0
AU       = 1000.D0
AD       = 1000.D0
NNLO (M) = 0
ON-SHELL = 0
ON-SH-WZ = 0
IPOLE    = 0
OFF-SUSY = 0
INDIDEC  = 0
NF-GG    = 5
IGOLD    = 0
MPLANCK  = 2.4D18
MGOLD    = 1.D-13
************** LAGRANGIAN 0 - chiral  1 - SILH  2 - MCHM4/5 **************
LAGPARAM = 1
**** Turn off (0) or on (1) the elw corrections for LAGPARAM = 1 or 2 ****
IELW     = {IELW}
******************* VARIATION OF HIGGS COUPLINGS *************************
CW       = 0D0
CZ       = 0D0
Ctau     = 0D0
Cmu      = 0D0
Ct       = 0D0
Cb       = 0D0
Cc       = 0D0
Cs       = 0D0
Cgaga    = 0D0
Cgg      = 0D0
CZga     = 0D0
CWW      = 0D0
CZZ      = 0D0
CWdW     = 0D0
CZdZ     = 0D0
**************************** SILH Lagrangian *****************************
CHbar    = {CHbar}
CTbar    = {CTbar}
Ctaubar  = {Ctaubar}
Cmubar   = {Cmubar}
Ctbar    = {Ctbar}
Cbbar    = {Cbbar}
Ccbar    = {Ccbar}
Csbar    = {Csbar}
CWbar    = {CWbar}
CBbar    = {CBbar}
CHWbar   = {CHWbar}
CHBbar   = {CHBbar}
Cgambar  = {Cgambar}
Cgbar    = {Cgbar}
******** MCHM4 (fermrepr=1), MCHM5 (fermrepr=2) parametrisation ********
fermrepr = 2
xi       = 0.D0
'''.format(**inp)

def ratios(inpt):
    '''
    Calculates ratio of the Higgs partial widths to the SM predictions  
    according to the numerical formulae given in A. Djouadi, J. Kalinowski, 
    M. Spira et al. accurate to linear order in the Wilson coefficients.
    '''
    # Determine rescaling factors from numerical formulae
    # Derive alpha2, s2w, c2w, aem
    aEM, Gf, MZ = 1./inpt['aEWM1'], inpt['Gf'], inpt['MZ'] # EW inputs 
    
    s2w = (1.- sqrt(1. - (4.*pi*aEM)/(sqrt(2.)*Gf*MZ**2)))/2. # sin^2(theta_W)
    a2 = aEM/s2w # SU(2) coupling squared

    t2w = s2w/(1.-s2w)
    
    cH, cT, cW = inpt['CHbar'], inpt['CTbar'], inpt['CWbar']
    cB, cHW, cHB = inpt['CBbar'], inpt['CHWbar'], inpt['CHBbar'] 
    cc, cs, cb, ct = inpt['Csbar'], inpt['Ccbar'], inpt['Cbbar'], inpt['Ctbar']
    cmu, ctau = inpt['Cmubar'], inpt['Ctaubar']
    cg, cgam = inpt['Cgbar'], inpt['Cgambar']
    
    rscl={}
    rscl[(5,-5)] = 1.- cH - 1.992*cb - 0.0085*ct
    
    rscl[(15,-15)] = 1.- cH - 2.*cmu
    
    rscl[(13,-13)] = 1.- cH - 2.*ctau
    
    rscl[(3,-3)] = 1.- cH - 1.971*cs - 0.029*ct
    
    rscl[(4,-4)] = 1.- cH - 1.985*cc - 0.015*ct
    
    rscl[(6,-6)] = 1.- cH - 2.*ct
    
    rscl[(21,21)] = (1.- cH - 2.12*ct + 0.024*cc + 0.1*cb 
                    + 22.2*cg*4.*pi/sqrt(aEM*a2))
    
    rscl[(22,22)] = (1.- cH + 0.54*ct - 0.003*cc - 0.007*cb - 0.007*ctau 
                    + 5.04*cW - 0.54*cgam*4.*pi/aEM)
    
    rscl[(23,22)] = (1.- cH + 0.12*ct - 5e-4*cc - 0.003*cb - 9e-5*ctau + 4.2*cW  
                    + 0.19*(cHW - cHB + 8.*cgam*s2w)*4.*pi/sqrt(a2*aEM) )
                 
    rscl[(24,-24)] = 1.- cH + 2.2*cW + 3.7*cHW

    rscl[(23,23)] = (1.- cH + 2.*cT+ 2.*(cW + t2w*cB) + 3.*(cHW + t2w*cHB) 
                     - 0.26*cgam)

    return rscl
    
def interpolated(inpt, electroweak=True, SM_BRs=None):
    '''
    Return the new Higgs BRs and Total widths based on the rescaling factors 
    computed in ratio(). SM BRs can be provided through the SM_BRs keyword 
    argument, otherwise they will be calculated using eHDECAY.
    '''
    # calculate SM Higgs Branching ratios
    if SM_BRs is None:
        BRs = SM_BR(inputs = inpt, electroweak=electroweak)
    else:
        BRs = SM_BRs
    
    # get partial width rescaling factors
    rscl = ratios(inpt)
    
    # rescale SM BRs
    SMwid, BRs['WTOT'] = BRs['WTOT'], 0.
    for k,fact in rscl.iteritems():
        PW = BRs[k]*SMwid*fact
        BRs[k] = PW
        BRs['WTOT'] += PW
    
    for k,fact in rscl.iteritems():
        BRs[k] /= BRs['WTOT']

    return BRs

    
    
    
    
    
    