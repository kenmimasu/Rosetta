import __init__ as Rosetta
from math import sqrt
from .. import SILHBasis as SB
import os
from tempfile import mkdtemp
import subprocess as sub
from collections import namedtuple
from settings import eHDECAY_dir
################################################################################
# required info
masses = {25,3,4,5,6,15,13,24,23} # H, c, b, t, tau, mu, Z, W masses
inputs = {1,2,3} # aEWM1, Gf, aS@MZ

# eHDECAY executable
executable = '{}/run'.format(eHDECAY_dir) 
################################################################################
__doc__='''
Interface with eHDECAY program (arXiv:1403.3381) to calculate new Higgs width 
and branching ratio to SM particles. Relies on the eixtence of a translation of 
ones basis to the SILH basis. The minimal set of SILH coefficients that should 
be translated to are:

    'sH','sT','sW','sB','sHW','sHB','sBB','sGG',
    'se33Re','se22Re','su33Re','su22Re','sd33Re','sd22Re'
    
The values of the SILH coefficients are then 
rescaled from the convention of the LHXSWG draft to match those of eHDECAY. 
Required SM inputs are: 

    'MH','aSMZ','MC','MB','MT','MTAU','MMU','aEWM1','Gf','MZ', 'MW'. 
    
The widths of the W and Z are also looked for in the Decay blocks of the silh 
instance and are set to default PDG values if not found. 
The absolute path to the local eHDECAY directory containing the executable 
should be specified in config.txt as:

eHDECAY_dir     /PATH/TO/eHDECAY
 
'''
################################################################################
def run(basis, electroweak=True):
    '''
    Run local installation of eHDECAY and return the resulting Higgs width 
    and branching fraction information.
    Keyword arguments:  
        electroweak - switch for electroweak corrections, IELW
    '''
    if not os.path.exists(executable):
        print ('Rosetta: could not find eHDECAY executable in {}'.format(
        eHDECAY_dir
        ))
    
    print ('########## eHDECAY ##########\n'
           'If you use this feature, please cite:\n'
           'R. Contino et al., Comput.Phys.Commun. 185 (2014) 3412\n'
           'A. Djouadi, J. Kalinowski, M. Spira et al., '
           'Comput.Phys.Commun. 108 (1998) 56 \n')
    
    # ensure required masses & inputs
    basis.check_masses(masses, message='eHDECAY interface')
    basis.check_sminputs(inputs, message='eHDECAY interface')
    
    # translate to silh instance
    thesilh = basis.translate(target='silh',verbose=False)

    thesilh.set_flavor(thesilh.flavor, 'general')
    inp = from_silh(thesilh, ew = electroweak) 
        
    # create temporary directory
    tmpdir = mkdtemp(prefix='eHDECAY_',dir = os.getcwd())
    
    # write out eHDECAY input file
    with open('{}/ehdecay.in'.format(tmpdir),'w') as infile:
        infile.write( create_input(inp) )
        
    process = sub.Popen(executable, stdout = sub.PIPE, 
                        stderr = sub.PIPE, cwd = tmpdir)
    out, err = process.communicate()
    
    if err: 
        raise RuntimeError('eHDECAY error: {}'.format(err))
    print 'eHDECAY output:\n{}'.format(out)
    
    # read BRs and total width
    result = read_output(tmpdir)
    
    # clean up temp directory
    sub.call(['cp','{}/ehdecay.in'.format(tmpdir),'.'])
    sub.call(['rm','-r',tmpdir])
    
    return result
    
def from_silh(silh_instance, ew=True):
    '''
    Rescaling of the relevant SILH parameters in the LHCXSWG convention to that 
    of the eHDECAY publication (arXiv:1403.3381).
    '''
    si = silh_instance
    s2w, c2w, ee2, gw2, gp2, MZ, vev, gs2 = silh_instance.calculate_inputs()
    g=sqrt(gw2)
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
    if si.card.has_decay(23):
        SILH['GAMZ'] = si.card.decays[23].total
    else:
        SILH['GAMZ'] = 2.4952 # default
    if si.card.has_decay(24):
        SILH['GAMW'] = si.card.decays[24].total
    else:
        SILH['GAMW'] = 2.085 # default
    # gauge & Higgs coefficients
    SILH['CHbar'] = 2.*si['sH']
    SILH['CTbar'] = 2.*si['sT']
    SILH['CWbar'] = gw2/4.*si['sW']
    SILH['CBbar'] = gw2/4.*si['sB']
    SILH['CHWbar'] = gw2/4.*si['sHW']
    SILH['CHBbar'] = gw2/4.*si['sHB']
    SILH['Cgambar'] = gw2/16.*si['sBB']
    SILH['Cgbar'] = gw2/16.*si['sGG']
    
    SILH['Ctaubar'] = sqrt(2.)*si['SBxe'][3,3].real
    SILH['Cmubar'] = sqrt(2.)*si['SBxe'][2,2].real
    SILH['Ctbar'] = sqrt(2.)*si['SBxu'][3,3].real
    SILH['Ccbar'] = sqrt(2.)*si['SBxu'][2,2].real
    SILH['Cbbar'] = sqrt(2.)*si['SBxd'][3,3].real
    SILH['Csbar'] = sqrt(2.)*si['SBxd'][2,2].real
    
    return SILH

def nonzero_mass(basis,PID):
    '''
    Force non-zero values for all relevant masses. eHDECAY gives nan otherwise.
    '''
    themass = basis.mass[PID]
    if themass == 0.:
        name = Rosetta.particle_names[PID]
        default = Rosetta.default_masses[PID]
        print ('eHDECAY requires nonzero mass for {}. '.format(name) +
                'Default value of {} GeV used.'.format(default))
        return default
    else:
        return themass
    
def read_output(workdir):
    '''
    Read eHDECAY output files br.eff1 and br.eff2
    '''
    with open('{}/br.eff1'.format(workdir),'r') as br1, open('{}/br.eff2'.format(workdir),'r') as br2:
        br1_dat, br2_dat = br1.readlines()[3], br2.readlines()[3]
        
    MH, BB, TATA, MUMU, SS, CC, TT   = tuple(map(float, br1_dat.split()))
    __, GG, AA,   ZA,   WW, ZZ, WTOT = tuple(map(float, br2_dat.split()))
    
    BR = { (5,5):BB, (15,15):TATA, (13,13):MUMU, 
           (3,3):SS, (4,4):CC, (6,6):TT , (21,21):GG,
           (22,22):AA, (23,22):ZA, (24,24):WW, (23,23):ZZ, 
           'WTOT':WTOT}
           
    return BR


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
