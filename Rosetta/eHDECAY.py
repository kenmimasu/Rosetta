from __init__ import eHDECAY_dir
import os
from tempfile import mkdtemp
import subprocess as sub
from collections import namedtuple
####################################################################################################
SM_inputs   = ['MH','aSMZ','MC','MB','MT','MTAU',
               'MMU','aEWM1','Gf','MZ','MW','IELW']
               
SILH_inputs = ['CHbar','CTbar','Ctaubar','Cmubar',
               'Ctbar','Cbbar','Ccbar','Csbar','CWbar',
               'CBbar','CHWbar','CHBbar','Cgambar','Cgbar']
               
executable = '{}/run'.format(eHDECAY_dir) # eHDECAY executable

SILH = namedtuple('SILH', SM_inputs+SILH_inputs) # Required inputs for eHDECAY    
####################################################################################################
__doc__='''
Interface with eHDECAY program (arXiv:1403.3381) to calculate new Higgs width and branching ratio
to SM particles. Currently takes inputs via a dictionary of parameter:value pairs to be returned by
the eHDECAY_inputs() function of a basis class.
These should contain the following SM inputs:
{}
They should also contain the values of the following coefficients of the SILH basis:
{}
This module will create a temporary directory to write out the input and output files, call the 
eHDECAY exectuable and store the output so that it may be written to the new parameter card. 
'''.format(', '.join(SM_inputs),', '.join(SILH_inputs))
####################################################################################################

def eHDECAY(basis):
    print 'Running eHDECAY'
    input_dict = basis.eHDECAY_inputs()
    print input_dict
    inp = SILH(**input_dict)
    # create temporary directory
    tmpdir = mkdtemp(prefix='eHDECAY_',dir = os.getcwd())
    # write out eHDECAY input file
    with open('{}/ehdecay.in'.format(tmpdir),'w') as infile:
        infile.write( create_input(inp) )
    process = sub.Popen(executable, stdout = sub.PIPE, stderr = sub.PIPE, cwd = tmpdir)
    out, err = process.communicate()
    if err: 
        raise RuntimeError('eHDECAY error: {}'.format(err))
    print 'eHDECAY output:\n{}'.format(out)
    # read BRs and total width
    result = read_output(tmpdir)
    # clean up temp directory
    sub.call(['rm','-r',tmpdir])
    return result

def read_output(workdir):
    with open('{}/br.eff1'.format(workdir),'r') as br1, open('{}/br.eff2'.format(workdir),'r') as br2:
        br1_dat, br2_dat = br1.readlines()[3], br2.readlines()[3]
    MH, BB, TATA, MUMU, SS, CC, TT   = tuple(map(float,br1_dat.split()))
    __, GG, AA,   ZA,   WW, ZZ, WTOT = tuple(map(float,br2_dat.split()))
    BR = { 'bb':BB, 'tata':TATA, 'mumu':MUMU, 
           'ss':SS, 'cc':CC, 'tt':TT , 'gg':GG,
           'aa':AA, 'Za':ZA, 'WW':WW, 'ZZ':ZZ, 
           'WTOT':WTOT}
    return BR

def create_input(inp):
    values = inp._asdict().values()
    return \
'''SLHAIN   = 0
SLHAOUT  = 0
COUPVAR  = 1
HIGGS    = 0
SM4      = 0
FERMPHOB = 0
MODEL    = 1
TGBET    = 1.D0
MABEG    = {}
MAEND    = 1000.D0
NMA      = 1
ALS(MZ)  = {}
MSBAR(2) = 0.100D0
MC       = {}
MB       = {}
MT       = {}
MTAU     = {}
MMUON    = {}
1/ALPHA  = {}
GF       = {}
GAMW     = 2.08856D0
GAMZ     = 2.49581D0
MZ       = {}
MW       = {}
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
IELW     = {}
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
CHbar    = {}
CTbar    = {}
Ctaubar  = {}
Cmubar   = {}
Ctbar    = {}
Cbbar    = {}
Ccbar    = {}
Csbar    = {}
CWbar    = {}
CBbar    = {}
CHWbar   = {}
CHBbar   = {}
Cgambar  = {}
Cgbar    = {}
******** MCHM4 (fermrepr=1), MCHM5 (fermrepr=2) parametrisation ********
fermrepr = 2
xi       = 0.D0
'''.format(*values)



