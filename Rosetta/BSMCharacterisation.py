import math
from math import sqrt

from internal import Basis
################################################################################
class BSMCharacterisation(Basis.Basis):
    '''
    The main output basis for the Rosetta translation module. Based entirely on 
    the Higgs basis parametrisation without assuming the SU(2)xU(1) preserving 
    relations. It is a general collection of most possible Lorentz and charge 
    conserving structures at dimension 6 excluding those not described in the 
    current incarnation of the Higgs basis. Most importantly, a FeynRules/UFO 
    implementation of the interactions in this model is associated with Rosetta 
    and can be found at https://feynrules.irmp.ucl.ac.be/wiki/EFTmassbasis.
    '''
    
    name ='bsmc'
    # kinetic terms [Eqn. (3.3)]
    BCxMASS = ['dM']
    # triple gauge couplings [Eqn. (3.6)]
    BCxTGC = ['dKa', 'dKz', 'dG1z', 'La', 'Lz', 'C3g', 
              'tKa', 'tKz', 'tLa', 'tLz', 'tC3g']
    # quartic gauge couplings [Eqn. (3.7)]
    BCxQGC = ['dGw4','dGw2z2','dGw2za', 
              'Lw4', 'Lw2z2', 'Lw2a2', 'Lw2az', 'Lw2za', 'C4g',
              'tLw4', 'tLw2z2', 'tLw2a2', 'tLw2az', 'tLw2za', 'tC4g']
    # single Higgs couplings [Eqn. (3.8)]
    BCxh = ['dCw', 'dCz', 'Cww', 'Cgg', 'Caa', 'Cza', 'Czz', 
            'Cwbx', 'Czbx', 'Cabx',
            'tCww', 'tCgg', 'tCaa', 'tCza', 'tCzz']
    # double Higgs couplings [Eqn. (3.13)]
    BCxhh = ['dCw2', 'dCz2', 'Cww2', 'Cgg2', 'Caa2', 'Cza2', 'Czz2', 
            'Cwbx2', 'Czbx2', 'Cabx2',
            'tCww2', 'tCgg2', 'tCaa2', 'tCza2', 'tCzz2']
    # Higgs self couplings [Eqn. (3.12)]
    BCxhself = ['dL3', 'dL4']
    # 4-fermion operators
    BCx4F = ['cll1122', 'cpuu3333', 'cll1221']
    
    blocks = {'BCxMASS':BCxMASS, 'BCxTGC':BCxTGC, 'BCxQGC':BCxQGC, 
              'BCxh':BCxh, 'BCxhh':BCxhh, 'BCxhself':BCxhself, 'BCx4F':BCx4F}
    
    flavored = {
        # Z Vertex Corrections [Eqn (3.4)]
        'BCxdGLze': {'kind':'hermitian', 'domain':'complex', 'cname':'dGLze'},
        'BCxdGRze': {'kind':'hermitian', 'domain':'complex', 'cname':'dGRze'},
        'BCxdGLzv': {'kind':'hermitian', 'domain':'complex', 'cname':'dGLzv'},
        'BCxdGLzu': {'kind':'hermitian', 'domain':'complex', 'cname':'dGLzu'},
        'BCxdGRzu': {'kind':'hermitian', 'domain':'complex', 'cname':'dGRzu'},
        'BCxdGLzd': {'kind':'hermitian', 'domain':'complex', 'cname':'dGLzd'},
        'BCxdGRzd': {'kind':'hermitian', 'domain':'complex', 'cname':'dGRzd'},
        # W Vertex Corrections [Eqn (3.4)]
        'BCxdGLwl': {'kind':'hermitian', 'domain':'complex', 'cname':'dGLwl'},
        'BCxdGLwq': {'kind':'general', 'domain':'complex', 'cname':'dGLwq'},
        'BCxdGRwq': {'kind':'general', 'domain':'complex', 'cname':'dGRwq'},
        # Dipole interactions with single gauge bosons [Eqn. (3.5)]
        'BCxdgu': {'kind':'hermitian', 'domain':'complex', 'cname':'dgu'},
        'BCxdgd': {'kind':'hermitian', 'domain':'complex', 'cname':'dgd'},
        'BCxdau': {'kind':'hermitian', 'domain':'complex', 'cname':'dau'},
        'BCxdad': {'kind':'hermitian', 'domain':'complex', 'cname':'dad'},
        'BCxdae': {'kind':'hermitian', 'domain':'complex', 'cname':'dae'},
        'BCxdzu': {'kind':'hermitian', 'domain':'complex', 'cname':'dzu'},
        'BCxdzd': {'kind':'hermitian', 'domain':'complex', 'cname':'dzd'},
        'BCxdze': {'kind':'hermitian', 'domain':'complex', 'cname':'dze'},
        'BCxdwu': {'kind':'general', 'domain':'complex', 'cname':'dwu'},
        'BCxdwd': {'kind':'general', 'domain':'complex', 'cname':'dwd'},
        'BCxdwe': {'kind':'general', 'domain':'complex', 'cname':'dwe'},
        'BCxtdgu': {'kind':'hermitian', 'domain':'complex', 'cname':'tdgu'},
        'BCxtdgd': {'kind':'hermitian', 'domain':'complex', 'cname':'tdgd'},
        'BCxtdau': {'kind':'hermitian', 'domain':'complex', 'cname':'tdau'},
        'BCxtdad': {'kind':'hermitian', 'domain':'complex', 'cname':'tdad'},
        'BCxtdae': {'kind':'hermitian', 'domain':'complex', 'cname':'tdae'},
        'BCxtdzu': {'kind':'hermitian', 'domain':'complex', 'cname':'tdzu'},
        'BCxtdzd': {'kind':'hermitian', 'domain':'complex', 'cname':'tdzd'},
        'BCxtdze': {'kind':'hermitian', 'domain':'complex', 'cname':'tdze'},
        # single Higgs couplings to fermions [Eqn. (3.8)]
        'BCxdYu': {'kind':'general', 'domain':'real', 'cname':'dYu'},
        'BCxdYd': {'kind':'general', 'domain':'real', 'cname':'dYd'},
        'BCxdYe': {'kind':'general', 'domain':'real', 'cname':'dYe'},
        'BCxSu': {'kind':'general', 'domain':'real', 'cname':'Su' },
        'BCxSd': {'kind':'general', 'domain':'real', 'cname':'Sd' },
        'BCxSe': {'kind':'general', 'domain':'real', 'cname':'Se' },
        # Higgs contact interactions HVff [Eqn. (3.10)]
        'BCxdGLhze': {'kind':'hermitian', 'domain':'complex', 'cname':'dGLhze'},
        'BCxdGRhze': {'kind':'hermitian', 'domain':'complex', 'cname':'dGRhze'},
        'BCxdGLhzv': {'kind':'hermitian', 'domain':'complex', 'cname':'dGLhzv'},
        'BCxdGLhzu': {'kind':'hermitian', 'domain':'complex', 'cname':'dGLhzu'},
        'BCxdGRhzu': {'kind':'hermitian', 'domain':'complex', 'cname':'dGRhzu'},
        'BCxdGLhzd': {'kind':'hermitian', 'domain':'complex', 'cname':'dGLhzd'},
        'BCxdGRhzd': {'kind':'hermitian', 'domain':'complex', 'cname':'dGRhzd'},
        'BCxdGLhwl': {'kind':'hermitian', 'domain':'complex', 'cname':'dGLhwl'},
        'BCxdGLhwq': {'kind':'general', 'domain':'complex', 'cname':'dGLhwq'},
        'BCxdGRhwq': {'kind':'general', 'domain':'complex', 'cname':'dGRhwq'},
        # Dipole interactions with single higgs and gauge boson [Eqn. (3.11)]
        'BCxdhgu': {'kind':'hermitian', 'domain':'complex', 'cname':'dhgu'},
        'BCxdhgd': {'kind':'hermitian', 'domain':'complex', 'cname':'dhgd'},
        'BCxdhau': {'kind':'hermitian', 'domain':'complex', 'cname':'dhau'},
        'BCxdhad': {'kind':'hermitian', 'domain':'complex', 'cname':'dhad'},
        'BCxdhae': {'kind':'hermitian', 'domain':'complex', 'cname':'dhae'},
        'BCxdhzu': {'kind':'hermitian', 'domain':'complex', 'cname':'dhzu'},
        'BCxdhzd': {'kind':'hermitian', 'domain':'complex', 'cname':'dhzd'},
        'BCxdhze': {'kind':'hermitian', 'domain':'complex', 'cname':'dhze'},
        'BCxdhwu': {'kind':'general', 'domain':'complex', 'cname':'dhwu'},
        'BCxdhwd': {'kind':'general', 'domain':'complex', 'cname':'dhwd'},
        'BCxdhwe': {'kind':'general', 'domain':'complex', 'cname':'dhwe'},
        'BCxtdhgu': {'kind':'hermitian', 'domain':'complex', 'cname':'tdhgu'},
        'BCxtdhgd': {'kind':'hermitian', 'domain':'complex', 'cname':'tdhgd'},
        'BCxtdhau': {'kind':'hermitian', 'domain':'complex', 'cname':'tdhau'},
        'BCxtdhad': {'kind':'hermitian', 'domain':'complex', 'cname':'tdhad'},
        'BCxtdhae': {'kind':'hermitian', 'domain':'complex', 'cname':'tdhae'},
        'BCxtdhzu': {'kind':'hermitian', 'domain':'complex', 'cname':'tdhzu'},
        'BCxtdhzd': {'kind':'hermitian', 'domain':'complex', 'cname':'tdhzd'},
        'BCxtdhze': {'kind':'hermitian', 'domain':'complex', 'cname':'tdhze'},
        # couplings of two Higgs bosons to fermions [Sec. 3.8]
        'BCxY2u': {'kind':'general', 'domain':'complex', 'cname':'Y2u'},
        'BCxY2d': {'kind':'general', 'domain':'complex', 'cname':'Y2d'},
        'BCxY2e': {'kind':'general', 'domain':'complex', 'cname':'Y2e'}
    }
    # All parameters independent
    independent = ( [c for v in blocks.values() for c in v] + 
                    [c for c in flavored.keys()] )

    required_inputs = {1, 2, 4}
    # all other undefined behaviour inherited from Basis.Basis by default
################################################################################
    def modify_inputs(self):
        '''
        W mass modification from dM.
        '''
        ee2 = 4.*math.pi/self.inputs['aEWM1'] # EM coupling squared
        Gf, MZ = self.inputs['Gf'], self.inputs['MZ']
        s2w = (1.- sqrt(1. - ee2/(sqrt(2.)*Gf*MZ**2)))/2. # sin^2(theta_W)
        c2w = (1.-s2w)
        MW = MZ*sqrt(c2w)

        if 24 in self.mass:
            self.mass[24] = MW + self['dM']
        else:
            self.mass.new_entry(24, MW + self['dM'], name = 'MW')
################################################################################
