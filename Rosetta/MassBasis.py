from internal import Basis
################################################################################
class MassBasis(Basis.Basis):
    '''
    The main output basis for the Rosetta translation module. Based entirely on 
    the Higgs basis parametrisation without assuming the SU(2)xU(1) preserving 
    relations. It is a general collection of most possible Lorentz and charge 
    conserving structures at dimension 6 excluding those not described in the 
    current incarnation of the Higgs basis. Most importantly, a FeynRules/UFO 
    implementation of the interactions in this model is associated with Rosetta 
    and can be found at https://feynrules.irmp.ucl.ac.be/wiki/EFTmassbasis.
    '''
    
    name ='mass'
    # kinetic terms [Eqn. (3.3)]
    MBxMASS = ['dM']
    # triple gauge couplings [Eqn. (3.6)]
    MBxTGC = ['dKa', 'tKa', 'dG1z', 'dKz', 'tKz', 
              'La', 'tLa', 'Lz', 'tLz', 'C3g', 'tC3g']
    # quartic gauge couplings [Eqn. (3.7)]
    MBxQGC = ['dGw4','dGw2z2','dGw2za', 
              'Lw4', 'Lw2z2', 'Lw2a2', 'Lw2az', 'Lw2za', 'C4g',
              'tLw4', 'tLw2z2', 'tLw2a2', 'tLw2az', 'tLw2za', 'tC4g']
    # single Higgs couplings [Eqn. (3.8)]
    MBxh = ['dCw', 'dCz', 'Cww', 'Cwbx', 
            'Cgg', 'Caa', 'Cza', 'Czz', 'Czbx', 'Cabx',
            'tCww', 'tCgg', 'tCaa', 'tCza', 'tCzz']
    # double Higgs couplings [Eqn. (3.13)]
    MBxhh =  ['dCw2', 'dCz2', 'Cww2', 'Cwbx2', 
            'Cgg2', 'Caa2', 'Cza2', 'Czz2', 'Czbx2', 'Cabx2',
            'tCww2', 'tCgg2', 'tCaa2', 'tCza2', 'tCzz2']
    # Higgs self couplings [Eqn. (3.12)]
    MBxhself = ['dL3', 'dL4']
    # 4-fermion operators
    MBx4F = ['cll1122', 'cpuu3333', 'cll1221']
    
    blocks = {'MBxMASS':MBxMASS, 'MBxTGC':MBxTGC, 'MBxQGC':MBxQGC, 
              'MBxh':MBxh, 'MBxhh':MBxhh, 'MBxhself':MBxhself, 'MBx4F':MBx4F}
    
    flavoured = {
        # Z Vertex Corrections [Eqn (3.4)]
        'MBxdGLze': {'kind':'hermitian', 'domain':'complex', 'cname':'dGLze'},
        'MBxdGRze': {'kind':'hermitian', 'domain':'complex', 'cname':'dGRze'},
        'MBxdGLzv': {'kind':'hermitian', 'domain':'complex', 'cname':'dGLzv'},
        'MBxdGLzu': {'kind':'hermitian', 'domain':'complex', 'cname':'dGLzu'},
        'MBxdGRzu': {'kind':'hermitian', 'domain':'complex', 'cname':'dGRzu'},
        'MBxdGLzd': {'kind':'hermitian', 'domain':'complex', 'cname':'dGLzd'},
        'MBxdGRzd': {'kind':'hermitian', 'domain':'complex', 'cname':'dGRzd'},
        # W Vertex Corrections [Eqn (3.4)]
        'MBxdGLwl': {'kind':'hermitian', 'domain':'complex', 'cname':'dGLwl'},
        'MBxdGLwq': {'kind':'general', 'domain':'complex', 'cname':'dGLwq'},
        'MBxdGRwq': {'kind':'general', 'domain':'complex', 'cname':'dGRwq'},
        # Dipole interactions with single gauge bosons [Eqn. (3.5)]
        'MBxdgu': {'kind':'hermitian', 'domain':'complex', 'cname':'dgu'},
        'MBxdgd': {'kind':'hermitian', 'domain':'complex', 'cname':'dgd'},
        'MBxdau': {'kind':'hermitian', 'domain':'complex', 'cname':'dau'},
        'MBxdad': {'kind':'hermitian', 'domain':'complex', 'cname':'dad'},
        'MBxdae': {'kind':'hermitian', 'domain':'complex', 'cname':'dae'},
        'MBxdzu': {'kind':'hermitian', 'domain':'complex', 'cname':'dzu'},
        'MBxdzd': {'kind':'hermitian', 'domain':'complex', 'cname':'dzd'},
        'MBxdze': {'kind':'hermitian', 'domain':'complex', 'cname':'dze'},
        'MBxdwq': {'kind':'general', 'domain':'complex', 'cname':'dwq'},
        'MBxdwl': {'kind':'general', 'domain':'complex', 'cname':'dwl'},
        'MBxtdgu': {'kind':'hermitian', 'domain':'complex', 'cname':'tdgu'},
        'MBxtdgd': {'kind':'hermitian', 'domain':'complex', 'cname':'tdgd'},
        'MBxtdau': {'kind':'hermitian', 'domain':'complex', 'cname':'tdau'},
        'MBxtdad': {'kind':'hermitian', 'domain':'complex', 'cname':'tdad'},
        'MBxtdae': {'kind':'hermitian', 'domain':'complex', 'cname':'tdae'},
        'MBxtdzu': {'kind':'hermitian', 'domain':'complex', 'cname':'tdzu'},
        'MBxtdzd': {'kind':'hermitian', 'domain':'complex', 'cname':'tdzd'},
        'MBxtdze': {'kind':'hermitian', 'domain':'complex', 'cname':'tdze'},
        'MBxtdwq': {'kind':'general', 'domain':'complex', 'cname':'tdwq'},
        'MBxtdwl': {'kind':'general', 'domain':'complex', 'cname':'tdwl'},
        # single Higgs couplings to fermions [Eqn. (3.8)]
        'MBxdYu': {'kind':'general', 'domain':'real', 'cname':'dYu'},
        'MBxdYd': {'kind':'general', 'domain':'real', 'cname':'dYd'},
        'MBxdYe': {'kind':'general', 'domain':'real', 'cname':'dYe'},
        'MBxSu': {'kind':'general', 'domain':'real', 'cname':'Su' },
        'MBxSd': {'kind':'general', 'domain':'real', 'cname':'Sd' },
        'MBxSe': {'kind':'general', 'domain':'real', 'cname':'Se' },
        # Higgs contact interactions HVff [Eqn. (3.10)]
        'MBxdGLhze': {'kind':'hermitian', 'domain':'complex', 'cname':'dGLhze'},
        'MBxdGRhze': {'kind':'hermitian', 'domain':'complex', 'cname':'dGRhze'},
        'MBxdGLhzv': {'kind':'hermitian', 'domain':'complex', 'cname':'dGLhzv'},
        'MBxdGLhzu': {'kind':'hermitian', 'domain':'complex', 'cname':'dGLhzu'},
        'MBxdGRhzu': {'kind':'hermitian', 'domain':'complex', 'cname':'dGRhzu'},
        'MBxdGLhzd': {'kind':'hermitian', 'domain':'complex', 'cname':'dGLhzd'},
        'MBxdGRhzd': {'kind':'hermitian', 'domain':'complex', 'cname':'dGRhzd'},
        'MBxdGLhwl': {'kind':'hermitian', 'domain':'complex', 'cname':'dGLhwl'},
        'MBxdGLhwq': {'kind':'general', 'domain':'complex', 'cname':'dGLhwq'},
        'MBxdGRhwq': {'kind':'general', 'domain':'complex', 'cname':'dGRhwq'},
        # Dipole interactions with single higgs and gauge boson [Eqn. (3.11)]
        'MBxdhgu': {'kind':'hermitian', 'domain':'complex', 'cname':'dhgu'},
        'MBxdhgd': {'kind':'hermitian', 'domain':'complex', 'cname':'dhgd'},
        'MBxdhau': {'kind':'hermitian', 'domain':'complex', 'cname':'dhau'},
        'MBxdhad': {'kind':'hermitian', 'domain':'complex', 'cname':'dhad'},
        'MBxdhae': {'kind':'hermitian', 'domain':'complex', 'cname':'dhae'},
        'MBxdhzu': {'kind':'hermitian', 'domain':'complex', 'cname':'dhzu'},
        'MBxdhzd': {'kind':'hermitian', 'domain':'complex', 'cname':'dhzd'},
        'MBxdhze': {'kind':'hermitian', 'domain':'complex', 'cname':'dhze'},
        'MBxdhwq': {'kind':'general', 'domain':'complex', 'cname':'dhwq'},
        'MBxdhwl': {'kind':'general', 'domain':'complex', 'cname':'dhwl'},
        'MBxtdhgu': {'kind':'hermitian', 'domain':'complex', 'cname':'tdhgu'},
        'MBxtdhgd': {'kind':'hermitian', 'domain':'complex', 'cname':'tdhgd'},
        'MBxtdhau': {'kind':'hermitian', 'domain':'complex', 'cname':'tdhau'},
        'MBxtdhad': {'kind':'hermitian', 'domain':'complex', 'cname':'tdhad'},
        'MBxtdhae': {'kind':'hermitian', 'domain':'complex', 'cname':'tdhae'},
        'MBxtdhzu': {'kind':'hermitian', 'domain':'complex', 'cname':'tdhzu'},
        'MBxtdhzd': {'kind':'hermitian', 'domain':'complex', 'cname':'tdhzd'},
        'MBxtdhze': {'kind':'hermitian', 'domain':'complex', 'cname':'tdhze'},
        'MBxtdhwq': {'kind':'general', 'domain':'complex', 'cname':'tdhwq'},
        'MBxtdhwl': {'kind':'general', 'domain':'complex', 'cname':'tdhwl'},
        # couplings of two Higgs bosons to fermions [Sec. 3.8]
        'MBxY2u': {'kind':'general', 'domain':'complex', 'cname':'Y2u'},
        'MBxY2d': {'kind':'general', 'domain':'complex', 'cname':'Y2d'},
        'MBxY2e': {'kind':'general', 'domain':'complex', 'cname':'Y2e'}
    }
    # All parameters independent
    independent = ( [c for v in blocks.values() for c in v] + 
                    [c for c in flavoured.keys()] )
    # all other undefined behavoiur inherited from Basis.Basis by default
################################################################################
