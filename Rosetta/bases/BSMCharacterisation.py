import math
from math import sqrt, pi
import re
# from warnings import warn

from . import basis
from ..internal import PID
from ..internal import matrix_mult, matrix_add, matrix_sub, matrix_eq
from ..internal.errors import TranslationWarning
from ..internal import session

################################################################################
class BSMCharacterisation(basis.Basis):
    '''
    The main output basis for the Rosetta translation module. Based entirely on 
    the Higgs basis parametrisation without assuming the SU(2)xU(1) preserving 
    relations. It is a general collection of most possible Lorentz and charge 
    conserving structures at dimension 6 excluding those not described in the 
    current incarnation of the Higgs basis. Most importantly, a FeynRules/UFO 
    implementation of the interactions in this model is associated with Rosetta 
    and can be found at 
    https://feynrules.irmp.ucl.ac.be/wiki/BSMCharacterisation .
    '''
    
    name ='bsmc'
    # kinetic terms 
    BCxMASS = ['dM']
    # triple gauge couplings
    BCxTGC = ['dKa', 'dKz', 'dG1z', 'La', 'Lz', 'C3g', 
              'tKa', 'tKz', 'tLa', 'tLz', 'tC3g']
    # quartic gauge couplings 
    BCxQGC = ['dGw4','dGw2z2','dGw2za', 
              'Lw4', 'Lw2z2', 'Lw2a2', 'Lw2az', 'Lw2za', 'C4g',
              'tLw4', 'tLw2z2', 'tLw2a2', 'tLw2az', 'tLw2za', 'tC4g']
    # single Higgs couplings 
    BCxh = ['dCw', 'dCz', 'Cww', 'Cgg', 'Caa', 'Cza', 'Czz', 
            'Cwbx', 'Czbx', 'Cabx',
            'tCww', 'tCgg', 'tCaa', 'tCza', 'tCzz']
    # double Higgs couplings 
    BCxhh = ['dCw2', 'dCz2', 'Cww2', 'Cgg2', 'Caa2', 'Cza2', 'Czz2', 
            'Cwbx2', 'Czbx2', 'Cabx2',
            'tCww2', 'tCgg2', 'tCaa2', 'tCza2', 'tCzz2']
    # Higgs self couplings 
    BCxhself = ['dL3', 'dL4']
    # 4-fermion operators
    BCx4F = ['cll1122', 'cpuu3333', 'cll1221']
    
    blocks = {'BCxMASS':BCxMASS, 'BCxTGC':BCxTGC, 'BCxQGC':BCxQGC, 
              'BCxh':BCxh, 'BCxhh':BCxhh, 'BCxhself':BCxhself, 'BCx4F':BCx4F}
    
    flavored = {
        # Z Vertex Corrections 
        'BCxdGLze': {'kind':'hermitian', 'domain':'complex', 'cname':'dGLze'},
        'BCxdGRze': {'kind':'hermitian', 'domain':'complex', 'cname':'dGRze'},
        'BCxdGLzv': {'kind':'hermitian', 'domain':'complex', 'cname':'dGLzv'},
        'BCxdGLzu': {'kind':'hermitian', 'domain':'complex', 'cname':'dGLzu'},
        'BCxdGRzu': {'kind':'hermitian', 'domain':'complex', 'cname':'dGRzu'},
        'BCxdGLzd': {'kind':'hermitian', 'domain':'complex', 'cname':'dGLzd'},
        'BCxdGRzd': {'kind':'hermitian', 'domain':'complex', 'cname':'dGRzd'},
        # W Vertex Corrections 
        'BCxdGLwl': {'kind':'hermitian', 'domain':'complex', 'cname':'dGLwl'},
        'BCxdGLwq': {'kind':'general', 'domain':'complex', 'cname':'dGLwq'},
        'BCxdGRwq': {'kind':'general', 'domain':'complex', 'cname':'dGRwq'},
        # Dipole interactions with single gauge bosons 
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
        'BCxdwl': {'kind':'general', 'domain':'complex', 'cname':'dwl'},
        'BCxtdgu': {'kind':'hermitian', 'domain':'complex', 'cname':'tdgu'},
        'BCxtdgd': {'kind':'hermitian', 'domain':'complex', 'cname':'tdgd'},
        'BCxtdau': {'kind':'hermitian', 'domain':'complex', 'cname':'tdau'},
        'BCxtdad': {'kind':'hermitian', 'domain':'complex', 'cname':'tdad'},
        'BCxtdae': {'kind':'hermitian', 'domain':'complex', 'cname':'tdae'},
        'BCxtdzu': {'kind':'hermitian', 'domain':'complex', 'cname':'tdzu'},
        'BCxtdzd': {'kind':'hermitian', 'domain':'complex', 'cname':'tdzd'},
        'BCxtdze': {'kind':'hermitian', 'domain':'complex', 'cname':'tdze'},
        # single Higgs couplings to fermions 
        'BCxdYu': {'kind':'general', 'domain':'real', 'cname':'dYu'},
        'BCxdYd': {'kind':'general', 'domain':'real', 'cname':'dYd'},
        'BCxdYe': {'kind':'general', 'domain':'real', 'cname':'dYe'},
        'BCxSu': {'kind':'general', 'domain':'real', 'cname':'Su' },
        'BCxSd': {'kind':'general', 'domain':'real', 'cname':'Sd' },
        'BCxSe': {'kind':'general', 'domain':'real', 'cname':'Se' },
        # Higgs contact interactions HVff 
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
        # Dipole interactions with single higgs and gauge boson 
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
        'BCxdhwl': {'kind':'general', 'domain':'complex', 'cname':'dhwl'},
        'BCxtdhgu': {'kind':'hermitian', 'domain':'complex', 'cname':'tdhgu'},
        'BCxtdhgd': {'kind':'hermitian', 'domain':'complex', 'cname':'tdhgd'},
        'BCxtdhau': {'kind':'hermitian', 'domain':'complex', 'cname':'tdhau'},
        'BCxtdhad': {'kind':'hermitian', 'domain':'complex', 'cname':'tdhad'},
        'BCxtdhae': {'kind':'hermitian', 'domain':'complex', 'cname':'tdhae'},
        'BCxtdhzu': {'kind':'hermitian', 'domain':'complex', 'cname':'tdhzu'},
        'BCxtdhzd': {'kind':'hermitian', 'domain':'complex', 'cname':'tdhzd'},
        'BCxtdhze': {'kind':'hermitian', 'domain':'complex', 'cname':'tdhze'},
        # couplings of two Higgs bosons to fermions
        'BCxY2u': {'kind':'general', 'domain':'complex', 'cname':'dYu2'},
        'BCxY2d': {'kind':'general', 'domain':'complex', 'cname':'dYd2'},
        'BCxY2e': {'kind':'general', 'domain':'complex', 'cname':'dYe2'}
    }
    # All parameters independent
    independent = ( [c for v in blocks.values() for c in v] + 
                    [c for c in flavored.keys()] )

    required_inputs = {1, 2, 3, 4}
    required_masses = {23, 24, 25}
    # all other undefined behaviour inherited from Basis.Basis by default

################################################################################
    def calculate_inputs(self): # calculate a few required EW params from aEWM1, Gf, MZ
        ee2 = 4.*math.pi/self.inputs['aEWM1'] # EM coupling squared
        gs2 = 4.*math.pi*self.inputs['aS'] # strong coupling squared
        Gf, MZ = self.inputs['Gf'], self.inputs['MZ']
        s2w = (1.- sqrt(1. - ee2/(sqrt(2.)*Gf*MZ**2)))/2. # sin^2(theta_W)
        c2w = (1.-s2w)
        gw2 = ee2/s2w # SU(2) coupling squared
        gp2 = gw2*s2w/c2w # Hypercharge coupling squared
        vev =  2.*MZ*sqrt(c2w/gw2)
        return s2w, c2w, ee2, gw2, gp2, MZ, vev, gs2


    @basis.translation('higgs')        
    def to_higgs(self, instance):
        msg = ('This translation maps an anomalous couplings description to a '
        'dimension-6 EFT and is not one-to-one in general. Relations for '
        'dependent parameters should be satisfied for a fully consistent '
        'translation. Inconsistent relations will be overwritten by the '
        'calculate_dependent() method of HiggsBasis class.')
        session.warnings.warn(msg, TranslationWarning)
        # trivial translation apart from dipole terms
        B = self
        B.fix_matrices()
        H = instance
        dipoles = ['Cdau','Cdzu']
        for k, v in self.iteritems():
            
            is_dipole = re.match(r'Ct{0,1}dh{0,1}[a,z,g][u,d,e]\dx\d', k)
            
            if not is_dipole:
                try:
                    H[k] = v
                except KeyError:
                    print '    ' + k + ' not found in Higgs Basis definition.'
        
        # splitting dipole coefficients
        ii = complex(0.,1.)
        for f in ('u','d','e'):

            glu, tglu = 'BCxdg'+f, 'BCxtdg'+f
            pho, tpho = 'BCxda'+f, 'BCxtda'+f
            zed, tzed = 'BCxdz'+f, 'BCxtdz'+f

            for i,j in B[zed].keys():
                if f in ('u','d'):
                    H['HBxdg'+f][i,j] = (B[glu][i,j]-ii*B[tglu][i,j])/2.
                    H['HBxdg'+f][j,i] = (B[glu][i,j]+ii*B[tglu][i,j])/2.

                H['HBxda'+f][i,j] = (B[pho][i,j]-ii*B[tpho][i,j])/2.
                H['HBxda'+f][j,i] = (B[pho][i,j]+ii*B[tpho][i,j])/2.

                H['HBxdz'+f][i,j] = (B[zed][i,j]-ii*B[tzed][i,j])/2.
                H['HBxdz'+f][j,i] = (B[zed][i,j]+ii*B[tzed][i,j])/2.
        
        
        # Explicitly look for bad values of dependent coefficients in the HB
        temp = H.__class__(dependent=True)
        for coeff in H.independent:
            temp[coeff] = H[coeff]
        temp.mass = B.mass
        temp.inputs = B.inputs
        temp.ckm = B.ckm
        temp.calculate_dependent()
        
        msg = ('dependent coefficient "{}" set to a value ({}) different from '
               'that obtained ({}) with calculate_dependent()').format
               
        for coeff, a in self.iteritems():
            is_dipole = re.match(r'Ct{0,1}dh{0,1}[a,z,g][u,d,e]\dx\d', coeff)
            
            if not is_dipole:
                try:
                    b = temp[coeff]
                    if abs(a-b) > 1e-10:
                        session.warnings.warn(msg(coeff,a,b), TranslationWarning)
                        # sys.exit()
                except KeyError:
                    print coeff
        return H
    
        
    @basis.translation('hc')        
    def to_hc(self, instance):

        s2w, c2w, ee2, gw2, gp2, MZ, vev, gs2 = self.calculate_inputs() 
        MWsq = MZ**2*c2w
        aS, aEM, Gf = self.inputs['aS'], 1./self.inputs['aEWM1'], self.inputs['Gf']
        B = self
        H = instance
        
        # issue warning for non compatible parameter values
        if B['dM']>1e-10:
            msg = ('"dM" parameter should be zero (dCw=dCz) to correctly map '
                   'from to BSMC to HC.')
            warnings.warn(msg, TranslationWarning)
        
        H['cosa'], H['Lambda'] = 1./sqrt(2.), vev
        Lam =  H['Lambda']
        cosa = H['cosa']
        sina = 1./sqrt(2.)
        
        # constants
        C = sqrt( aEM*Gf*MZ**2/(8.*sqrt(2.)*pi) )
        gHZZ, gHWW = 2.*MZ**2/vev, 2.*MWsq/vev
        gHaa, gAaa = 47.*aEM/(18.*pi)/vev, 4.*aEM/(3.*pi)/vev
        gHza, gAza = C*(94*c2w-13.)/(9.*pi)/vev, 2.*C*(8.*c2w-5.)/(3.*pi)/vev
        gHgg, gAgg = -aS/(3.*pi)/vev, aS/(2.*pi)/vev

        # Gauge        
        H['kSM'] = (1. + B['dCz'])*vev*(gw2+gp2)/(2.*cosa*gHZZ)
        # B['dCz'] = 2.*cosa*gHZZ*H['kSM']/(vev*(gw2+gp2))-1.
        
        H['kHaa'] = -B['Caa']*ee2/(cosa*gHaa*vev)
        # B['Caa'] = -cosa*gHaa*H['kHaa']/ee2*vev

        H['kAaa'] = -B['tCaa']*ee2/(sina*gAaa*vev)
        # B['tCaa'] = -sina*gAaa*H['kAaa']/ee2*vev
        
        H['kHza'] = -B['Cza']*sqrt(ee2*(gw2+gp2))/(cosa*gHza*vev)
        # B['Cza'] = -cosa*gHza*H['kHza']*vev/sqrt(ee2*(gw2+gp2))

        H['kAza'] = -B['tCza']*sqrt(ee2*(gw2+gp2))/(sina*gAza*vev)
        # B['tCza'] = -sina*gAza*H['kAza']*vev/sqrt(ee2*(gw2+gp2))

        H['kHgg'] = -B['Cgg']*gs2/(2.*cosa*gHgg*vev)
        # B['Cgg'] = -2.*cosa*gHgg*H['kHgg']/gs2

        H['kAgg'] = -B['tCgg']*gs2/(2.*sina*gAgg*vev)
        # B['tCgg'] = -2.*sina*gAgg*H['kAgg']/gs2

        H['kHzz'] = -B['Czz']*(gw2+gp2)*Lam/(cosa*vev)
        # B['Czz'] = -cosa*H['kHzz']/(gw2+gp2)*(vev/Lam)

        H['kAzz'] = -B['tCzz']*(gw2+gp2)*Lam/(sina*vev)
        # B['tCzz'] = -sina*H['kAzz']/(gw2+gp2)*(vev/Lam)

        H['kHww'] = -B['Cww']*gw2*Lam/(cosa*vev)        
        # B['Cww'] = -cosa*H['kHww']/gw2*(vev/Lam)

        H['kAww'] = -B['tCww']*gw2*Lam/(sina*vev)        
        # B['tCww'] = -sina*H['kAww']/gw2*(vev/Lam)

        H['kHda'] = B['Cabx']*sqrt(gw2*gp2)*Lam/(cosa*vev)
        # B['Cabx'] = cosa*H['kHda']/sqrt(gw2*gp2)*(vev/Lam)

        H['kHdz'] = B['Czbx']*gw2*Lam/(cosa*vev)
        # B['Czbx'] = cosa*H['kHdz']/gw2*(vev/Lam)

        H['kHdwR'] = B['Cwbx']*gw2*Lam/(cosa*vev)
        # B['Cwbx'] = cosa*H['kHdwR']/gw2*(vev/Lam)

        # Yukawa
        H['kHtt'] = (1. - B['BCxdYu'][3,3]*sqrt(1. - B['BCxSu'][3,3]**2))/cosa
        H['kAtt'] = B['BCxdYu'][3,3]*B['BCxSu'][3,3]/sina
        # B['BCxdYu'][3,3] = sqrt((cosa*H['kHtt']-1.)**2 + (H['kAtt']*sina)**2 )
        # B['BCxSu'][3,3] = H['kAtt']*sina/B['BCxdYu'][3,3]
        
        H['kHcc'] = (1. - B['BCxdYu'][2,2]*sqrt(1. - B['BCxSu'][2,2]**2))/cosa
        H['kAcc'] = B['BCxdYu'][2,2]*B['BCxSu'][2,2]/sina
        # B['BCxdYu'][2,2] = sqrt((cosa*H['kHcc']-1.)**2 + (H['kAcc']*sina)**2 )
        # B['BCxSu'][2,2] = H['kAcc']*sina/B['BCxdYu'][2,2]

        H['kHbb'] = (1. - B['BCxdYd'][3,3]*sqrt(1. - B['BCxSd'][3,3]**2))/cosa
        H['kAbb'] = B['BCxdYd'][3,3]*B['BCxSd'][3,3]/sina
        # B['BCxdYd'][3,3] = sqrt((cosa*H['kHbb']-1.)**2 + (H['kAbb']*sina)**2 )
        # B['BCxSd'][3,3] = H['kAbb']*sina/B['BCxdYd'][3,3]
        
        H['kHll'] = (1. - B['BCxdYe'][3,3]*sqrt(1. - B['BCxSe'][3,3]**2))/cosa
        H['kAll'] = B['BCxdYe'][3,3]*B['BCxSe'][3,3]/sina
        # B['BCxdYe'][3,3] = sqrt((cosa*H['kHll']-1.)**2 + (H['kAll']*sina)**2 )
        # B['BCxSe'][3,3] = H['kAll']*sina/B['BCxdYe'][3,3]

        return H

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
