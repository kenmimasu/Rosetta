from internal import Basis

import math, re
from math import sqrt
from itertools import combinations_with_replacement as comb
from itertools import product
from internal import PID
########################################################################
flavmat = Basis.flavour_matrix
# Higgs basis class
class HiggsBasis(Basis.Basis):
    ##### declare blocks
    # Kinetic terms
    KIN_ind =['dM']
    HBKIN = KIN_ind
    
    # Z Vertex Corrections
    dGLze = flavmat('dGLze', kind='hermitian', domain='complex')
    dGRze = flavmat('dGRze', kind='hermitian', domain='complex')
    dGLzv = flavmat('dGLzv', kind='hermitian', domain='complex')
    dGLzu = flavmat('dGLzu', kind='hermitian', domain='complex')
    dGRzu = flavmat('dGRzu', kind='hermitian', domain='complex')
    dGLzd = flavmat('dGLzd', kind='hermitian', domain='complex')
    dGRzd = flavmat('dGRzd', kind='hermitian', domain='complex')

    # W Vertex Corrections
    dGLwl = flavmat('dGLwl', kind='hermitian', domain='complex')
    dGLwq = flavmat('dGLwq', kind='hermitian', domain='complex')
    dGRwq = flavmat('dGRwq', kind='general'  , domain='complex')

    VERTEX_ind = dGLze + dGRze+ dGLzu+ dGRzu+ dGLzd+ dGRzd + dGLwl + dGRwq
    VERTEX_dep = dGLzv + dGLwq
    HBVERTEX = VERTEX_ind + VERTEX_dep
    
    # Single Higgs couplings to gauge bosons
    HVV_ind = ['dCz','Cgg','Czz','Caa','Cza','Czbx',
                       'CTgg','CTzz','CTaa','CTza']
    HVV_dep = ['dCw','Cww','CTww','Cwbx','Cabx']
    HBHVV = HVV_ind + HVV_dep
    
    # Single Higgs couplings to fermions
    dYu = flavmat('dYu', kind='general', domain='real')
    dYd = flavmat('dYd', kind='general', domain='real')
    dYe = flavmat('dYe', kind='general', domain='real')
    Su  = flavmat('Su', kind='general', domain='real')
    Sd  = flavmat('Sd', kind='general', domain='real')
    Se  = flavmat('Se', kind='general', domain='real')
    
    HFF_ind = dYu + dYd + dYe + Su + Sd + Se
    HBHFF = HFF_ind
    
    # Higgs contact interactions HVff
    CLze = flavmat('CLze', kind='hermitian', domain='complex')
    CRze = flavmat('CRze', kind='hermitian', domain='complex')
    CLzv = flavmat('CLzv', kind='hermitian', domain='complex')
    CLzu = flavmat('CLzu', kind='hermitian', domain='complex')
    CRzu = flavmat('CRzu', kind='hermitian', domain='complex')
    CLzd = flavmat('CLzd', kind='hermitian', domain='complex')
    CRzd = flavmat('CRzd', kind='hermitian', domain='complex')
    CLwl = flavmat('CLwl', kind='hermitian', domain='complex')
    CLwq = flavmat('CLwq', kind='hermitian', domain='complex')
    CRwq = flavmat('CRwq', kind='general'  , domain='complex')
    
    HVFF_dep = (CLze + CRze + CLzv + CLzu + CRzu 
              + CLzd + CRzd + CLwl + CLwq + CRwq)
    HBHVFF = HVFF_dep
    
    # Triple and quartic gauge couplings [Sec. 3.7]
    
    V3_ind = [ 'Lz','C3g','LTz','CT3g' ]
    V3_dep = [ 'dG1z','dKa','dKz','La','KTa','KTz','LTa' ]
    HBV3 = V3_ind + V3_dep
    
    V4_dep = ['dGw4','dGw2z2','dGw2za']
    HBV4 = V4_dep
    
    D2V4_dep = ['Ldw4', 'Ldzdwzw', 'Ldzdwaw', 
                'Ldadwaw', 'Ldadwzw', 'Ldgg3',
                'LTdw4', 'LTdzdwzw', 'LTdzdwaw',
                'LTdadwaw', 'LTdadwzw', 'LTdgg3']
    HBD2V4 = D2V4_dep
    
    # Couplings of two Higgs bosons [Sec. 3.8]
    H3_ind = ['dL3']
    HBH3 = H3_ind
    
    HHVV_dep = ['Cgg2','CTgg2']
    HBHHVV = HHVV_dep
    
    Y2u = flavmat('Y2u', kind='symmetric', domain='complex')
    Y2d = flavmat('Y2d', kind='symmetric', domain='complex')
    Y2e = flavmat('Y2e', kind='symmetric', domain='complex')
    HHFF_dep = Y2u + Y2d + Y2e
    HBHHFF = HHFF_dep
    
    # 4-fermion operators
    fourfermi_ind = ['cll1122','cpuu3333'] # needed for SILH <-> Warsaw
    fourfermi_dep = ['cll1221'] # affects Gf input
    HB4F = fourfermi_ind + fourfermi_dep
    
    # block structure
    blocks = {'HBKIN':HBKIN,'HBVERTEX':HBVERTEX, 'HBHVV':HBHVV, 'HBHFF':HBHFF,
              'HBHVFF':HBHVFF, 'HBV3':HBV3, 'HBV4':HBV4,'HBD2V4':HBD2V4,
              'HBH3':HBH3, 'HBHHVV':HBHHVV, 'HBHHFF':HBHHFF, 'HB4F':HB4F}  
                
    # Full set of independent coefficients
    independent = (KIN_ind + VERTEX_ind  + HVV_ind + HFF_ind + V3_ind
                  + H3_ind + fourfermi_ind )
                
    # Required inputs/masses             
    required_masses = {1, 2, 3, 4, 5, 6, 11, 12, 13, 14, 15, 16}
    required_inputs = {1, 2, 3, 4} # aEWM1, Gf, aS, MZ
    
    translate_to={'mass'}
    
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

    def calculate_dependent(self):
        self.newname = 'Mass'
        s2w, c2w, ee2, gw2, gp2, MZ, vev, gs2 = self.calculate_inputs() 
        A = self
        MW = MZ*sqrt(c2w)
        # Higgs and EW gauge bosons [Sec 3.4] [eqn (3.11)]
        A['dCw']  = A['dCz']  + A['dM']*4. 
        A['Cww']  = A['Czz']  + A['Cza']*2.*s2w  + A['Caa'] *s2w**2
        A['CTww'] = A['CTzz'] + A['CTza']*2.*s2w + A['CTaa']*s2w**2 
        A['Cwbx'] = (A['Czbx']*gw2 + A['Czz']*gp2 - A['Caa']*ee2*s2w 
                    - A['Cza']*(gw2-gp2)*s2w )/(gw2-gp2)
        A['Cabx'] = (A['Czbx']*2.*gw2 + A['Czz']*(gw2+gp2) 
                    - A['Caa']*ee2 - A['Cza']*(gw2-gp2))/(gw2-gp2)
        
        # Gauge-current and Higgs-gauge-current contact interactions [Sec 3.6]
        for i,j in comb((1,2,3),2):# dependent dgV coeffs [eqn (3.5)]
            indx = '{}{}'.format(i,j)
            if i==j:
                A['dGLzv'+indx] = A['dGLze'+indx] + A['dGLwl'+indx] 
                A['dGLwq'+indx] = A['dGLzu'+indx] - A['dGLzd'+indx]
            else:
                for part in ('Re', 'Im'):
                    tail = indx + part
                    A['dGLzv'+tail] = A['dGLze'+tail] + A['dGLwl'+tail] 
                    A['dGLwq'+tail] = A['dGLzu'+tail] - A['dGLzd'+tail]
     
        # list of all z/w vertex corrections
        for dG in self.HBVERTEX: # 4-point coeffs [eqn (3.18)]
            cvff = dG.replace('dG','C')
            A[cvff] = A[dG]
                
        # Triple gauge couplings [Sec 3.7] [eqn (3.21)] 
        A['dG1z'] = (A['Caa']*ee2*gp2 + A['Cza']*(gw2-gp2)*gp2 
                    - A['Czz']*(gw2+gp2)*gp2 - A['Czbx']*(gw2+gp2)*gw2 
                    )/2./(gw2-gp2)
        A['dKa'] = - (A['Caa']*ee2  + A['Cza']*(gw2-gp2) 
                    - A['Czz']*(gw2+gp2) )*gw2/2./(gw2+gp2)
        A['KTa'] = - ( A['CTaa']*ee2 + A['CTza']*(gw2-gp2) 
                    - A['CTzz']*(gw2+gp2))*gw2/2./(gw2+gp2)
        A['dKz'] = A['dG1z'] - gp2/gw2*A['dKa']
        A['KTz'] = - A['KTa']*gp2/gw2
        A['La'] = A['Lz']
        A['LTa'] = A['LTz']
        
        # Quartic gauge couplings [Sec 3.7] [eqn (3.23)] 
        A['dGw4'] = 2.*c2w*A['dG1z']
        A['dGw2z2'] = 2.*A['dG1z']
        A['dGw2za'] = A['dG1z']
        # two derivative quartic gauge couplings [Sec 3.7] [eqn (3.24)] 
        A['Ldw4'] = -gw2/2./MW**2*A['Lz']
        A['LTdw4'] = -gw2/2./MW**2*A['LTz']
        A['Ldzdw_zw'] = -gw2*c2w/MW**2*A['Lz']
        A['LTdzdw_zw'] = -gw2*c2w/MW**2*A['LTz']
        A['Ldzdw_aw'] = -ee2/MW**2*A['Lz']
        A['LTdzdw_aw'] = -ee2/MW**2*A['LTz']
        A['Ldadw_aw'] = -sqrt(ee2*gw2*c2w)/MW**2*A['Lz']
        A['LTdadw_aw'] = -sqrt(ee2*gw2*c2w)/MW**2*A['LTz']
        A['Ldadw_zw'] = -sqrt(ee2*gw2*c2w)/MW**2*A['Lz']
        A['LTdadw_zw'] = -sqrt(ee2*gw2*c2w)/MW**2*A['LTz']
        A['Ldg_g3'] = 3.*sqrt(gs2)**3/vev**2*A['C3g']
        A['LTdg_g3'] = 3.*sqrt(gs2)**3/vev**2*A['CT3g']
        
        # Couplings of two Higgs bosons [Sec 3.8] [eqn (3.27)]
        def delta(i,j):
            return 1. if i==j else 0.
        A['Cgg2'], A['CTgg2'] = A['Cgg'], A['CTgg']
        for i,j in comb((1,2,3),2):
            for f in ('u','d','e'):
                name = '{}{}{}'.format(f,i,j)
                Yij   = A['dY' + name]
                sinij = A['S' + name] 
                cosij = sqrt(1. - sinij**2)
                A['Y2{}Re'.format(name)] = (3.*Yij*cosij 
                                               - A['dCz']*delta(i,j))
                A['Y2{}Im'.format(name)] = 3.*Yij*sinij
        # 4-fermion operators [Sec. 3.9]
        # [eqn (3.32)]
        A['cll1221'] = 2.*(A['dGLwl11'] + A['dGLwl22'] - 2.*A['dM']) 
        
        self.mass[24] = MW + A['dM']


########################################################################
