from Basis import Basis, flavour_matrix
import math, re
from math import sqrt
from itertools import combinations_with_replacement as comb
from itertools import product
from __init__ import PID
########################################################################
# Higgs basis class
class HiggsBasis(Basis):
    # Kinetic terms
    kinetic_ind=['dM']
    
    # Z Vertex Corrections
    dGLze = flavour_matrix('dGLze', kind='hermitian', domain='complex')
    dGRze = flavour_matrix('dGRze', kind='hermitian', domain='complex')
    dGLzv = flavour_matrix('dGLzv', kind='hermitian', domain='complex')
    dGLzu = flavour_matrix('dGLzu', kind='hermitian', domain='complex')
    dGRzu = flavour_matrix('dGRzu', kind='hermitian', domain='complex')
    dGLzd = flavour_matrix('dGLzd', kind='hermitian', domain='complex')
    dGRzd = flavour_matrix('dGRzd', kind='hermitian', domain='complex')
    
    zvertex_ind = dGLze + dGRze+ dGLzu+ dGRzu+ dGLzd+ dGRzd
    zvertex_dep = dGLzv
    
    # W Vertex Corrections
    dGLwl = flavour_matrix('dGLwl', kind='hermitian', domain='complex')
    dGLwq = flavour_matrix('dGLwq', kind='hermitian', domain='complex')
    dGRwq = flavour_matrix('dGRwq', kind='general'  , domain='complex')

    wvertex_ind = dGLwl + dGRwq
    wvertex_dep = dGLwq
    
    # Single Higgs couplings to gauge bosons
    higgs_gauge_ind = ['dCz','Cgg','Czz','Caa','Cza','Czbx',
                       'CTgg','CTzz','CTaa','CTza']
    higgs_gauge_dep = ['dCw','Cww','CTww','Cwbx','Cabx']
    
    # Single Higgs couplings to fermions
    dYu = flavour_matrix('dYu', kind='general', domain='real')
    dYd = flavour_matrix('dYd', kind='general', domain='real')
    dYe = flavour_matrix('dYe', kind='general', domain='real')
    Su  = flavour_matrix('Su', kind='general', domain='real')
    Sd  = flavour_matrix('Sd', kind='general', domain='real')
    Se  = flavour_matrix('Se', kind='general', domain='real')
    
    higgs_ferm_ind = dYu + dYd + dYe + Su + Sd + Se

    # Higgs contact interactions HVff
    CLze = flavour_matrix('CLze', kind='hermitian', domain='complex')
    CRze = flavour_matrix('CRze', kind='hermitian', domain='complex')
    CLzv = flavour_matrix('CLzv', kind='hermitian', domain='complex')
    CLzu = flavour_matrix('CLzu', kind='hermitian', domain='complex')
    CRzu = flavour_matrix('CRzu', kind='hermitian', domain='complex')
    CLzd = flavour_matrix('CLzd', kind='hermitian', domain='complex')
    CRzd = flavour_matrix('CRzd', kind='hermitian', domain='complex')
    CLwl = flavour_matrix('CLwl', kind='hermitian', domain='complex')
    CLwq = flavour_matrix('CLwq', kind='hermitian', domain='complex')
    CRwq = flavour_matrix('CRwq', kind='general'  , domain='complex')
    
    higgs_contact_dep = (CLze + CRze + CLzv + CLzu + CRzu 
                       + CLzd + CRzd + CLwl + CLwq + CRwq)
    
    # Triple and quartic gauge couplings [Sec. 3.7]
    
    triple_quartic_ind = [ 'Lz','C3G','LTz','CT3G' ]
    triple_quartic_dep = [ 'dG1z','dKa','dKz','La','KTa','KTz','LTa' ]
    
    # Couplings of two Higgs bosons [Sec. 3.8]
    two_higgs_ind = ['dL3']
    gg2 = ['Cgg2','CTgg2']
    Y2u = flavour_matrix('Y2u', kind='symmetric', domain='complex')
    Y2d = flavour_matrix('Y2d', kind='symmetric', domain='complex')
    Y2e = flavour_matrix('Y2e', kind='symmetric', domain='complex')
    two_higgs_dep = gg2 + Y2u + Y2d + Y2e
    
    # 4-fermion operators
    fourfermi_ind = ['cuu3333','cpuu3333'] # needed for SILH <-> Warsaw
    fourfermi_dep = ['cll1221'] # affects Gf input
    
    # Full set of independent and dependent coefficients
    independent = (kinetic_ind + zvertex_ind + wvertex_ind 
                 + higgs_gauge_ind + higgs_ferm_ind + triple_quartic_ind 
                 + two_higgs_ind + fourfermi_ind )
                
    dependent   = (zvertex_dep + wvertex_dep + higgs_gauge_dep 
                 + higgs_contact_dep + triple_quartic_dep 
                 + two_higgs_dep + fourfermi_dep)
                 
    # Required inputs/masses             
    required_masses = {1, 2, 3, 4, 5, 6, 11, 12, 13, 14, 15, 16}
    required_inputs = {1, 2, 4} # aEWM1, Gf, MZ

    def calculate_inputs(self): # calculate a few required EW params from aEWM1, Gf, MZ
        ee2 = 4.*math.pi/self.input['aEWM1'] # EM coupling squared
        Gf, MZ = self.input['Gf'], self.input['MZ']
        s2w = (1.- sqrt(1. - ee2/(sqrt(2.)*Gf*MZ**2)))/2.#sin^2(theta_W)
        c2w = (1.-s2w)
        gw2 = ee2/s2w # SU(2) coupling squared
        gp2 = gw2*s2w/c2w # Hypercharge coupling squared
        self.input['s2w'],self.input['ee2'] = s2w, ee2
        self.input['gw2'],self.input['gp2'] = gw2, gp2
        return s2w, ee2, gw2, gp2
         
    def calculate_dependent(self):
        p = self.par_dict
        s2w, ee2, gw2, gp2 = self.calculate_inputs() # EW parameters
        
        # Higgs and EW gauge bosons [Sec 3.4] [eqn (3.11)]
        p['dCw']  = p['dCz']  + p['dM']*4. 
        p['Cww']  = p['Czz']  + p['Cza']*2.*s2w  + p['Caa'] *s2w**2
        p['CTww'] = p['CTzz'] + p['CTza']*2.*s2w + p['CTaa']*s2w**2 
        p['Cwbx'] = (p['Czbx']*gw2 + p['Czz']*gp2 - p['Caa']*ee2*s2w 
                    - p['Cza']*(gw2-gp2)*s2w )/(gw2-gp2)
        p['Cabx'] = (p['Czbx']*2.*gw2 + p['Czz']*(gw2+gp2) 
                    - p['Caa']*ee2 - p['Cza']*(gw2-gp2))/(gw2-gp2)
        
        # Gauge-current and Higgs-gauge-current contact interactions [Sec 3.6]
        for i,j in comb((1,2,3),2):# dependent dgV coeffs [eqn (3.5)]
            ind = '{}{}'.format(i,j)
            if i==j:
                p['dGLzv'+ind] = p['dGLze'+ind] + p['dGLwl'+ind] 
                p['dGLwq'+ind] = p['dGLzu'+ind] - p['dGLzd'+ind]
            else:
                for part in ('_Re', '_Im'):
                    tail = ind + part
                    p['dGLzv'+tail] = p['dGLze'+tail] + p['dGLwl'+tail] 
                    p['dGLwq'+tail] = p['dGLzu'+tail] - p['dGLzd'+tail]
     
        # list of all z/w vertex corrections
        zw_vertex_coeffs = (self.zvertex_ind + self.zvertex_dep 
                          + self.wvertex_ind + self.wvertex_dep )
        for dG in zw_vertex_coeffs: # 4-point coeffs [eqn (3.18)]
            cvff = dG.replace('dG','C')
            p[cvff] = p[dG]
                
        # Triple gauge couplings [Sec 3.7] [eqn (3.21)] 
        p['dG1z'] = (p['Caa']*ee2*gp2 + p['Cza']*(gw2-gp2)*gp2 
                    - p['Czz']*(gw2+gp2)*gp2 - p['Czbx']*(gw2+gp2)*gw2 
                    )/2./(gw2-gp2)
        p['dKa'] = - (p['Caa']*ee2  + p['Cza']*(gw2-gp2) 
                    - p['Czz']*(gw2+gp2) )*gw2/2./(gw2+gp2)
        p['KTa'] = - ( p['CTaa']*ee2 + p['CTza']*(gw2-gp2) 
                    - p['CTzz']*(gw2+gp2))*gw2/2./(gw2+gp2)
        p['dKz'] = p['dG1z'] - gp2/gw2*p['dKa']
        p['KTz'] = - p['KTa']*gp2/gw2
        p['La'] = p['Lz']
        p['LTa'] = p['LTz']
        
        # Couplings of two Higgs bosons [Sec 3.8]
        # [eqn (3.27)]
        def delta(i,j):
            return 1. if i==j else 0.
        p['Cgg2'], p['CTgg2'] = p['Cgg'], p['CTgg']
        for i,j in comb((1,2,3),2):
            for f in ('u','d','e'):
                name = '{}{}{}'.format(f,i,j)
                Yij   = p['dY' + name]
                sinij = p['S' + name] 
                cosij = sqrt(1. - sinij**2)
                p['Y2{}_Re'.format(name)] = (3.*Yij*cosij - p['dCz']*delta(i,j))
                p['Y2{}_Im'.format(name)] = 3.*Yij*sinij
        # 4-fermion operators [Sec. 3.9]
        # [eqn (3.32)]
        p['cll1221'] = 2.*(p['dGLwl11'] + p['dGLwl22'] - 2.*p['dM']) 
            
########################################################################
