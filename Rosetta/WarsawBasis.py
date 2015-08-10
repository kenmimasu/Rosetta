from internal import Basis
import math
from math import sqrt
from itertools import combinations_with_replacement as comb
from itertools import product
from internal import PID
from internal.matrices import matrix_mult, matrix_add, matrix_sub, matrix_eq
################################################################################
class WarsawBasis(Basis.Basis):
    '''
    Basis class for Rosetta based on [Grzadkowski et al., JHEP 1010 (2010) 085]. 
    Part of the three intrinsic basis implementations in Rosetta along with the 
    Higgs and SILH bases. The exact list of operators included as well as the 
    equations for the translation to the Higgs basis can be found in the HXSWG 
    note, for which all references to tables and equations are in this 
    implementation. Table 1 of the note shows full list of operators. The block 
    structure of the basis implementation maps to this table.
    '''
    
    name = 'warsaw'
    ###### 
    # [Tab. 1]
    WBV2H2 = ['cGG','tcGG','cWW','tcWW','cBB','tcBB','cWB','tcWB']
    
    WBH4D2 = ['cH','cT']
    
    WBH6 = ['c6H']
    
    WBV3D3 = ['c3W','c3G','tc3W','tc3G']
    
    # affects Gf input, Warsaw <-> SILH translation
    WB4F = ['cll1221','cll1122','cpuu3333'] 
    ######
    blocks = {'WBxV2H2':WBV2H2, 'WBxH4D2':WBH4D2, 'WBxH6':WBH6, 
              'WBxV3D3':WBV3D3, 'WBx4F':WB4F}
              
    flavoured = {
        'WBxu': {'cname':'cu', 'kind':'general', 'domain':'complex'},
        'WBxd': {'cname':'cd', 'kind':'general', 'domain':'complex'},
        'WBxe': {'cname':'ce', 'kind':'general', 'domain':'complex'},
        'WBxHl': {'cname':'cHl' , 'kind':'hermitian', 'domain':'complex'},
        'WBxHpl': {'cname':'cHpl', 'kind':'hermitian', 'domain':'complex'},
        'WBxHe': {'cname':'cHe' , 'kind':'hermitian', 'domain':'complex'},
        'WBxHq': {'cname':'cHq' , 'kind':'hermitian', 'domain':'complex'},
        'WBxHpq': {'cname':'cHpq', 'kind':'hermitian', 'domain':'complex'},
        'WBxHu': {'cname':'cHu' , 'kind':'hermitian', 'domain':'complex'},
        'WBxHd' : {'cname':'cHd' , 'kind':'hermitian', 'domain':'complex'},
        'WBxHud' : {'cname':'cHud', 'kind':'general', 'domain':'complex'}
    }
    
    independent = [c for c in blocks.keys()] + [c for c in flavoured.keys()]
    
    required_masses = set([y for x in PID.values() for y in x.values()])
    required_inputs = {1, 2, 3, 4, 25} # aEWM1, Gf, MZ, MH, aS(MZ)
        
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

    @Basis.translation('mass')        
    def to_mass(self, instance):
        s2w, c2w, ee2, gw2, gp2, MZ, vev, gs2 = self.calculate_inputs() 
        MH = self.mass[25]
        MW = MZ*sqrt(c2w)
        
        W = self
        M = instance
        
        def delta(i,j):
            return 1. if i==j else 0.
                    
        def f(T3,Q,i,j): # [eqn (4.11)]
                return delta(i,j)*( -Q*W['cWB']*gw2*gp2/(gw2-gp2)
                                   + (W['cT']-dv)*(T3 + Q*gp2/(gw2-gp2)))

        # Higgs vev shift [eqn (4.8)]
        dv = (W['WBxHpl'][1,1].real + W['WBxHpl'][2,2].real)/2.-W['cll1221']/4.
        
        # W mass shift [eqn (4.9)]
        M['dM'] = ( gw2*W['cT'] - gp2*gw2*W['cWB']-gp2*dv )/(gw2-gp2)
        
        # W/Z chiral coupling deviations
        for i,j in W['WBxHpl'].keys():
            M['MBxdGLzv'][i,j] = (1./2.*W['WBxHpl'][i,j] 
                                - 1./2.*W['WBxHl'][i,j] + f(1./2.,0.,i,j))
            M['MBxdGLze'][i,j] = (-1./2.*W['WBxHpl'][i,j] 
                                - 1./2.*W['WBxHl'][i,j] + f(-1./2.,-1.,i,j))
            M['MBxdGRze'][i,j] = (- 1./2.*W['WBxHe'][i,j] + f(0.,-1.,i,j))
            M['MBxdGLzu'][i,j] = (1./2.*W['WBxHpq'][i,j] 
                                - 1./2.*W['WBxHq'][i,j] + f(1./2.,2./3.,i,j))
            M['MBxdGLzd'][i,j] = (-1./2.*W['WBxHpq'][i,j] 
                                - 1./2.*W['WBxHq'][i,j] + f(-1./2.,-1./3.,i,j))
            M['MBxdGRzu'][i,j] = (- 1./2.*W['WBxHu'][i,j] + f(0.,2./3.,i,j))
            M['MBxdGRzd'][i,j] = (- 1./2.*W['WBxHd'][i,j] + f(0.,-1./3.,i,j))
            M['MBxdGLwl'][i,j] = (W['WBxHpl'][i,j] + f(1./2.,0.,i,j) 
                                - f(-1./2.,-1.,i,j))

        for k,v in W['WBxHud'].iteritems():
            M['MBxdGRwq'][k] = -v/2.

        matrix_sub(matrix_mult(M['MBxdGLzu'], W.ckm),
                   matrix_mult(W.ckm, M['MBxdGLzd']),
                   M['MBxdGLwq'])
        
        vertex = ['MBxdGLze', 'MBxdGRze', 'MBxdGLzv', 'MBxdGLzu', 'MBxdGRzu', 
                  'MBxdGLzd', 'MBxdGRzd', 'MBxdGLwl', 'MBxdGLwq', 'MBxdGRwq']

        # HVFF coeffs [Eqns (3.10)]
        for dG in vertex: 
            dGh = dG[:-2]+'h'+dG[-2:]
            matrix_eq(M[dG], M[dGh])

        # Higgs couplings to W/Z [eqn (4.14)]
        M['dCw'] = ( -W['cH']*(gw2-gp2) - W['cWB']*4.*gw2*gp2 
                    + W['cT']*4.*gw2 - dv*(3.*gw2+gp2) )/(gw2-gp2)
        M['dCz'] = -W['cH'] - 3.*dv 
        
        # Two derivative field strength interactions [eqn (4.15)]
        M['Cgg']  = W['cGG'] 
        M['Caa']  = W['cWW'] + W['cBB'] - 4.*W['cWB']
        M['Czz']  = ( gw2**2*W['cWW'] + gp2**2*W['cBB'] 
                     + 4.*gw2*gp2*W['cWB'] )/(gw2+gp2)**2
        M['Czbx'] =  -(2./gw2)*(W['cT'] - dv) 
        M['Cza']  = ( gw2*W['cWW'] - gp2*W['cBB'] 
                      - 2.*(gw2-gp2)*W['cWB'] )/(gw2+gp2)
        M['Cabx'] =  2./(gw2-gp2)*((gw2+gp2)*W['cWB'] - 2.*W['cT'] + 2.*dv) 
        M['Cww']  =  W['cWW']
        M['Cwbx'] =  2./(gw2-gp2)*(gp2*W['cWB'] - W['cT'] + dv) 
        M['tCgg'] = W['tcGG'] 
        M['tCaa'] =  W['tcWW'] + W['tcBB'] - 4.*W['tcWB']
        M['tCzz'] = ( gw2**2*W['tcWW'] + gp2**2*W['tcBB'] 
                      + 4.*gw2*gp2*W['tcWB'] )/(gw2+gp2)**2
        M['tCza'] = ( gw2*W['tcWW'] - gp2*W['tcBB'] 
                      - 2.*(gw2-gp2)*W['tcWB'] )/(gw2+gp2)
        M['tCww'] =  W['tcWW']
        
        
        # solution for  [eqn (4.16)]
        # dy*cos(phi) == X
        # dy*sin(phi) == Y
        def dy_sf(X,Y): 
            R = sqrt(X**2+Y**2)
            if R==0: 
                return 0., 0.
            elif Y==0.:
                return X, 0.
            else:
                signY = Y/abs(Y)
                if X*Y > 0.:
                    return R*signY, abs(Y)/R
                else:
                    return -R*signY, -abs(Y)/R

        # Yukawa type interaction coefficients [eqns. (5.16) & (4.17)]
        for f in ('u','d','e'):
            matrix = 'WBx' + f
            for i,j in W['WBx'+f].keys(): 
                diag = delta(i,j)*(W['cH']+dv)
                diag2 = 2.*delta(i,j)*W['cH']
                mi, mj = self.mass[ PID[f][i] ],self.mass[ PID[f][j] ] 
                if mi and mj:
                    s_Re, s_Im = W[matrix][i,j].real, W[matrix][i,j].imag
                    dy_cosphi = vev*s_Re/sqrt(2.*mi*mj) - diag
                    dy_sinphi = vev*s_Im/sqrt(2.*mi*mj)
                    
                    M['MBxdY'+f][i,j], M['MBxS'+f][i,j] = dy_sf(dy_cosphi, 
                                                                dy_sinphi)
                    
                    re = 3.*vev*s_Re/sqrt(2.*mi*mj) - diag2
                    im = 3.*vev*s_Im/sqrt(2.*mi*mj)
                    M['MBxY2'+f][i,j] = complex(re, im)
                    
        # Triple gauge couplings [eqn. (4.18)]
        M['dG1z'] = (gw2+gp2)/(gw2-gp2)*( -W['cWB']*gp2 + W['cT'] - dv )
        M['dKa']  = W['cWB']*gw2
        M['dKz']  = (-W['cWB']*2.*gw2*gp2 + (gw2+gp2)*(W['cT'] - dv))/(gw2-gp2)
        M['La']   = -W['c3W']*3./2.*gw2**2
        M['Lz']   = M['La']
        M['tKa']  = W['tcWB']*gw2
        M['tKz']  = -W['tcWB']*gp2
        M['tLa']  = -W['tc3W']*3./2.*gw2**2
        M['tLz']  = M['tLa']
        M['C3g']  = W['c3G']
        M['tC3g']  = W['tc3G']
        
        # Quartic gauge couplings [Sec 3.7] [eqn (3.23)] 
        M['dGw4'] = 2.*c2w*M['dG1z']
        M['dGw2z2'] = 2.*M['dG1z']
        M['dGw2za'] = M['dG1z']
        
        # two derivative quartic gauge couplings [Sec 3.7] [eqn (3.24)] 
        M['Lw4'] = -gw2/2./MW**2*M['Lz']
        M['tLw4'] = -gw2/2./MW**2*M['tLz']
        M['Lw2z2'] = -gw2*c2w/MW**2*M['Lz']
        M['tLw2z2'] = -gw2*c2w/MW**2*M['tLz']
        M['Lw2za'] = -ee2/MW**2*M['Lz']
        M['tLw2za'] = -ee2/MW**2*M['tLz']
        M['Lw2a2'] = -sqrt(ee2*gw2*c2w)/MW**2*M['Lz']
        M['tLw2a2'] = -sqrt(ee2*gw2*c2w)/MW**2*M['tLz']
        M['Lw2az'] = -sqrt(ee2*gw2*c2w)/MW**2*M['Lz']
        M['tLw2az'] = -sqrt(ee2*gw2*c2w)/MW**2*M['tLz']
        M['C4g'] = 3.*sqrt(gs2)**3/vev**2*M['C3g']
        M['tC4g'] = 3.*sqrt(gs2)**3/vev**2*M['tC3g']
        
        # Higgs cubic interaction [eqn. (4.19)]
        M['dL3']  =  -MH**2/(2.*vev**2) * (3.*W['cH'] + dv) - W['c6H']
        M['dL4'] = 3./2.*M['dL3'] - MH**2/vev**2/6.*M['dCz']
        
        # Couplings of two Higgs bosons to gluons [Sec 3.8]
        # [eqn (3.27)] copied from HiggsBasis implemetation
        M['Cgg2'], M['tCgg2'] = M['Cgg'], M['tCgg']
        
        M['dCz2'] = M['dCz']
        M['dCw2'] = M['dCw'] + 3.*M['dM']
                
        hvv = ['Cww', 'Cwbx', 'Cgg', 'Caa', 'Cza', 'Czz', 'Czbx', 'Cabx',
               'tCww', 'tCgg', 'tCaa', 'tCza', 'tCzz']
        
        for cvv in hvv:
            cvv2 = cvv +'2'
            M[cvv2] = M[cvv]
        
        # 4-fermion 
        M['cll1122'] = W['cll1122']
        M['cpuu3333'] = W['cpuu3333']
        M['cll1221'] = W['cll1221']
        
        # print B.card.blocks
        # print B.card.blocks['mbvertex']
        # W mass shift
        self.mass[24] = MW + M['dM']
        return M
    
    @Basis.translation('higgs')        
    def to_higgs(self, instance):
        s2w, c2w, ee2, gw2, gp2, MZ, vev, gs2 = self.calculate_inputs() 
        MH = self.mass[25]
        MW = MZ*sqrt(c2w)
        
        W = self
        H = instance
        
        def delta(i,j):
            return 1. if i==j else 0.
                    
        def f(T3,Q,i,j): # [eqn (4.11)]
                return delta(i,j)*( -Q*W['cWB']*gw2*gp2/(gw2-gp2)
                                   + (W['cT']-dv)*(T3 + Q*gp2/(gw2-gp2)))

        # Higgs vev shift [eqn (4.8)]
        dv = (W['WBxHpl'][1,1].real + W['WBxHpl'][2,2].real)/2.-W['cll1221']/4.
        
        # W mass shift [eqn (4.9)]
        H['dM'] = ( gw2*W['cT'] - gp2*gw2*W['cWB']-gp2*dv )/(gw2-gp2)
        
        # W/Z chiral coupling deviations
        for i,j in W['WBxHpl'].keys():
            H['HBxdGLzv'][i,j] = (1./2.*W['WBxHpl'][i,j] 
                                - 1./2.*W['WBxHl'][i,j] + f(1./2.,0.,i,j))
            H['HBxdGLze'][i,j] = (-1./2.*W['WBxHpl'][i,j] 
                                - 1./2.*W['WBxHl'][i,j] + f(-1./2.,-1.,i,j))
            H['HBxdGRze'][i,j] = (- 1./2.*W['WBxHe'][i,j] + f(0.,-1.,i,j))
            H['HBxdGLzu'][i,j] = (1./2.*W['WBxHpq'][i,j] 
                                - 1./2.*W['WBxHq'][i,j] + f(1./2.,2./3.,i,j))
            H['HBxdGLzd'][i,j] = (-1./2.*W['WBxHpq'][i,j] 
                                - 1./2.*W['WBxHq'][i,j] + f(-1./2.,-1./3.,i,j))
            H['HBxdGRzu'][i,j] = (- 1./2.*W['WBxHu'][i,j] + f(0.,2./3.,i,j))
            H['HBxdGRzd'][i,j] = (- 1./2.*W['WBxHd'][i,j] + f(0.,-1./3.,i,j))
            H['HBxdGLwl'][i,j] = (W['WBxHpl'][i,j] + f(1./2.,0.,i,j) 
                                - f(-1./2.,-1.,i,j))

        for k,v in W['WBxHud'].iteritems():
            H['HBxdGRwq'][k] = -v/2.
        
        matrix_sub(matrix_mult(H['HBxdGLzu'], W.ckm),
                   matrix_mult(W.ckm, H['HBxdGLzd']),
                   H['HBxdGLwq'])
        
        vertex = ['HBxdGLze', 'HBxdGRze', 'HBxdGLzv', 'HBxdGLzu', 'HBxdGRzu', 
                  'HBxdGLzd', 'HBxdGRzd', 'HBxdGLwl', 'HBxdGLwq', 'HBxdGRwq']

        # HVFF coeffs [Eqns (3.10)]
        for dG in vertex: 
            dGh = dG[:-2]+'h'+dG[-2:]
            matrix_eq(H[dG], H[dGh])

        # Higgs couplings to W/Z [eqn (4.14)]
        H['dCw'] = ( -W['cH']*(gw2-gp2) - W['cWB']*4.*gw2*gp2 
                    + W['cT']*4.*gw2 - dv*(3.*gw2+gp2) )/(gw2-gp2)
        H['dCz'] = -W['cH'] - 3.*dv 
        
        # Two derivative field strength interactions [eqn (4.15)]
        H['Cgg']  = W['cGG'] 
        H['Caa']  = W['cWW'] + W['cBB'] - 4.*W['cWB']
        H['Czz']  = ( gw2**2*W['cWW'] + gp2**2*W['cBB'] 
                     + 4.*gw2*gp2*W['cWB'] )/(gw2+gp2)**2
        H['Czbx'] =  -(2./gw2)*(W['cT'] - dv) 
        H['Cza']  = ( gw2*W['cWW'] - gp2*W['cBB'] 
                      - 2.*(gw2-gp2)*W['cWB'] )/(gw2+gp2)
        H['Cabx'] =  2./(gw2-gp2)*((gw2+gp2)*W['cWB'] - 2.*W['cT'] + 2.*dv) 
        H['Cww']  =  W['cWW']
        H['Cwbx'] =  2./(gw2-gp2)*(gp2*W['cWB'] - W['cT'] + dv) 
        H['tCgg'] = W['tcGG'] 
        H['tCaa'] =  W['tcWW'] + W['tcBB'] - 4.*W['tcWB']
        H['tCzz'] = ( gw2**2*W['tcWW'] + gp2**2*W['tcBB'] 
                      + 4.*gw2*gp2*W['tcWB'] )/(gw2+gp2)**2
        H['tCza'] = ( gw2*W['tcWW'] - gp2*W['tcBB'] 
                      - 2.*(gw2-gp2)*W['tcWB'] )/(gw2+gp2)
        H['tCww'] =  W['tcWW']
        
        
        # solution for  [eqn (4.16)]
        # dy*cos(phi) == X;  dy*sin(phi) == Y
        def dy_sf(X,Y): 
            R = sqrt(X**2+Y**2)
            if R==0: 
                return 0., 0.
            elif Y==0.:
                return X, 0.
            else:
                signY = Y/abs(Y)
                if X*Y > 0.:
                    return R*signY, abs(Y)/R
                else:
                    return -R*signY, -abs(Y)/R

        # Yukawa type interaction coefficients [eqns. (5.16) & (4.17)]
        for f in ('u','d','e'):
            matrix = 'WBx' + f
            for i,j in W['WBx'+f].keys(): 
                diag = delta(i,j)*(W['cH']+dv)
                diag2 = 2.*delta(i,j)*W['cH']
                mi, mj = self.mass[ PID[f][i] ],self.mass[ PID[f][j] ] 
                if mi and mj:
                    s_Re, s_Im = W[matrix][i,j].real, W[matrix][i,j].imag
                    dy_cosphi = vev*s_Re/sqrt(2.*mi*mj) - diag
                    dy_sinphi = vev*s_Im/sqrt(2.*mi*mj)
                    
                    H['HBxdY'+f][i,j], H['HBxS'+f][i,j] = dy_sf(dy_cosphi, 
                                                                dy_sinphi)
                    
                    re = 3.*vev*s_Re/sqrt(2.*mi*mj) - diag2
                    im = 3.*vev*s_Im/sqrt(2.*mi*mj)
                    H['HBxY2'+f][i,j] = complex(re, im)
                    
        # Triple gauge couplings [eqn. (4.18)]
        H['dG1z'] = (gw2+gp2)/(gw2-gp2)*( -W['cWB']*gp2 + W['cT'] - dv )
        H['dKa']  = W['cWB']*gw2
        H['dKz']  = (-W['cWB']*2.*gw2*gp2 + (gw2+gp2)*(W['cT'] - dv))/(gw2-gp2)
        H['La']   = -W['c3W']*3./2.*gw2**2
        H['Lz']   = H['La']
        H['tKa']  = W['tcWB']*gw2
        H['tKz']  = -W['tcWB']*gp2
        H['tLa']  = -W['tc3W']*3./2.*gw2**2
        H['tLz']  = H['tLa']
        H['C3g']  = W['c3G']
        H['tC3g']  = W['tc3G']
        
        # Quartic gauge couplings [Sec 3.7] [eqn (3.23)] 
        H['dGw4'] = 2.*c2w*H['dG1z']
        H['dGw2z2'] = 2.*H['dG1z']
        H['dGw2za'] = H['dG1z']
        
        # two derivative quartic gauge couplings [Sec 3.7] [eqn (3.24)] 
        H['Lw4'] = -gw2/2./MW**2*H['Lz']
        H['tLw4'] = -gw2/2./MW**2*H['tLz']
        H['Lw2z2'] = -gw2*c2w/MW**2*H['Lz']
        H['tLw2z2'] = -gw2*c2w/MW**2*H['tLz']
        H['Lw2za'] = -ee2/MW**2*H['Lz']
        H['tLw2za'] = -ee2/MW**2*H['tLz']
        H['Lw2a2'] = -sqrt(ee2*gw2*c2w)/MW**2*H['Lz']
        H['tLw2a2'] = -sqrt(ee2*gw2*c2w)/MW**2*H['tLz']
        H['Lw2az'] = -sqrt(ee2*gw2*c2w)/MW**2*H['Lz']
        H['tLw2az'] = -sqrt(ee2*gw2*c2w)/MW**2*H['tLz']
        H['C4g'] = 3.*sqrt(gs2)**3/vev**2*H['C3g']
        H['tC4g'] = 3.*sqrt(gs2)**3/vev**2*H['tC3g']
        
        # Higgs cubic interaction [eqn. (4.19)]
        H['dL3']  =  -MH**2/(2.*vev**2) * (3.*W['cH'] + dv) - W['c6H']
        H['dL4'] = 3./2.*H['dL3'] - MH**2/vev**2/6.*H['dCz']
        
        # Couplings of two Higgs bosons to gluons [Sec 3.8]
        # [eqn (3.27)] copied from HiggsBasis implemetation
        H['Cgg2'], H['tCgg2'] = H['Cgg'], H['tCgg']
        
        H['dCz2'] = H['dCz']
        H['dCw2'] = H['dCw'] + 3.*H['dM']
                
        hvv = ['Cww', 'Cwbx', 'Cgg', 'Caa', 'Cza', 'Czz', 'Czbx', 'Cabx',
               'tCww', 'tCgg', 'tCaa', 'tCza', 'tCzz']
        
        for cvv in hvv:
            cvv2 = cvv +'2'
            H[cvv2] = H[cvv]
        
        # 4-fermion 
        H['cll1122'] = W['cll1122']
        H['cpuu3333'] = W['cpuu3333']
        H['cll1221'] = W['cll1221']
        
        self.mass[24] = MW + H['dM']
        return H
        
    @Basis.translation('silh')
    def to_silh(self, instance):
        s2w, c2w, ee2, gw2, gp2, MZ, vev, gs2 = self.calculate_inputs() 
        MH = self.mass[25]
        lam = -MH**2/(2.*vev**2) # Higgs self-coupling
        
        W = self
        S = instance
        
        def delta(i,j):
            return 1. if i==j else 0.
        
        S['sBB'] = W['cBB'] - 4.*W['cWB'] + W['cWW']
        S['tsBB'] = W['tcBB'] - 4.*W['tcWB'] + W['tcWW']
        S['s2W'] = W['cll1221']/gw2
        S['s2B'] = (2.*W['cll1122'] + W['cll1221'])/gp2
        S['s2G'] = 3.*W['cpuu3333']/gs2
        
        S['sB'] = (4.*W['cWB'] - W['cWW']  - (4.*W['WBxHl'][1,1].real 
                 + 4.*W['cll1122'] + 2.*W['cll1221'])/gp2)
        S['sW'] = W['cWW'] - (2.*W['cll1221'] - 4.*W['WBxHpl'][1,1].real)/gw2 
        S['sHW'] = - W['cWW']
        S['tsHW'] = - W['tcWW']
        S['sHB'] = W['cWW'] - 4.*W['cWB']
        S['tsHB'] = W['tcWW'] - 4.*W['tcWB']
        S['sH'] = W['cH'] - 3./4.*W['cll1221'] + 3.*W['WBxHpl'][1,1,].real
        S['sT'] = (W['cT'] - W['WBxHl'][1,1,].real - W['cll1122']/2. 
                 - W['cll1221']/4.)
    
        S['s6H'] = W['c6H'] - 2.*lam*(W['cll1221'] - 4.*W['WBxHpl'][1,1,].real)
        
        def sHf(wc, Yf, i, j):
                return W[wc][i,j] + delta(i,j)*2.*Yf*W['WBxHl'][1,1,].real 

        def spHf(wc, i, j):
                return W[wc][i,j] - delta(i,j)*W['WBxHpl'][1,1,].real
        
        for i,j in W['WBxHl'].keys():
            S['SBxHl'][i,j] = sHf('WBxHl', -1./2., i, j)
            S['SBxHe'][i,j] = sHf('WBxHe', -1., i, j)
            S['SBxHq'][i,j] = sHf('WBxHq', 1./6., i, j)
            S['SBxHu'][i,j] = sHf('WBxHu', 2./3., i, j)
            S['SBxHd'][i,j] = sHf('WBxHd', -1./3., i, j)
            S['SBxHpl'][i,j] = spHf('WBxHpl', i, j)
            S['SBxHpq'][i,j] = spHf('WBxHpq', i, j)
            if i==j==1:
                S['SBxHl'][i,j] = 0.
                S['SBxHpl'][i,j] = 0.
        for k in W['WBxHud'].keys():
            S['SBxHud'][k] = W['WBxHud'][k]            

        def sf(wc, i, j):
            mass = self.mass[ PID[f][i] ]
            yuk = mass/vev/sqrt(2.)
            return W[wc][i,j] - delta(i,j)*yuk*(W['cll1221'] 
                  - 4.*W['WBxHpl'][1,1,].real)

            
        # [eqn (5.9)]
        for f in ('u','d','e'): # fermion loop
            wmatrix, smatrix = 'WBx'+f, 'SBx'+f
            for i,j in W[wmatrix].keys(): # flavour loop
                S[smatrix][i,j] = sf(wmatrix, i, j)

        # trivial translation, sX==cX
        trivial = ['sGG','tsGG','s3W','ts3W','s3G','ts3G']
        for coeff in trivial:
            wcoeff = coeff.replace('s','c')
            S[coeff] = W[wcoeff]
        
        return S

################################################################################