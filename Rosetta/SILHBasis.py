from internal import Basis
import math
from math import sqrt
from itertools import combinations_with_replacement as comb
from itertools import product
from internal import PID
################################################################################
flavmat = Basis.flavour_matrix
################################################################################
class SILHBasis(Basis.Basis):
    '''
    Basis class for Rosetta based on [Giudice et al., JHEP 0706 (2007) 045 ] 
    and a number of succesive publications. Part of the three intrinsic basis 
    implementations in Rosetta along with the Higgs and Warsaw bases. The exact 
    list of operators included as well as the equations for the translation to 
    the Higgs and Warsaw bases can be found in the HXSWG note, for which all 
    references to tables and equations are in this implementation. The basis 
    differs slightly to the versions defined in the original publications. 
    Please refer to the HXSWG note for further details. Table 1 contains the 
    majority of the operators included, with the remaining structure described 
    in section 5.
    '''

    name = 'silh'
    ##########################
    # declare coefficients
    
    # operators present in SILH but not Warsaw, see [Sec. 5].
    # cWW, ctWW, cWB, ctWB, cHl11, cpHl11, cll1221, cll1122, cpuu3333 
    # absent compared to Warsaw.
    SBNOTWARSAW = ['sW','sB','sHW','sHB','tsHW','tsHB','s2W','s2B','s2G']
    # [Tab. 1]
    SBV2H2 = ['sGG','tsGG','sBB', 'tsBB']
    
    SBH4D2 = ['sH','sT']
    
    SBH6 = ['s6H']
    
    SBV3D3 = ['s3W','s3G','ts3W','ts3G']

    ##########################
    # block structure
    blocks = {'SBxV2H2':SBV2H2, 'SBxH4D2':SBH4D2, 'SBxH6':SBH6,
              'SBxV3D3':SBV3D3, 'SBxNOTWARSAW':SBNOTWARSAW}    
    
    flavoured={
        'SBxu': {'cname':'su', 'kind':'general', 'domain':'complex'},
        'SBxd': {'cname':'sd', 'kind':'general', 'domain':'complex'},
        'SBxe': {'cname':'se', 'kind':'general', 'domain':'complex'},
        'SBxHl': {'cname':'sHl' , 'kind':'hermitian', 'domain':'complex'},
        'SBxHpl': {'cname':'sHpl', 'kind':'hermitian', 'domain':'complex'},
        'SBxHe': {'cname':'sHe' , 'kind':'hermitian', 'domain':'complex'},
        'SBxHq': {'cname':'sHq' , 'kind':'hermitian', 'domain':'complex'},
        'SBxHpq': {'cname':'sHpq', 'kind':'hermitian', 'domain':'complex'},
        'SBxHu': {'cname':'sHu' , 'kind':'hermitian', 'domain':'complex'},
        'SBxHd' : {'cname':'sHd' , 'kind':'hermitian', 'domain':'complex'},
        'SBxHud' : {'cname':'sHud', 'kind':'general', 'domain':'complex'}
    }
    #
    # independent = ( [c for v in blocks.values() for c in v] +
    #                 [c for c in flavoured.keys()] )
                    
    independent = [c for c in blocks.keys()] + [c for c in flavoured.keys()]
                  
    # two coefficients are zero by construction (set in calculate_dependent())
    dependent = ['RsHl1x1', 'RsHpl1x1']
    
    required_masses = set([y for x in PID.values() for y in x.values()])
    required_inputs = {1, 2, 3, 4, 25} # aEWM1, Gf, aS, MZ, MH
    ##########################
    def calculate_dependent(self):
        # These coefficients are implictly set to zero
        self['RsHl1x1'], self['RsHl1x1']=0., 0.
        
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
    def to_mass(self,instance):
        s2w, c2w, ee2, gw2, gp2, MZ, vev, gs2 = self.calculate_inputs() 
        MH = self.mass[25]
        MW = MZ*sqrt(c2w)
                
        A = self
        B = instance
        
        # W mass shift [eqn (5.11)]
        B['dM'] = - gw2*gp2/(4.*(gw2-gp2))*(A['sW'] + A['sB'] + A['s2W'] 
                                            + A['s2B'] - 4./gp2*A['sT'] 
                                            + 2./gw2*A['RsHpl2x2'])
        def f(T3,Q,i,j): # [eqn (5.12)]
            Qcoeff = gp2/4./(gw2-gp2)*( -(2.*gw2-gp2)*A['s2B'] 
                       - gw2*(A['s2W'] + A['sW'] + A['sB'] ) 
                       + 4.*A['sT'] - 2.*A['RsHpl2x2'])
            T3coeff = ( gw2*A['s2W'] + gp2*A['s2B'] 
                        + 4.*A['sT']-2.*A['RsHpl2x2'])/4.
            return  delta(i,j)*(Q*Qcoeff + T3*T3coeff)
            
        def delta(i,j):
            return 1. if i==j else 0.
            
        # W/Z chiral coupling deviations
        for px in ('','IM'):
            for i,j in A[px+'SBxHpl'].keys():
                B[px+'MBxdGLzv'][i,j] = (1./2.*A[px+'SBxHpl'][i,j] 
                                       - 1./2.*A[px+'SBxHl'][i,j] 
                                       + f(1./2.,0.,i,j))
                B[px+'MBxdGLze'][i,j] = (-1./2.*A[px+'SBxHpl'][i,j] 
                                       - 1./2.*A[px+'SBxHl'][i,j] 
                                       + f(-1./2.,-1.,i,j))
                B[px+'MBxdGRze'][i,j] = (- 1./2.*A[px+'SBxHe'][i,j] 
                                       + f(0.,-1.,i,j))
                B[px+'MBxdGLzu'][i,j] = (1./2.*A[px+'SBxHpq'][i,j] 
                                       - 1./2.*A[px+'SBxHq'][i,j] 
                                       + f(1./2.,2./3.,i,j))
                B[px+'MBxdGLzd'][i,j] = (-1./2.*A[px+'SBxHpq'][i,j] 
                                       - 1./2.*A[px+'SBxHq'][i,j] 
                                       + f(-1./2.,-1./3.,i,j))
                B[px+'MBxdGRzu'][i,j] = (- 1./2.*A[px+'SBxHu'][i,j] 
                                       + f(0.,2./3.,i,j))
                B[px+'MBxdGRzd'][i,j] = (- 1./2.*A[px+'SBxHd'][i,j]
                                       + f(0.,-1./3.,i,j))
                B[px+'MBxdGLwl'][i,j] = (A[px+'SBxHpl'][i,j] + f(1./2.,0.,i,j) 
                                       - f(-1./2.,-1.,i,j))
                B[px+'MBxdGLwq'][i,j] = (A[px+'SBxHpq'][i,j] 
                                       + f(1./2.,2./3.,i,j)
                                       - f(-1./2.,-1./3.,i,j))
            for k,v in A[px+'SBxHud'].iteritems():
                B[px+'MBxdGRwq'][k] = -v/2.

        
        vertex = ['MBxdGLze', 'MBxdGRze', 'MBxdGLzv', 'MBxdGLzu', 'MBxdGRzu', 
                  'MBxdGLzd', 'MBxdGRzd', 'MBxdGLwl', 'MBxdGLwq', 'MBxdGRwq']

        # HVFF coeffs [Eqns (3.10)]
        for dG in vertex: 
            dGh = dG[:-2]+'h'+dG[-2:]
            for k,v in B[dG].iteritems():
                B[dGh][k] = v
            for k,v in B['IM'+dG].iteritems():
                B['IM'+dGh][k] = v
        
        # Higgs couplings to W/Z [eqn (4.14)]
        B['dCw'] =  -A['sH'] - gw2*gp2/(gw2-gp2)*( A['sW'] + A['sB'] 
                         + A['s2W'] + A['s2B'] - 4./gp2*A['sT'] 
                         + (3.*gw2 + gp2)/(2.*gw2*gp2)*A['RsHpl2x2'])
        B['dCz'] = -A['sH'] - 3./2.*A['RsHpl2x2'] 
        
        # Two derivative field strength interactions [eqn (4.15)]
        B['Cgg']  = A['sGG']
        B['Caa']  = A['sBB']
        B['Czz']  = -1./(gw2+gp2)*( gw2*A['sHW'] + gp2*A['sHB'] 
                                     - gp2*s2w*A['sBB'] )
        B['Czbx'] =  1./(2.*gw2)*( gw2*(A['sW'] + A['sHW'] +A ['s2W'])
                                  + gp2*(A['sB'] + A['sHB'] +A ['s2B'])
                                  - 4.*A['sT'] + 2.*A['RsHpl2x2'] ) 
        B['Cza']  = ( A['sHB'] - A['sHW'])/2. - s2w*A['sBB']
        
        B['Cabx'] =  ( A['sHW'] - A['sHB'])/2.+ 1./(gw2-gp2)*(
                           gw2*(A['sW']+A['s2W']) + gp2*(A['sB']+A['s2B'])
                         - 4.*A['sT'] + 2.*A['RsHpl2x2'] )
        B['Cww']  =  -A['sHW']
        # factor 2 wrong in note here 
        B['Cwbx'] =  A['sHW']/2. + 1./2./(gw2-gp2)*(
                           gw2*(A['sW']+A['s2W']) + gp2*(A['sB']+A['s2B'])
                         - 4.*A['sT'] + 2.*A['RsHpl2x2'] )
        B['tCgg'] = A['tsGG'] 
        B['tCaa'] = A['tsBB']
        B['tCzz'] = -1./(gw2+gp2)*( gw2*A['tsHW'] + gp2*A['tsHB'] 
                                     - gp2*s2w*A['tsBB'] )
        B['tCza'] = ( A['tsHB'] - A['tsHW'])/2. - s2w*A['tsBB']
        B['tCww'] =  -A['tsHW']
        
        
        # solution for  [eqn (5.16)]
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
            matrix = 'SBx' + f
            for i,j in A['SBx'+f].keys(): 
                diag = delta(i,j)*(A['sH'] + A['RsHpl2x2']/2.)
                diag2 = 2.*delta(i,j)*(A['sH'])
                mi, mj = self.mass[ PID[f][i] ],self.mass[ PID[f][j] ] 
                if mi and mj:
                    s_Re, s_Im = A[matrix][i,j], A['IM'+matrix][i,j]
                    dy_cosphi = vev*s_Re/sqrt(2.*mi*mj) - diag
                    dy_sinphi = vev*s_Im/sqrt(2.*mi*mj)
                    
                    B['MBxdY'+f][i,j], B['MBxS'+f][i,j] = dy_sf(dy_cosphi, 
                                                                dy_sinphi)
                    
                    B['MBxY2'+f][i,j] = 3.*vev*s_Re/sqrt(2.*mi*mj) - diag2
                    B['IMMBxY2'+f][i,j] = 3.*vev*s_Im/sqrt(2.*mi*mj)

        # Triple gauge couplings [eqn. (5.18)]
        B['dG1z'] = -(gw2+gp2)/(gw2-gp2)/4.*( (gw2-gp2)*A['sHW'] 
                        + gw2*(A['sW'] + A['s2W']) + gp2*(A['sB'] + A['s2B']) 
                        - 4.*A['sT'] + 2.*A['RsHpl2x2'] )
        B['dKa'] = - gw2/4.*(A['sHW']+A['sHB'])
        B['dKz'] = ( -1./4.*(gw2*A['sHW'] - gp2*A['sHB'])
                      -(gw2+gp2)/(gw2-gp2)/4.*(
                          gw2*(A['sW'] + A['s2W'])
                          + gp2*(A['sB'] + A['s2B']) 
                          - 4.*A['sT'] + 2.*A['RsHpl2x2'] 
                          )
                      )
        B['La'] = -A['s3W']*3./2.*gw2**2
        B['Lz'] =  B['La']
        B['tKa'] =  - gw2/4.*(A['tsHW']+A['tsHB'])
        B['tKz'] = gp2/4.*(A['tsHW']+A['tsHB'])
        B['tLa'] = -A['ts3W']*3./2.*gw2**2
        B['tLz'] =  B['tLa']
        B['C3g']  = A['s3G']
        B['tC3g']  = A['ts3G']

        # Quartic gauge couplings [Sec 3.7] [eqn (3.23)] 
        B['dGw4'] = 2.*c2w*B['dG1z']
        B['dGw2z2'] = 2.*B['dG1z']
        B['dGw2za'] = B['dG1z']
        
        # two derivative quartic gauge couplings [Sec 3.7] [eqn (3.24)] 
        B['Lw4'] = -gw2/2./MW**2*B['Lz']
        B['tLw4'] = -gw2/2./MW**2*B['tLz']
        B['Lw2z2'] = -gw2*c2w/MW**2*B['Lz']
        B['tLw2z2'] = -gw2*c2w/MW**2*B['tLz']
        B['Lw2za'] = -ee2/MW**2*B['Lz']
        B['tLw2za'] = -ee2/MW**2*B['tLz']
        B['Lw2a2'] = -sqrt(ee2*gw2*c2w)/MW**2*B['Lz']
        B['tLw2a2'] = -sqrt(ee2*gw2*c2w)/MW**2*B['tLz']
        B['Lw2az'] = -sqrt(ee2*gw2*c2w)/MW**2*B['Lz']
        B['tLw2az'] = -sqrt(ee2*gw2*c2w)/MW**2*B['tLz']
        B['C4g'] = 3.*sqrt(gs2)**3/vev**2*B['C3g']
        B['tC4g'] = 3.*sqrt(gs2)**3/vev**2*B['tC3g']

        # Higgs cubic interaction [eqn. (4.19)]
        B['dL3']  =  -MH**2/(2.*vev**2)*(3.*A['sH'] + A['RsHpl2x2']/2.) - A['s6H']
        B['dL4'] = 3./2.*B['dL3'] - MH**2/vev**2/6.*B['dCz']
        
        # Couplings of two Higgs bosons to gluons [Sec 3.8]
        # [eqn (3.27)] copied from HiggsBasis implemetation
        B['Cgg2'], B['tCgg2'] = B['Cgg'], B['tCgg']
        
        B['dCz2'] = B['dCz']
        B['dCw2'] = B['dCw'] + 3.*B['dM']
        
        hvv = ['Cww', 'Cwbx', 'Cgg', 'Caa', 'Cza', 'Czz', 'Czbx', 'Cabx',
               'tCww', 'tCgg', 'tCaa', 'tCza', 'tCzz']
        
        for cvv in hvv:
            cvv2 = cvv +'2'
            B[cvv2] = B[cvv]
        
        B['cll1122'] = (gp2*A['s2B'] - gw2*A['s2W'])/2.
        B['cll1221'] = gw2*A['s2W']
        B['cpuu3333'] = (1./3.)*gs2*A['s2G']
        
        self.mass[24] = MW + B['dM']

        return B
    
    @Basis.translation('higgs')        
    def to_higgs(self, instance):
        s2w, c2w, ee2, gw2, gp2, MZ, vev, gs2 = self.calculate_inputs() 
        MH = self.mass[25]
        MW = MZ*sqrt(c2w)
                
        W = self
        H = instance
        
        # W mass shift [eqn (5.11)]
        H['dM'] = - gw2*gp2/(4.*(gw2-gp2))*(W['sW'] + W['sB'] + W['s2W'] 
                                            + W['s2B'] - 4./gp2*W['sT'] 
                                            + 2./gw2*W['RsHpl2x2'])
        def f(T3,Q,i,j): # [eqn (5.12)]
            Qcoeff = gp2/4./(gw2-gp2)*( -(2.*gw2-gp2)*W['s2B'] 
                       - gw2*(W['s2W'] + W['sW'] + W['sB'] ) 
                       + 4.*W['sT'] - 2.*W['RsHpl2x2'])
            T3coeff = ( gw2*W['s2W'] + gp2*W['s2B'] 
                        + 4.*W['sT']-2.*W['RsHpl2x2'])/4.
            return  delta(i,j)*(Q*Qcoeff + T3*T3coeff)
            
        def delta(i,j):
            return 1. if i==j else 0.
            
        # W/Z chiral coupling deviations
        for px in ('','IM'):
            for i,j in W[px+'SBxHpl'].keys():
                H[px+'HBxdGLzv'][i,j] = (1./2.*W[px+'SBxHpl'][i,j] 
                                       - 1./2.*W[px+'SBxHl'][i,j] 
                                       + f(1./2.,0.,i,j))
                H[px+'HBxdGLze'][i,j] = (-1./2.*W[px+'SBxHpl'][i,j] 
                                       - 1./2.*W[px+'SBxHl'][i,j] 
                                       + f(-1./2.,-1.,i,j))
                H[px+'HBxdGRze'][i,j] = (- 1./2.*W[px+'SBxHe'][i,j] 
                                       + f(0.,-1.,i,j))
                H[px+'HBxdGLzu'][i,j] = (1./2.*W[px+'SBxHpq'][i,j] 
                                       - 1./2.*W[px+'SBxHq'][i,j] 
                                       + f(1./2.,2./3.,i,j))
                H[px+'HBxdGLzd'][i,j] = (-1./2.*W[px+'SBxHpq'][i,j] 
                                       - 1./2.*W[px+'SBxHq'][i,j] 
                                       + f(-1./2.,-1./3.,i,j))
                H[px+'HBxdGRzu'][i,j] = (- 1./2.*W[px+'SBxHu'][i,j] 
                                       + f(0.,2./3.,i,j))
                H[px+'HBxdGRzd'][i,j] = (- 1./2.*W[px+'SBxHd'][i,j]
                                       + f(0.,-1./3.,i,j))
                H[px+'HBxdGLwl'][i,j] = (W[px+'SBxHpl'][i,j] + f(1./2.,0.,i,j) 
                                       - f(-1./2.,-1.,i,j))
                H[px+'HBxdGLwq'][i,j] = (W[px+'SBxHpq'][i,j] 
                                       + f(1./2.,2./3.,i,j)
                                       - f(-1./2.,-1./3.,i,j))
            for k,v in W[px+'SBxHud'].iteritems():
                H[px+'HBxdGRwq'][k] = -v/2.

        
        vertex = ['HBxdGLze', 'HBxdGRze', 'HBxdGLzv', 'HBxdGLzu', 'HBxdGRzu', 
                  'HBxdGLzd', 'HBxdGRzd', 'HBxdGLwl', 'HBxdGLwq', 'HBxdGRwq']

        # HVFF coeffs [Eqns (3.10)]
        for dG in vertex: 
            dGh = dG[:-2]+'h'+dG[-2:]
            for k,v in H[dG].iteritems():
                H[dGh][k] = v
            for k,v in H['IM'+dG].iteritems():
                H['IM'+dGh][k] = v

        # Higgs couplings to W/Z [eqn (4.14)]
        H['dCw'] =  -W['sH'] - gw2*gp2/(gw2-gp2)*( W['sW'] + W['sB'] 
                         + W['s2W'] + W['s2B'] - 4./gp2*W['sT'] 
                         + (3.*gw2 + gp2)/(2.*gw2*gp2)*W['RsHpl2x2'])
        H['dCz'] = -W['sH'] - 3./2.*W['RsHpl2x2'] 
        
        # Two derivative field strength interactions [eqn (4.15)]
        H['Cgg']  = W['sGG']
        H['Caa']  = W['sBB']
        H['Czz']  = -1./(gw2+gp2)*( gw2*W['sHW'] + gp2*W['sHB'] 
                                     - gp2*s2w*W['sBB'] )
        H['Czbx'] =  1./(2.*gw2)*( gw2*(W['sW'] + W['sHW'] +W['s2W'])
                                  + gp2*(W['sB'] + W['sHB'] +W['s2B'])
                                  - 4.*W['sT'] + 2.*W['RsHpl2x2'] ) 
        H['Cza']  = ( W['sHB'] - W['sHW'])/2. - s2w*W['sBB']
        
        H['Cabx'] =  ( W['sHW'] - W['sHB'])/2.+ 1./(gw2-gp2)*(
                           gw2*(W['sW']+W['s2W']) + gp2*(W['sB']+W['s2B'])
                         - 4.*W['sT'] + 2.*W['RsHpl2x2'] )
        H['Cww']  =  -W['sHW']
        # factor 2 wrong in note here 
        H['Cwbx'] =  W['sHW']/2. + 1./2./(gw2-gp2)*(
                           gw2*(W['sW']+W['s2W']) + gp2*(W['sB']+W['s2B'])
                         - 4.*W['sT'] + 2.*W['RsHpl2x2'] )
        H['tCgg'] = W['tsGG'] 
        H['tCaa'] = W['tsBB']
        H['tCzz'] = -1./(gw2+gp2)*( gw2*W['tsHW'] + gp2*W['tsHB'] 
                                     - gp2*s2w*W['tsBB'] )
        H['tCza'] = ( W['tsHB'] - W['tsHW'])/2. - s2w*W['tsBB']
        H['tCww'] =  -W['tsHW']
        
        
        # solution for  [eqn (5.16)]
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
            matrix = 'SBx' + f
            for i,j in W['SBx'+f].keys(): 
                diag = delta(i,j)*(W['sH'] + W['RsHpl2x2']/2.)
                diag2 = 2.*delta(i,j)*(W['sH'])
                mi, mj = self.mass[ PID[f][i] ],self.mass[ PID[f][j] ] 
                if mi and mj:
                    s_Re, s_Im = W[matrix][i,j], W['IM'+matrix][i,j]
                    dy_cosphi = vev*s_Re/sqrt(2.*mi*mj) - diag
                    dy_sinphi = vev*s_Im/sqrt(2.*mi*mj)
                    
                    H['HBxdY'+f][i,j], H['HBxS'+f][i,j] = dy_sf(dy_cosphi, 
                                                                dy_sinphi)
                    
                    H['HBxY2'+f][i,j] = 3.*vev*s_Re/sqrt(2.*mi*mj) - diag2
                    H['IMHBxY2'+f][i,j] = 3.*vev*s_Im/sqrt(2.*mi*mj)

        # Triple gauge couplings [eqn. (5.18)]
        H['dG1z'] = -(gw2+gp2)/(gw2-gp2)/4.*( (gw2-gp2)*W['sHW'] 
                        + gw2*(W['sW'] + W['s2W']) + gp2*(W['sB'] + W['s2B']) 
                        - 4.*W['sT'] + 2.*W['RsHpl2x2'] )
        H['dKa'] = - gw2/4.*(W['sHW']+W['sHB'])
        H['dKz'] = ( -1./4.*(gw2*W['sHW'] - gp2*W['sHB'])
                      -(gw2+gp2)/(gw2-gp2)/4.*(
                          gw2*(W['sW'] + W['s2W'])
                          + gp2*(W['sB'] + W['s2B']) 
                          - 4.*W['sT'] + 2.*W['RsHpl2x2'] 
                          )
                      )
        H['La'] = -W['s3W']*3./2.*gw2**2
        H['Lz'] =  H['La']
        H['tKa'] =  - gw2/4.*(W['tsHW']+W['tsHB'])
        H['tKz'] = gp2/4.*(W['tsHW']+W['tsHB'])
        H['tLa'] = -W['ts3W']*3./2.*gw2**2
        H['tLz'] =  H['tLa']
        H['C3g']  = W['s3G']
        H['tC3g']  = W['ts3G']

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
        H['dL3']  =  -MH**2/(2.*vev**2)*(3.*W['sH'] + W['RsHpl2x2']/2.) - W['s6H']
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
        
        H['cll1122'] = (gp2*W['s2B'] - gw2*W['s2W'])/2.
        H['cll1221'] = gw2*W['s2W']
        H['cpuu3333'] = (1./3.)*gs2*W['s2G']
        
        self.mass[24] = MW + H['dM']

        return H
    
    @Basis.translation('warsaw')
    def to_warsaw(self, instance):

        s2w, c2w, ee2, gw2, gp2, MZ, vev, gs2 = self.calculate_inputs() 
        MH = self.mass[25]
        lam = -MH**2/(2.*vev**2) # Higgs self-coupling     
           
        S = self
        W = instance
        
        def delta(i,j):
            return 1. if i==j else 0.
        
        # [eqn (5.7)]
        W['cH'] = S['sH'] - 3./4.*gw2*(S['sW'] + S['sHW'] + S['s2W'])
        W['cT'] = S['sT'] - 1./4.*gp2*(S['sB'] + S['sHB'] + S['s2B'])
        W['c6H'] = S['s6H'] - 2.*lam*gw2*(S['sW'] + S['sHW'] + S['s2W'])
        W['cWB'] = - 1./4.*(S['sHB'] + S['sHW'])
        W['cBB'] = S['sBB'] - S['sHB']
        W['cWW'] = -S['sHW']
        W['tcWB'] = - 1./4.*(S['tsHB'] + S['tsHW'])
        W['tcBB'] = S['tsBB'] - S['tsHB']
        W['tcWW'] = -S['tsHW']
        
        def cHf(sc, Yf, i, j):
            if 'IM'!=sc[:2]:
                return S[sc][i,j] + delta(i,j)*gp2*Yf/2.*(S['sB'] 
                                                    + S['sHB'] + 2.*S['s2B'])
            else:
                return S[sc][i,j]
            
        def cpHf(sc, i, j):
            if 'IM'!=sc[:2]:
                return S[sc][i,j] + delta(i,j)*gw2/4.*(S['sW'] 
                                                     + S['sHW'] + 2.*S['s2W'])
            else:
                return S[sc][i,j]
        # [eqn (5.8)]
        
        for px in ('','IM'):
            for i,j in S[px+'SBxHl'].keys():
                W[px+'WBxHl'][i,j] = cHf(px+'SBxHl', -1./2., i, j)
                W[px+'WBxHe'][i,j] = cHf(px+'SBxHe', -1., i, j)
                W[px+'WBxHq'][i,j] = cHf(px+'SBxHq', 1./6., i, j)
                W[px+'WBxHu'][i,j] = cHf(px+'SBxHu', 2./3., i, j)
                W[px+'WBxHd'][i,j] = cHf(px+'SBxHd', -1./3., i, j)
                W[px+'WBxHpl'][i,j] = cpHf(px+'sBxHpl', i, j)
                W[px+'WBxHpq'][i,j] = cpHf(px+'sBxHpq', i, j)
            for k in S[px+'SBxHud'].keys():
                W[px+'WBxHud'][k] = S[px+'SBxHud'][k]
                
        def cf(sc, i, j):
            mass = self.mass[ PID[f][i] ]
            yuk = sqrt(2.)*mass/vev
            if 'IM'!=sc[:2]:
                return S[sc][i,j] - delta(i,j)*gw2*yuk*(
                                      S['sW'] + S['sHW'] + S['s2W'])/2.
            else:
                return S[sc][i,j]
                
        # [eqn (5.9)]
        for f in ('u','d','e'): # fermion loop
            for px in ('','IM'):
                wmatrix, smatrix = px+'WBx'+f, px+'SBx'+f
                for i,j in S[smatrix].keys(): # flavour loop
                    W[wmatrix][i,j] = cf(smatrix, i, j)
                    
        # [eqn (5.10)]
        W['cll1221'] = gw2*S['s2W'] # cll1221==0 in SILH
        # derived from [eqn (5.4)]
        W['cll1122']=( gp2*S['s2B'] - gw2*S['s2W'] )/2.
        W['cpuu3333'] = (1./3.)*gs2*S['s2G']
        
        # trivial translation, cX==sX
        trivial = ['sGG','tsGG','s3W','ts3W','s3G','ts3G']
        for coeff in trivial:
            wcoeff = coeff.replace('s','c')
            W[wcoeff] = S[coeff]
        
        return W
        
################################################################################