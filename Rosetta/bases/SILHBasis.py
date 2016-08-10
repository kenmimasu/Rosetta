import math
from math import sqrt

from ..internal import basis
from ..internal import PID
from ..internal import matrix_mult, matrix_add, matrix_sub, matrix_eq
################################################################################
class SILHBasis(basis.Basis):
    '''
    Basis class for Rosetta based on [Giudice et al., JHEP 0706 (2007) 045 ] 
    and a number of successive publications. Part of the three intrinsic basis 
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
    # [Tab. 1]
    SBV2H2 = ['sg','sa','sW','sB','sHW','sHB','s2W','s2B','s2G',
              'tsg','tsa','tsHW','tsHB']
    
    SBH4D2 = ['sH','sT']
    
    SBH6 = ['s6']
    
    SBV3D3 = ['s3W','s3G','ts3W','ts3G']

    ##########################
    # block structure
    blocks = {'SBxV2H2':SBV2H2, 'SBxH4D2':SBH4D2, 
              'SBxH6':SBH6, 'SBxV3D3':SBV3D3} 
    
    flavored={
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
        'SBxHud' : {'cname':'sHud', 'kind':'general', 'domain':'complex'},
        'SBxeW' : {'cname':'seW' , 'kind':'general', 'domain':'complex'},
        'SBxeB' : {'cname':'seB' , 'kind':'general', 'domain':'complex'},
        'SBxuG' : {'cname':'suG' , 'kind':'general', 'domain':'complex'},
        'SBxuW' : {'cname':'suW' , 'kind':'general', 'domain':'complex'},
        'SBxuB' : {'cname':'suB' , 'kind':'general', 'domain':'complex'},
        'SBxdG' : {'cname':'sdG' , 'kind':'general', 'domain':'complex'},
        'SBxdW' : {'cname':'sdW' , 'kind':'general', 'domain':'complex'},
        'SBxdB' : {'cname':'sdB' , 'kind':'general', 'domain':'complex'}
    }

    independent = blocks.keys() + flavored.keys()
                  
    # two coefficients are zero by construction (set in calculate_dependent())
    dependent = ['CsHl1x1', 'CsHpl1x1']
    
    # required_masses = set([y for x in PID.values() for y in x.values()])
    required_masses = {25, 24, 1, 2, 3, 4, 5, 6, 11, 13, 15}  # Higgs & W masses
    required_inputs = {1, 2, 3, 4, 25} # aEWM1, Gf, aS, MZ, MH

    ##########################
    def calculate_dependent(self):
        # These coefficients are implictly set to zero
        self['SBxHl'][1,1], self['SBxHpl'][1,1] = 0., 0.
        
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
        
        
    @basis.translation('warsaw')
    def to_warsaw(self, instance):

        s2w, c2w, ee2, gw2, gp2, MZ, vev, gs2 = self.calculate_inputs() 
        
        gw, gp, gs = sqrt(gw2), sqrt(gp2), sqrt(gs2)
        
        MH = self.mass[25]
        lam = MH**2/(2.*vev**2) # Higgs self-coupling     
           
        S = self
        W = instance
        
        def delta(i,j):
            return 1. if i==j else 0.
        
        W['cH'] = ((S['s2W'] + S['sW'] + S['sHW'])*8. - S['s6'])*lam
        
        W['cHbx'] = ( (S['s2W'] + S['sW'] + S['sHW'])*3.
                     + (S['s2B'] + S['sB'] + S['sHB'])*gp2/gw2
                     - (S['sH'] + S['sT'])/2.)
        
        W['cHD'] = (S['s2B'] + S['sB'] + S['sHB'])*4.*gp2/gw2 - S['sT']*2.
        
        W['cHG'] = 4.*gs2/gw2*S['sg']
        
        W['cHW'] =  -S['sHW']
        
        W['cHB'] = (4.*S['sa'] - S['sHB'])*gp2/gw2
        
        W['cHWB'] = -(S['sHW'] + S['sHB'])*gp/gw
    
        W['tcHG'] = 4.*gs2/gw2*S['tsg']
        
        W['tcHW'] =  -S['tsHW']
        
        W['tcHB'] = (4.*S['tsa'] - S['tsHB'])*gp2/gw2
        
        W['tcHWB'] = -(S['tsHW'] + S['tsHB'])*gp/gw
        
        W['cG'] = 4.*gs**3/gw2*S['s3G']
        
        W['cW'] = 4.*gw*S['s3W']
    
        W['tcG'] = 4.*gs**3/gw2*S['ts3G']
    
        W['tcW'] = 4.*gw*S['ts3W']
        
        # Gauge current operators
        def cH1f(sc, Yf, i, j):
            try: # special case for SBxHl[1,1] which is dependent in SILH
                ss = S[sc][i,j]
            except KeyError:
                ss = 0.
            coeff = delta(i,j)*2.*Yf*gp2/gw2
            return  ss + (2.*S['s2B']+ S['sB'] + S['sHB'])*coeff
            
        def cH3f(sc, i, j):
            try: # special case for SBxHpl[1,1] which is dependent in SILH
                ss = S[sc][i,j]
            except KeyError:
                ss = 0.
            return ss + delta(i,j)*(2.*S['s2W'] + S['sW'] + S['sHW'])
        
        for i,j in S['SBxHq'].keys():
            W['WBxH1l'][i,j] = cH1f('SBxHl', -1./2., i, j)
            W['WBxHe'][i,j] = cH1f('SBxHe', -1., i, j)
            W['WBxH1q'][i,j] = cH1f('SBxHq', 1./6., i, j)
            W['WBxHu'][i,j] = cH1f('SBxHu', 2./3., i, j)
            W['WBxHd'][i,j] = cH1f('SBxHd', -1./3., i, j)
            W['WBxH3l'][i,j] = cH3f('SBxHpl', i, j)
            W['WBxH3q'][i,j] = cH3f('SBxHpq', i, j)
            
        for k in S['SBxHud'].keys():
            W['WBxHud'][k] = S['SBxHud'][k]

        # Yukawa operators
        for f in ('u','d','e'): # fermion loop
            wmatrix, smatrix = 'WBx'+f+'H', 'SBx'+f
            for i,j in S[smatrix].keys(): # flavor loop
            
                mi, mj = self.mass[PID[f][i]], self.mass[PID[f][j]]
                
                diag = 2.*delta(i,j)*(S['s2W'] + S['sW'] + S['sHW'])
                
                W[wmatrix][i,j] = sqrt(2.*mi*mj)/vev*(S[smatrix][i,j] + diag)

        # Dipole operators
        for f in ('u','d','e'):
            eta = 1 if f=='u' else -1
            for i,j in S['SBx'+f+'W'].keys():
                
                mi, mj = self.mass[PID[f][i]], self.mass[PID[f][j]]
                
                MFVnorm = gw2*vev/sqrt(2.*mi*mj)/4.
                # gluon
                if f in ('u','d'):
                    W['WBx'+f+'G'][i,j] = S['SBx'+f+'G'][i,j]/MFVnorm
                # Weak
                W['WBx'+f+'W'][i,j] = S['SBx'+f+'W'][i,j]/MFVnorm
                # Hypercharge
                W['WBx'+f+'B'][i,j] = S['SBx'+f+'B'][i,j]/MFVnorm

        # cll1221==0 in SILH
        W['cll1221'] = 4.*S['s2W'] 

        W['cll1122']=  (S['s2B']*gp2/gw2 - S['s2W'])*2.
        
        W['cpuu3333'] = (1./3.)*gs2*S['s2G']
        
        return W
        
        
    @basis.translation('higgs')
    def to_higgs(self, instance, target='bsmc'):
        '''
        Translation function to  Higgs basis.
        '''
        # Higgs Basis prefix
        XB = 'HB'
        
        s2w, c2w, ee2, gw2, gp2, MZ, vev, gs2 = self.calculate_inputs() 
        
        MH = self.mass[25]
        lam = MH**2/(2.*vev**2) # Higgs self-coupling     
        
        MW = self.mass[24]
        # MW = MZ*sqrt(c2w)
                
        S = self
        M = instance
        

        M['dM'] = -( S['sW'] + S['sB'] + S['s2W'] + S['s2B'] 
                   - S['sT']*gw2/(2.*gp2)
                   + S['SBxHpl'][2,2].real/2. )*gp2/(gw2-gp2)
        
        def delta(i,j):
            return 1. if i==j else 0.
            
        def f(T3,Q,i,j):
            
            Qcoeff = ( S['s2B']*(2.*gw2-gp2)/gw2 + S['s2W'] + S['sW'] + S['sB'] 
                      - S['sT']/2. + S['SBxHpl'][2,2].real/2. )*gp2/(gw2-gp2)
                      
            T3coeff = ( S['s2W'] + S['s2B']*gp2/gw2 + S['sT']/2. 
                      - S['SBxHpl'][2,2].real/2. )
                      
            return  delta(i,j)*( T3*T3coeff - Q*Qcoeff )
        
        # W/Z chiral coupling deviations
        for i,j in S['SBxHpl'].keys():

            M[XB+'xdGLze'][i,j] = (-S['SBxHpl'][i,j]/2. - S['SBxHl'][i,j]/2. 
                                  + f(-1./2.,-1.,i,j) )
                                
            M[XB+'xdGRze'][i,j] = -S['SBxHe'][i,j]/2. + f(0.,-1.,i,j)
            
            M[XB+'xdGLzu'][i,j] = ( S['SBxHpq'][i,j]/2.- S['SBxHq'][i,j]/2.
                                  + f(1./2.,2./3.,i,j))
            
            M[XB+'xdGLzd'][i,j] = -1./2.*matrix_mult( S.ckm.dag(),
                                          matrix_mult(
                                             matrix_add(S['SBxHpq'],S['SBxHq']), 
                                             S.ckm)
                                         )[i,j] + f(-1./2.,-1./3.,i,j)
            
            M[XB+'xdGRzu'][i,j] = (- 1./2.*S['SBxHu'][i,j] + f(0.,2./3.,i,j))
            
            M[XB+'xdGRzd'][i,j] = (- 1./2.*S['SBxHd'][i,j] + f(0.,-1./3.,i,j))
            
            M[XB+'xdGLwl'][i,j] = (S['SBxHpl'][i,j] + f(1./2.,0.,i,j) 
                                - f(-1./2.,-1.,i,j))

        for k,v in S['SBxHud'].iteritems():
            M[XB+'xdGRwq'][k] = -v/2.
        
        
        # Dipole interactions
        ii = complex(0.,1.)
        for f in ('u','d','e'):
            eta = 1 if f=='u' else -1
            
            glu, pho, zed = XB+'xdg'+f,  XB+'xda'+f, XB+'xdz'+f
            
            for i,j in M[zed].keys():
                if f in ('u','d'):
                    M[glu][i,j] = -16./gw2*(S['SBx'+f+'G'][i,j])

                M[pho][i,j] = -16./gw2*(eta*S['SBx'+f+'W'][i,j] 
                                       + S['SBx'+f+'B'][i,j])

                M[zed][i,j] = -16./gw2*(eta*c2w*S['SBx'+f+'W'][i,j] 
                                       - s2w*S['SBx'+f+'B'][i,j])
        
        # Triple gauge couplings [eqn. (5.18)]
        M['Lz'] = -S['s3W']*6.*gw2
        M['tLz'] = -S['ts3W']*6.*gw2
        M['C3g']  = S['s3G']*4./gw2
        M['tC3g']  = S['ts3G']*4./gw2
        
        def dy_sf(X,Y): 
            '''
            Return solution for:
                dy*cos(phi) == X
                dy*sin(phi) == Y
            as ( dY, sin(phi) )
            '''
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
        
        # Yukawa type interaction
        for f in ('u','d','e'):
            matrix = 'SBx' + f
            for i,j in S['SBx'+f].keys(): 
                diag = delta(i,j)*(S['sH'] + S['SBxHpl'][2,2].real)/2.

                s_Re, s_Im = S[matrix][i,j].real, S[matrix][i,j].imag
                
                dy_cosphi = -s_Re - diag
                dy_sinphi = s_Im
                
                M[XB+'xdY'+f][i,j], M[XB+'xS'+f][i,j] = dy_sf(dy_cosphi, 
                                                              dy_sinphi)
        
        # Higgs couplings to Z 
        M['dCz'] = -S['sH']/2. - 3./2.*S['SBxHpl'][2,2].real
        
        # Two derivative field strength interactions
        M['Cgg']  = S['sg']*(16./gw2)
        
        M['Caa']  = S['sa']*(16./gw2)
        
        M['Czz']  = -4./(gw2+gp2)*( S['sHW'] + S['sHB']*gp2/gw2 
                                  - S['sa']*4.*gp2/gw2*s2w )
                                     
        M['Czbx'] =  2./gw2*( (S['sW'] + S['sHW'] +S['s2W'])
                            + (S['sB'] + S['sHB'] +S['s2B'])*gp2/gw2
                            - S['sT']/2. + S['SBxHpl'][2,2].real/2. ) 
                            
        M['Cza']  = 2./gw2*( S['sHB'] - S['sHW'] - S['sa']*8.*s2w)
        
        M['tCgg']  = S['tsg']*(16./gw2)
        
        M['tCaa']  = S['tsa']*(16./gw2)
        
        M['tCzz']  = -4./(gw2+gp2)*( S['tsHW'] + S['tsHB']*gp2/gw2 
                                  - S['tsa']*4.*gp2/gw2*s2w )
                            
        M['tCza']  = 2./gw2*( S['tsHB'] - S['tsHW'] - S['tsa']*8.*s2w)
        
        
        # # Higgs cubic interaction
        M['dL3']  =  lam*( S['s6'] - 3./2.*S['sH'] - S['SBxHpl'][2,2].real/2.)

        # Four fermion interactions
        M['cll1122'] = (S['s2B']*gp2/gw2 - S['s2W'])*2.
        
        M['cpuu3333'] = S['s2G']*gs2/3.
        
        return M
        
################################################################################