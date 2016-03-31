import math
from math import sqrt

from ..internal import basis
from ..internal import PID
from ..internal import matrix_mult, matrix_add, matrix_sub, matrix_eq
################################################################################
class mSILHBasis(basis.Basis):
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

    name = 'm-silh'
    ##########################
    # declare coefficients
    # [Tab. 1]
    SBV2H2 = ['sGG','sBB','sW','sB','sHW','sHB','s2W','s2B','s2G',
              'tsGG','tsBB','tsHW','tsHB']
    
    SBH4D2 = ['sH','sT']
    
    SBH6 = ['s6H']
    
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
    required_masses = {25, 24}  # Higgs & W masses
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
        
    @basis.translation('higgs')
    def to_higgs(self, instance, target='bsmc'):
        '''
        Translation function to  Higgs basis.
        '''
        # Higgs Basis prefix
        XB = 'HB'
        
        s2w, c2w, ee2, gw2, gp2, MZ, vev, gs2 = self.calculate_inputs() 
        MH = self.mass[25]
        MW = self.mass[24]
        # MW = MZ*sqrt(c2w)
                
        S = self
        M = instance
        
        # W mass shift [eqn (5.11)]
        M['dM'] = - gw2*gp2/(4.*(gw2-gp2))*(S['sW'] + S['sB'] + S['s2W'] 
                                            + S['s2B'] - 4./gp2*S['sT'] 
                                            + 2./gw2*S['SBxHpl'][2,2].real)
        def f(T3,Q,i,j): # [eqn (5.12)]
            Qcoeff = gp2/4./(gw2-gp2)*( -(2.*gw2-gp2)*S['s2B'] 
                       - gw2*(S['s2W'] + S['sW'] + S['sB'] ) 
                       + 4.*S['sT'] - 2.*S['SBxHpl'][2,2].real)
            T3coeff = ( gw2*S['s2W'] + gp2*S['s2B'] 
                        + 4.*S['sT']-2.*S['SBxHpl'][2,2].real)/4.
            return  delta(i,j)*(Q*Qcoeff + T3*T3coeff)
            
        def delta(i,j):
            return 1. if i==j else 0.
            
        # W/Z chiral coupling deviations
        for i,j in S['SBxHpl'].keys():
            # M[XB+'xdGLzv'][i,j] = (1./2.*S['SBxHpl'][i,j]
            #                     - 1./2.*S['SBxHl'][i,j] + f(1./2.,0.,i,j))
            M[XB+'xdGLze'][i,j] = (-1./2.*S['SBxHpl'][i,j] 
                                - 1./2.*S['SBxHl'][i,j] + f(-1./2.,-1.,i,j))
            M[XB+'xdGRze'][i,j] = (- 1./2.*S['SBxHe'][i,j] + f(0.,-1.,i,j))
            M[XB+'xdGLzu'][i,j] = (1./2.*S['SBxHpq'][i,j] 
                                - 1./2.*S['SBxHq'][i,j] + f(1./2.,2./3.,i,j))
            
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
        
        onames = ['xdGLze', 'xdGRze', 'xdGLzv', 'xdGLzu', 'xdGRzu', 
                  'xdGLzd', 'xdGRzd', 'xdGLwl', 'xdGLwq', 'xdGRwq']
                  
        vertex = [XB + o for o in onames]
        
        # Higgs couplings to Z [eqn (4.14)]
        M['dCz'] = -S['sH'] - 3./2.*S['SBxHpl'][2,2].real 
        
        # Two derivative field strength interactions [eqn (4.15)]
        M['Cgg']  = S['sGG']
        M['Caa']  = S['sBB']
        M['Czz']  = -1./(gw2+gp2)*( gw2*S['sHW'] + gp2*S['sHB'] 
                                     - gp2*s2w*S['sBB'] )
        M['Czbx'] =  1./(2.*gw2)*( gw2*(S['sW'] + S['sHW'] +S['s2W'])
                                  + gp2*(S['sB'] + S['sHB'] +S['s2B'])
                                  - 4.*S['sT'] + 2.*S['SBxHpl'][2,2].real ) 
        M['Cza']  = ( S['sHB'] - S['sHW'])/2. - s2w*S['sBB']
        
        M['tCgg'] = S['tsGG'] 
        M['tCaa'] = S['tsBB']
        M['tCzz'] = -1./(gw2+gp2)*( gw2*S['tsHW'] + gp2*S['tsHB'] 
                                     - gp2*s2w*S['tsBB'] )
        M['tCza'] = ( S['tsHB'] - S['tsHW'])/2. - s2w*S['tsBB']        
        
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
            for i,j in S['SBx'+f].keys(): 
                diag = delta(i,j)*(S['sH'] + S['SBxHpl'][2,2].real/2.)
                diag2 = 2.*delta(i,j)*(S['sH'])
                s_Re, s_Im = S[matrix][i,j].real, S[matrix][i,j].imag
                dy_cosphi = s_Re/sqrt(2.) - diag
                dy_sinphi = s_Im/sqrt(2.)
                
                M[XB+'xdY'+f][i,j], M[XB+'xS'+f][i,j] = dy_sf(dy_cosphi, 
                                                              dy_sinphi)

        # # Dipole interactions
        # ii = complex(0.,1.)
        # for f in ('u','d','e'):
        #     eta = 1 if f=='u' else -1
        #
        #     glu, tglu = XB+'xdg'+f, XB+'xtdg'+f
        #     pho, tpho = XB+'xda'+f, XB+'xtda'+f
        #     zed, tzed = XB+'xdz'+f, XB+'xtdz'+f
        #     dub = XB+'xdwl' if f=='e' else XB+'xdw'+f
        #
        #     for i,j in M[zed].keys():
        #         if f in ('u','d'):
        #             Gij, Gji = S['SBx'+f+'G'][i,j], S['SBx'+f+'G'][j,i]
        #             M[glu][i,j] = -sqrt(2.)*(Gij + Gji.conjugate())
        #             M[tglu][i,j] =  -ii*sqrt(2.)*(Gij - Gji.conjugate())
        #
        #         Aij = eta*S['SBx'+f+'W'][i,j] + S['SBx'+f+'B'][i,j]
        #         Aji = eta*S['SBx'+f+'W'][j,i] + S['SBx'+f+'B'][j,i]
        #         M[pho][i,j] = -sqrt(2.)*(Aij + Aji.conjugate())
        #         M[tpho][i,j] =  -ii*sqrt(2.)*(Aij - Aji.conjugate())
        #
        #         Zij = gw2*eta*S['SBx'+f+'W'][i,j] - gp2*S['SBx'+f+'B'][i,j]
        #         Zji = gw2*eta*S['SBx'+f+'W'][j,i] - gp2*S['SBx'+f+'B'][j,i]
        #         M[zed][i,j] = -sqrt(2.)/(gw2+gp2)*(Zij + Zji.conjugate())
        #         M[tzed][i,j] =  -ii*sqrt(2.)/(gw2+gp2)*(Zij - Zji.conjugate())
        
        # Dipole interactions
        ii = complex(0.,1.)
        for f in ('u','d','e'):
            eta = 1 if f=='u' else -1
            
            glu, pho, zed = XB+'xdg'+f,  XB+'xda'+f, XB+'xdz'+f
            
            for i,j in M[zed].keys():
                if f in ('u','d'):
                    M[glu][i,j] = -sqrt(2.)*(S['SBx'+f+'G'][i,j])

                M[pho][i,j] = -(eta*S['SBx'+f+'W'][i,j] + S['SBx'+f+'B'][i,j])

                M[zed][i,j] = -(eta*c2w*S['SBx'+f+'W'][i,j] 
                               - s2w*S['SBx'+f+'B'][i,j])

        # Triple gauge couplings [eqn. (5.18)]
        M['Lz'] = -S['s3W']*3./2.*gw2**2
        M['tLz'] = -S['ts3W']*3./2.*gw2**2
        M['C3g']  = S['s3G']
        M['tC3g']  = S['ts3G']
        
        # # Higgs cubic interaction [eqn. (4.19)]
        M['dL3']  =  -MH**2/(2.*vev**2)*(3.*S['sH'] + S['SBxHpl'][2,2].real/2.) - S['s6H']
        M['dL4'] = 3./2.*M['dL3'] - MH**2/vev**2/6.*M['dCz']

        # Four fermion interactions
        M['cll1122'] = (gp2*S['s2B'] - gw2*S['s2W'])/2.
        M['cpuu3333'] = (1./3.)*gs2*S['s2G']
        
        return M
    
    @basis.translation('m-warsaw')
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
            try: # special case for SBxHl[1,1] which is dependent in SILH
                ss = S[sc][i,j]
            except KeyError:
                ss = 0.
            return  ss + delta(i,j)*gp2*Yf/2.*(S['sB'] + S['sHB'] + 2.*S['s2B'])
            
        def cpHf(sc, i, j):
            try: # special case for SBxHpl[1,1] which is dependent in SILH
                ss = S[sc][i,j]
            except KeyError:
                ss = 0.
            return ss + delta(i,j)*gw2/4.*(S['sW'] + S['sHW'] + 2.*S['s2W'])
        
        for i,j in S['SBxHq'].keys():
            W['WBxHl'][i,j] = cHf('SBxHl', -1./2., i, j)
            W['WBxHe'][i,j] = cHf('SBxHe', -1., i, j)
            W['WBxHq'][i,j] = cHf('SBxHq', 1./6., i, j)
            W['WBxHu'][i,j] = cHf('SBxHu', 2./3., i, j)
            W['WBxHd'][i,j] = cHf('SBxHd', -1./3., i, j)
            W['WBxHpl'][i,j] = cpHf('sBxHpl', i, j)
            W['WBxHpq'][i,j] = cpHf('sBxHpq', i, j)
        for k in S['SBxHud'].keys():
            W['WBxHud'][k] = S['SBxHud'][k]
        
        def cf(sc, i, j):
            return S[sc][i,j] - delta(i,j)*gw2*(S['sW'] + 
                                               S['sHW'] + S['s2W'])/sqrt(2.)
                
        # [eqn (5.9)]
        for f in ('u','d','e'): # fermion loop
            wmatrix, smatrix = 'WBx'+f, 'SBx'+f
            for i,j in S[smatrix].keys(): # flavor loop
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
        
        trivial2 = ['xeW','xeB','xuG','xuW','xuB','xdG','xdW','xdB']
        for coeff in trivial2:
            matrix_eq(S['SB'+coeff], W['WB'+coeff])
        
        return W
        
################################################################################