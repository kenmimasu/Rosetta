import math
from math import sqrt

from ..internal import basis
from ..internal import PID
from ..internal import matrix_mult, matrix_add, matrix_sub, matrix_eq
################################################################################
class mWarsawBasis(basis.Basis):
    '''
    Basis class for Rosetta based on [Grzadkowski et al., JHEP 1010 (2010) 085]. 
    Part of the three intrinsic basis implementations in Rosetta along with the 
    Higgs and modified-SILH bases. The exact list of operators included as well as the 
    equations for the translation to the Higgs basis can be found in the HXSWG 
    note, for which all references to tables and equations are in this 
    implementation. Table 1 of the note shows full list of operators. The block 
    structure of the basis implementation maps to this table. 
    '''
    
    name = 'm-warsaw'
    ###### 
    # [Tab. 1]
    WBV2H2 = ['cGG','cWW','cBB','cWB','tcGG','tcWW','tcBB','tcWB']
    
    WBH4D2 = ['cH','cT']
    
    WBH6 = ['c6H']
    
    WBV3D3 = ['c3W','c3G','tc3W','tc3G']
    
    # affects Gf input, Warsaw <-> SILH translation
    WB4F = ['cll1221','cll1122','cpuu3333'] 
    ######
    blocks = {'WBxV2H2':WBV2H2, 'WBxH4D2':WBH4D2, 'WBxH6':WBH6, 
              'WBxV3D3':WBV3D3, 'WBx4F':WB4F}
              
    flavored = {
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
        'WBxHud' : {'cname':'cHud', 'kind':'general', 'domain':'complex'},
        'WBxeW' : {'cname':'ceW' , 'kind':'general', 'domain':'complex'},
        'WBxeB' : {'cname':'ceB' , 'kind':'general', 'domain':'complex'},
        'WBxuG' : {'cname':'cuG' , 'kind':'general', 'domain':'complex'},
        'WBxuW' : {'cname':'cuW' , 'kind':'general', 'domain':'complex'},
        'WBxuB' : {'cname':'cuB' , 'kind':'general', 'domain':'complex'},
        'WBxdG' : {'cname':'cdG' , 'kind':'general', 'domain':'complex'},
        'WBxdW' : {'cname':'cdW' , 'kind':'general', 'domain':'complex'},
        'WBxdB' : {'cname':'cdB' , 'kind':'general', 'domain':'complex'}   
    }
    
    independent = [c for c in blocks.keys()] + [c for c in flavored.keys()]
    
    # required_masses = set([y for x in PID.values() for y in x.values()])
    required_masses = {25, 24} # Higgs & W masses
    required_inputs = {1, 2, 3, 4, 25} # aEWM1, Gf, MZ, MH, aS(MZ)
        
    def calculate_inputs(self): # calculate a few required EW params from aEWM1, Gf, MZ
        ee2 = 4.*math.pi/self.inputs['aEWM1'] # EM coupling squared at MZ
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
        '''
        Translation function to Mass basis or Higgs basis, which differs only in 
        the prefix of the flavored blocks.
        '''
        # NOTE all assignments to dependent coefficients have been removed, 
        # must now be taken care of in the target basis' calculate_dependents() 
        # function
        # Higgs Basis prefix 
        XB = 'HB'
            
        s2w, c2w, ee2, gw2, gp2, MZ, vev, gs2 = self.calculate_inputs() 
        MH = self.mass[25]
        MW = self.mass[24]
        
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
            M[XB+'xdGLze'][i,j] = (-1./2.*W['WBxHpl'][i,j] 
                                - 1./2.*W['WBxHl'][i,j] + f(-1./2.,-1.,i,j))
            M[XB+'xdGRze'][i,j] = (- 1./2.*W['WBxHe'][i,j] + f(0.,-1.,i,j))
            M[XB+'xdGLzu'][i,j] = (1./2.*W['WBxHpq'][i,j] 
                                   - 1./2.*W['WBxHq'][i,j] + f(1./2.,2./3.,i,j))
            M[XB+'xdGLzd'][i,j] = -1./2.*matrix_mult( W.ckm.dag(),
                                          matrix_mult(
                                             matrix_add(W['WBxHpq'],W['WBxHq']), 
                                             W.ckm)
                                         )[i,j] + f(-1./2.,-1./3.,i,j)
            
            M[XB+'xdGRzu'][i,j] = (- 1./2.*W['WBxHu'][i,j] + f(0.,2./3.,i,j))
            M[XB+'xdGRzd'][i,j] = (- 1./2.*W['WBxHd'][i,j] + f(0.,-1./3.,i,j))
            M[XB+'xdGLwl'][i,j] = (W['WBxHpl'][i,j] + f(1./2.,0.,i,j) 
                                - f(-1./2.,-1.,i,j))
        
        
        for k,v in W['WBxHud'].iteritems():
            M[XB+'xdGRwq'][k] = -v/2.

        M['dCz'] = -W['cH'] - 3.*dv 
        
        # Two derivative field strength interactions [eqn (4.15)]
        M['Cgg']  = W['cGG'] 
        M['Caa']  = W['cWW'] + W['cBB'] - 4.*W['cWB']
        M['Czz']  = ( gw2**2*W['cWW'] + gp2**2*W['cBB'] 
                     + 4.*gw2*gp2*W['cWB'] )/(gw2+gp2)**2
        M['Czbx'] =  -(2./gw2)*(W['cT'] - dv) 
        M['Cza']  = ( gw2*W['cWW'] - gp2*W['cBB'] 
                      - 2.*(gw2-gp2)*W['cWB'] )/(gw2+gp2)

        M['tCgg'] = W['tcGG'] 
        M['tCaa'] =  W['tcWW'] + W['tcBB'] - 4.*W['tcWB']
        M['tCzz'] = ( gw2**2*W['tcWW'] + gp2**2*W['tcBB'] 
                      + 4.*gw2*gp2*W['tcWB'] )/(gw2+gp2)**2
        M['tCza'] = ( gw2*W['tcWW'] - gp2*W['tcBB'] 
                      - 2.*(gw2-gp2)*W['tcWB'] )/(gw2+gp2)
        
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
            for i,j in W[matrix].keys(): 
                diag = delta(i,j)*(W['cH']+dv)
                diag2 = 2.*delta(i,j)*W['cH']
                s_Re, s_Im = W[matrix][i,j].real, W[matrix][i,j].imag
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
        #             Gij, Gji = W['WBx'+f+'G'][i,j], W['WBx'+f+'G'][j,i]
        #             M[glu][i,j] = -sqrt(2.)*(Gij + Gji.conjugate())
        #             M[tglu][i,j] =  -ii*sqrt(2.)*(Gij - Gji.conjugate())
        #
        #         Aij = eta*W['WBx'+f+'W'][i,j] + W['WBx'+f+'B'][i,j]
        #         Aji = eta*W['WBx'+f+'W'][j,i] + W['WBx'+f+'B'][j,i]
        #         M[pho][i,j] = -sqrt(2.)*(Aij + Aji.conjugate())
        #         M[tpho][i,j] =  -ii*sqrt(2.)*(Aij - Aji.conjugate())
        #
        #         Zij = gw2*eta*W['WBx'+f+'W'][i,j] - gp2*W['WBx'+f+'B'][i,j]
        #         Zji = gw2*eta*W['WBx'+f+'W'][j,i] - gp2*W['WBx'+f+'B'][j,i]
        #         M[zed][i,j] = -sqrt(2.)/(gw2+gp2)*(Zij + Zji.conjugate())
        #         M[tzed][i,j] =  -ii*sqrt(2.)/(gw2+gp2)*(Zij - Zji.conjugate())
        #
        ii = complex(0.,1.)
        for f in ('u','d','e'):
            eta = 1 if f=='u' else -1
            
            glu, pho, zed = XB+'xdg'+f,  XB+'xda'+f, XB+'xdz'+f
            
            for i,j in M[zed].keys():
                if f in ('u','d'):
                    M[glu][i,j] = -sqrt(2.)*(W['WBx'+f+'G'][i,j])

                M[pho][i,j] = -(eta*W['WBx'+f+'W'][i,j] + W['WBx'+f+'B'][i,j])

                M[zed][i,j] = -(eta*c2w*W['WBx'+f+'W'][i,j] 
                               - s2w*W['WBx'+f+'B'][i,j])
        
        # TGC's [eqn. (5.18)]
        M['Lz']   = -W['c3W']*3./2.*gw2**2

        M['tLz']  = -W['tc3W']*3./2.*gw2**2

        M['C3g']  = W['c3G']
        
        M['tC3g']  = W['tc3G']
        
        # Higgs cubic interaction [eqn. (4.19)]
        M['dL3']  =  -MH**2/(2.*vev**2) * (3.*W['cH'] + dv) - W['c6H']

        # Four fermion interactions
        M['cll1122'] = W['cll1122']
        M['cpuu3333'] = W['cpuu3333']
        
        return M
    
    @basis.translation('m-silh')
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
            # avoid dependent coefficients
            if (i,j) != (1,1):
                S['SBxHl'][i,j] = sHf('WBxHl', -1./2., i, j)
                S['SBxHpl'][i,j] = spHf('WBxHpl', i, j)
            S['SBxHe'][i,j] = sHf('WBxHe', -1., i, j)
            S['SBxHq'][i,j] = sHf('WBxHq', 1./6., i, j)
            S['SBxHu'][i,j] = sHf('WBxHu', 2./3., i, j)
            S['SBxHd'][i,j] = sHf('WBxHd', -1./3., i, j)
            S['SBxHpq'][i,j] = spHf('WBxHpq', i, j)
        for k in W['WBxHud'].keys():
            S['SBxHud'][k] = W['WBxHud'][k]            

        def sf(wc, i, j):
            return W[wc][i,j] - delta(i,j)*(W['cll1221'] 
                                           - 4.*W['WBxHpl'][1,1,].real)/sqrt(2.)
            
        # [eqn (5.9)]
        for f in ('u','d','e'): # fermion loop
            wmatrix, smatrix = 'WBx'+f, 'SBx'+f
            for i,j in W[wmatrix].keys(): # flavor loop
                S[smatrix][i,j] = sf(wmatrix, i, j)

        # trivial translation, sX==cX
        trivial = ['sGG','tsGG','s3W','ts3W','s3G','ts3G']
        for coeff in trivial:
            wcoeff = coeff.replace('s','c')
            S[coeff] = W[wcoeff]
            
        trivial2 = ['xeW','xeB','xuG','xuW','xuB','xdG','xdW','xdB']
        
        for coeff in trivial2:
            matrix_eq(W['WB'+coeff], S['SB'+coeff])

        return S
        