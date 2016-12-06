import math
from math import sqrt

from ..internal import basis
from ..internal import PID
from ..internal import matrix_mult, matrix_add, matrix_sub, matrix_eq
################################################################################
class WarsawBasis(basis.Basis):
    '''
    Basis class for Rosetta based on [Grzadkowski et al., JHEP 1010 (2010) 085]. 
    Part of the three intrinsic basis implementations in Rosetta along with the 
    Higgs and modified-SILH bases. The exact list of operators included as well as the 
    equations for the translation to the Higgs basis can be found in the HXSWG 
    note, for which all references to tables and equations are in this 
    implementation.
    '''
    ###### 
    name = 'warsaw'
    ###### 
    # [Tab. 1]
    WBX2H2 = ['cHG','cHW','cHB','cHWB','tcHG','tcHW','tcHB','tcHWB']
    
    WBH6 = ['cH']

    WBH4D2 = ['cHbx', 'cHD']
    
    WBX3 = ['cW','cG','tcW','tcG']
    
    # cll1221 affects Gf input
    # cll1221, cll1122 & cpuu3333 needed for Warsaw <-> SILH translation
    WB4F = ['cll1111','cll1122','cll1221','cll1133','cll1331','cll2332',
            'cle1111','cle1122','cle2211','cle1133','cle3311',
            'cee1111','cee1122','cee1133',
            'cpuu3333']
    
    ######
    blocks = {'WBxX2H2':WBX2H2, 'WBxH4D2':WBH4D2, 'WBxH6':WBH6, 
              'WBxX3':WBX3, 'WBx4F':WB4F}
              
    flavored = {
        # Yukawa operators
        'WBxuH': {'cname':'cuH', 'kind':'general', 'domain':'complex'},
        'WBxdH': {'cname':'cdH', 'kind':'general', 'domain':'complex'},
        'WBxeH': {'cname':'ceH', 'kind':'general', 'domain':'complex'},
        # Gauge current operators
        'WBxH1l': {'cname':'cH1l', 'kind':'hermitian', 'domain':'complex'},
        'WBxH3l': {'cname':'cH3l', 'kind':'hermitian', 'domain':'complex'},
        'WBxHe': {'cname':'cHe' , 'kind':'hermitian', 'domain':'complex'},
        'WBxH1q': {'cname':'cH1q' , 'kind':'hermitian', 'domain':'complex'},
        'WBxH3q': {'cname':'cH3q', 'kind':'hermitian', 'domain':'complex'},
        'WBxHu': {'cname':'cHu' , 'kind':'hermitian', 'domain':'complex'},
        'WBxHd' : {'cname':'cHd' , 'kind':'hermitian', 'domain':'complex'},
        'WBxHud' : {'cname':'cHud', 'kind':'general', 'domain':'complex'},
        # Dipole operators
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
    
    # Higgs, W & fermion  masses
    required_masses = {25, 24, 1, 2, 3, 4, 5, 6, 11, 13 ,15} 
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

    @basis.translation('silh')
    def to_silh(self, instance):
        s2w, c2w, ee2, gw2, gp2, MZ, vev, gs2 = self.calculate_inputs()
        gw, gp, gs = sqrt(gw2), sqrt(gp2), sqrt(gs2)
        
        MH = self.mass[25]
        lam = MH**2/(2.*vev**2) # Higgs self-coupling

        W = self
        S = instance

        def delta(i,j):
            return 1. if i==j else 0.

        S['sg'] = W['cHG']*gw2/(4.*gs2)

        S['sa'] = (W['cHW'] + W['cHB']*gw2/gp2 - W['cHWB']*gw/gp)/4.
        
        S['tsg'] = W['tcHG']*gw2/(4.*gs2)

        S['tsa'] = (W['tcHW'] + W['tcHB']*gw2/gp2 - W['tcHWB']*gw/gp)/4.

        S['sW'] = W['cHW'] + W['WBxH3l'][1,1].real - W['cll1221']/2.

        S['sB'] = (- W['cHW'] + W['cHWB']*gw/gp - W['WBxH1l'][1,1].real*gw2/gp2
                   - W['cll1221']*gw2/gp2/2. - W['cll1122']*gw2/gp2)

        S['sH'] = ( W['WBxH3l'][1,1].real*6. - W['cHbx']*2. + W['cHD']/2. 
                  - W['cll1221']*3./2. )

        S['sT'] = (- W['cHD']/2. - W['WBxH1l'][1,1].real*2.
                   - W['cll1221']/2. - W['cll1122'])

        S['s6'] = - W['cH']/lam + W['WBxH3l'][1,1,].real*8. - W['cll1221']*2.

        S['sHW'] = - W['cHW']

        S['sHB'] = W['cHW'] - W['cHWB']*gw/gp
        
        S['tsHW'] = - W['tcHW']

        S['tsHB'] = W['tcHW'] - W['tcHWB']*gw/gp

        S['s2W'] = W['cll1221']/4.
        
        S['s2B'] = ( W['cll1221']/2. + W['cll1122'] )*gw2/gp2/2.

        S['s2G'] = W['cpuu3333']*3./gs2

        S['s3W'] = W['cW']/gw/4.
        
        S['s3G'] = W['cG']*gw2/gs**3/4.
        
        S['ts3W'] = W['tcW']/gw/4.
        
        S['ts3G'] = W['tcG']*gw2/gs**3/4.

        # Gauge current operators
        def sHf(wc, Yf, i, j):
                return W[wc][i,j] + delta(i,j)*2.*Yf*W['WBxH1l'][1,1,].real

        def spHf(wc, i, j):
                return W[wc][i,j] - delta(i,j)*W['WBxH3l'][1,1,].real

        for i,j in W['WBxH1l'].keys():
            # avoid dependent coefficients that are zero in SILH
            if (i,j) != (1,1):
                S['SBxHl'][i,j] = sHf('WBxH1l', -1./2., i, j)
                S['SBxHpl'][i,j] = spHf('WBxH3l', i, j)
            S['SBxHe'][i,j] = sHf('WBxHe', -1., i, j)
            S['SBxHq'][i,j] = sHf('WBxH1q', 1./6., i, j)
            S['SBxHu'][i,j] = sHf('WBxHu', 2./3., i, j)
            S['SBxHd'][i,j] = sHf('WBxHd', -1./3., i, j)
            S['SBxHpq'][i,j] = spHf('WBxH3q', i, j)
            
        for k in W['WBxHud'].keys():
            S['SBxHud'][k] = W['WBxHud'][k]

        # Yukawa operators
        for f in ('u','d','e'): # loop over fermions
            wmatrix, smatrix = 'WBx'+f+'H', 'SBx'+f
            for i,j in W[wmatrix].keys(): # flavor loop
                mi, mj = self.mass[PID[f][i]], self.mass[PID[f][j]]
                
                diag = delta(i,j)*(W['cll1221']/2. - W['WBxH3l'][1,1,].real*2.)
                
                S[smatrix][i,j] = vev/sqrt(2.*mi*mj)*W[wmatrix][i,j] + diag
        
        # Dipole operators
        for f in ('u','d','e'):
            eta = 1 if f=='u' else -1
            for i,j in S['SBx'+f+'W'].keys():
                
                mi, mj = self.mass[PID[f][i]], self.mass[PID[f][j]]
                
                MFVnorm = gw2*vev/sqrt(2.*mi*mj)/4.
                # gluon
                if f in ('u','d'):
                    S['SBx'+f+'G'][i,j] = W['WBx'+f+'G'][i,j]*MFVnorm
                # Weak
                S['SBx'+f+'W'][i,j] = W['WBx'+f+'W'][i,j]*MFVnorm
                # Hypercharge
                S['SBx'+f+'B'][i,j] = W['WBx'+f+'B'][i,j]*MFVnorm

        # Four fermion operators
        sll = ['sll1111','sll1133','sll1331','sll2332']
        for coeff in sll:
            i,j,k,l = [ind for ind in coeff[-4:]]
            if i==j==k==l:
                S[coeff] = W['c'+coeff[1:]] - (S['s2B']*gp2/gw2 + S['s2W'])
            elif i==j and k==l:
                S[coeff] = W['c'+coeff[1:]] - (S['s2B']*gp2/gw2 - S['s2W'])*2.
            elif i==l and j==k:
                S[coeff] = W['c'+coeff[1:]] - 4.*S['s2W']
                
        trivial_4f = ['sle1111','sle1122','sle2211','sle1133','sle3311',
                      'see1111','see1122','see1133']

        for c in trivial_4f:
            S[c] = W['c'+c[1:]]

        return S


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
        gw, gp, gs = sqrt(gw2), sqrt(gp2), sqrt(gs2)
        
        MH = self.mass[25]
        lam = MH**2/(2.*vev**2) # Higgs self-coupling
        
        MW = self.mass[24]
        
        W = self
        M = instance
        
        M['dM'] = ( -W['cHWB']*gw*gp - W['cHD']*gw2/4.
                  + ( W['cll1221'] - 2.*W['WBxH3l'][1,1].real 
                    - 2.*W['WBxH3l'][2,2].real)*gp2/4.
                  )/(gw2-gp2)
        
        
        def delta(i,j):
            return 1. if i==j else 0.
                    

        def f(T3,Q,i,j): # [eqn (4.11)]
            if i==j:
                Acoeff = - gw*gp/(gw2-gp2)*W['cHWB']
                Zcoeff = (W['cll1221']/4.- W['WBxH3l'][1,1].real/2. 
                         - W['WBxH3l'][2,2].real/2. - W['cHD']/4.)
                return Acoeff*Q + Zcoeff*(T3 + Q*gp2/(gw2-gp2))
            else: 
                return 0.

        # W/Z chiral coupling deviations
        for i,j in W['WBxH3l'].keys():
            M[XB+'xdGLwl'][i,j] = (W['WBxH3l'][i,j] + f(1./2.,0.,i,j) 
                                  - f(-1./2.,-1.,i,j))
            
            M[XB+'xdGLze'][i,j] = (-1./2.*W['WBxH3l'][i,j] 
                                - 1./2.*W['WBxH1l'][i,j] + f(-1./2.,-1.,i,j))
                                
                                
            M[XB+'xdGRze'][i,j] = (- 1./2.*W['WBxHe'][i,j] + f(0.,-1.,i,j))
            
            M[XB+'xdGLzu'][i,j] = (1./2.*W['WBxH3q'][i,j] 
                                  - 1./2.*W['WBxH1q'][i,j] + f(1./2.,2./3.,i,j))
                                  
            M[XB+'xdGLzd'][i,j] = -1./2.*matrix_mult( W.ckm.dag(),
                                          matrix_mult(
                                             matrix_add(W['WBxH3q'],W['WBxH1q']), 
                                             W.ckm)
                                         )[i,j] + f(-1./2.,-1./3.,i,j)
            
            M[XB+'xdGRzu'][i,j] = (- 1./2.*W['WBxHu'][i,j] + f(0.,2./3.,i,j))

            M[XB+'xdGRzd'][i,j] = (- 1./2.*W['WBxHd'][i,j] + f(0.,-1./3.,i,j))

        for k,cHud_ij in W['WBxHud'].iteritems():
            M[XB+'xdGRwq'][k] = -cHud_ij/2.
        
        
        # Dipole interactions
        for f in ('u','d','e'):
            eta = 1 if f=='u' else -1
            
            glu, pho, zed = XB+'xdg'+f,  XB+'xda'+f, XB+'xdz'+f
            
            for i,j in M[zed].keys():
                mi, mj = self.mass[PID[f][i]], self.mass[PID[f][j]]
                
                MFVnorm = 2.*sqrt(2.)*vev/sqrt(mi*mj)
                
                if f in ('u','d'):
                    M[glu][i,j] = -MFVnorm*(W['WBx'+f+'G'][i,j])

                M[pho][i,j] = -MFVnorm*(eta*W['WBx'+f+'W'][i,j] 
                                       + W['WBx'+f+'B'][i,j])

                M[zed][i,j] = -MFVnorm/(gw2+gp2)*(eta*gw2*W['WBx'+f+'W'][i,j] 
                                                 - gp2*W['WBx'+f+'B'][i,j])
        
        # TGCs
        M['Lz']   = -3./2.*gw*W['cW']

        M['tLz']  = -3./2.*gw*W['tcW']

        M['C3g']  = W['cG']/gs**3
        
        M['tC3g']  = W['tcG']/gs**3
        
        # Yukawa type interaction
 
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

        for f in ('u','d','e'):
            matrix = 'WBx' + f +'H'
            for i,j in W[matrix].keys(): 
                
                diag = delta(i,j)*( W['cll1221']/4. + W['cHbx'] - W['cHD']/4.
                        - W['WBxH3l'][1,1].real/2. - W['WBxH3l'][2,2].real/2. )
                s_Re, s_Im = W[matrix][i,j].real, W[matrix][i,j].imag
                
                mi, mj = self.mass[PID[f][i]], self.mass[PID[f][j]]
                
                MFVnorm = vev/sqrt(2.*mi*mj)
                
                dy_cosphi = -MFVnorm*s_Re + diag
                dy_sinphi = MFVnorm*s_Im
                
                M[XB+'xdY'+f][i,j], M[XB+'xS'+f][i,j] = dy_sf(dy_cosphi, 
                                                              dy_sinphi)
        
        M['dCz'] = (W['cHbx'] - W['cHD']/4. + 3./4.*W['cll1221']
                   - 3./2.*W['WBxH3l'][1,1].real - 3./2.*W['WBxH3l'][2,2].real) 
        
        # Two derivative field strength interactions 
        M['Czbx'] = (-W['cll1221']/2. + W['cHD']/2.
                    + W['WBxH3l'][1,1].real + W['WBxH3l'][2,2].real)/gw2
        
        M['Cgg']  = (4./gs2)*W['cHG'] 
        
        M['Caa']  = 4.*(W['cHW']/gw2 + W['cHB']/gp2 - W['cHWB']/gw/gp)
        
        M['Czz']  = 4.*( gw2*W['cHW'] + gp2*W['cHB'] 
                       + gw*gp*W['cHWB'] )/(gw2+gp2)**2
                     
        M['Cza']  = (4.*W['cHW'] - 4.*W['cHB'] 
                    - 2.*(gw2-gp2)/(gw*gp)*W['cHWB'])/(gw2+gp2)

        M['tCgg']  = (4./gs2)*W['tcHG'] 
        
        M['tCaa']  = 4.*(W['tcHW']/gw2 + W['tcHB']/gp2 - W['tcHWB']/gw/gp)
        
        M['tCzz']  = 4.*( gw2*W['tcHW'] + gp2*W['tcHB'] 
                       + gw*gp*W['tcHWB'] )/(gw2+gp2)**2
                     
        M['tCza']  = (4.*W['tcHW'] - 4.*W['tcHB'] 
                    - 2.*(gw2-gp2)/(gw*gp)*W['tcHWB'])/(gw2+gp2)
        
        # Higgs cubic interaction 
        M['dL3']  =  (lam*( 3.*W['cHbx'] - 3./4.*W['cHD'] + 1./4.*W['cll1221']
                         - W['WBxH3l'][1,1].real/2. - W['WBxH3l'][2,2].real/2.)
                         - W['cH'])
        
        # Four fermion interactions
        trivial_4f =  ['cll1111','cll1122','cll1133','cll1331','cll2332',
                       'cle1111','cle1122','cle2211','cle1133','cle3311',
                       'cee1111','cee1122','cee1133',
                       'cpuu3333']
        
        for c in trivial_4f:
            M[c] = W[c]
        
        return M
    
