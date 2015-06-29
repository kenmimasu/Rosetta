from Basis import Basis, flavour_matrix
import MassBasis as MB
import WarsawBasis as WB
import HiggsBasis as HB
import math
from math import sqrt
from itertools import combinations_with_replacement as comb
from itertools import product
from __init__ import PID
################################################################################
# SILH basis class
class SILHBasis(Basis):
    SBNOTWARSAW = ['sW','sB','sHW','sHB','stHW','stHB','s2W','s2B','s2G']
    # [Tab. 1]
    # cWW, ctWW, sWB, stWB absent compared to Warsaw
    SBV2H2 = ['sGG','stGG','sBB', 'stBB']
    
    SBH4D2 = ['sH','sT']
    
    SBH6 = ['s6H']
    
    SBV3D3 = ['s3W','s3G','st3W','st3G']
    
    cu = flavour_matrix('su',kind='general',domain='complex')
    cd = flavour_matrix('sd',kind='general',domain='complex')
    ce = flavour_matrix('se',kind='general',domain='complex')
    SBF2H3 = cu + cd + ce
    
    cHl  = flavour_matrix('sHl' ,kind='hermitian',domain='complex')
    cpHl = flavour_matrix('spHl',kind='hermitian',domain='complex')
    cHe  = flavour_matrix('sHe' ,kind='hermitian',domain='complex')
    cHq  = flavour_matrix('sHq' ,kind='hermitian',domain='complex')
    cpHq = flavour_matrix('spHq',kind='hermitian',domain='complex')
    cHu  = flavour_matrix('sHu' ,kind='hermitian',domain='complex')
    cHd  = flavour_matrix('sHd' ,kind='hermitian',domain='complex')
    cHud = flavour_matrix('sHud',kind='general',domain='complex')
    
    # Two vertex operators absent also
    cHl.remove('sHl11')
    cpHl.remove('spHl11')
    
    SBF2H2D = cHl + cpHl + cHe + cHq + cpHq + cHu + cHd + cHud
        
    independent = ( SBH4D2 + SBH6 + SBV3D3 + SBF2H3 
                  + SBV2H2 + SBF2H2D + SBNOTWARSAW )
    
    dependent = ['sHl11', 'spHl11']
    
    required_masses = set([y for x in PID.values() for y in x.values()])
    required_inputs = {1, 2, 3, 4, 8} # aEWM1, Gf, aS, MZ, MH

    blocks = {'SBV2H2':SBV2H2, 'SBH4D2':SBH4D2, 'SBH6':SBH6,
              'SBV3D3':SBV3D3, 'SBF2H3':SBF2H3, 
              'SBF2H2D':SBF2H2D+['sHl11', 'spHl11'], 'SBNOTWARSAW':SBNOTWARSAW}

    translate_to={'mass','warsaw','higgs'}
    
    def calculate_dependent(self):
        # These coefficients are implictly set to zero
        self['sHl11'], self['spHl11']=0.,0.
        
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
        
    def translate(self):
        if self.target_basis=='mass': 
            self.newname='Mass'
            instance = MB.MassBasis()
            return self.translate_to_mass(instance)
        elif self.target_basis=='higgs':
            self.newname='Higgs'
            instance = HB.HiggsBasis() 
            return self.translate_to_mass(instance)
        elif self.target_basis=='warsaw':
            self.newname='Warsaw'
            instance = WB.WarsawBasis()
            return self.translate_to_warsaw(instance)
        else: 
            raise NotImplementedError
                
    def translate_to_mass(self,instance):
        s2w, c2w, ee2, gw2, gp2, MZ, vev, gs2 = self.calculate_inputs() 
        MH = self.mass[25]
        MW = MZ*sqrt(c2w)
                
        A = self
        B = instance
        
        # W mass shift [eqn (5.11)]
        B['dM'] = - gw2*gp2/(4.*(gw2-gp2))*(A['sW'] + A['sB'] + A['s2W'] 
                                            + A['s2B'] - 4./gp2*A['sT'] 
                                            + 2./gw2*A['spHl22'])
        def f(T3,Q): # [eqn (5.12)]
            Qcoeff = gp2/4./(gw2-gp2)*( -(2.*gw2-gp2)*A['s2B'] 
                       - gw2*(A['s2W'] + A['sW'] + A['sB'] ) 
                       + 4.*A['sT'] - 2.*A['spHl22'])
            T3coeff = ( gw2*A['s2W'] + gp2*A['s2B'] 
                        + 4.*A['sT']-2.*A['spHl22'])/4.
            return  delta(i,j)*(Q*Qcoeff + T3*T3coeff)
            
        def delta(i,j):
            return 1. if i==j else 0.
            
        # W/Z chiral coupling deviations
        for i,j in comb((1,2,3),2):
            ind = '{}{}'.format(i,j)
            tail = [ind] if i==j else [ind+'_Re', ind+'_Im']
            for t in tail:
                # [eqn (5.13)]
                B['dGLzv%s' % t] = (1./2.*A['spHl%s' % t]  
                                  - 1./2.*A['sHl%s' % t] + f(1./2.,0.))
                B['dGLze%s' % t] = (-1./2.*A['spHl%s' % t] 
                                   - 1./2.*A['sHl%s' % t] + f(-1./2.,-1.))
                B['dGRze%s' % t] = - 1./2.*A['sHe%s' % t] + f(0.,-1.)
                B['dGLzu%s' % t] = (1./2.*A['spHq%s' % t] 
                                  - 1./2.*A['sHq%s' % t] + f(1./2.,2./3.))
                B['dGLzd%s' % t] = (-1./2.*A['spHq%s' % t] 
                                  - 1./2.*A['sHq%s' % t] + f(-1./2.,-1./3.))
                B['dGRzu%s' % t] = - 1./2.*A['sHu%s' % t] + f(0.,2./3.)
                B['dGRzd%s' % t] = - 1./2.*A['sHd%s' % t] + f(0.,-1./3.)
                B['dGLwl%s' % t] = (A['spHl%s' % t] + f(1./2.,0.) 
                                                    - f(-1./2.,-1.))
                B['dGLwq%s' % t] = (A['spHq%s' % t] + f(1./2.,2./3.)
                                                    - f(-1./2.,-1./3.))
                # [eqn (5.14)]
                B['CLwl%s' % t] = B['dGLwl%s' % t]
                B['CLzv%s' % t] = B['dGLzv%s' % t]
                B['CLze%s' % t] = B['dGLze%s' % t]
                B['CRze%s' % t] = B['dGRze%s' % t]
                B['CLwq%s' % t] = B['dGLwq%s' % t]
                B['CLzu%s' % t] = B['dGLzu%s' % t]
                B['CLzd%s' % t] = B['dGLzd%s' % t]
                B['CRzu%s' % t] = B['dGRzu%s' % t]
                B['CRzd%s' % t] = B['dGRzd%s' % t]
        
        # Treat dGRwq separately as it has more flavour components
        dGRwq = flavour_matrix('dGRwq', kind='general', domain='complex')
        for coeffW, coeffM in zip(self.cHud, dGRwq):
            B[coeffM] = -1./2.*A[coeffW]
            cvff = coeffM.replace('dG','C')
            B[cvff] = B[coeffM]
        
        # Higgs couplings to W/Z [eqn (4.14)]
        B['dCw'] =  -A['sH'] - gw2*gp2/(gw2-gp2)*( A['sW'] + A['sB'] 
                         + A['s2W'] + A['s2B'] - 4./gp2*A['sT'] 
                         + (3.*gw2 + gp2)/(2.*gw2*gp2)*A['spHl22'])
        B['dCz'] = -A['sH'] - 3./2.*A['spHl22'] 
        
        # Two derivative field strength interactions [eqn (4.15)]
        B['Cgg']  = A['sGG']
        B['Caa']  = A['sBB']
        B['Czz']  = -1./(gw2+gp2)*( gw2*A['sHW'] + gp2*A['sHB'] 
                                     - gp2*s2w*A['sBB'] )
        B['Czbx'] =  1./(2.*gw2)*( gw2*(A['sW'] + A['sHW'] +A ['s2W'])
                                  + gp2*(A['sB'] + A['sHB'] +A ['s2B'])
                                  - 4.*A['sT'] + 2.*A['spHl22'] ) 
        B['Cza']  = ( A['sHB'] - A['sHW'])/2. - s2w*A['sBB']
        
        B['Cabx'] =  ( A['sHW'] - A['sHB'])/2.+ 1./(gw2-gp2)*(
                           gw2*(A['sW']+A['s2W']) + gp2*(A['sB']+A['s2B'])
                         - 4.*A['sT'] + 2.*A['spHl22'] )
        B['Cww']  =  -A['sHW']
        # factor 2 wrong in note here 
        B['Cwbx'] =  A['sHW']/2. + 1./2./(gw2-gp2)*(
                           gw2*(A['sW']+A['s2W']) + gp2*(A['sB']+A['s2B'])
                         - 4.*A['sT'] + 2.*A['spHl22'] ) 
        B['CTgg'] = A['stGG'] 
        B['CTaa'] = A['stBB']
        B['CTzz'] = -1./(gw2+gp2)*( gw2*A['stHW'] + gp2*A['stHB'] 
                                     - gp2*s2w*A['stBB'] )
        B['CTza'] = ( A['stHB'] - A['stHW'])/2. - s2w*A['stBB']
        B['CTww'] =  -A['stHW']
        
        
        # solution for  [eqn (5.16)]
        # dy*cos(phi) == X
        # dy*sin(phi) == Y
        def dy_sf(X,Y): 
            R = math.sqrt(X**2+Y**2)
            sf = -Y/R # solution is +-Y/R 
            dy = (-Y**2 + X*abs(X))/R
            return dy,sf
        
        # Yukawa type interaction coefficients [eqn. (4.16)]
        for i,j in comb((1,2,3),2): 
            for f in ('u','d','e'):
                mi, mj = self.mass[ PID[f][i] ],self.mass[ PID[f][j] ] 
                name = '{}{}{}'.format(f,i,j)
                if mi and mj:
                    dy_cosphi = (vev*A['s'+name+'_Re']/sqrt(2.*mi*mj) -
                                 delta(i,j)*(A['sH'] + A['spHl22']/2.))
                    dy_sinphi = vev*A['s'+name+'_Im']/sqrt(2.*mi*mj)
                    if (dy_cosphi == 0.) and (dy_sinphi == 0.): # check this consistency
                        B['dY'+name], B['S'+name] = 0., 0. 
                    else:
                        B['dY'+name], B['S'+name]  = dy_sf(dy_cosphi, dy_sinphi)
        
        # Double Higgs Yukawa type interaction coefficients [eqn. (4.17)]
        for i,j in comb((1,2,3),2):
            for f in ('u','d','e'):
                name = '{}{}{}'.format(f,i,j)
                Yij   = B['dY' + name]
                sinij = B['S' + name] 
                cosij = sqrt(1. - sinij**2)
                B['Y2{}_Re'.format(name)] = (3.*Yij*cosij + delta(i,j)*(
                                                  A['sH'] + 3./2.*A['spHl22']))
                B['Y2{}_Im'.format(name)] =  3.*Yij*sinij
        
        # Triple gauge couplings [eqn. (4.18)]
        B['dG1z'] = -(gw2+gp2)/(gw2-gp2)/4.*( (gw2-gp2)*A['sHW'] 
                        + gw2*(A['sW'] + A['s2W']) + gp2*(A['sB'] + A['s2B']) 
                        - 4.*A['sT'] + 2.*A['spHl22'] )
        B['dKa'] = - gw2/4.*(A['sHW']+A['sHB'])
        B['dKz'] = ( -1./4.*(gw2*A['sHW'] - gp2*A['sHB'])
                      -(gw2+gp2)/(gw2-gp2)/4.*(
                          gw2*(A['sW'] + A['s2W'])
                          + gp2*(A['sB'] + A['s2B']) 
                          - 4.*A['sT'] + 2.*A['spHl22'] 
                          )
                      )
        B['La'] = -A['s3W']*3./2.*gw2**2
        B['Lz'] =  B['La']
        B['KTa'] =  - gw2/4.*(A['stHW']+A['stHB'])
        B['KTz'] = gp2/4.*(A['stHW']+A['stHB'])
        B['LTa'] = -A['st3W']*3./2.*gw2**2
        B['LTz'] =  B['LTa']

        # Quartic gauge couplings [Sec 3.7] [eqn (3.23)] 
        B['dGw4'] = 2.*c2w*B['dG1z']
        B['dGw2z2'] = 2.*B['dG1z']
        B['dGw2za'] = B['dG1z']
        
        # two derivative quartic gauge couplings [Sec 3.7] [eqn (3.24)] 
        B['Ldw4'] = -gw2/2./MW**2*B['Lz']
        B['LTdw4'] = -gw2/2./MW**2*B['LTz']
        B['Ldzdw_zw'] = -gw2*c2w/MW**2*B['Lz']
        B['LTdzdw_zw'] = -gw2*c2w/MW**2*B['LTz']
        B['Ldzdw_aw'] = -ee2/MW**2*B['Lz']
        B['LTdzdw_aw'] = -ee2/MW**2*B['LTz']
        B['Ldadw_aw'] = -sqrt(ee2*gw2*c2w)/MW**2*B['Lz']
        B['LTdadw_aw'] = -sqrt(ee2*gw2*c2w)/MW**2*B['LTz']
        B['Ldadw_zw'] = -sqrt(ee2*gw2*c2w)/MW**2*B['Lz']
        B['LTdadw_zw'] = -sqrt(ee2*gw2*c2w)/MW**2*B['LTz']
        B['Ldg_g3'] = 3.*sqrt(gs2)**3/vev**2*B['C3G']
        B['LTdg_g3'] = 3.*sqrt(gs2)**3/vev**2*B['CT3G']

        # Higgs cubic interaction [eqn. (4.19)]
        B['dL3']  =  -MH**2/(2.*vev**2)*(3.*A['sH'] + A['spHl22']/2.) - A['s6H']
        
        # Couplings of two Higgs bosons to gluons [Sec 3.8]
        # [eqn (3.27)] copied from HiggsBasis implemetation
        B['Cgg2'], B['CTgg2'] = B['Cgg'], B['CTgg']
        
        B['cll1122']=(gw2*A['s2W']+gp2*A['s2B'])/2.
        B['cll1221']=gw2*A['s2W']
        B['cpuu3333'] = (1./3.)*gs2*A['s2G']
        
        self.mass[24] = MW + B['dM']

        return B
        
    def translate_to_warsaw(self, instance):
        s2w, c2w, ee2, gw2, gp2, MZ, vev, gs2 = self.calculate_inputs() 
        MH = self.mass[25]
        lam = -MH**2/(2.*vev**2) # Higgs self-coupling     
           
        A = self
        B = instance
        
        def delta(i,j):
            return 1. if i==j else 0.
        
        # [eqn (5.7)]
        B['cH'] = A['sH'] - 3./4.*gw2*(A['sW'] + A['sHW'] + A['s2W'])
        B['cT'] = A['sT'] - 1./4.*gp2*(A['sB'] + A['sHB'] + A['s2B'])
        B['c6H'] = A['s6H'] - 2.*lam*gw2*(A['sW'] + A['sHW'] + A['s2W'])
        B['cWB'] = - 1./4.*(A['sHB'] + A['sHW'])
        B['cBB'] = A['sBB'] - A['sHB']
        B['cWW'] = -A['sHW']
        B['ctWB'] = - 1./4.*(A['stHB'] + A['stHW'])
        B['ctBB'] = A['stBB'] - A['stHB']
        B['ctWW'] = -A['stHW']
        
        def cHf(sc, Yf, i, j):
            if '_Im' not in sc: 
                return A[sc] + delta(i,j)*gp2*Yf/2.*(A['sB'] + 
                                                     A['sHB'] + 2.*A['s2B'])
            else:
                return A[sc]
            
        def cpHf(sc, i, j):
            if '_Im' not in sc:
                return A[sc] + delta(i,j)*gw2/4.*(A['sW'] +
                                                  A['sHW'] + 2.*A['s2W'])
            else:
                return A[sc]
            
        # [eqn (5.8)]
        for i,j in comb((1,2,3),2):
            ind = '{}{}'.format(i,j)
            tail = [ind] if i==j else [ind+'_Re', ind+'_Im']
            for t in tail:
                B['cHl'+t] = cHf('sHl'+t, -1./2., i, j)
                B['cHe'+t] = cHf('sHe'+t, -1., i, j)
                B['cHq'+t] = cHf('sHq'+t, 1./6., i, j)
                B['cHu'+t] = cHf('sHu'+t, 2./3., i, j)
                B['cHd'+t] = cHf('sHd'+t, -1./3., i, j)
                B['cpHl'+t] = cpHf('spHl'+t, i, j)
                B['cpHq'+t] = cpHf('spHq'+t, i, j)

        # [eqn (5.9)]
        # self.mass[ PID[f][i] ],self.mass[ PID[f][j] ]
        for i,j in comb((1,2,3),2):
            ind = ['{}{}_{}'.format(i,j,cplx) for cplx in ('Re','Im')]
            for t in ind:
                for f in ('u','d','e'):
                    mass = self.mass[ PID[f][i] ]
                    yuk = sqrt(2.)*mass/vev
                    wcoeff, scoeff = 'c{}{}'.format(f,t),'s{}{}'.format(f,t)
                    if '_Im' not in scoeff:
                        B[wcoeff] = A[scoeff] - delta(i,j)*gw2*yuk*(
                                              A['sW'] + A['sHW'] + A['s2W'])/2.
                    else:
                        B[wcoeff] = A[scoeff] 

        # [eqn (5.10)]
        B['cll1221'] = gw2*A['s2W'] # cll1221==0 in SILH
        # derived from [eqn (5.4)]
        B['cll1122']=(gw2*A['s2W']+gp2*A['s2B'])/2.
        B['cpuu3333'] = (1./3.)*gs2*A['s2G']
        
        # trivial translation, cX==sX
        cHud = flavour_matrix('cHud', kind='general', domain='complex')
        others = ['cGG','ctGG','c3W','ct3W','c3G','ct3G']
        trivial = cHud+others
        for coeff in trivial:
            scoeff = 's'+coeff[1:]
            B[coeff] = A[scoeff]
        
        return B
        
################################################################################