from Basis import Basis, flavour_matrix
from MassBasis import MassBasis
import math
from math import sqrt
from itertools import combinations_with_replacement as comb
from itertools import product
from __init__ import PID
################################################################################
# SILH basis class
class SILHBasis(Basis):
    not_warsaw = ['sW','sB','sHW','sHB','stHW','stHB','s2W','s2B','s2G']
    # [Tab. 1]
    # cWW, ctWW, sWB, stWB absent compared to Warsaw
    V2H2 = ['sGG','stGG','sBB', 'stBB']
    
    H4D2 = ['sH','sT']
    
    H6 = ['s6H']
    
    V3D3 = ['s3W','s3G','st3W','st3G']
    
    cu = flavour_matrix('su',kind='symmetric',domain='complex')
    cd = flavour_matrix('sd',kind='symmetric',domain='complex')
    ce = flavour_matrix('se',kind='symmetric',domain='complex')
    f2H3 = cu + cd + ce
    
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
    
    f2H2D = cHl + cpHl + cHe + cHq + cpHq + cHu + cHd + cHud
    
    fourfermi = []
    
    independent = (H4D2 + H6 + V3D3 + f2H3 + V2H2 + f2H2D 
                +  not_warsaw + fourfermi)
    
    required_masses = set([y for x in PID.values() for y in x.values()])
    required_inputs = {1, 2, 4, 8} # aEWM1, Gf, MZ, MH
    
    translate_to={'mass','warsaw'}
        
    def calculate_inputs(self): # calculate a few required EW params from aEWM1, Gf, MZ
        ee2 = 4.*math.pi/self.input['aEWM1'] # EM coupling squared
        Gf, MZ = self.input['Gf'], self.input['MZ']
        s2w = (1.- sqrt(1. - ee2/(sqrt(2.)*Gf*MZ**2)))/2.#sin^2(theta_W)
        c2w = 1. - s2w
        gw2 = ee2/s2w # SU(2) coupling squared
        gp2 = gw2*s2w/c2w # Hypercharge coupling squared
        vev =  2.*MZ*sqrt(c2w/gw2)
        self.input['s2w'],self.input['ee2'] = s2w, ee2
        self.input['gw2'],self.input['gp2'] = gw2, gp2
        self.input['vev'] = vev
        return s2w, ee2, gw2, gp2, vev
        
    def translate(self):
        if self.target_basis=='mass': 
            self.translate_to_mass()
        else: 
            raise NotImplementedError
    # def translate_to_warsaw(self):
    #     self.newname='Warsaw'
    #     s2w, ee2, gw2, gp2, vev = self.calculate_inputs()
            
    def translate_to_mass(self):
        self.newname='Mass'
        s2w, ee2, gw2, gp2, vev = self.calculate_inputs()
        MH = self.input['MH']

        A = self.coeffs._asdict()
        B = MassBasis().coeffs._asdict()
        # These coefficients are implictly set to zero
        A['sHl11'], A['spHl11']=0.,0.
        
        
        # W mass shift [eqn (5.11)]
        B['dM'] = - gw2*gp2/(4.*(gw2-gp2))*(A['sW'] + A['sB'] + A['s2W'] + A['s2B']
                                      - 4./gp2*A['sT'] + 2./gw2*A['spHl22'])
        def f(T3,Q): # [eqn (5.12)]
            Qcoeff = gp2/4./(gw2-gp2)*( (2.*gw2-gp2)*A['s2B'] 
                       - gw2*(A['s2W'] + A['sW'] + A['sB'] ) 
                       + 4.*A['sT'] - 2.*A['spHl22'])
            T3coeff = ( gw2*A['s2W'] + gp2*A['s2B'] 
                        + 4.*A['sT']-2.*A['spHl22'])/4.
            return  Q*Qcoeff + T3*T3coeff
            
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
        B['Czbx'] =  -1./(2.*gw2)*( gw2*(A['sW'] + A['sHW'] +A ['s2W'])
                                  + gp2*(A['sB'] + A['sHB'] +A ['s2B'])
                                  - 4.*A['sT'] + 2.*A['spHl22'] ) 
        B['Cza']  = ( A['sHB'] - A['sHW'])/2. - s2w*A['sBB']
        
        B['Cabx'] =  ( A['sHW'] - A['sHB'])/2.+ 1./(gw2-gp2)*(
                           gw2*(A['sW']+A['s2W']) + gp2*(A['sB']+A['s2B'])
                         - 4.*A['sT'] + 2.*A['spHl22'] )
        B['Cww']  =  -A['sHW']
        B['Cwbx'] =  A['sHW']/2. + 1./(gw2-gp2)*(
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
                                 delta(i,j)*(A['sH']
                                     + 3./4.*gw2*(A['sW'] + A['sHW'] + A['s2W'])
                                     + A['spHl22']/2.))
                    dy_sinphi = vev*A['s'+name+'_Im']/sqrt(2.*mi*mj) 
                    B['dY'+name], B['S'+name]  = dy_sf(dy_cosphi, dy_sinphi)
        
        # Double Higgs Yukawa type interaction coefficients [eqn. (4.17)]
        for i,j in comb((1,2,3),2):
            for f in ('u','d','e'):
                name = '{}{}{}'.format(f,i,j)
                Yij   = B['dY' + name]
                sinij = B['S' + name] 
                cosij = sqrt(1. - sinij**2)
                B['Y2{}_Re'.format(name)] = (3.*Yij*cosij 
                                     + delta(i,j)*(A['sH']
                                     + 3./4.*gw2*(A['sW'] + A['sHW'] + A['s2W']) 
                                     + 3./2.* A['spHl22']))
                B['Y2{}_Im'.format(name)] = 3.*Yij*sinij
        
        # Triple gauge couplings [eqn. (4.18)]
        B['dG1z'] = -(gw2+gp2)/(gw2-gp2)/4.*( (gw2-gp2)*A['sHW'] 
                        + gw2*(A['sW'] + A['s2W']) + gp2*(A['sB'] + A['s2B']) 
                        - 4.*A['sT'] + 2.*A['spHl22'] )
        B['dKa']  = - gw2/4.*(A['sHW']+A['sHB'])
        B['dKz']  = (-1./4.*(gw2*A['sHW'] + gp2*A['sHB'])
                    -(gw2+gp2)/(gw2-gp2)/4.*(gw2*(A['sW'] + A['s2W'])
                         + gp2*(A['sB'] + A['s2B']) 
                         - 4.*A['sT'] + 2.*A['spHl22'] ))
        B['La']   = -A['s3W']*3./2.*gw2**2
        B['Lz']   =  B['La']
        B['KTa']  =  - gw2/4.*(A['stHW']+A['stHB'])
        B['KTz']  = gp2/4.*(A['stHW']+A['stHB'])
        B['LTa']  = -A['st3W']*3./2.*gw2**2
        B['LTz']  =  B['LTa']
        
        # Higgs cubic interaction [eqn. (4.19)]
        B['dL3']  =  -MH**2/(2.*vev**2) * (3.*A['sH'] + A['spHl22']/2.) - A['s6H']
        
        # Couplings of two Higgs bosons to gluons [Sec 3.8]
        # [eqn (3.27)] copied from HiggsBasis implemetation
        B['Cgg2'], B['CTgg2'] = B['Cgg'], B['CTgg']
        
        B['cll1221']=gw2*A['s2W']
        self.newpar = B
        
        self.newmass[24] = self.mass[24]+self.newpar['dM'] # W mass shift
        
################################################################################