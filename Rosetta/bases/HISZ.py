import math
from math import sqrt

from ..internal import basis
from ..internal import PID
from ..internal import matrix_mult, matrix_add, matrix_sub, matrix_eq
################################################################################
class HISZ(basis.Basis):
    '''
    Basis class for Rosetta based on [ Hagiwara et al., 
    Phys.Rev. D48 (1993) 2182-2203 ]. The exact list of operators included as 
    well as the equations for the translation to and from the Higgs basis can 
    be found in the HXSWG note, for which all references to equations are in 
    this implementation.
    '''

    name = 'hisz'
    ##########################
    # declare coefficients
    # [Tab. 1]
    LAMBDA = ['Lam']
    EVEN = ['fGG', 'fBB', 'fWW', 'fB', 'fW', 'fWWW', 'fH2']
    ODD = ['tfGG', 'tfBB', 'tfWW', 'tfW', 'tfWWW']

    ##########################
    # block structure
    blocks = {'LAMBDA':LAMBDA, 'EVEN':EVEN, 'ODD':ODD} 
    
    flavored={
        'HZxu': {'cname':'fu', 'kind':'general', 'domain':'complex'},
        'HZxd': {'cname':'fd', 'kind':'general', 'domain':'complex'},
        'HZxe': {'cname':'fe', 'kind':'general', 'domain':'complex'}
    }

    independent = blocks.keys() + flavored.keys()

    required_masses = {25, 24}  # Higgs & W masses
    required_inputs = {1, 2, 3, 9} # aEWM1, Gf, aS, MW
    ##########################
        
    def calculate_inputs(self): # calculate a few required EW params from aEWM1, Gf, MW
        ee2 = 4.*math.pi/self.inputs['aEWM1'] # EM coupling squared
        gs2 = 4.*math.pi*self.inputs['aS'] # strong coupling squared
        Gf, MW = self.inputs['Gf'], self.inputs['MW']
        # s2w = (1.- sqrt(1. - ee2/(sqrt(2.)*Gf*MZ**2)))/2. # sin^2(theta_W)
        s2w = ee2/(4.*sqrt(2.)*Gf*MW**2) # sin^2(theta_W)
        c2w = (1.-s2w)
        gw2 = ee2/s2w # SU(2) coupling squared
        gp2 = gw2*s2w/c2w # Hypercharge coupling squared
        vev =  2.*MW/sqrt(gw2)
        return s2w, c2w, ee2, gw2, gp2, MW, vev, gs2
    
    @basis.translation('higgs')
    def to_higgs(self, instance):
        '''
        Translation function to Higgs basis.
        '''
        # Higgs basis prefix
        XB = 'HB'
        
        s2w, c2w, ee2, gw2, gp2, MW, vev, gs2 = self.calculate_inputs() 

        H = self
        M = instance
        
        vLsq = vev**2/H['Lam']**2
        
        # Two derivative field strength interactions
        # [eqn (A.8)]
        # CP even
        M['Cgg']  = -H['fGG']*(1./(8.*math.pi**2))*vLsq
        
        M['dCz'] = -H['fH2']/2.*vLsq 
        
        M['Caa']  = -(H['fWW'] + H['fBB'])*vLsq 
        
        M['Cza']  = (H['fW']/4. - H['fB']/4. 
                    - c2w*H['fWW'] + s2w*H['fBB'])*vLsq
                    
        M['Czz']  = (c2w*H['fW']/2. + s2w*H['fB']/2. 
                    - c2w**2*H['fWW'] - s2w**2*H['fBB'])*vLsq
        
        M['Czbx'] =  -(H['fW'] + s2w/c2w*H['fB'])/4.*vLsq 

        # [eqn (A.9)]
        # CP odd                       
        M['tCgg'] = -H['tfGG']*(1./(8.*math.pi**2))*vLsq
                
        M['tCaa'] = -(H['tfWW'] + H['tfBB'])*vLsq 
        
        # M['tCza'] = (H['tfW']/4. - c2w*H['tfWW'] + s2w*H['tfBB'])*vLsq
        M['tCza'] = (H['tfW']/4. - c2w*H['tfWW'] + s2w*H['tfBB'])*vLsq
        
        # M['tCzz'] = (c2w*H['tfW']/2. - c2w**2*H['tfWW'] - s2w**2*H['tfBB'])*vLsq
        M['tCzz'] = (c2w*H['tfW']/2. - c2w**2*H['tfWW'] - s2w**2*H['tfBB'])*vLsq
        
        # Triple gauge couplings [eqn. (A.10)]
        M['Lz'] = 3.*gw2**2/8.*H['fWWW']
                
        M['tLz'] = 3.*gw2**2/8.*H['tfWWW']
        
        # solution for Higgs-fermion couplings
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
        
        def delta(i,j):
            return 1. if i==j else 0.
            
        # Yukawa type interaction coefficients [eqns. (A.8) & (A.9)]
        for f in ('u','d','e'):
            matrix = 'HZx' + f
            for i,j in H[matrix].keys(): 
                diag = delta(i,j)*H['fH2']/2.
                f_Re, f_Im = H[matrix][i,j].real, H[matrix][i,j].imag
                dy_cosphi = -(f_Re/sqrt(2.) + delta(i,j)*H['fH2']/2.)*vLsq
                dy_sinphi = f_Im/sqrt(2.)*vLsq
                
                M[XB+'xdY'+f][i,j], M[XB+'xS'+f][i,j] = dy_sf(dy_cosphi, 
                                                              dy_sinphi)
        
        # Provide the right Z mass for input...
        MZ = MW/sqrt(c2w)
        H.mass[23]= MZ
        try:
            H.inputs.new_entry(4, MZ, name='MZ')
        except KeyError:
            H.inputs[4] = MZ
        return M

        
################################################################################