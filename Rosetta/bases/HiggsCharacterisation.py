import math
from math import sqrt, pi

from . import basis
from ..internal import PID
from ..internal import matrix_mult, matrix_add, matrix_sub, matrix_eq
################################################################################
class HiggsCharacterisation(basis.Basis):
    '''
    Implementation of anomalous Higgs couplings in FeynRules upon which the 
    BSMCharacterisation implementation is based. 
    See https://feynrules.irmp.ucl.ac.be/wiki/HiggsCharacterisation.
    '''
    
    name ='hc'
    # kinetic terms [Eqn. (3.3)]
    FRBlock = [
    'Lambda','cosa','kSM',
    'kHtt','kAtt',
    'kHbb','kAbb',
    'kHll','kAll',
    'kHaa','kAaa',
    'kHza','kAza',
    'kHgg','kAgg',
    'kHzz','kAzz',
    'kHww','kAww',
    'kHda','kHdz',
    'kHdwR','kHdwI',
    'kHcc','kAcc'
    ]
    
    blocks = {'FRBlock':FRBlock}

    numbers = {'kHcc':100, 'kAcc':101}
    
    # All parameters independent
    # dependent = ['cosa', 'Lambda']
    
    # independent = [c for c in FRBlock if c not in dependent ]
    independent = FRBlock

    required_masses = {25, 24, 1, 2, 3, 4, 5, 6, 11, 13, 15} 

    required_inputs = {1, 2, 3, 4}
    # all other undefined behaviour inherited from Basis.Basis by default
################################################################################
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
        
    @basis.translation('bsmc')        
    def to_bsmc(self, instance):
        
        s2w, c2w, ee2, gw2, gp2, MZ, vev, gs2 = self.calculate_inputs() 
        MWsq = MZ**2*c2w
        aS, aEM, Gf = self.inputs['aS'], self.inputs['aEWM1'], self.inputs['Gf']
        H = self
        B = instance
        
        cosa, Lam = H['cosa'], H['Lambda']
        sina = sqrt(1.-cosa**2)
        
        # constants
        C = sqrt( aEM*Gf*MZ**2/(8.*sqrt(2.)*pi) )
        gHZZ, gHWW = 2.*MZ**2/vev, 2.*MWsq/vev
        gHaa, gAaa = 47.*aEM/(18.*pi)/vev, 4.*aEM/(3.*pi)/vev
        gHza, gAza = C*(94*c2w-13.)/(9.*pi)/vev, 2.*C*(8.*c2w-5.)/(3.*pi)/vev
        gHgg, gAgg = -aS/(3.*pi)/vev, aS/(2.*pi)/vev

        # Gauge        
        B['dCz'] = 2.*cosa*gHZZ*H['kSM']/(vev*(gw2+gp2))-1.

        B['Caa'] = -cosa*gHaa*H['kHaa']/ee2

        B['tCaa'] = -sina*gAaa*H['kAaa']/ee2
        
        B['Cza'] = -cosa*gHza*H['kHza']*vev/sqrt(ee2*(gw2+gp2))

        B['tCza'] = -sina*gAza*H['kAza']*vev/sqrt(ee2*(gw2+gp2))

        B['Cgg'] = -2.*cosa*gHgg*H['kHgg']/gs2

        B['tCgg'] = -2.*sina*gAgg*H['kAgg']/gs2

        B['Czz'] = -cosa*H['kHzz']/(gw2+gp2)*(vev/Lam)

        B['tCzz'] = -sina*H['kAzz']/(gw2+gp2)*(vev/Lam)
        
        B['Cww'] = -cosa*H['kHww']/gw2*(vev/Lam)

        B['tCww'] = -sina*H['kAww']/gw2*(vev/Lam)

        B['Cabx'] = -cosa*H['kHda']/sqrt(gw2*gp2)*(vev/Lam)

        B['Czbx'] = -cosa*H['kHdz']/gw2*(vev/Lam)

        B['Cwbx'] = -cosa*H['kHdwR']/gw2*(vev/Lam)

        # Yukawa
        B['BCxdYu'][3,3] = sqrt( (cosa*H['kHtt']-1.)**2 + (H['kAtt']*sina)**2 )
        
        B['BCxSu'][3,3] = H['kAtt']*sina/B['BCxdYu'][3,3]
        
        B['BCxdYu'][2,2] = sqrt( (cosa*H['kHcc']-1.)**2 + (H['kAcc']*sina)**2 )
        
        B['BCxSu'][2,2] = H['kAcc']*sina/B['BCxdYu'][2,2]

        B['BCxdYd'][3,3] = sqrt( (cosa*H['kHbb']-1.)**2 + (H['kAbb']*sina)**2 )
        
        B['BCxSd'][3,3] = H['kAbb']*sina/B['BCxdYd'][3,3]
        
        B['BCxdYe'][3,3] = sqrt( (cosa*H['kHll']-1.)**2 + (H['kAll']*sina)**2 )
        
        B['BCxSe'][3,3] = H['kAll']*sina/B['BCxdYe'][3,3]

        return B

    # def modify_inputs(self):
    #     '''
    #     W mass modification from dM.
    #     '''
    #     ee2 = 4.*math.pi/self.inputs['aEWM1'] # EM coupling squared
    #     Gf, MZ = self.inputs['Gf'], self.inputs['MZ']
    #     s2w = (1.- sqrt(1. - ee2/(sqrt(2.)*Gf*MZ**2)))/2. # sin^2(theta_W)
    #     c2w = (1.-s2w)
    #     MW = MZ*sqrt(c2w)
    #
    #     if 24 in self.mass:
    #         self.mass[24] = MW + self['dM']
    #     else:
    #         self.mass.new_entry(24, MW + self['dM'], name = 'MW')
################################################################################
