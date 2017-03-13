import math
from math import sqrt, pi

from . import basis
from ..internal import PID
from ..internal import matrix_mult, matrix_add, matrix_sub, matrix_eq
################################################################################
class HiggsPO(basis.Basis):
    '''
    Implementation of anomalous Higgs pseudo-observables in FeynRules. 
    See http://www.physik.uzh.ch/data/HiggsPO .
    '''
    
    name ='higgspo'

    HPO2F = ['kb','kc','ktau','kmu','lb','lc','ltau','lmu']
    HPO4F = ['kZZ','kWW','kAA','kZA',
             'eZZ','eWW','lAACP','lZACP','eZZCP','eWWCP',
             'eZeL','eZmuL','eZtauL','eZeR','eZmuR','eZtauR','eZv',
             'eWe','eWmu','eWtau','phiWe','phiWmu','phiWtau']
    HPOSM = ['eAASM','eZASM','ybeff','yceff','ytaueff','ymueff']
    WZPOLE = ['gZeL','gZmuL','gZtauL','gZeR','gZmuR','gZtauR','gZv',
              'gWe','gWmu','gWtau']
    
    blocks = {'HPO2F':HPO2F, 'HPO4F':HPO4F, 'HPOSM':HPOSM, 'WZPOLE':WZPOLE}

    numbers = {'lb':11, 'lc':12,'ltau':13,'lmu':14}
    
    # All parameters independent
    independent = [x for v in blocks.values() for x in v ]

    # MH, MZ, MW, Ms, Mc, Mb, Mt, Mtau
    required_masses = {25, 23, 24, 3, 4, 5, 6, 15} 

    # aEWM1, vF, aS
    required_inputs = {1, 10, 3} # 
    
    inputs_blockname = 'smparam'
    # all other undefined behaviour inherited from Basis.Basis by default
################################################################################
    # calculate a few required EW params from aEWM1, vF, MZ, aS
    def calculate_inputs(self): 
        ee2 = 4.*math.pi/self.inputs['aEWM1'] # EM coupling squared
        gs2 = 4.*math.pi*self.inputs['aS'] # strong coupling squared
        vev, MZ = self.inputs['vF'], self.mass[23]
        s2w = (1.- sqrt(1. - ee2*vev**2/MZ**2))/2. # sin^2(theta_W)
        c2w = (1.-s2w)
        gw2 = ee2/s2w # SU(2) coupling squared
        gp2 = gw2*s2w/c2w # Hypercharge coupling squared
        return s2w, c2w, ee2, gw2, gp2, MZ, vev, gs2
        
    @basis.derived_input
    def Gf(self):
        return 1./sqrt(2.)/self.inputs['vF']**2
    
    @basis.derived_input
    def vF(instance):
        '''
        Derive the value of vF input from the Gf input parameter 
        of a basis translating to HiggsPO.
        '''
        return sqrt( 1./sqrt(2.)/instance.inputs['Gf'] )
        
    @basis.translation('bsmc')
    def to_bsmc(self, instance):
        s2w, c2w, ee2, gw2, gp2, MZ, vev, gs2 = self.calculate_inputs()
        # MWsq = MZ**2*c2w
        # aS, aEM, Gf = self.inputs['aS'], self.inputs['aEWM1'], self.inputs['Gf']
        H = self
        B = instance
        
        # Gauge-Higgs
        B['dCz'] = H['kZZ']-1.
        
        B['dCw'] = H['kWW']-1.
        
        B['Czz'] = -2.*H['eZZ']/(gw2+gp2)
        
        B['tCzz'] = -2.*H['eZZCP']/(gw2+gp2)
        
        B['Cza'] = -2.*H['kZA']/sqrt(ee2*(gw2+gp2))
        
        B['tCza'] = -2.*H['lZACP']/sqrt(ee2*(gw2+gp2))
        
        B['Caa'] = -2.*H['kAA']/ee2
        
        B['tCaa'] = -2.*H['lAACP']/ee2
        
        B['Cww'] = -2.*H['eWW']/gw2
        
        B['tCww'] = -2.*H['eWWCP']/gw2

        B['Cabx'] = 0.

        B['Czbx'] = 0.

        B['Cwbx'] = 0.
        
        # Higgs-Z-f-f contact   
        B['BCxdGLhze'][1,1]=H['eZeL']/sqrt(gw2+gp2)
        B['BCxdGLhze'][2,2]=H['eZmuL']/sqrt(gw2+gp2)
        B['BCxdGLhze'][3,3]=H['eZtauL']/sqrt(gw2+gp2)
        
        B['BCxdGRhze'][1,1]=H['eZeR']/sqrt(gw2+gp2)
        B['BCxdGRhze'][2,2]=H['eZmuR']/sqrt(gw2+gp2)
        B['BCxdGRhze'][3,3]=H['eZtauR']/sqrt(gw2+gp2)
        
        for i in (1,2,3):
            B['BCxdGLhzv'][i,i]=H['eZv']/sqrt(gw2+gp2)

        # Higgs-W-f-f contact   
        B['BCxdGLhwl'][1,1]=H['eWe']*sqrt(2./gw2)
        B['BCxdGLhwl'][2,2]=H['eWmu']*sqrt(2./gw2)
        B['BCxdGLhwl'][3,3]=H['eWtau']*sqrt(2./gw2)
        # complex phase parameters do not enter diagonals, dGLhwl is Hermitian.
        
        # Z-f-f coupling deviation
        B['BCxdGLze'][1,1]=H['gZeL'] - (s2w - 1./2.)
        B['BCxdGLze'][2,2]=H['gZmuL'] - (s2w - 1./2.)
        B['BCxdGLze'][3,3]=H['gZtauL'] - (s2w - 1./2.)
        
        B['BCxdGRze'][1,1]=H['gZeR'] - s2w  
        B['BCxdGRze'][2,2]=H['gZmuR'] - s2w 
        B['BCxdGRze'][3,3]=H['gZtauR'] - s2w  
        
        for i in (1,2,3):
            B['BCxdGLzv'][i,i]=H['gZv'] - 1./2.
        
        # W-f-f coupling deviation
        B['BCxdGLwl'][1,1]=H['gWe']-1.
        B['BCxdGLwl'][2,2]=H['gWmu']-1.
        B['BCxdGLwl'][3,3]=H['gWtau']-1.
        
        # NOTE: current UFO version does not contain any of the couplings
        # involving quarks analogous to the above.
        
        # Yukawa
        B['BCxdYd'][3,3] = sqrt( (H['kb']-1.)**2 + (H['lb'])**2 )
        
        B['BCxSd'][3,3] = 0. if B['BCxdYd'][3,3]==0. else H['lb']/B['BCxdYd'][3,3]
                
        B['BCxdYu'][2,2] = sqrt( (H['kc']-1.)**2 + (H['lc'])**2 )
        
        B['BCxSu'][2,2] = 0. if B['BCxdYu'][2,2]==0. else H['lc']/B['BCxdYu'][2,2]

        B['BCxdYe'][3,3] = sqrt( (H['ktau']-1.)**2 + (H['ltau'])**2 )
        
        B['BCxSe'][3,3] = 0. if B['BCxdYe'][3,3]==0. else H['ltau']/B['BCxdYe'][3,3]
        
        B['BCxdYe'][2,2] = sqrt( (H['kmu']-1.)**2 + (H['lmu'])**2 )
        
        B['BCxSe'][2,2] = 0. if B['BCxdYe'][3,3]==0. else H['lmu']/B['BCxdYe'][3,3]
        
        return B

    @basis.translation_from('bsmc')
    def from_bsmc(instance, self):
        
        s2w, c2w, ee2, gw2, gp2, MZ, vev, gs2 = instance.calculate_inputs()
        
        H = self
        B = instance
            
        # Gauge-Higgs
        H['kZZ'] = 1. + B['dCz'] + B['Czbx']*gw2
        
        H['kWW'] = 1. + B['dCw'] + B['Cwbx']*gw2
        
        H['eZZ'] = -B['Czz']*(gw2+gp2)/2.
        
        H['eZZCP'] = -B['tCzz']*(gw2+gp2)/2.
        
        H['kZA'] = -B['Cza']*sqrt(ee2*(gw2+gp2))/2.
        
        H['lZACP'] = -B['tCza']*sqrt(ee2*(gw2+gp2))/2.
        
        H['kAA'] = -B['Caa']*ee2/2.
        
        H['lAACP'] = -B['tCaa']*ee2/2.
        
        H['eWW'] = -B['Cww']*gw2/2.
        
        H['eWWCP'] = -B['tCww']*gw2/2.
        
        # Higgs-Z-f-f contact
        CLze = ( sqrt(gw2)**3/sqrt(c2w)/2./sqrt(2)*(s2w - 1./2.)*B['Czbx']
                - sqrt(gw2*gp2*ee2)/2.*B['Cabx'] )
        
        H['eZeL'] = CLze + sqrt(gw2+gp2)*B['BCxdGLhze'][1,1].real
        
        H['eZmuL'] = CLze + sqrt(gw2+gp2)*B['BCxdGLhze'][2,2].real
        
        H['eZtauL'] = CLze + sqrt(gw2+gp2)*B['BCxdGLhze'][3,3].real
        
        CRze = ( sqrt(gw2)**3/sqrt(c2w)/2./sqrt(2)*s2w*B['Czbx']
                - sqrt(gw2*gp2*ee2)/2.*B['Cabx'] )
        
        H['eZeR'] = CRze + sqrt(gw2+gp2)*B['BCxdGRhze'][1,1].real
        
        H['eZmuR'] = CRze + sqrt(gw2+gp2)*B['BCxdGRhze'][1,1].real
        
        H['eZtauR'] = CRze + sqrt(gw2+gp2)*B['BCxdGRhze'][1,1].real
        
        H['eZv'] = ( sqrt(gw2)**3/sqrt(c2w)/2./sqrt(2)*(1./2.)*B['Czbx'] 
                    + sqrt(gw2+gp2)*B['BCxdGLhzv'][1,1].real )
        
        # Higgs-W-f-f contact
        H['eWe'] = ( sqrt(gw2)**3/2./sqrt(2)*(1. + B['BCxdGLwl'][1,1].real)*B['Cwbx']
                    + sqrt(gw2/2.)* B['BCxdGLhwl'][1,1].real)
        
        H['eWmu'] = ( sqrt(gw2)**3/2./sqrt(2)*(1. + B['BCxdGLwl'][2,2].real)*B['Cwbx']
                    + sqrt(gw2/2.)* B['BCxdGLhwl'][2,2].real)
        
        H['eWtau'] = ( sqrt(gw2)**3/2./sqrt(2)*(1. + B['BCxdGLwl'][3,3].real)*B['Cwbx']
                    + sqrt(gw2/2.)* B['BCxdGLhwl'][3,3].real)
        
        
        # Z-f-f coupling
        
        H['gZeL'] = (- 1./2. + s2w + B['BCxdGLze'][1,1].real)
        
        H['gZmuL'] = (- 1./2. + s2w + B['BCxdGLze'][2,2].real)
        
        H['gZtauL'] =( - 1./2. + s2w + B['BCxdGLze'][3,3].real)
        
        H['gZeR'] = (s2w + B['BCxdGRze'][1,1].real)
        
        H['gZmuR'] = (s2w + B['BCxdGRze'][2,2].real)
        
        H['gZtauR'] = (s2w + B['BCxdGRze'][3,3].real)
        
        H['gZv'] = (1./2. + B['BCxdGLze'][1,1].real)
        
        # W-f-f coupling
        H['gWe'] = (1. + B['BCxdGLwl'][1,1].real)
        
        H['gWmu'] = (1. + B['BCxdGLwl'][2,2].real)
        
        H['gWtau'] = (1. + B['BCxdGLwl'][3,3].real)
        
        # Yukawa
        
        H['kb'] = 1.+ B['BCxdYd'][3,3]*sqrt(1.-B['BCxSd'][3,3]**2)
        H['lb'] = B['BCxdYd'][3,3]*B['BCxSd'][3,3]
        H['kc'] = 1.+ B['BCxdYu'][2,2]*sqrt(1.-B['BCxSu'][2,2]**2)
        H['lc'] = B['BCxdYu'][2,2]*B['BCxSu'][2,2]
        H['ktau'] = 1.+ B['BCxdYe'][3,3]*sqrt(1.-B['BCxSe'][3,3]**2)
        H['ltau'] = B['BCxdYe'][3,3]*B['BCxSe'][3,3]
        H['kmu'] = 1.+ B['BCxdYe'][2,2]*sqrt(1.-B['BCxSe'][2,2]**2)
        H['lmu'] = B['BCxdYe'][2,2]*B['BCxSe'][2,2]
        
        return H

################################################################################
