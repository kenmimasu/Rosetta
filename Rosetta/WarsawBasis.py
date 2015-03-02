from Basis import Basis
from MassBasis import MassBasis
import math
from itertools import combinations_with_replacement as comb
from itertools import product
from Rosetta import PID
####################################################################################################
# Warsaw basis class        
class WarsawBasis(Basis):
    independent = ['cH','cT','cGG','cWW','cBB','cWB','ctGG','ctWW','ctBB','ctWB','cpll',
              'cHl11','cHl12','cHl13','cHl22','cHl23','cHl33',
              'cpHl11','cpHl12','cpHl13','cpHl22','cpHl23','cpHl33',
              'cHe11','cHe12','cHe13','cHe22','cHe23','cHe33',
              'cHq11','cHq12','cHq13','cHq22','cHq23','cHq33',
              'cpHq11','cpHq12','cpHq13','cpHq22','cpHq23','cpHq33',
              'cHu11','cHu12','cHu13','cHu22','cHu23','cHu33',
              'cHd11','cHd12','cHd13','cHd22','cHd23','cHd33',
              'cHud11','cHud12','cHud13','cHud22','cHud23','cHud33',
              'cu11Re','cu12Re','cu13Re','cu22Re','cu23Re','cu33Re',
              'cu11Im','cu12Im','cu13Im','cu22Im','cu23Im','cu33Im',
              'cd11Re','cd12Re','cd13Re','cd22Re','cd23Re','cd33Re',
              'cd11Im','cd12Im','cd13Im','cd22Im','cd23Im','cd33Im',
              'ce11Re','ce12Re','ce13Re','ce22Re','ce23Re','ce33Re',
              'ce11Im','ce12Im','ce13Im','ce22Im','ce23Im','ce33Im']
    required_masses = set([y for x in PID.values() for y in x.values()])
    required_inputs = {'aEWM1','Gf','MZ','MH'}
    translate_to={'mass'}
        
    def calculate_inputs(self): # calculate a few required EW params from aEWM1, Gf, MZ
        self.input['ee2'] = 4.*math.pi/self.input['aEWM1'] # EM coupling squared
        self.input['s2w'] = (1.- math.sqrt(1. - self.input['ee2']/(math.sqrt(2.)*self.input['Gf']*self.input['MZ']**2)))/2. #sin^2(theta_W)
        self.input['gw2'] = self.input['ee2']/self.input['s2w'] # SU(2) coupling squared
        self.input['gp2'] = self.input['gw2']*self.input['s2w']/(1.-self.input['s2w']) # Hypercharge coupling squared
        self.input['vev'] =  2.*self.input['MZ']*math.sqrt(1-self.input['s2w'])/math.sqrt(self.input['gw2'])
        
    def translate(self):
        if self.target_basis=='mass': 
            self.translate_to_mass()
        else: 
            raise NotImplementedError
            
    def translate_to_mass(self):
        self.newname='Mass'
        self.calculate_inputs()
        s2w, ee2, gw2, gp2 = tuple([self.input[x] for x in ('s2w', 'ee2', 'gw2', 'gp2')]) # get EW params
        A = self.coeffs._asdict()
        B = MassBasis().coeffs._asdict()
        
        def f(T3,Q,i,j): # [eqn (4.8)]
            if i==j:
                return -Q*A['cWB']*gw2*gp2/(gw2-gp2) + (A['cT']-dv)*(T3 + Q*gp2/(gw2-gp2))
            else:
                return 0.
        def dy_sf(X,Y): # solution for dy,sin(phi) of [eqn (4.12)]
            R = math.sqrt(X**2+Y**2)
            sf = -Y/R # solution is +-Y/R 
            dy = ( X*abs(X)+Y*abs(Y) )/R
            return dy,sf
        # Higgs vev shift [eqn (4.5)]
        dv = (A['cpHl11']+A['cpHl22'])/2.-A['cpll'] 
        # W mass shift [eqn (4.6)]
        B['dM'] = 1./(gw2-gp2)*(gw2*A['cT']-gp2*gw2*A['cWB']-gp2*dv)
        # W/Z chiral coupling deviations
        for i,j in comb((1,2,3),2):
            B['dGLwl{}{}'.format(i,j)] = A['cpHl{}{}'.format(i,j)] + f(1./2.,0.,i,j) - f(-1./2.,-1.,i,j) # [eqn (4.7)]
            B['dGLzv{}{}'.format(i,j)] = 1./2.*A['cpHl{}{}'.format(i,j)] - 1./2.*A['cHl{}{}'.format(i,j)] + f(1./2.,0.,i,j)
            B['dGLze{}{}'.format(i,j)] = -1./2.*A['cpHl{}{}'.format(i,j)] - 1./2.*A['cHl{}{}'.format(i,j)] + f(-1./2.,-1.,i,j)
            B['dGRze{}{}'.format(i,j)] = - 1./2.*A['cHe{}{}'.format(i,j)] + f(0.,-1.,i,j)
            B['dGLwq{}{}'.format(i,j)] = A['cpHq{}{}'.format(i,j)] + f(1./2.,2./3.,i,j) - f(-1./2.,-1./3.,i,j) # [eqn (4.9)]
            B['dGRwq{}{}'.format(i,j)] = -1./2.*A['cHud{}{}'.format(i,j)] 
            B['dGLzu{}{}'.format(i,j)] = 1./2.*A['cpHq{}{}'.format(i,j)] - 1./2.*A['cHq{}{}'.format(i,j)] + f(1./2.,2./3.,i,j)
            B['dGLzd{}{}'.format(i,j)] = -1./2.*A['cpHq{}{}'.format(i,j)] - 1./2.*A['cHq{}{}'.format(i,j)] + f(-1./2.,-1./3.,i,j)
            B['dGRzu{}{}'.format(i,j)] = - 1./2.*A['cHu{}{}'.format(i,j)] + f(0.,2./3.,i,j)
            B['dGRzd{}{}'.format(i,j)] = - 1./2.*A['cHd{}{}'.format(i,j)] + f(0.,-1./3.,i,j)
        # Higgs couplings to W/Z [eqn (4.10)]
        B['dCw'] = -A['cH'] - 2.*A['cWB']*gw2*gp2/(gw2-gp2) + 2.*A['cT']*gw2/(gw2-gp2) - dv*(gw2+gp2)/(gw2-gp2)
        B['dCz'] = -A['cH'] - 2.*A['cT'] - dv 
        # Two derivative field strength interactions [eqn (4.11)]
        B['Cgg'] = A['cGG'] 
        B['Caa'] =  A['cWW'] + A['cBB'] - 4.*A['cWB']
        B['Czz'] = (gw2**2*A['cWW'] + gp2**2*A['cBB'] + 4.*gw2*gp2*A['cWB'])/(gw2+gp2)**2
        B['Cza'] = (gw2*A['cWW'] - gp2*A['cBB'] - 2.*(gw2-gp2)*A['cWB'])/(gw2+gp2)
        B['Cww'] =  A['cWW']
        B['CTgg'] = A['ctGG'] 
        B['CTaa'] =  A['ctWW'] + A['ctBB'] - 4.*A['ctWB']
        B['CTzz'] = (gw2**2*A['ctWW'] + gp2**2*A['ctBB'] + 4.*gw2*gp2*A['ctWB'])/(gw2+gp2)**2
        B['CTza'] = (gw2*A['ctWW'] - gp2*A['ctBB'] - 2.*(gw2-gp2)*A['ctWB'])/(gw2+gp2)
        B['CTww'] =  A['ctWW']
        # Yukawa type interaction coefficients [eqn. (4.12)]
        for i,j in comb((1,2,3),2): 
            for f in ('u','d','e'):
                mi, mj = self.mass[ PID[f][i] ],self.mass[ PID[f][j] ] 
                if mi and mj:
                    dy_cosphi = self.input['vev']*(A['c{}{}{}Re'.format(f,i,j)])
                    dy_sinphi = self.input['vev']*(A['c{}{}{}Im'.format(f,i,j)])
                    if i==j: dy_cosphi -= mi*(dv+A['cH'])
                    Rdysf = math.sqrt( dy_sinphi**2 + dy_cosphi**2)
                    B['S{}{}{}'.format(f,i,j)], B['dY{}{}{}'.format(f,i,j)] = dy_sf( dy_cosphi, dy_sinphi )
        # Four-point hVff interactions
        for i,j in comb((1,2,3),2): 
            B['CLwl{}{}'.format(i,j)] = A['cpHl{}{}'.format(i,j)] # [eqn (4.13)]
            B['CLwq{}{}'.format(i,j)] = A['cpHq{}{}'.format(i,j)]
            B['CRwq{}{}'.format(i,j)] = -A['cHud{}{}'.format(i,j)]/2.
            B['CLzv{}{}'.format(i,j)] = (A['cpHl{}{}'.format(i,j)]-A['cHl{}{}'.format(i,j)])/2. # [eqn (4.14)]
            B['CLze{}{}'.format(i,j)] = (-A['cpHl{}{}'.format(i,j)]-A['cHl{}{}'.format(i,j)])/2.
            B['CRze{}{}'.format(i,j)] = -A['cHe{}{}'.format(i,j)]/2.
            B['CLzu{}{}'.format(i,j)] = (A['cpHq{}{}'.format(i,j)]-A['cHq{}{}'.format(i,j)])/2. # [eqn (4.15)]
            B['CLzd{}{}'.format(i,j)] = (-A['cpHq{}{}'.format(i,j)]-A['cHq{}{}'.format(i,j)])/2.
            B['CRzu{}{}'.format(i,j)] = -A['cHu{}{}'.format(i,j)]/2.
            B['CRzd{}{}'.format(i,j)] = -A['cHd{}{}'.format(i,j)]/2.
        self.newpar = B
        
        self.newmass[24] = self.mass[24]+self.newpar['dM'] # W mass shift
        
        

####################################################################################################