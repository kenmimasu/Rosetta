from Basis import Basis
import math, re
from itertools import combinations_with_replacement as comb
from itertools import product
from __init__ import PID
####################################################################################################
# Higgs basis class
class HiggsBasis(Basis):
    independent = ['dCw','dCz','Cgg','Czz','Caa','Cza','CTgg','CTzz','CTaa','CTza',
                   'dYu11','dYu12','dYu13','dYu22','dYu23','dYu33',
                   'dYd11','dYd12','dYd13','dYd22','dYd23','dYd33',
                   'Su11','Su12','Su13','Su22','Su23','Su33',
                   'Sd11','Sd12','Sd13','Sd22','Sd23','Sd33',
                   'dYe11','dYe12','dYe13','dYe22','dYe23','dYe33',
                   'Se11','Se12','Se13','Se22','Se23','Se33',
                   'dGLze11','dGLze12','dGLze13','dGLze22','dGLze23','dGLze33',
                   'dGRze11','dGRze12','dGRze13','dGRze22','dGRze23','dGRze33',
                   'dGLzu11','dGLzu12','dGLzu13','dGLzu22','dGLzu23','dGLzu33',
                   'dGLzd11','dGLzd12','dGLzd13','dGLzd22','dGLzd23','dGLzd33',
                   'dGRzu11','dGRzu12','dGRzu13','dGRzu22','dGRzu23','dGRzu33',
                   'dGRzd11','dGRzd12','dGRzd13','dGRzd22','dGRzd23','dGRzd33',
                   'dGLwl11','dGLwl12','dGLwl13','dGLwl22','dGLwl23','dGLwl33',
                   'dGRwq11','dGRwq12','dGRwq13','dGRwq22','dGRwq23','dGRwq33','dM']
    dependent   = ['Cww','CTww',
                   'dGLzv11','dGLzv12','dGLzv13','dGLzv22','dGLzv23','dGLzv33',
                   'dGLwq11','dGLwq12','dGLwq13','dGLwq22','dGLwq23','dGLwq33',
                   'CLwq11','CLwq12','CLwq13','CLwq22','CLwq23','CLwq33',
                   'CRwq11','CRwq12','CRwq13','CRwq22','CRwq23','CRwq33',
                   'CLwl11','CLwl12','CLwl13','CLwl22','CLwl23','CLwl33',
                   'CLze11','CLze12','CLze13','CLze22','CLze23','CLze33',
                   'CLzv11','CLzv12','CLzv13','CLzv22','CLzv23','CLzv33',
                   'CRze11','CRze12','CRze13','CRze22','CRze23','CRze33',
                   'CLzu11','CLzu12','CLzu13','CLzu22','CLzu23','CLzu33',
                   'CLzd11','CLzd12','CLzd13','CLzd22','CLzd23','CLzd33',
                   'CRzu11','CRzu12','CRzu13','CRzu22','CRzu23','CRzu33',
                   'CRzd11','CRzd12','CRzd13','CRzd22','CRzd23','CRzd33']
                   
    required_masses = {1, 2, 3, 4, 5, 6, 11, 12, 13, 14, 15, 16}
    required_inputs = {1, 2, 4} # aEWM1, Gf, MZ

    def calculate_inputs(self): # calculate a few required EW params from aEWM1, Gf, MZ
        self.input['ee2'] = 4.*math.pi/self.input['aEWM1'] # EM coupling squared
        self.input['s2w'] = (1.- math.sqrt(1. - self.input['ee2']/(math.sqrt(2.)*self.input['Gf']*self.input['MZ']**2)))/2. #sin^2(theta_W)
        self.input['gw2'] = self.input['ee2']/self.input['s2w'] # SU(2) coupling squared
        self.input['gp2'] = self.input['gw2']*self.input['s2w']/(1.-self.input['s2w']) # Hypercharge coupling squared
        
    def calculate_dependent(self): # calculate dependent parameters
        p = self.par_dict
        self.calculate_inputs() # set a few useful EW parameters
        s2w, ee2, gw2, gp2 = tuple([self.input[x] for x in ('s2w', 'ee2', 'gw2', 'gp2')]) # get EW params
        
        def hVff(cname): # calculate hVff coeffs [eqn (3.11)]
            def T3(f,chi): # return t3 
                t3 = {'u':0.5,'d':-0.5,'v':0.5,'e':-0.5}
                return 0. if chi=='R' else t3[f]
            def Y(f,chi): # return hypercharge
                y = {'u':{'L':1./6.,'R':2./3.},'d':{'L':1./6.,'R':-1./3.},'v':{'L':-0.5,'R':0.},'e':{'L':-0.5,'R':-1.}}
                return y[f][chi]
            
            match = re.match(r'C(\w)(\w)(\w)(\d)(\d)',cname) # parse coefficient name
            chi, V, f, i,j = tuple([match.group(i) for i in range(1,6)])
            dgV  = cname.replace('C','dG') # corresponding V coupling deviation
            cVff = p[dgV]
            if i==j: # diagonal part
                if V=='z':
                    cVff += ( (p['dCw'] - p['dCz'] - 4*p['dM'])*( T3(f,chi)*gw2/gp2 + Y(f,chi) ) \
                    + T3(f,chi)*gw2*( p['Caa']*ee2/(gw2+gp2) + p['Cza']*(gw2-gp2)/(gw2+gp2) -p['Czz'] ))/2.
                elif V=='z' and chi=='L':
                    cVff +=( (p['dCw'] - p['dCz'] - 4*p['dM'])*gw2/gp2 \
                    + gw2*( p['Caa']*ee2/(gw2+gp2) + p['Cza']*(gw2-gp2)/(gw2+gp2) - p['Czz'] ))/2.
            return cVff
            
        p['Cww'] = p['Czz'] + 2.*s2w*p['Cza'] + s2w**2*p['Caa'] # Cww [eqn (3.9)]
        p['CTww'] = p['CTzz'] + 2.*s2w*p['CTza'] + s2w**2*p['CTaa'] # Ctww [eqn (3.9)]
        for i,j in comb((1,2,3),2):# dependent dgV coeffs [eqn (3.7)]
            p['dGLzv{}{}'.format(i,j)] = p['dGLze{}{}'.format(i,j)] + p['dGLwl{}{}'.format(i,j)] 
            p['dGLwq{}{}'.format(i,j)] = p['dGLzu{}{}'.format(i,j)] - p['dGLzd{}{}'.format(i,j)]
            for f,chi in product(('u','d','e','v'),('L','R')): # hVff 4-point coeffs [eqn (3.11)]
                if f=='v' and chi=='R': continue
                cvff = 'C{}z{}{}{}'.format(chi,f,i,j)
                p[cvff] = hVff(cvff)
            for f,chi in (('q','L'),('q','R'),('l','L')):
                cvff = 'C{}w{}{}{}'.format(chi,f,i,j)
                p[cvff] = hVff(cvff)

####################################################################################################
