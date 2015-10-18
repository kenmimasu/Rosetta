from internal import Basis
from internal.matrices import matrix_mult, matrix_add, matrix_sub, matrix_eq
import MassBasis as MB
import sys
import math, re
from math import sqrt
from itertools import combinations_with_replacement as comb
from itertools import product
from internal.constants import PID
################################################################################
MassBasis = MB.MassBasis
################################################################################
class HiggsBasis(Basis.Basis):
    '''
    Main basis class for Rosetta, based on the reccomendation for the 
    parametrisation of Higgs Effective Field theory in the LHC Higgs cross 
    section working group note (LHCHXSWG-INT-2015-001). Part of the three 
    intrinsic basis implementations in Rosetta along with the Warsaw and SILH 
    bases. This implementation includes almost all operators listed in the note 
    apart from the dipole-type Lorentz structures (\sigma^{\mu\nu}). The 
    implementation only includes the minimal set (3) of four-fermion operators 
    required to consistently map between the three intrinsic bases as defined 
    in the note. A number of structures involving more than two Higgs fields 
    (apart from the triple Higgs vertex) or interactions between a single higgs 
    field and 3 or more gauge bosons are not currently included. Besides this, 
    the Higgs Basis encodes relations between certain parameters to ensure 
    SU(2)xU(1) invariance such that it is consistent with a dimension six 
    operator basis for an effective field theory with linearly realised 
    electroweak symmetry breaking (in unitary gauge) and a general flavor 
    structure. These relations are implemented in calculate_dependent().
    '''
    
    name='higgs'
    ##########################
    # declare coefficients
    # kinetic terms [Eqn. (3.3)]
    HBxMASS = ['dM']
    # triple gauge couplings [Eqn. (3.6)]
    HBxTGC = ['Lz', 'tLz',  'C3g', 'tC3g', 'dKa', 'tKa', 
              'dG1z', 'dKz', 'tKz', 'La', 'tLa']
    # quartic gauge couplings [Eqn. (3.7)]
    HBxQGC = ['dGw4','dGw2z2','dGw2za', 
              'Lw4', 'Lw2z2', 'Lw2a2', 'Lw2az', 'Lw2za', 'C4g',
              'tLw4', 'tLw2z2', 'tLw2a2', 'tLw2az', 'tLw2za', 'tC4g']
    # single Higgs couplings [Eqn. (3.8)]
    HBxh = ['dCz', 
            'Cgg', 'Caa', 'Cza', 'Czz', 'Czbx', 'Cabx',
            'tCgg', 'tCaa', 'tCza', 'tCzz',
            'dCw', 'Cww', 'Cwbx', 'tCww' ]
    # double Higgs couplings [Eqn. (3.13)]
    HBxhh =  ['dCw2', 'dCz2', 'Cww2', 'Cwbx2', 
              'Cgg2', 'Caa2', 'Cza2', 'Czz2', 'Czbx2', 'Cabx2',
              'tCww2', 'tCgg2', 'tCaa2', 'tCza2', 'tCzz2']
    # Higgs self couplings [Eqn. (3.12)]
    HBxhself = ['dL3', 'dL4']
    # 4-fermion operators
    HBx4F = ['cll1122', 'cpuu3333', 'cll1221']
    ##########################
    # block structure
    blocks = {'HBxMASS':HBxMASS, 'HBxTGC':HBxTGC, 'HBxQGC':HBxQGC, 
              'HBxh':HBxh, 'HBxhh':HBxhh, 'HBxhself':HBxhself, 'HBx4F':HBx4F}
              
    # copy flavored block structure from MassBasis.MassBasis          
    flavored = {k.replace('MB','HB'):v for k,v in MassBasis.flavored.iteritems()} 

    # independent coefficients
    independent = [
    # [Eqn. (5.1)]
    'dM', 
    'HBxdGLze', 'HBxdGRze', 'HBxdGLwl', 'HBxdGLzu', 
    'HBxdGRzu', 'HBxdGLzd', 'HBxdGRzd', 'HBxdGRwq', 
    'HBxdgu', 'HBxdgd', 'HBxdau', 'HBxdad', 'HBxdae', 
    'HBxdzu', 'HBxdzd', 'HBxdze', 'HBxdwq', 'HBxdwl', 
    'HBxtdgu', 'HBxtdgd', 'HBxtdau', 'HBxtdad', 'HBxtdae', 
    'HBxtdzu', 'HBxtdzd', 'HBxtdze', 'HBxtdwq', 'HBxtdwl',
    # [Eqn. (5.2)]
    'Cgg', 'dCz', 'Caa', 'Cza', 'Czz', 'Czbx', 'tCgg', 'tCaa', 'tCza', 'tCzz', 
    'HBxdYu', 'HBxdYd', 'HBxdYe', 'HBxSu', 'HBxSd', 'HBxSe', 'dL3',
    # [Eqn. (5.3)]
    'Lz', 'tLz', 'C3g', 'tC3g',
    'cll1122','cpuu3333'
    ]
                
    # Required inputs/masses             
    required_masses = {25}
    required_inputs = {1, 2, 3, 4} # aEWM1, Gf, aS, MZ
    ########################## 
    
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

    def calculate_dependent(self):
        s2w, c2w, ee2, gw2, gp2, MZ, vev, gs2 = self.calculate_inputs() 
        A = self
        MW = MZ*sqrt(c2w)
        MH = self.mass[25]
        
        def delta(i,j):
            return 1. if i==j else 0.
        
        # Higgs and EW gauge bosons [Sec 5.2] [eqn (5.5)]
        A['dCw']  = A['dCz']  + A['dM']*4. 
        A['Cww']  = A['Czz']  + A['Cza']*2.*s2w  + A['Caa'] *s2w**2
        A['tCww'] = A['tCzz'] + A['tCza']*2.*s2w + A['tCaa']*s2w**2 
        A['Cwbx'] = (A['Czbx']*gw2 + A['Czz']*gp2 - A['Caa']*ee2*s2w 
                    - A['Cza']*(gw2-gp2)*s2w )/(gw2-gp2)
        A['Cabx'] = (A['Czbx']*2.*gw2 + A['Czz']*(gw2+gp2) 
                    - A['Caa']*ee2 - A['Cza']*(gw2-gp2))/(gw2-gp2)
        
        # Higgs self couplings [eqn (5.7)]
        A['dL4'] = 3./2.*A['dL3'] - MH**2/vev**2/6.*A['dCz']
        
        # dependent dgV coeffs [eqn (5.9)]
        matrix_add(A['HBxdGLze'], A['HBxdGLwl'], A['HBxdGLzv'])

        # dGLwq = dGLzu.VCKM - VCKM.dGLzd [eqn (5.9)]
        matrix_sub(matrix_mult(A['HBxdGLzu'], A.ckm),
                   matrix_mult(A.ckm, A['HBxdGLzd']),
                   A['HBxdGLwq'])

        # list of all z/w vertex correction blocks
        vertex = ['HBxdGLze', 'HBxdGRze', 'HBxdGLzv', 'HBxdGLzu', 'HBxdGRzu', 
                  'HBxdGLzd', 'HBxdGRzd', 'HBxdGLwl', 'HBxdGLwq', 'HBxdGRwq']
        # list of all z/w/a dipole interaction blocks
        dipole = ['HBxdgu', 'HBxdgd', 'HBxdau', 'HBxdad', 'HBxdae', 
                  'HBxdzu', 'HBxdzd', 'HBxdze', 'HBxdwq', 'HBxdwl', 
                  'HBxtdgu', 'HBxtdgd', 'HBxtdau', 'HBxtdad', 'HBxtdae', 
                  'HBxtdzu', 'HBxtdzd', 'HBxtdze', 'HBxtdwq', 'HBxtdwl']
                  
        # HVFF coeffs and dipole-like Higgs couplings [Eqns. (3.10) & (5.6)]
        for dG in vertex + dipole: 
            dGh = dG[:-2]+'h'+dG[-2:]
            matrix_eq(A[dG], A[dGh])
               
        # Triple gauge couplings [Sec 3.7] [eqn (3.21)] 
        A['dG1z'] = (A['Caa']*ee2*gp2 + A['Cza']*(gw2-gp2)*gp2 
                    - A['Czz']*(gw2+gp2)*gp2 - A['Czbx']*(gw2+gp2)*gw2 
                    )/2./(gw2-gp2)
        A['dKa'] = - (A['Caa']*ee2  + A['Cza']*(gw2-gp2) 
                    - A['Czz']*(gw2+gp2) )*gw2/2./(gw2+gp2)
        A['tKa'] = - ( A['tCaa']*ee2 + A['tCza']*(gw2-gp2) 
                    - A['tCzz']*(gw2+gp2))*gw2/2./(gw2+gp2)
        A['dKz'] = A['dG1z'] - gp2/gw2*A['dKa']
        A['tKz'] = - A['tKa']*gp2/gw2
        A['La'] = A['Lz']
        A['tLa'] = A['tLz']
        
        # Quartic gauge couplings [Sec 3.7] [eqn (3.23)] 
        A['dGw4'] = 2.*c2w*A['dG1z']
        A['dGw2z2'] = 2.*A['dG1z']
        A['dGw2za'] = A['dG1z']
        
        # two derivative quartic gauge couplings [Sec 3.7] [eqn (3.24)] 
        A['Lw4'] = -gw2/2./MW**2*A['Lz']
        A['tLw4'] = -gw2/2./MW**2*A['tLz']
        A['Lw2z2'] = -gw2*c2w/MW**2*A['Lz']
        A['tLw2z2'] = -gw2*c2w/MW**2*A['tLz']
        A['Lw2za'] = -ee2/MW**2*A['Lz']
        A['tLw2za'] = -ee2/MW**2*A['tLz']
        A['Lw2a2'] = -sqrt(ee2*gw2*c2w)/MW**2*A['Lz']
        A['tLw2a2'] = -sqrt(ee2*gw2*c2w)/MW**2*A['tLz']
        A['Lw2az'] = -sqrt(ee2*gw2*c2w)/MW**2*A['Lz']
        A['tLw2az'] = -sqrt(ee2*gw2*c2w)/MW**2*A['tLz']
        A['C4g'] = 3.*sqrt(gs2)**3/vev**2*A['C3g']
        A['tC4g'] = 3.*sqrt(gs2)**3/vev**2*A['tC3g']
        
        # Couplings of two Higgs bosons [Sec 3.8] [eqn (3.27)]
        A['Cgg2'], A['tCgg2'] = A['Cgg'], A['tCgg']
        A['dCz2'] = A['dCz']
        A['dCw2'] = A['dCw'] + 3.*A['dM']
        
        hvv = ['Cww', 'Cwbx', 'Cgg', 'Caa', 'Cza', 'Czz', 'Czbx', 'Cabx',
               'tCww', 'tCgg', 'tCaa', 'tCza', 'tCzz']
        for cvv in hvv:
            cvv2 = cvv + '2'
            A[cvv2] = A[cvv]
        
        for f in ('u','d','e'):
            yuk = 'HBxdY' + f
            sin = 'HBxS' + f
            yuk2 = 'HBxY2' + f
            for i,j in A[yuk].keys():
                Yij = A[yuk][i,j]
                sinij = A[sin][i,j]
                cosij = sqrt(1. - sinij**2)
                re, im = (3.*Yij*cosij - A['dCz']*delta(i,j)), 3.*Yij*sinij
                A[yuk2][i,j] = complex(re, im)

        # 4-fermion operators [Sec. 3.9]
        # [eqn (3.32)]
        A['cll1221'] = (A['HBxdGLwl'][1,1].real 
                      + A['HBxdGLwl'][2,2].real - 2.*A['dM'])*2.
        
        
    @Basis.translation('mass')        
    def to_mass(self, instance):
        # trivial translation
        for k, v in self.iteritems():
            instance[k] = v
        return instance
        
    @Basis.translation('warsaw')
    def to_warsaw(self, wbinstance):

        def delta(i,j):
            return 1. if i==j else 0.
        
        H = self
        W = wbinstance

        s2w, c2w, ee2, gw2, gp2, MZ, vev, gs2 = self.calculate_inputs()
        MH = self.mass[25]
        dM, dCz = H['dM'], H['dCz']
        
        W['cll1221'] = (H['HBxdGLwl'][1,1].real 
                      + H['HBxdGLwl'][2,2].real - 2.*dM)*2.
        
        for w,h in [('',''),('t','T')]: # CP even and odd sum
            Cgg, Caa = H['%sCgg' % h], H['%sCaa' % h]
            Czz, Cza = H['%sCzz' % h], H['%sCza' % h]
            W['%scGG' % w] = Cgg
        
            W['%scWW' % w] = Czz + ( Caa*gp2 + 2.*Cza*(gw2 + gp2) 
                                   )*gp2/(gw2 + gp2)**2
                            
            W['%scBB' % w] = Czz + ( Caa*gw2 - 2.*Cza*(gw2 + gp2) 
                                   )*gw2/(gw2 + gp2)**2
        
            W['%scWB' % w] = Czz/2. - ( Caa*gw2*gp2 + Cza*(gw2**2-gp2**2) 
                                      )/2./(gw2 + gp2)**2
        Czz, Caa = H['Czz'], H['Caa']
        Czbx, Cza = H['Czbx'], H['Cza']
        
        W['cH'] = ( ( Caa*gw2**2*gp2**2 
                    + Cza*gw2*gp2*(gw2**2 - gp2**2) 
                    - Czz*gw2*gp2*(gp2 + gw2)**2
                    - Czbx*gw2**2*(gp2 + gw2)**2
                    )*3./(2.*(gw2 - gp2)*(gw2 + gp2)**2)
                   - 3.*dM - dCz)
        cH = W['cH']
        
        W['cT'] = ((- Caa*gp2*gw2 
                    - Cza*(gw2**2-gp2**2) 
                    + (Czbx+Czz)*(gp2 + gw2)**2
                   )*gp2*gw2/(2.*(gw2 - gp2)*(gw2 + gp2)**2)
                    + dM)
        cT = W['cT']
        for i,j in H['HBxdGLzu'].keys():
            fac = (dM-cT)*delta(i,j)
            facp =  -(dM + (cH + dCz)/3.)*delta(i,j)
            W['WBxHq'][i,j] = (fac/3. - H['HBxdGLzu'][i,j] - H['HBxdGLzd'][i,j])
            W['WBxHu'][i,j] = 4./3.*fac - 2.*H['HBxdGRzu'][i,j]
            W['WBxHd'][i,j] = -2./3.*fac - 2.*H['HBxdGRzd'][i,j] 
            W['WBxHpq'][i,j] = (facp + H['HBxdGLzu'][i,j] - H['HBxdGLzd'][i,j])
            W['WBxHl'][i,j] = -(fac + 2.*H['HBxdGLze'][i,j] + H['HBxdGLwl'][i,j]) 
            W['WBxHpl'][i,j] = facp + H['HBxdGLwl'][i,j] 
            W['WBxHe'][i,j] = -2.*fac - 2.*H['HBxdGRze'][i,j] 
        for k,v in H['HBxdGRwq'].iteritems():
            W['WBxHud'][k] = -2.*v
                
        W['c6H'] = (3.*dCz + 8.*dM + 
                    4*( - Caa*gw2**2*gp2**2 
                        - Cza*gw2*gp2*(gw2**2-gp2**2) 
                        + Czz*gw2*gp2*(gw2+gp2)**2
                        + Czbx*gw2**2*(gw2+gp2)**2
                      )/(gw2 - gp2)/(gw2+gp2)**2
                    )*MH**2/(2.*vev**2) - H['dL3']

        for f in ('u','d','e'):
            for i,j in H['HBxdY'+f].keys(): 
                diag = delta(i,j)*(dCz-2.*cH)/3.
                yuk = H['HBxdY'+f][i,j]
                sin = H['HBxS'+f][i,j]
                cos = sqrt(1.-sin**2)
                re = ( yuk*cos - diag )*sqrt(2.)
                im = yuk*sin*sqrt(2.)
                W['WBx'+f][i,j] = complex(re, im)
        
        W['c3G'], W['tc3G'] = H['C3g'], H['tC3g']
        W['c3W'], W['tc3W'] = -2./3./gw2**2*H['Lz'], -2./3./gw2**2*H['tLz']
        W['cll1122'], W['cpuu3333'] = H['cll1122'], H['cpuu3333']

        return W
        
    @Basis.translation('silh')
    def to_silh(self,instance):
        
        H = self
        S = instance
        
        s2w, c2w, ee2, gw2, gp2, MZ, vev, gs2 = self.calculate_inputs()
        MH = self.mass[25]
        dM, dCz = H['dM'], H['dCz']
        
        def delta(i,j):
            return 1. if i==j else 0.
        
        for s,h in [('',''),('t','T')]: # CP even and odd sum
            Cgg, Caa = H['%sCgg' % h], H['%sCaa' % h]
            Czz, Cza = H['%sCzz' % h], H['%sCza' % h]
            S['%ssGG' % s] = Cgg
        
            S['%ssHW' % s] = - ( gp2**2*Caa + 2.*gp2*(gp2 + gw2)*Cza
                               )/(gw2 + gp2)**2 - Czz
            
            S['%ssHB' % s] = ( gp2*(gp2 + 2.*gw2)*Caa + 2.*gw2*(gp2 + gw2)*Cza
                             )/(gw2 + gp2)**2 - Czz
            
            S['%ssBB' % s] = Caa
        
        Czz, Caa = H['Czz'], H['Caa']
        Czbx, Cza = H['Czbx'], H['Cza']
        
        S['sBB'] = Caa
        S['s2W'] = (H['HBxdGLwl'][1,1].real 
                  + H['HBxdGLwl'][2,2].real  - 2.*dM)*2./gw2
        S['s2B'] = (H['cll1122'] + H['HBxdGLwl'][1,1].real  
                   + H['HBxdGLwl'][2,2].real - 2.*dM)*2./gp2
                   
        S['s2G'] = 3.*H['cpuu3333']/gs2

        S['sB'] = (( gp2**2*Caa 
                   - 2.*gw2*(gp2+gw2)*Czbx
                   - (gp2+gw2)**2*Czz
                   )/(gw2**2-gp2**2) - 
                   4.*( H['cll1122'] - 2.*H['HBxdGLze'][1,1,].real 
                      + H['HBxdGLwl'][2,2].real - 2.*dM
                      )/gp2)
                      
        S['sW'] = (( -gp2**2*Caa 
                   + 2.*gw2*(gp2+gw2)*Czbx
                   + (gp2+gw2)**2*Czz
                   )/(gw2**2-gp2**2) - 
                   4.*(  H['HBxdGLwl'][2,2].real - 2.*dM
                      )/gw2)
        
        S['sHW'] = -( gp2**2*Caa 
                    + 2.*gp2*(gp2 + gw2)*Cza
                    )/(gw2 + gp2)**2 - Czz
                    
        S['sHB'] = ( gp2*(gp2 + 2.*gw2)*Caa 
                   + 2.*gw2*(gp2 + gw2)*Cza
                   )/(gw2 + gp2)**2 - Czz
        
        S['sH'] = (H['HBxdGLwl'][1,1].real - H['HBxdGLwl'][2,2].real)*3./2.-dCz
        S['sT'] = (H['HBxdGLwl'][1,1].real - H['HBxdGLwl'][2,2].real 
                 - H['cll1122'] + 4.*dM + 4.*H['HBxdGLze'][1,1].real )/2. 
        
        S['s3G'], S['ts3G'] = H['C3g'], H['tC3g']
        S['s3W'], S['ts3W'] = -2./3./gw2**2*H['Lz'], -2./3./gw2**2*H['tLz']
        
        S['s6H'] = (3.*dCz - 4.*H['HBxdGLwl'][1,1].real 
                  + 4.*H['HBxdGLwl'][2,2].real )*MH**2/(2.*vev**2) - H['dL3']
        
        for i, j in H['HBxdGLzu'].keys():
            fac = -delta(i,j)*(4.*H['HBxdGLze'][1,1].real 
                              + 2.*H['HBxdGLwl'][1,1].real)
            facp = -delta(i,j)*H['HBxdGLwl'][1,1].real
        
            S['SBxHq'][i,j] = (fac/6. - H['HBxdGLzu'][i,j]
                                - H['HBxdGLzd'][i,j]) 
            S['SBxHu'][i,j] = 2./3.*fac - 2.*H['HBxdGRzu'][i,j]
            S['SBxHd'][i,j] = -1./3.*fac - 2.*H['HBxdGRzd'][i,j] 
            S['SBxHpq'][i,j] = (facp + H['HBxdGLzu'][i,j] 
                                 - H['HBxdGLzd'][i,j])
            S['SBxHl'][i,j] = (- fac/2. - 2.*H['HBxdGLze'][i,j] 
                                - H['HBxdGLwl'][i,j]) 
            S['SBxHpl'][i,j] = facp + H['HBxdGLwl'][i,j] 
            S['SBxHe'][i,j] = - fac - 2.*H['HBxdGRze'][i,j]
            
        for k, v in H['HBxdGRwq'].iteritems():
            S['SBxHud'][k] = -2.*v            
        
        for f in ('u','d','e'):
            for i,j in H['HBxdY'+f].keys(): 
                diag = delta(i,j)*(dCz - H['HBxdGLwl'][1,1].real 
                                  + H['HBxdGLwl'][2,2].real)                
                yuk = H['HBxdY'+f][i,j]
                sin = H['HBxS'+f][i,j]
                cos = sqrt(1.-sin**2)
                
                re = ( yuk*cos - diag )*sqrt(2.)
                im = yuk*sin*sqrt(2.)
                S['SBx'+f][i,j] = complex(re, im)

        return S
    
    def modify_inputs(self):
        '''
        W mass modification from dM.
        '''
        s2w, c2w, ee2, gw2, gp2, MZ, vev, gs2 = self.calculate_inputs() 
        MW = MZ*sqrt(c2w)
        if 24 in self.mass:
            self.mass[24] = MW + self['dM']
        else:
            self.mass.new_entry(24, MW + self['dM'], name = 'MW')
        
        
################################################################################ 