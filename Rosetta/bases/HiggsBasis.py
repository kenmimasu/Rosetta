# import sys
import math
import re
from math import sqrt
#
from ..internal import basis
from ..internal import PID
from ..internal import matrix_mult, matrix_add, matrix_sub, matrix_eq
from . import BSMCharacterisation as BSMC
################################################################################
BSMCharacterisation = BSMC.BSMCharacterisation
################################################################################
class HiggsBasis(basis.Basis):
    '''
    Main basis class for Rosetta, based on the recommendation for the 
    parametrisation of Higgs Effective Field theory in the LHC Higgs cross 
    section working group note, LHCHXSWG-INT-2015-001 
    
    https://cds.cern.ch/record/2001958
    
    Part of the three intrinsic basis implementations in Rosetta along with the 
    Warsaw and SILH bases. The implementation only includes the minimal set (3) 
    of four-fermion operators required to consistently map between the three 
    intrinsic bases as defined in the note. A number of structures involving 
    more than two Higgs fields (apart from the triple Higgs vertex) or 
    interactions between a single higgs field and 3 or more gauge bosons are 
    not currently included. 
    The Higgs Basis encodes relations between certain parameters to ensure 
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
    # 4-fermion operators. cll1221 defined as dependent
    HBx4F = ['cll1111','cll1122','cll1221','cll1133','cll1331','cll2332',
             'cle1111','cle1122','cle2211','cle1133','cle3311',
             'cee1111','cee1122','cee1133',
             'cpuu3333']
    
    ##########################
    # block structure
    blocks = {'HBxMASS':HBxMASS, 'HBxTGC':HBxTGC, 'HBxQGC':HBxQGC, 
              'HBxh':HBxh, 'HBxhh':HBxhh, 'HBxhself':HBxhself, 'HBx4F':HBx4F}
              
    # same flavored block structure as BSMCharacterisation except dipole 
    # operators, which are general complex matrices instead of a pair of 
    # hermitian ones for the CP conserving and CP-violating interactions
    flavored = {
        # Z Vertex Corrections [Eqn (3.4)]
        'HBxdGLze': {'kind':'hermitian', 'domain':'complex', 'cname':'dGLze'},
        'HBxdGRze': {'kind':'hermitian', 'domain':'complex', 'cname':'dGRze'},
        'HBxdGLzv': {'kind':'hermitian', 'domain':'complex', 'cname':'dGLzv'},
        'HBxdGLzu': {'kind':'hermitian', 'domain':'complex', 'cname':'dGLzu'},
        'HBxdGRzu': {'kind':'hermitian', 'domain':'complex', 'cname':'dGRzu'},
        'HBxdGLzd': {'kind':'hermitian', 'domain':'complex', 'cname':'dGLzd'},
        'HBxdGRzd': {'kind':'hermitian', 'domain':'complex', 'cname':'dGRzd'},
        # W Vertex Corrections [Eqn (3.4)]
        'HBxdGLwl': {'kind':'hermitian', 'domain':'complex', 'cname':'dGLwl'},
        'HBxdGLwq': {'kind':'general', 'domain':'complex', 'cname':'dGLwq'},
        'HBxdGRwq': {'kind':'general', 'domain':'complex', 'cname':'dGRwq'},
        # Dipole interactions with single gauge bosons [Eqn. (3.5)]
        'HBxdgu': {'kind':'general', 'domain':'complex', 'cname':'dgu'},
        'HBxdgd': {'kind':'general', 'domain':'complex', 'cname':'dgd'},
        'HBxdau': {'kind':'general', 'domain':'complex', 'cname':'dau'},
        'HBxdad': {'kind':'general', 'domain':'complex', 'cname':'dad'},
        'HBxdae': {'kind':'general', 'domain':'complex', 'cname':'dae'},
        'HBxdzu': {'kind':'general', 'domain':'complex', 'cname':'dzu'},
        'HBxdzd': {'kind':'general', 'domain':'complex', 'cname':'dzd'},
        'HBxdze': {'kind':'general', 'domain':'complex', 'cname':'dze'},
        'HBxdwu': {'kind':'general', 'domain':'complex', 'cname':'dwu'},
        'HBxdwd': {'kind':'general', 'domain':'complex', 'cname':'dwd'},
        'HBxdwl': {'kind':'general', 'domain':'complex', 'cname':'dwl'},
        # single Higgs couplings to fermions [Eqn. (3.8)]
        'HBxdYu': {'kind':'general', 'domain':'real', 'cname':'dYu'},
        'HBxdYd': {'kind':'general', 'domain':'real', 'cname':'dYd'},
        'HBxdYe': {'kind':'general', 'domain':'real', 'cname':'dYe'},
        'HBxSu': {'kind':'general', 'domain':'real', 'cname':'Su' },
        'HBxSd': {'kind':'general', 'domain':'real', 'cname':'Sd' },
        'HBxSe': {'kind':'general', 'domain':'real', 'cname':'Se' },
        # Higgs contact interactions HVff [Eqn. (3.10)]
        'HBxdGLhze': {'kind':'hermitian', 'domain':'complex', 'cname':'dGLhze'},
        'HBxdGRhze': {'kind':'hermitian', 'domain':'complex', 'cname':'dGRhze'},
        'HBxdGLhzv': {'kind':'hermitian', 'domain':'complex', 'cname':'dGLhzv'},
        'HBxdGLhzu': {'kind':'hermitian', 'domain':'complex', 'cname':'dGLhzu'},
        'HBxdGRhzu': {'kind':'hermitian', 'domain':'complex', 'cname':'dGRhzu'},
        'HBxdGLhzd': {'kind':'hermitian', 'domain':'complex', 'cname':'dGLhzd'},
        'HBxdGRhzd': {'kind':'hermitian', 'domain':'complex', 'cname':'dGRhzd'},
        'HBxdGLhwl': {'kind':'hermitian', 'domain':'complex', 'cname':'dGLhwl'},
        'HBxdGLhwq': {'kind':'general', 'domain':'complex', 'cname':'dGLhwq'},
        'HBxdGRhwq': {'kind':'general', 'domain':'complex', 'cname':'dGRhwq'},
        # Dipole interactions with single higgs and gauge boson [Eqn. (3.11)]
        'HBxdhgu': {'kind':'general', 'domain':'complex', 'cname':'dhgu'},
        'HBxdhgd': {'kind':'general', 'domain':'complex', 'cname':'dhgd'},
        'HBxdhau': {'kind':'general', 'domain':'complex', 'cname':'dhau'},
        'HBxdhad': {'kind':'general', 'domain':'complex', 'cname':'dhad'},
        'HBxdhae': {'kind':'general', 'domain':'complex', 'cname':'dhae'},
        'HBxdhzu': {'kind':'general', 'domain':'complex', 'cname':'dhzu'},
        'HBxdhzd': {'kind':'general', 'domain':'complex', 'cname':'dhzd'},
        'HBxdhze': {'kind':'general', 'domain':'complex', 'cname':'dhze'},
        'HBxdhwu': {'kind':'general', 'domain':'complex', 'cname':'dhwu'},
        'HBxdhwd': {'kind':'general', 'domain':'complex', 'cname':'dhwd'},
        'HBxdhwl': {'kind':'general', 'domain':'complex', 'cname':'dhwl'},
        # couplings of two Higgs bosons to fermions [Sec. 3.8]
        'HBxY2u': {'kind':'general', 'domain':'complex', 'cname':'Y2u'},
        'HBxY2d': {'kind':'general', 'domain':'complex', 'cname':'Y2d'},
        'HBxY2e': {'kind':'general', 'domain':'complex', 'cname':'Y2e'}
    }
    
    # independent coefficients
    independent = [
    # [Eqn. (5.1)]
    'dM', 
    'HBxdGLze', 'HBxdGRze', 'HBxdGLwl', 'HBxdGLzu', 
    'HBxdGRzu', 'HBxdGLzd', 'HBxdGRzd', 'HBxdGRwq', 
    'HBxdgu', 'HBxdgd', 'HBxdau', 'HBxdad', 'HBxdae', 
    'HBxdzu', 'HBxdzd', 'HBxdze',
    # [Eqn. (5.2)]
    'Cgg', 'dCz', 'Caa', 'Cza', 'Czz', 'Czbx', 'tCgg', 'tCaa', 'tCza', 'tCzz', 
    'HBxdYu', 'HBxdYd', 'HBxdYe', 'HBxSu', 'HBxSd', 'HBxSe', 'dL3',
    # [Eqn. (5.3)]
    'Lz', 'tLz', 'C3g', 'tC3g',
    'cll1111','cll1122','cll1133','cll1331','cll2332',
    'cle1111','cle1122','cle2211','cle1133','cle3311',
    'cee1111','cee1122','cee1133',
    'cpuu3333'
    ]
                
    # Required inputs/masses
    # Higgs, W & fermion masses
    required_masses = {25, 24, 1, 2, 3, 4, 5, 6, 11, 13, 15} 
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

        MH = self.mass[25]
        MW = self.mass[24]
               
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
        
        # W dipole interaction
        for f in ('u','d','e'):
            eta = 1. if f=='u' else -1.
            
            if f=='e':
                wdip = 'HBxdwl'
            else:
                wdip = 'HBxdw'+f
                
            adip, zdip = 'HBxda'+f, 'HBxdz'+f

            for i,j in A[wdip].keys():
                A[wdip][i,j] = eta*( A[zdip][i,j] + s2w*A[adip][i,j] )

        # list of all z/w vertex correction blocks
        vertex = ['HBxdGLze', 'HBxdGRze', 'HBxdGLzv', 'HBxdGLzu', 'HBxdGRzu', 
                  'HBxdGLzd', 'HBxdGRzd', 'HBxdGLwl', 'HBxdGLwq', 'HBxdGRwq']
        # list of all z/w/a dipole interaction blocks
        dipole = ['HBxdgu', 'HBxdgd', 'HBxdau', 'HBxdad', 'HBxdae', 
                  'HBxdzu', 'HBxdzd', 'HBxdze', 'HBxdwu', 'HBxdwd', 'HBxdwl', 
                  # 'HBxtdgu', 'HBxtdgd', 'HBxtdau', 'HBxtdad', 'HBxtdae',
                  # 'HBxtdzu', 'HBxtdzd', 'HBxtdze'
              ]
                  
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
        # Gauge
        A['Cgg2'], A['tCgg2'] = A['Cgg'], A['tCgg']
        A['dCz2'] = A['dCz']
        A['dCw2'] = A['dCw'] + 3.*A['dM']
        
        hvv = ['Cww', 'Cwbx', 'Cgg', 'Caa', 'Cza', 'Czz', 'Czbx', 'Cabx',
               'tCww', 'tCgg', 'tCaa', 'tCza', 'tCzz']
        for cvv in hvv:
            cvv2 = cvv + '2'
            A[cvv2] = A[cvv]
        
        # Yukawas
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
        A['cll1221'] = ( A['HBxdGLwl'][1,1].real 
                       + A['HBxdGLwl'][2,2].real - 2.*A['dM'])*2.
        
        
    @basis.translation('bsmc')        
    def to_bsmc(self, instance):
        # trivial translation apart from dipole terms
        H = self
        B = instance
        dipoles = ['Cdau','Cdzu']
        for k, v in self.iteritems():
            
            is_dipole = re.match(r'Cdh{0,1}[a,z,g][u,d,e]\dx\d', k)
            
            if not is_dipole:
                try:
                    B[k] = v
                except KeyError:
                    print '    ' + k + ' not found in BSMC definition.'
        
        # splitting dipole coefficients
        ii = complex(0.,1.)
        for f in ('u','d','e'):

            glu, tglu = 'BCxdg'+f, 'BCxtdg'+f
            pho, tpho = 'BCxda'+f, 'BCxtda'+f
            zed, tzed = 'BCxdz'+f, 'BCxtdz'+f

            for i,j in B[zed].keys():
                if f in ('u','d'):
                    Gij, Gji = H['HBxdg'+f][i,j], H['HBxdg'+f][j,i]
                    B[glu][i,j] = -(Gij + Gji.conjugate())
                    B[tglu][i,j] =  -ii*(Gij - Gji.conjugate())

                Aij, Aji = H['HBxda'+f][i,j], H['HBxda'+f][j,i]
                B[pho][i,j] = -(Aij + Aji.conjugate())
                B[tpho][i,j] =  -ii*(Aij - Aji.conjugate())

                Zij, Zji = H['HBxdz'+f][i,j], H['HBxdz'+f][j,i]
                B[zed][i,j] = -(Zij + Zji.conjugate())
                B[tzed][i,j] =  -ii*(Zij - Zji.conjugate())
        
        return instance


    @basis.translation('warsaw')
    def to_warsaw(self, wbinstance):

        def delta(i,j):
            return 1. if i==j else 0.
        
        H = self
        W = wbinstance

        s2w, c2w, ee2, gw2, gp2, MZ, vev, gs2 = self.calculate_inputs()
        
        gw, gp, gs = sqrt(gw2), sqrt(gp2), sqrt(gs2)
        
        MH = self.mass[25]
        
        lam = MH**2/(2.*vev**2) # Higgs self-coupling     

        dM, dCz = H['dM'], H['dCz']
        
        W['cH'] = ( 3.*lam*dCz + 8.*lam*dM - H['dL3']
                  + ( H['Czbx']*gw2 + H['Czz']*gp2
                    - H['Cza']*gp2*(c2w-s2w) - H['Caa']*gw2*s2w**2
                    )*4.*gw2*lam/(gw2-gp2) 
                  )
        
        W['cHbx'] = ( dCz + 2.*dM
                    + ( H['Czbx']*(3.*gw2-gp2) + 2.*H['Czz']*gp2 
                      - 2.*H['Cza']*gp2*(c2w-s2w) - 2.*H['Caa']*gw2*s2w**2
                      )*gw2/(gw2-gp2)/2.
                    )
        
        W['cHD'] = ( - ( H['Czbx'] + H['Czz'] - H['Cza']*(c2w-s2w) 
                       - H['Caa']*s2w*c2w
                       )*2.*gw2*gp2/(gw2-gp2) - 4.*dM
                   )
                  
        W['cHG'] = H['Cgg']*gs2/4.
        
        W['cHW'] = (H['Czz'] + H['Cza']*2.*s2w + H['Caa']*s2w**2)*gw2/4.
        
        W['cHB'] = (H['Czz'] - H['Cza']*2.*c2w + H['Caa']*c2w**2)*gp2/4.

        W['cHWB'] = (H['Czz'] - H['Cza']*(c2w-s2w) - H['Caa']*c2w*s2w)*gw*gp/2.
        
        W['tcHG'] = H['tCgg']*gs2/4.
        
        W['tcHW'] = (H['tCzz'] + H['tCza']*2.*s2w + H['tCaa']*s2w**2)*gw2/4.
        
        W['tcHB'] = (H['tCzz'] - H['tCza']*2.*c2w + H['tCaa']*c2w**2)*gp2/4.

        W['tcHWB'] = (H['tCzz'] - H['tCza']*(c2w-s2w) - H['tCaa']*c2w*s2w)*gw*gp/2.

        W['cW'] = -H['Lz']*2./gw/3. # one or two powers?

        W['cG'] = H['C3G']*gs**3
        
        W['tcW'] = -H['tLz']*2./gw/3.

        W['tcG'] = H['tC3G']*gs**3
        
        C1 = ( H['Czbx']*gw2 + H['Czz']*gp2 - H['Cza']*gp2*(c2w-s2w) 
                - H['Caa']*gw2*s2w**2 )*gw2/(gw2-gp2) 
                
        C2 = ( H['Czbx'] + H['Czz'] - H['Cza']*(c2w-s2w) - H['Caa']*s2w*c2w
                )*gw2*gp2/(gw2-gp2)
        
        # Yukawa
        for f in ('u','d','e'):
            for i,j in H['HBxdY'+f].keys(): 
                
                mi, mj = self.mass[PID[f][i]], self.mass[PID[f][j]]
                
                diag = ( dCz + 2.*dM + C1 )*delta(i,j)
               
                yuk = H['HBxdY'+f][i,j]
                sin = H['HBxS'+f][i,j]
                cos = sqrt(1.-sin**2)
                
                re = -sqrt(2.*mi*mj)/vev*(yuk*cos - diag)
                im = sqrt(2.*mi*mj)/vev*yuk*sin
                W['WBx'+f+'H'][i,j] = complex(re, im)
        
        RdGLzd = matrix_mult(matrix_mult(H.ckm, H['HBxdGLzd']),H.ckm.dag()) 
        
        for i,j in H['HBxdGLzu'].keys():
            
            facp = delta(i,j)*C1
            fac = delta(i,j)*C2
            
            W['WBxH3l'][i,j] = H['HBxdGLwl'][i,j] + facp/2.
            
            W['WBxH3q'][i,j] = H['HBxdGLzu'][i,j] - RdGLzd[i,j] + facp/2.
            
            W['WBxH1l'][i,j] = -2.*H['HBxdGLze'][i,j] - H['HBxdGLwl'][i,j] + fac/2.
            
            W['WBxH1q'][i,j] = -H['HBxdGLzu'][i,j] - RdGLzd[i,j] - fac/6.

            W['WBxHe'][i,j] = - 2.*H['HBxdGRze'][i,j] + fac
            
            W['WBxHu'][i,j] = - 2.*H['HBxdGRzu'][i,j] - 2./3.*fac
            
            W['WBxHd'][i,j] = - 2.*H['HBxdGRzd'][i,j] + 1./3.*fac

        for k,v in H['HBxdGRwq'].iteritems():
            W['WBxHud'][k] = -2.*v
        
        # Dipole operators
        for f in ('u','d','e'):
            eta = 1 if f=='u' else -1
            for i,j in W['WBx'+f+'W'].keys():
                
                mi, mj = self.mass[PID[f][i]], self.mass[PID[f][j]]
                
                MFVnorm = sqrt(mi*mj)/(2.*sqrt(2.)*vev)
                
                # gluon
                if f in ('u','d'):
                    W['WBx'+f+'G'][i,j] = - H['HBxdg'+f][i,j]*MFVnorm
                # Weak
                W['WBx'+f+'W'][i,j] = - eta*(H['HBxdz'+f][i,j] 
                                            + s2w*H['HBxda'+f][i,j])*MFVnorm
                # Hypercharge
                W['WBx'+f+'B'][i,j] = (H['HBxdz'+f][i,j] 
                                      - c2w*H['HBxda'+f][i,j])*MFVnorm
                                      
        # Some four fermion operators
        W['cll1221'] = 2.*( H['HBxdGLwl'][1,1].real + H['HBxdGLwl'][2,2].real) - 4.*dM
        
        trivial_4f =  ['cll1111','cll1122','cll1133','cll1331','cll2332',
                       'cle1111','cle1122','cle2211','cle1133','cle3311',
                       'cee1111','cee1122','cee1133',
                       'cpuu3333']
        
        for c in trivial_4f:
            W[c] = H[c]
        # W['cll1122'], W['cpuu3333'] = H['cll1122'], H['cpuu3333']

        return W
        
        
    @basis.translation('silh')
    def to_silh(self, instance):
        
        H = self
        S = instance
        
        s2w, c2w, ee2, gw2, gp2, MZ, vev, gs2 = self.calculate_inputs()
        sw, cw = sqrt(s2w), sqrt(c2w)
        
        MH = self.mass[25]
        lam = MH**2/(2.*vev**2) # Higgs self-coupling     
        
        
        def delta(i,j):
            return 1. if i==j else 0.
            
        
        S['sg'] = gw2/16.*H['Cgg']

        S['tsg'] = gw2/16.*H['tCgg']

        S['sa'] = gw2/16.*H['Caa']

        S['tsa'] = gw2/16.*H['tCaa']

        S['sW'] = ( ( 2.*H['Czbx']*gw2 + H['Czz']*(gw2+gp2) 
                    - H['Caa']*s2w*gp2 )*gw2/(gw2-gp2)/4. 
                  - H['HBxdGLwl'][2,2].real + 2.*H['dM'] )
        
        S['sB'] = ( (-2.*H['Czbx']*gw2 - H['Czz']*(gw2+gp2)
                    + H['Caa']*s2w*gp2 )*gw2/(gw2-gp2)/4. 
                  + (-H['HBxdGLwl'][2,2].real + 2.*H['HBxdGLze'][1,1].real
                    + 2.*H['dM'] - H['cll1122'])*gw2/gp2 )
        
        S['sH'] = (-2.*H['dCz'] + 3.*H['HBxdGLwl'][1,1].real
                                - 3.*H['HBxdGLwl'][2,2].real )
        
        S['sT'] = ( H['HBxdGLwl'][1,1].real - H['HBxdGLwl'][2,2].real 
                  + 4.*H['HBxdGLze'][1,1].real + 4.*H['dM'] - H['cll1122'] )
        
        S['s6'] = ( H['dL3']/lam - 3.*H['dCz'] + 4.*H['HBxdGLwl'][1,1].real 
                                               - 4.*H['HBxdGLwl'][2,2].real )
        
        S['sHW'] = (-H['Czz'] - 2.*H['Cza']*s2w - H['Caa']*s2w**2 )*gw2/4.

        S['tsHW'] = (-H['tCzz'] - 2.*H['tCza']*s2w - H['tCaa']*s2w**2)*gw2/4.
        
        S['sHB'] = (-H['Czz'] + 2.*H['Cza']*c2w + H['Caa']*s2w*(1.+c2w))*gw2/4.

        S['tsHB'] = (-H['tCzz'] + 2.*H['tCza']*c2w + H['tCaa']*s2w*(1.+c2w))*gw2/4.
        
        
        S['s2W'] = H['HBxdGLwl'][1,1].real/2. + H['HBxdGLwl'][2,2].real/2. - H['dM'] 

        S['s2B'] = ( H['HBxdGLwl'][1,1].real + H['HBxdGLwl'][2,2].real 
                   - 2.*H['dM'] + H['cll1122'] )*gw2/(2.*gp2)

        S['s2G']= H['cpuu3333']*3./gs2 

        S['s3W'] = -H['Lz']/gw2/6.

        S['ts3W'] = -H['tLz']/gw2/6.

        S['s3G'] = H['C3G']*gw2/4.

        S['ts3G'] = H['tC3G']*gw2/4.

        # Yukawa interactions
        for f in ('u','d','e'):
            for i,j in H['HBxdY'+f].keys(): 
                diag = delta(i,j)*( H['dCz'] - H['HBxdGLwl'][1,1].real 
                                  + H['HBxdGLwl'][2,2].real )     
                                             
                yuk = H['HBxdY'+f][i,j]
                sin = H['HBxS'+f][i,j]
                cos = sqrt(1.-sin**2)
                
                re = -yuk*cos + diag 
                im = yuk*sin
                S['SBx'+f][i,j] = complex(re, im)


        # Gauge boson
        RdGLzd = matrix_mult(matrix_mult(H.ckm, H['HBxdGLzd']),H.ckm.dag()) 
        for i, j in H['HBxdGLzu'].keys():
            
            diagp = delta(i,j)*H['HBxdGLwl'][1,1]
            
            S['SBxHpl'][i,j] = H['HBxdGLwl'][i,j] - diagp

            S['SBxHpq'][i,j] = H['HBxdGLzu'][i,j] - RdGLzd[i,j] - diagp 

            diag = delta(i,j)*(H['HBxdGLwl'][1,1] + 2.*H['HBxdGLze'][1,1])

            S['SBxHl'][i,j] = -H['HBxdGLwl'][i,j] - 2.*H['HBxdGLze'][i,j] + diag

            S['SBxHq'][i,j] = -H['HBxdGLzu'][i,j] - RdGLzd[i,j] - diag/3.

            S['SBxHe'][i,j] = - 2.*H['HBxdGRze'][i,j] + 2.*diag

            S['SBxHu'][i,j] = - 2.*H['HBxdGRzu'][i,j] - 4./3.*diag

            S['SBxHd'][i,j] = - 2.*H['HBxdGRzd'][i,j] + 2./3.*diag
            
        for k, v in H['HBxdGRwq'].iteritems():
            S['SBxHud'][k] = -2.*v              
        
        
        # Dipole interactions
        ii = complex(0.,1.)
        for f in ('u','d','e'):
            eta = 1 if f=='u' else -1
            for i,j in S['SBx'+f+'W'].keys():
                # gluon
                if f in ('u','d'):
                    S['SBx'+f+'G'][i,j] = - H['HBxdg'+f][i,j]*gw2/16.
                # Weak
                S['SBx'+f+'W'][i,j] = - eta*(H['HBxdz'+f][i,j] 
                                            + s2w*H['HBxda'+f][i,j])*gw2/16.
                # Hypercharge
                S['SBx'+f+'B'][i,j] = (H['HBxdz'+f][i,j] 
                                      - c2w*H['HBxda'+f][i,j])*gw2/16.
        
        # Some four fermion operators
        trivial_4f = ['cll1111','cll1133','cll1331','cll2332',
                      'cle1111','cle1122','cle2211','cle1133','cle3311',
                      'cee1111','cee1122','cee1133']

        for c in trivial_4f:
            S['s'+c[1:]] = H[c]
        
        return S


    @basis.translation('hisz')
    def to_hisz(self, instance):
        
        s2w, c2w, ee2, gw2, gp2, MZ, vev, gs2 = self.calculate_inputs() 

        dg = gw2 - gp2
        dg_inv = 1./dg
        
        H = self
        Z = instance
        
        Z['Lam'] = vev
        
        Z['fGG'] = -8.*math.pi**2*H['Cgg']
    
        Z['fH2'] = -2.*H['dCz']

        Z['fW'] = -4.*dg_inv*(gw2*H['Czbx'] + gp2*H['Czz'] 
                             - s2w*ee2*H['Caa'] - s2w*dg*H['Cza'])
        
        Z['fB'] = 4.*dg_inv*(gw2*H['Czbx'] + gw2*H['Czz'] 
                            - c2w*ee2*H['Caa'] - c2w*dg*H['Cza'])

        Z['fWW'] = -dg_inv*(2.*gw2*H['Czbx'] + (gw2 + gp2)*H['Czz'] 
                           - s2w*gp2*H['Caa'])
                           
        Z['fBB'] = dg_inv*(2.*gw2*H['Czbx'] + (gw2 + gp2)*H['Czz'] 
                          - c2w*gw2*H['Caa'])
                           
        Z['fWWW'] = 8./(3.*gw2**2)*H['Lz']
    
        Z['tfGG'] = -8.*math.pi**2*H['tCgg']

        Z['tfW'] = 4.*(H['tCzz']- (c2w-s2w)*H['tCza']  - s2w*c2w*H['tCaa'] )

        Z['tfWW'] = (H['tCzz']- 2.*c2w*H['tCza'] - s2w*(c2w + 1.)*H['tCaa'])
                           
        Z['tfBB'] = -H['tCzz']+ 2.*c2w*H['tCza'] - c2w**2*H['tCaa']
        
        Z['tfWWW'] = 8./(3.*gw2**2)*H['tLz']
        
        
        def delta(i,j):
            return 1. if i==j else 0.
        
        # Yukawa interactions
        for f in ('u','d','e'):
            for i,j in H['HBxdY'+f].keys(): 
               
                yuk = H['HBxdY'+f][i,j]
                sin = H['HBxS'+f][i,j]
                cos = sqrt(1.-sin**2)
                
                re = ( delta(i,j)*H['dCz'] - yuk*cos )*sqrt(2.)
                im = yuk*sin*sqrt(2.)
                Z['HZx'+f][i,j] = complex(re, im)
        
        # Provide the right W mass for input...
        MW = MZ*sqrt(c2w)
        Z.mass[24]= MW
        try:
            Z.inputs.new_entry(9, MW, name='MW')
        except KeyError:
            Z.inputs[9] = MW
        
        return Z


################################################################################
    def modify_inputs(self):
        '''
        W mass modification from dM.
        '''
        try:
            self.mass[24] += self['dM']
        except KeyError:
            self.mass.new_entry(24, MW + self['dM'], name = 'MW')
        
################################################################################ 