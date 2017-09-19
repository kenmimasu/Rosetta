from math import sqrt
from SLHA import NamedMatrix, CNamedMatrix
from matrices import CTwoDMatrix
from itertools import product
################################################################################
__doc__ = '''
Some useful lookup dictionaries and default values for particle masses, SM  
inputs and the CKM matrix.
'''
# PDG ID dictionary
PID = {'u':{1:1 ,2:4, 3:6},'d':{1:2, 2:3, 3:5},
       'e':{1:11, 2:13, 3:15},'v':{1:12, 2:14, 3:16}}
       
# PID:name dictionary for particles
particle_names = {1:'u', 2:'d', 3:'s', 4:'c', 5:'b', 6:'t', 
                  11:'e', 12:'ve', 13:'mu', 14:'vmu', 15:'ta', 16:'vta', 
                  21:'g', 22:'a', 23:'Z', 24:'W', 25:'H'}
                  
# ID:name dictionary for SLHA inputs
input_names = {1:'aEWM1', 2:'Gf', 3:'aS', 4:'MZ', 
               5:'MB', 6:'MT', 7:'MTAU', 
               9:'MW', 25:'MH', 10:'vF'}
               
# ID:PID dictionary for SLHA inputs that are particle masses
input_to_PID = {4:23, 5:5, 6:6, 7:15, 25:25, 9:24}
PID_to_input = {v:k for k,v in input_to_PID.items()}

# PID:value dictionary for default particle masses when undefined
default_masses = {1:0.0048, 2:0.0023, 3:0.095, 4:1.42, 5:4.7, 6:173., 
                  11:0.000511, 12:0., 13:0.105658367, 14:0., 15:1.7768, 16:0.,
                  23:9.118760e+01, 24:8.0385e+01, 25:125.}
                  
# ID:value dictionary for default SHLA inputs when undefined
default_inputs = {1: 1.27916e+02, 2: 1.166379e-05, 3: 1.184e-01, 
                  4: default_masses[23], 5: default_masses[5], 
                  6: default_masses[6], 7: default_masses[15],
                  25:default_masses[25], 9:default_masses[24],
                  10:2.462200e+02 }

# Z, W & top decay widths for default decay blocks
GammaZ, GammaW, Gammat= 2.4952, 2.085, 1.4915

# Construct default CKM matrix
# Wolfenstein parameterisation (PDG best fit)
lam, A, rho, eta = 0.22535, 0.811, 0.131, 0.345
# CKM mixing angles and CP phase
s12, s23 = lam, A*lam**2
s13cd, s13sd = A*lam**3*rho, A*lam**3*eta
s13 = sqrt(s13sd**2+s13cd**2)
c12, c23, c13 = map(lambda x: sqrt(1.-x**2), (s12, s23, s13 ))
# real and imaginary parts
VCKM = {1:{1:c12*c13,                2:s12*c13,                3:s13cd  },
        2:{1:-s12*c23-c12*s23*s13cd, 2:c12*c23-s12*s23*s13cd,  3:s23*c13},
        3:{1:s12*s23-c12*c23*s13cd,  2:-c12*s23-s12*c23*s13cd, 3:c23*c13}}
IMVCKM = {1:{1:0.,             2:0.,             3:-s13sd},
          2:{1:-c12*s23*s13sd, 2:-s12*s23*s13sd, 3:0.   },
          3:{1:-c12*c23*s13sd, 2:-s12*c23*s13sd, 3:0.   }}
VCKMele = {(1,1):'ud',(1,2):'us',(1,3):'ub',(2,1):'cd',
           (2,2):'cs',(2,3):'cb',(3,1):'td',(3,2):'ts',(3,3):'tb'}

preamble = ('\n###################################\n'
        + '## CKM INFORMATION\n'
        + '###################################\n')
        
ckm = NamedMatrix(name='VCKM', preamble=preamble)
for i,j in product((1,2,3),(1,2,3)):
    cname = 'RV{}{}x{}'.format(VCKMele[(i,j)], i, j)
    ckm.new_entry((i,j), VCKM[i][j], name=cname)
    
imckm = NamedMatrix(name='IMVCKM')
for i,j in product((1,2,3),(1,2,3)):
    cname = 'IV{}{}x{}'.format(VCKMele[(i,j)], i, j)
    imckm.new_entry((i,j), IMVCKM[i][j], name=cname) 
    
# vckm = SLHA.CNamedMatrix(ckm, imckm)
default_ckm = CTwoDMatrix(CNamedMatrix(ckm, imckm))

################################################################################
           