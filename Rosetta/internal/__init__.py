################################################################################
# PDG ID dictionary
PID = {'u':{1:1 ,2:4, 3:6},'d':{1:2, 2:3, 3:5},
       'e':{1:11, 2:13, 3:15},'v':{1:12, 2:14, 3:16}}
# PID:name dictionary for particles
particle_names = {1:'u', 2:'d', 3:'s', 4:'c', 5:'b', 6:'t', 
                  11:'e', 12:'ve', 13:'mu', 14:'vmu', 15:'ta', 16:'vta', 
                  21:'a', 22:'g', 23:'Z', 24:'W', 25:'H'}
# ID:name dictionary for SLHA inputs
input_names = {1:'aEWM1', 2:'Gf', 3:'aS', 4:'MZ', 
               5:'MB', 6:'MT', 7:'MTAU', 25:'MH'}
# ID:PID dictionary for SLHA inputs that are particle masses
input_to_PID = {4:23, 5:5, 6:6, 7:15, 25:25}
PID_to_input = {v:k for k,v in input_to_PID.items()}
# PID:value dictionary for default particle masses when undefined
default_masses = {1:0.0048, 2:0.0023, 3:0.095, 4:1.42, 5:4.7, 6:173., 
                  11:0.000511, 12:0., 13:0.105658367, 14:0., 15:1.7768, 16:0.,
                  23:9.118800e+01, 24:7.982400e+01, 25:125.}
# ID:value dictionary for default SHLA inputs when undefined
default_inputs = {1: 1.325070e+02, 2: 1.166390e-05, 3: 1.180000e-01, 
                  4: default_masses[23], 5: default_masses[5], 
                  6: default_masses[6], 7: default_masses[15],
                  25:default_masses[25], 9:default_masses[24] }
################################################################################
class RosettaError(Exception):
    '''Fancy error name.'''
    pass
################################################################################
