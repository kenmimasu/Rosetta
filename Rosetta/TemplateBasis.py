from Basis import Basis
from MassBasis import MassBasis
####################################################################################################
# Template basis class
class TemplateBasis(Basis):
    independent = ['a','b','c']
    dependent=['d']
    required_masses = {1,2,3,4,5,6}
    required_inputs = {1, 4} # aEWM1, MZ
    
    def __init__(self,*args,**kwargs): # This can be overridden if you want
        super(TemplateBasis, self).__init__(*args,**kwargs) # ensure a call to the base class constructor is made as so
        # additional instructions can live here
        
    def calculate_dependent(self):
        '''
        Calculate dependent parameters here by assigning values to self.par_dict 
        corresponding to the keys in self.dependent.
        '''
        p = self.par_dict
        p['d'] = p['a']+ p['b']*p['c']
        
    def translate(self):
        '''
        Translate to the mass basis by creating an empty MassBasis and modifying its par_dict or coeffs._asdict()
        set self.newpar with resulting par_dict
        '''
        A = self.coeffs._asdict()
        B = MassBasis().coeffs._asdict()
        for k in B.keys(): # set all values of mass basis coeffs according to nonsense formula coeff_a*m_top/a_EW
            B[k] = self.myfunc( A['d'], self.mass[6], self.input['aEWM1'])
        self.newmass[23]=91.19 # MZ in newmass
        self.newinput[8]=126. # MH in newinput
        self.newpar = B
        self.newname = 'Mass'
        
    @staticmethod # pure function, no access to class instance
    def myfunc(x,y,z):
        return x*y/z
    
    def myfunc2(self,x,y,z): # access to class instance
        my_one = self.SLHA_sminputs[4]/self.SLHA_sminputs[4]
        return x*y/z*my_one

        
####################################################################################################
