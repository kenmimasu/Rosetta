from internal import Basis
from internal import default_masses, default_inputs, PID
from internal.matrices import matrix_mult, matrix_add, matrix_sub, matrix_eq

################################################################################
class TemplateBasis(Basis.Basis):
    '''
    Toy implementation of a user defined basis class. Makes use of the internal 
    function flavour_matrix to generate indexed coefficients e11,e12,...,e33. A 
    template input card can be fond in the Cards directory or can be generated 
    with write_template_card().
    '''
    
    name = 'template'
    ##########################
    # declare coefficients
    # flavoured = flavmat('e', kind='symmetric', domain='real')
    ########################## 
    independent = ['a','b','c','AA']
    
    required_masses = {1,2,3,4,5,6} # quark masses rquired
    
    required_inputs = {1, 4} # aEWM1, MZ
    
    blocks = {'newcoup':['a','b','c'],
              'dep':['d']}
    
    flavoured = {'AA':{'domain':'complex', 'kind':'general', 'cname':'cAA'},
                 'BB':{'domain':'complex', 'kind':'general', 'cname':'cBB'},
                 'CC':{'domain':'complex', 'kind':'general', 'cname':'cCC'}}
              
    ##########################
    # This can be overridden if you want but is typically unneccesary. Exclude 
    # it entirely and the class will inherit the  constructor from the parent 
    # class, Basis.Basis.
    def __init__(self,*args,**kwargs): 
        # If you want to add additional behaviour by redefining the 
        # constructor, ensure a call to the base class constructor is made 
        # like so:
        super(TemplateBasis, self).__init__(*args,**kwargs) 

        # additional instructions can live here
        
    def calculate_dependent(self):
        '''
        Calculate dependent parameters here by assigning values to parameters 
        not listed in self.independent.
        '''
        p = self
        p['d'] = p['a']+ p['b']*p['c']
        # element-wise asssignment for matrix BB: BB[i,j] = -A[i,j]/4.
        for k, v in p['AA'].iteritems():
            p['BB'][k] = -v/4.
        # matrix multiplication function assigning CC -> AA.BB
        matrix_mult(p['AA'], p['BB'], p['CC'])
        
    @Basis.translation('silh')
    def to_silh(self, instance):
        '''
        Toy translation to the SILH basis setting all coefficients to 1e-3, and 
        modifying the Z and Higgs masses.
        '''
        A = self
        B = instance
        # set all values of silh basis coeffs to 0.01:
        # coeff_d*m_top*a_EW
        for k in B: 
            B[k] = 0.001
        self.mass[23]=91.19 # MZ in newmass
        self.inputs[8]=126. # MH in newinput
        
        return B
        
    @Basis.translation('mass')        
    def to_mass(self, instance):
        '''
        Toy translation to the Mass basis setting all coefficients according to 
        a nonsense formula, myfunc, defined below. Also modifies the Z and 
        Higgs masses.
        '''
        A = self
        B = instance
        
        # set all values of mass basis coeffs according to nonsense formula 
        # y = d*m_top*a_EW
        for k in B: 
            B[k] = self.myfunc( A['d'], self.mass[6], self.inputs['aEWM1'] )
        self.mass[24]=91.19 # MZ in newmass
        self.inputs[8]=126. # MH in newinput
        
        return B
        
    @staticmethod # pure function, no access to class instance
    def myfunc(x,y,z):
        return x*y/z
    
    def myfunc2(self,x,y,z): # access to class instance, "self"
        my_one = self.inputs[4]/self.inputs[4]
        return x*y/z*my_one

################################################################################
if __name__=='__main__':
    # Testing area which can be executed when running `python TemplateBasis.py` 
    instance = TemplateBasis(param_card='../Cards/param_card_TemplateBasis.dat')
    instance.write_param_card('test_template.dat')
    # inputs = instance.eHDECAY_inputs()
    # print eHDECAY(inputs)
    