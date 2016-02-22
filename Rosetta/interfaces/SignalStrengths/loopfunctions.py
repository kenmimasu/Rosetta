from math import asin, sqrt, log, pi

def f(tau):
    if tau <= 1.:
        return complex(asin( sqrt(tau) )**2, 0.)
    else:
        return -complex(log((1. + sqrt(1. - 1./tau))/
                            (1. - sqrt(1. - 1./tau))),
                        -pi)**2/4.
        
def Af(tau):
    '''
    Fermion triangle loop function.
    '''
    return ( f(tau)*(tau - 1.) + tau )/(2.*tau**2)

def AW(tau):
    '''
    W boson triangle loop function.
    '''
    return -( 2.*tau**2 + 3.*tau + 3.*f(tau)*(2.*tau - 1.) )/tau**2