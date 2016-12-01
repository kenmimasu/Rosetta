from math import sqrt, pi


from ...internal.basis import checkers as check
from .production_xs import xsfb

masses = {25,6} # H, b, t masses


def get_xs_and_br(basis, sqrts = 13):
    from ..SignalStrengths.decay import decay

    decays = decay(basis)
    
    kl, kt, c2, cg, c2g = get_dihiggs_params(basis)
    
    xs, err = xsfb(sqrts, kl, kt, c2, cg, c2g, error=True)
    
    return xs, err, decays


def get_dihiggs_params(basis):
    
    bsmc = basis.translate(target='bsmc')
    
    check.masses(bsmc, masses, message='DiHiggs interface')
     
    MH = bsmc.mass[25] # Higgs mass
    vev = sqrt(1./sqrt(2)/bsmc.inputs['Gf']) # Higgs vev
    lam = MH**2/(2.*vev**2) # Higgs self-coupling
        
    kl = (1. + bsmc['dL3']/lam)
    
    kt = 1. + bsmc['BCxdYu'][3,3]
    
    c2 = bsmc['BCxY2u'][3,3].real/2.
    
    cg = bsmc['Cgg']*12.*pi**2
    
    c2g = -bsmc['Cgg2']*12.*pi**2
    
    return (kl, kt, c2, cg, c2g)