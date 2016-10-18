from errors import SqrtsError

# SM predictions for di-higgs production with various relative uncertainties at
# different pp collider energies
# errors in normalization, NNLL for mh = 125 GeV
# https://twiki.cern.ch/twiki/bin/view/LHCPhysics/LHCHXSWGHH#Current_recommendations_for_di_H
SM_predictions = {
    7  :{ 'xs':7.718, 'scalep':4.0, 'scalem' :-5.7, 'PDF':3.4, 'alphas':2.8},
    8  :{ 'xs':11.18, 'scalep':4.1, 'scalem' :-5.7, 'PDF':3.1, 'alphas':2.6},
    13 :{ 'xs':37.95, 'scalep':4.3, 'scalem' :-6.0, 'PDF':2.1, 'alphas':2.3},
    14 :{ 'xs':45.05, 'scalep':4.4, 'scalem' :-6.0, 'PDF':2.1, 'alphas':2.2},
    100:{ 'xs':1749., 'scalep':5.1, 'scalem' :-6.6, 'PDF':1.7, 'alphas':2.1}
}

# Polynomial coefficients for ratio at different pp collider energies
pp_coeffs = {
    7  :[2.20968, 9.82091, 0.332842, 0.120743, 1.13516, 
         -8.76709, -1.54253, 3.09385, 1.64789, -5.14831, 
         -0.790689, 2.12522, 0.385807, -0.952469, -0.618337],
    8  :[2.17938, 9.88152, 0.31969, 0.115609, 1.16772, 
         -8.69692, -1.49906, 3.02278, 1.59905, -5.09201, 
         -0.761032, 2.06131, 0.369, -0.922398, -0.604222],
    13 :[2.09078, 10.1517, 0.282307, 0.101205, 1.33191, 
         -8.51168, -1.37309, 2.82636, 1.45767, -4.91761, 
         -0.675197, 1.86189, 0.321422, -0.836276, -0.568156],
    14 :[2.07992, 10.2036, 0.277868, 0.0995436, 1.36558, 
         -8.492, -1.35778, 2.80127, 1.44117, -4.89626, 
         -0.664721, 1.83596, 0.315808, -0.826019, -0.564388],
    100:[2.17938, 9.88152, 0.31969, 0.115609, 1.16772, 
         -8.69692, -1.49906, 3.02278, 1.59905, -5.09201, 
         -0.761032, 2.06131, 0.369, -0.922398, -0.604222]
}

def xsfb(E, kl, kt, c2, cg, c2g, error=False):
    try:
        A = pp_coeffs[E]
        SM = SM_predictions[E]
    except KeyError:
        raise SqrtsError(sqrt(sqrts))
    
    R_SM = ( A[0]*kt**4 + A[1]*c2**2 + (A[2]*kt**2 + A[3]*cg**2)*kl**2 
           + A[4]*c2g**2 + ( A[5]*c2 + A[6]*kt*kl )*kt**2 
           + ( A[7]*kt*kl + A[8]*cg*kl )*c2 + A[9]*c2*c2g  
           + ( A[10]*cg*kl + A[11]*c2g )*kt**2 
           + ( A[12]*kl*cg + A[13]*c2g )*kt*kl + A[14]*cg*c2g*kl)
    
    if not error:
        return SM['xs']*R_SM
    else:
        return SM['xs']*R_SM, {k:v*R_SM for k,v in SM.iteritems() if k!='xs'}

# print "RHH = %f" % f(kl,kt,c2,cg,c2g)

########################

#print(f(1,1,0,0,0)*xs)
# print "XS =  %f fb, scale (+ %f, %f) fb, PDF +- %f fb, alphaS +- %f fb, top mass +- %f fb, RHH +- 0 " % (f(kl,kt,c2,cg,c2g)*xs, f(kl,kt,c2,cg,c2g)*scalep, f(kl,kt,c2,cg,c2g)*scalem, f(kl,kt,c2,cg,c2g)*PDF, f(kl,kt,c2,cg,c2g)*alphas, f(kl,kt,c2,cg,c2g)*0.1 )





