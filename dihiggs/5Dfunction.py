#!/usr/bin/env python
# RUN: python 5Dfunction.py -n 7 --kl 1 --kt 1 --c2 0 --cg 0 --c2g 0
from array import array
from optparse import OptionParser
parser = OptionParser()

'''A simple function returning a value, used in an expression'''

parser.add_option("-n", type="int", dest="num", help="LHC CM energy in TeV")
parser.add_option("--kl", type="float", dest="kll", help="LHC CM energy in TeV")
parser.add_option("--kt", type="float", dest="ktt", help="LHC CM energy in TeV")
parser.add_option("--c2", type="float", dest="c22", help="LHC CM energy in TeV")
parser.add_option("--cg", type="float", dest="cgg", help="LHC CM energy in TeV")
parser.add_option("--c2g", type="float", dest="c2gg", help="LHC CM energy in TeV")

(options, args) = parser.parse_args()
print "Let's talk about LHC at %s TeV" % options.num
number = options.num
kl = options.kll
kt = options.ktt
c2 = options.c22
cg = options.cgg
c2g = options.c2gg

A7tev = [2.20968, 9.82091, 0.332842, 0.120743, 1.13516, -8.76709, -1.54253, 3.09385, 1.64789, -5.14831, -0.790689, 2.12522, 0.385807, -0.952469, -0.618337]
A8tev = [2.17938, 9.88152, 0.31969, 0.115609, 1.16772, -8.69692, -1.49906, 3.02278, 1.59905, -5.09201, -0.761032, 2.06131, 0.369, -0.922398, -0.604222]
A13tev = [2.09078, 10.1517, 0.282307, 0.101205, 1.33191, -8.51168, -1.37309, 2.82636, 1.45767, -4.91761, -0.675197, 1.86189, 0.321422, -0.836276, -0.568156]
A14tev = [2.07992, 10.2036, 0.277868, 0.0995436, 1.36558, -8.492, -1.35778, 2.80127, 1.44117, -4.89626, -0.664721, 1.83596, 0.315808, -0.826019, -0.564388]
A100tev = [2.17938, 9.88152, 0.31969, 0.115609, 1.16772, -8.69692, -1.49906, 3.02278, 1.59905, -5.09201, -0.761032, 2.06131, 0.369, -0.922398, -0.604222]

if number == 7 : A = A7tev 
elif number == 8 : A = A8tev
elif number == 13 : A = A13tev
elif number == 14 : A = A14tev
elif number == 100 : A = A100tev
else : print ("invalid LHC energy")

def f(kl,kt,c2,cg,c2g):
    return A[0]*kt**4 + A[1]*c2**2 + (A[2]*kt**2 + A[3]*cg**2)*kl**2  + A[4]*c2g**2 + ( A[5]*c2 + A[6]*kt*kl )*kt**2  + (A[7]*kt*kl + A[0]*cg*kl )*c2 + A[9]*c2*c2g  + (A[10]*cg*kl + A[11]*c2g)*kt**2+ (A[12]*kl*cg + A[13]*c2g )*kt*kl + A[14]*cg*c2g*kl


print "RHH = %f" % f(kl,kt,c2,cg,c2g)

########################
# errors in normalization, NNLL
# https://twiki.cern.ch/twiki/bin/view/LHCPhysics/LHCHXSWGHH#Current_recommendations_for_di_H
# mh 125 GeV 

if number == 7 : 
   xs = 7.718
   scalep= 4.0
   scalem = -5.7
   PDF = 3.4
   alphas = 2.8
elif number == 8 : 
   xs = 11.18
   scalep= 4.1 
   scalem = -5.7
   PDF = 3.1
   alphas = 2.6
elif number == 13 : 
   xs = 37.95
   scalep= 4.3 
   scalem = -6.0
   PDF = 2.1
   alphas = 2.3
elif number == 14 : 
   xs = 45.05
   scalep= 4.4 
   scalem = -6.0
   PDF = 2.1
   alphas = 2.2
elif number ==100 : 
   xs = 1749
   scalep = 5.1 
   scalem = -6.6
   PDF = 1.7
   alphas = 2.1
else : print ("invalid LHC energy")

#print(f(1,1,0,0,0)*xs)
# by now we ignore errors in the RHH errors
print "XS =  %f fb, scale (+ %f, %f) fb, PDF +- %f fb, alphaS +- %f fb, top mass +- %f fb, RHH +- 0 " % (f(kl,kt,c2,cg,c2g)*xs, f(kl,kt,c2,cg,c2g)*scalep, f(kl,kt,c2,cg,c2g)*scalem, f(kl,kt,c2,cg,c2g)*PDF, f(kl,kt,c2,cg,c2g)*alphas, f(kl,kt,c2,cg,c2g)*0.1 )

# Have the XS in the H basis - with error expansion




