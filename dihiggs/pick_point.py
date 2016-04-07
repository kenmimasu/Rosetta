#!/usr/bin/env python
# RUN: python pick_point.py -n 12 --LHC 13 --kl 1 --kt 1 --c2 0 --cg 0 --c2g 0
from array import array
import math
import bisect
from optparse import OptionParser
import numpy as np
parser = OptionParser()

'''A simple function returning a value, used in an expression'''

parser.add_option("-n", type="int", dest="num", help="Number of clusters")
parser.add_option("--LHC", type="int", dest="lhc", help="LHC CM energy in TeV")
parser.add_option("--kl", type="float", dest="kll", help="LHC CM energy in TeV")
parser.add_option("--kt", type="float", dest="ktt", help="LHC CM energy in TeV")
parser.add_option("--c2", type="float", dest="c22", help="LHC CM energy in TeV")
parser.add_option("--cg", type="float", dest="cgg", help="LHC CM energy in TeV")
parser.add_option("--c2g", type="float", dest="c2gg", help="LHC CM energy in TeV")

(options, args) = parser.parse_args()
print "Let's talk about %s clusters" % options.num
print "LHC @ %s TeV" % options.lhc
print " "

number = options.num
CM =options.lhc 
kl = options.kll
kt = options.ktt
c2 = options.c22
cg = options.cgg
c2g = options.c2gg

# print "RHH = %f" % f(kl,kt,c2,cg,c2g)

########################
# 
filne = "Translation_1507points/clustering_nev20k_Nclu"+str(number)+"_50_5.asc"
f = open(filne, 'r+')
lines = f.readlines() # get all lines as a list (array)
head = []
clusters = [[]]
nclu=0
# Iterate over each line, printing each line and then move to the next
for line in lines:
    #print line
    clusters.append([])
    l = []
    counter =0
    tokens = line.split()
    for token in tokens:
        num = ""
        for char in token: 
            if char.isdigit() or (char in num_char): num = num + char
        try: 
           l.append(int(num))
           clusters[nclu].append(int(num))
        except ValueError: pass
        #print(l[counter])
        if counter==0 : head.append(int(num))
        counter += 1
    nclu += 1
f.close()
#print "heads"
#for x in range(0, 12): print(head[x])
#for x in range(0, len(clusters[0]) ) : print(clusters[0][x])
#
########################

########################
#
fipoints = "list_all_translation_1507.txt"
fp = open(fipoints, 'r+')
Vkl = []
Vkt = []
Vc2 = []
Vcg = []
Vc2g = []
linesp = fp.readlines() # get all lines as a list (array)
for linep in linesp:
    #print linep
    ll = []
    counterp =0
    tokens = linep.split()
    for token in tokens:
        num = ""
        num_char="."
        num_char2="-"
        for char in token: 
            if char.isdigit() or (char in num_char) or (char in num_char2): num = num + char
        try: ll.append(float(num))
        except ValueError: pass
        #print(ll[counterp])
        if counterp == 0 : Vkl.append(float(num))
        if counterp == 1 : Vkt.append(float(num))
        if counterp == 2 : Vc2.append(float(num))
        if counterp == 3 : Vcg.append(float(num))
        if counterp == 4 : Vc2g.append(float(num))
        counterp += 1
f.close()
#for x in range(0, 1507): print(Vkl[x], Vkt[x],Vc2[x],Vcg[x],Vc2g[x])
#
########################

########################
#
# just to check XS range
XS = [[]]
for x in range(0, 12 ): XS.append([])

A7tev = [2.20968, 9.82091, 0.332842, 0.120743, 1.13516, -8.76709, -1.54253, 3.09385, 1.64789, -5.14831, -0.790689, 2.12522, 0.385807, -0.952469, -0.618337]
A8tev = [2.17938, 9.88152, 0.31969, 0.115609, 1.16772, -8.69692, -1.49906, 3.02278, 1.59905, -5.09201, -0.761032, 2.06131, 0.369, -0.922398, -0.604222]
A13tev = [2.09078, 10.1517, 0.282307, 0.101205, 1.33191, -8.51168, -1.37309, 2.82636, 1.45767, -4.91761, -0.675197, 1.86189, 0.321422, -0.836276, -0.568156]
A14tev = [2.07992, 10.2036, 0.277868, 0.0995436, 1.36558, -8.492, -1.35778, 2.80127, 1.44117, -4.89626, -0.664721, 1.83596, 0.315808, -0.826019, -0.564388]
A100tev = [2.17938, 9.88152, 0.31969, 0.115609, 1.16772, -8.69692, -1.49906, 3.02278, 1.59905, -5.09201, -0.761032, 2.06131, 0.369, -0.922398, -0.604222]

if CM == 7 : A = A7tev 
elif CM == 8 : A = A8tev
elif CM == 13 : A = A13tev
elif CM == 14 : A = A14tev
elif CM == 100 : A = A100tev
else : print ("invalid LHC energy")

def f(kl,kt,c2,cg,c2g):
    return A[0]*kt**4 + A[1]*c2**2 + (A[2]*kt**2 + A[3]*cg**2)*kl**2  + A[4]*c2g**2 + ( A[5]*c2 + A[6]*kt*kl )*kt**2  + (A[7]*kt*kl + A[0]*cg*kl )*c2 + A[9]*c2*c2g  + (A[10]*cg*kl + A[11]*c2g)*kt**2+ (A[12]*kl*cg + A[13]*c2g )*kt*kl + A[14]*cg*c2g*kl


print "RHH = %f" % f(kl,kt,c2,cg,c2g)

########################
# errors in normalization, NNLL
# https://twiki.cern.ch/twiki/bin/view/LHCPhysics/LHCHXSWGHH#Current_recommendations_for_di_H
# mh 125 GeV 

if CM == 7 : 
    xs = 7.718
    scalep= 4.0
    scalem = -5.7
    PDF = 3.4
    alphas = 2.8
elif CM == 8 : 
    xs = 11.18
    scalep= 4.1 
    scalem = -5.7
    PDF = 3.1
    alphas = 2.6
elif CM == 13 : 
    xs = 37.95
    scalep= 4.3 
    scalem = -6.0
    PDF = 2.1
    alphas = 2.3
elif CM == 14 : 
    xs = 45.05
    scalep= 4.4 
    scalem = -6.0
    PDF = 2.1
    alphas = 2.2
elif CM ==100 : 
    xs = 1749
    scalep = 5.1 
    scalem = -6.6
    PDF = 1.7
    alphas = 2.1
else : print ("invalid LHC energy")
#
#######################


########################
#
# make the distances
distance =[]
for x in range(0, 1507): 
   distance.append(math.sqrt( ((Vkl[x] - kl)/2.5)**2 + ((Vkt[x] - kt)/0.5)**2 + (Vc2[x] - c2)**2 + ((Vcg[x] - cg)/0.1)**2 +((Vc2g[x] - c2g)/0.1)**2 ))
min_index = np.amin(distance)
min_value = np.argmin(distance)
list4min = np.argsort(distance)[:5]
#print (min_index,min_value)
#print("Center again: ", Vkl[list4min[0]], Vkt[list4min[0]],Vc2[list4min[0]],Vcg[list4min[0]],Vc2g[list4min[0]])

#
# check the cluster
#

######################
# enter the range
counterf=0
for x in range(0, len(clusters) ) :
    for y in range(0, len(clusters[x]) ) : 
       if clusters[x][y] == min_value : 
           print("Center: ",Vkl[min_value], Vkt[min_value],Vc2[min_value],Vcg[min_value],Vc2g[min_value], "Cluster ", x+1, " sample ", clusters[x][y], " distance: ", distance[min_value])
           print " "
           counterf+=1
       if clusters[x][y] == list4min[1] : 
           print("Neigbours: ", Vkl[list4min[1]], Vkt[list4min[1]],Vc2[list4min[1]],Vcg[list4min[1]],Vc2g[list4min[1]], "Cluster ", x+1, " sample ", clusters[x][y], " distance: ", distance[list4min[1]])
           print " "
       if clusters[x][y] == list4min[2] : 
           print("Neigbours: ", Vkl[list4min[2]], Vkt[list4min[2]],Vc2[list4min[2]],Vcg[list4min[2]],Vc2g[list4min[2]], "Cluster ", x+1, " sample ", clusters[x][y], " distance: ", distance[list4min[2]])
           print " "
       if clusters[x][y] == list4min[3] : 
           print("Neigbours: ", Vkl[list4min[3]], Vkt[list4min[3]],Vc2[list4min[3]],Vcg[list4min[3]],Vc2g[list4min[3]], "Cluster ", x+1, " sample ", clusters[x][y], " distance: ", distance[list4min[3]])
           print " "
       if clusters[x][y] == list4min[4] : 
           print("Neigbours: ", Vkl[list4min[4]], Vkt[list4min[4]],Vc2[list4min[4]],Vcg[list4min[4]],Vc2g[list4min[4]], "Cluster ", x+1, " sample ", clusters[x][y], " distance: ", distance[list4min[4]])
           print " "
       XS[x].append(f(Vkl[clusters[x][y]], Vkt[clusters[x][y]],Vc2[clusters[x][y]],Vcg[clusters[x][y]],Vc2g[clusters[x][y]])*xs)

#print XS[11]
print "SigmaHH = %f" % (f(kl,kt,c2,cg,c2g)*xs)

       #else :  continue
       # check neighbors
#for x in range(0, len(clusters) ) :
#    for y in range(0, len(clusters[x]) ) : 
#       if Vkl[clusters[x][y]] == Vkl[min_value] + 5 and Vkt[clusters[x][y]] ==kt and Vc2[clusters[x][y]] ==c2 and Vcg[clusters[x][y]] ==cg and Vc2g[clusters[x][y]] == c2g and counterf==1 :
#          print ("one kl more, Cluster:", clusters[x][y] )
#          print("one kl more: ",Vkl[clusters[x][y]], Vkt[clusters[x][y]],Vc2[clusters[x][y]],Vcg[clusters[x][y]],Vc2g[clusters[x][y]], "Cluster ",x+1)
#       #
#       if Vkl[clusters[x][y]] == Vkl[min_value] - 5 and Vkt[clusters[x][y]] ==kt and Vc2[clusters[x][y]] ==c2 and Vcg[clusters[x][y]] ==cg and Vc2g[clusters[x][y]] == c2g and counterf==1 :
#          print ("one kl less, Cluster:", clusters[x][y] )
#          print("one kl less: ",Vkl[clusters[x][y]], Vkt[clusters[x][y]],Vc2[clusters[x][y]],Vcg[clusters[x][y]],Vc2g[clusters[x][y]], "Cluster ",x+1)
#print("one less: ",Vkl[min_value], Vkt[min_value],Vc2[min_value],Vcg[min_value],Vc2g[min_value])
#if counterf ==1 :  continue

