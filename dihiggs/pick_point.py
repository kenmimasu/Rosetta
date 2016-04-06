#!/usr/bin/env python
# RUN: python pick_point.py -n 7 --kl 1 --kt 1 --c2 0 --cg 0 --c2g 0
from array import array
import math
import bisect
from optparse import OptionParser
import numpy as np
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

# print "RHH = %f" % f(kl,kt,c2,cg,c2g)

########################
# 
filne = "Translation_1507points/clustering_nev20k_Nclu12_50_5.asc"
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
# make the distances
distance =[]
for x in range(0, 1507): 
   distance.append(math.sqrt( (Vkl[x] - kl)**2 + (Vkt[x] - kt)**2 + (Vc2[x] - c2)**2 + (Vcg[x] - cg)**2 +(Vc2g[x] - c2g)**2 ))
min_index = np.amin(distance)
min_value = np.argmin(distance)
print (min_index,min_value)
print(Vkl[min_value], Vkt[min_value],Vc2[min_value],Vcg[min_value],Vc2g[min_value])
#
# check the cluster
#
for x in range(0, len(clusters) ) :
    for y in range(0, len(clusters[x]) ) : 
       if clusters[x][y] == min_value : print ("Cluster ", x+1)



