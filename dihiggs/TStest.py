#!/usr/bin/env python

import numpy as np
import os

from python.functions import reweighter_from_histogram_and_file
from python.functions import AnalyticalReweighter

#os.system("python python/functions.py")
#scriptpath = "python/functions.py"
# Add the directory containing your module to the Python path (wants absolute paths)
#import sys
#sys.path.append(os.path.abspath("python/functions.py"))
#import AnalyticalReweighter

import argparse
parser = argparse.ArgumentParser(prog='TStest', description='Return the closest shape benchmark topology of a EFT point')

# arguments to chose the (B)SM point in training and application
parser.add_argument('--kl', type=float,    default=1.0, help='Benchmark to calculate the limit')
parser.add_argument('--kt', type=float,    default=1.0, help='Benchmark to calculate the limit')
parser.add_argument('--c2', type=float,    default=0.0, help='Benchmark to calculate the limit')
parser.add_argument('--cg', type=float,    default=0.0, help='Benchmark to calculate the limit')
parser.add_argument('--c2g', type=float,    default=0.0, help='Benchmark to calculate the limit')
args = parser.parse_args()

ar = reweighter_from_histogram_and_file()

# First test example: find which benchmark point is the closest to the point passed as argument
BM, TS = ar.TS_test(args.kl, args.kt, args.c2, args.cg, args.c2g)
print("closest BM is # {} with TS {}".format(BM, TS))

# Second test example: find point in the kl scan that is the closest to each BM point
for ibm, bm in enumerate(ar.JHEP_BM):
    bm_kl, bm_kt, bm_c2, bm_cg, bm_c2g = bm
    kl_scan = []
    for ikl in xrange(-200, 201, 1):
        kl = ikl / 10.
        kt = 1.
        kl_scan.append((kl, kt, 0., 0., 0.))
    BM, TS = ar.find_closest_points(bm_kl, bm_kt, bm_c2, bm_cg, bm_c2g, kl_scan)
    coordinates = [kl_scan[x] for x in BM]
    print("closest to BM # {} {} are points {} with TS {}".format(ibm, bm, coordinates, TS[0]))
