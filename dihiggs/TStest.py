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
parser = argparse.ArgumentParser(prog='TStest', description='Return the closest shape benchmark topoligy of a EFT point')

# arguments to chose the (B)SM point in training and application
parser.add_argument('--kl', type=float,    default=1.0, help='Benchmark to calculate the limit')
parser.add_argument('--kt', type=float,    default=1.0, help='Benchmark to calculate the limit')
parser.add_argument('--c2', type=float,    default=0.0, help='Benchmark to calculate the limit')
parser.add_argument('--cg', type=float,    default=0.0, help='Benchmark to calculate the limit')
parser.add_argument('--c2g', type=float,    default=0.0, help='Benchmark to calculate the limit')
args = parser.parse_args()

ar = reweighter_from_histogram_and_file()
BM = ar.TS_test(args.kl, args.kt, args.c2, args.cg, args.c2g)
