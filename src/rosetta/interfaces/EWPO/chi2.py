import numpy as np
from scipy.stats import chi2
import os
################################################################################
# read in likelihood data
import importlib.resources as resources

module = 'rosetta.interfaces.EWPO'

with resources.open_text(module, 'likelihood/c0.dat') as dat:
    c0 = np.loadtxt(dat)
with resources.open_text(module, 'likelihood/sigmainv2.dat') as dat:
    sigmainv2 = np.loadtxt(dat)
with resources.open_text(module, 'likelihood/c0_MFV.dat') as dat:
    c0_MFV = np.loadtxt(dat)
with resources.open_text(module, 'likelihood/sigmainv2_MFV.dat') as dat:
    sigmainv2_MFV = np.loadtxt(dat)
################################################################################
def chisq(c):
    deltac = np.array(c) - c0
    return float(np.mat(deltac)*np.mat(sigmainv2)*np.mat(deltac).T)

c = np.full(36,1e-3)
x2 = chisq(c)
print(x2)
thirtysix = chi2(36)
print(thirtysix.sf(x2))
