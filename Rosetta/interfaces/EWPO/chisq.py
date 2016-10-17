import numpy as np
from scipy.stats import chi2
import os
################################################################################
# read in likelihood data
data_dir = os.path.dirname(__file__) + '/likelihood/'
with open(data_dir+'c0.dat') as dat:
    c0 = np.loadtxt(dat)
with open(data_dir+'sigmainv2.dat') as dat:
    sigmainv2 = np.loadtxt(dat)
with open(data_dir+'c0_MFV.dat') as dat:
    c0_MFV = np.loadtxt(dat)
with open(data_dir+'sigmainv2_MFV.dat') as dat:
    sigmainv2_MFV = np.loadtxt(dat)

chisq_36 = chi2(36)
chisq_23 = chi2(23)
################################################################################
def chisq(c, flavor='general'):
    if flavor=='universal':
        deltac = np.array(c) - c0_MFV
        sinv2 = sigmainv2_MFV
    else:
        deltac = np.array(c) - c0
        sinv2 = sigmainv2
    return float(np.mat(deltac)*np.mat(sinv2)*np.mat(deltac).T)

def pvalue(csq, flavor='general'):
    if flavor=='universal':
        return chisq_23.sf(csq)
    else:
        return chisq_36.sf(csq)
        
    
def chisq_and_pvalue(c, flavor='general'):
    chisq_val = chisq(c, flavor=flavor)
    return chisq_val, pvalue(chisq_val, flavor='flavor')
    
# c = np.full(36,1e-3)
# x2 = chisq(c)
# print x2
# thirtysix = chi2(36)
# print thirtysix.sf(x2)
