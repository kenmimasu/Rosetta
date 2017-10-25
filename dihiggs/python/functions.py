from __future__ import print_function

import numpy as np
#import array as array
from array import array
import math
from scipy.special import gammaln
#from root_numpy import hist2array


class AnalyticalReweighter(object):
    def __init__(self, bins, ref_effs, den_effs, bin_coefs, xs_coefs): 
        self.bins = bins
        self.ref_effs = ref_effs
        self.den_effs = den_effs
        self.bin_coefs = bin_coefs 
        self.xs_coefs = xs_coefs 
        
        # load benchmarks
        self.klJHEP=[1.0,  7.5,  1.0,  1.0,  -3.5, 1.0, 2.4, 5.0, 15.0, 1.0, 10.0, 2.4, 15.0]
        self.ktJHEP=[1.0,  1.0,  1.0,  1.0,  1.5,  1.0, 1.0, 1.0, 1.0,  1.0, 1.5,  1.0, 1.0]
        self.c2JHEP=[0.0,  -1.0, 0.5, -1.5, -3.0,  0.0, 0.0, 0.0, 0.0,  1.0, -1.0, 0.0, 1.0]
        self.cgJHEP=[0.0,  0.0, -0.8,  0.0, 0.0,   0.8, 0.2, 0.2, -1.0, -0.6, 0.0, 1.0, 0.0]
        self.c2gJHEP=[0.0, 0.0, 0.6, -0.8, 0.0, -1.0, -0.2,-0.2,  1.0,  0.6, 0.0, -1.0, 0.0] 
        # same thing, different format for easier argument passing
        self.JHEP_BM = []
        for i in xrange(len(self.klJHEP)):
            self.JHEP_BM.append((self.klJHEP[i], self.ktJHEP[i], self.c2JHEP[i], self.cgJHEP[i], self.c2gJHEP[i],))
        print ("initialize")
        
        """ 
    Class to do event reweighting based on a fitted analytical formula.

    Args:
        bins (list): list of array-like with bin edges used
        ref_effs (array-like): multi-dim array with reference efficiency
            for each bin (e.g. SM).
        den_effs (array-like): multi-dim array with denominator efficiency
            for each bin (e.g. SM).
        bin_coefs (array-like): multi-dim array with the analytical formula
            coeficients per bin
        xs_coefs (array-like): 1D array with the coeficients for the cross
            sections
        """
    
    @staticmethod
    def analytical_formula(kl, kt, c2, cg, c2g, A): 
        return A[0]*kt**4 + A[1]*c2**2 + (A[2]*kt**2 + A[3]*cg**2)*kl**2 \
               + A[4]*c2g**2 + ( A[5]*c2 + A[6]*kt*kl )*kt**2  + (A[7]*kt*kl \
               + A[8]*cg*kl )*c2 + A[9]*c2*c2g  + (A[10]*cg*kl \
               + A[11]*c2g)*kt**2+ (A[12]*kl*cg + A[13]*c2g )*kt*kl \
               + A[14]*cg*c2g*kl
        """ 
        Returns the R value according to the fitted analytical formula.
        Returns:    
            An array with the evaluation of the formula for those parameters,
            might have to be transposed to recover the original shape.
        """
    
    def parameter_point(self, kl, kt, c2, cg, c2g): 
        #precompute formula results for all bins and xs
        # transpose to have coefs in external index
        self.bin_results = self.analytical_formula(kl, kt, c2, cg, c2g,
                self.bin_coefs.T).T
        self.xs_results = self.analytical_formula(kl, kt, c2, cg, c2g,
                self.xs_coefs.T).T
        # compute numerator efficiency and total normalization
        self.num_effs = self.ref_effs*self.bin_results/self.xs_results
        self.Cnorm = self.num_effs.sum()
        # compute weight per bin (normalized by sum)
        self.bin_weights =  (self.num_effs/self.den_effs)/self.Cnorm
        bin_weights_ret =  (self.ref_effs*self.bin_results/self.xs_results)/self.Cnorm
        return bin_weights_ret
        """ Precompute all required information for a certain EFT point.

        Note:
            It has to be called before find_bin and weight methods.

        Args:
            kl (float): EFT parameter
            kt (float): EFT parameter
            c2 (float): EFT parameter
            cg (float): EFT parameter
            c2g (float): EFT parameter
        """


    def compute_TS(self, kl_1, kt_1, c2_1, cg_1, c2g_1, kl_2, kt_2, c2_2, cg_2, c2g_2):
        normEv = 1200000000
        TS = np.zeros(1)
        EvByBin_1 = np.absolute(np.rint(normEv * self.parameter_point(kl_1, kt_1, c2_1, cg_1, c2g_1)))
        logPoint_1 = gammaln(EvByBin_1)
        EvByBin_2 = np.absolute(np.rint(normEv * self.parameter_point(kl_2, kt_2, c2_2, cg_2, c2g_2)))
        log_2 = gammaln(EvByBin_2) # last bins are giving negative = increase fit precision
        NSumInt = np.rint((EvByBin_1 + EvByBin_2) / 2)
        logSum = gammaln(NSumInt)
        test = np.float64(-2 * (logPoint_1 + log_2 - 2 * logSum ))
        TS = test.sum()
        return TS


    def TS_test(self, kl, kt, c2, cg, c2g, verbose=True):
        DEBUG = True
        if verbose:
            print ("Calculating TS")
        TS = np.zeros(13) 
        for bench in xrange(13):
            TS[bench] = self.compute_TS(kl, kt, c2, cg, c2g, *self.JHEP_BM[bench])
        closestBM = np.argmin(TS)
        if DEBUG:
            print(np.where(TS == TS.min()))
            print(TS)
        if verbose:
            print ("Closest benchmark is: {} with TS {}".format(closestBM, TS[closestBM]))
        return int(closestBM)


    def find_bin(self, variables):  
        bin_ids = np.empty_like(variables, dtype=np.int64)
        for i_v, bins_arr in enumerate(self.bins):
            # substract 1 so bins start in 0
            bin_ids[:,i_v] = np.digitize(variables[:,i_v], bins=bins_arr) - 1
        return bin_ids
        """ Finds the bins corresponding to each event.

        Args:
            variable (np.array): 2-D float np.array of shape (n_events, n_rw_vars) 
        Returns:
            A 2D int64 array of shape (n_events, n_rw_vars) with the indexes
            for each event obtained using np.digitize for each variable.
        """


    def weight(self, variables):
        bin_ids =  self.find_bin(variables)
        weights = self.bin_weights[[ bin_ids[:,i_v] for i_v in range(bin_ids.shape[1])]]
        """ Finds the weight corresponding to each event.

        Args:
            variable (np.array): 2-D float np.array of shape (n_events, n_rw_vars) 
        Returns:
            A 1D float array of shape (n_events,) with the corresponding weight
            per event.
        """
        return weights


def reweighter_from_histogram_and_file():

        bins = [([   250.,    260.,    270.,    280.,    290.,    300.,    310.,
          320.,    330.,    340.,    350.,    360.,    370.,    380.,
          390.,    400.,    410.,    420.,    430.,    440.,    450.,
          460.,    470.,    480.,    490.,    500.,    510.,    520.,
          530.,    540.,    550.,    560.,    570.,    580.,    590.,
          600.,    610.,    620.,    630.,    640.,    650.,    660.,
          670.,    680.,    690.,    700.,    750.,    800.,    850.,
          900.,    950.,   1000.,   1100.,   1200.,   1300.,   1400.,
         1500.,   1750.,   2000.,  50000.]), ([ 0.,  0.40000001,  0.60000002,  0.80000001,  1.  ])]
        n_ev_V0_sum=np.loadtxt('data/n_ev_V0_sum_59_4.out')
        den_effs=n_ev_V0_sum/n_ev_V0_sum.sum()

        column_names = ["n_samples"]
        var_center_names = [ "mhh_center", "costhst_center"]
        column_names += var_center_names
        n_ev_names = ["n_ev_SM", "n_ev_JHEP_sum"]
        column_names += n_ev_names 
        coef_names = ["A_{}".format(i) for i in range(0,15)]
        column_names += coef_names
        coef_err_names = ["A_err_{}".format(i) for i in range(0,15)]
        column_names += coef_err_names

        data_from_file = np.genfromtxt("data/coefficientsByBin_extended_3M_costHHSim_59-4.txt", names = column_names)

        # asumming consistent ordering in text file (could digitize otherwise)
        n_ev_SM = data_from_file[n_ev_names[0]].reshape(den_effs.shape)
        ref_effs = n_ev_SM/n_ev_SM.sum() 

        bin_coefs = data_from_file[coef_names].view(np.float)\
            .reshape(den_effs.shape+(len(coef_names),))

        # hardcoded coeficients for total cross section
        xs_coefs = np.array([2.09078, 10.1517, 0.282307, 0.101205, 1.33191,
                         -8.51168, -1.37309, 2.82636, 1.45767, -4.91761,
                         -0.675197, 1.86189, 0.321422, -0.836276,
                         -0.568156])

        analytical_reweighter = AnalyticalReweighter(bins, ref_effs, den_effs, bin_coefs, xs_coefs)

        return analytical_reweighter 

