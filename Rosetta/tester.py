#!/usr/bin/env python
from HiggsBasis import HiggsBasis
from WarsawBasis import WarsawBasis
from SILHBasis import SILHBasis
from MassBasis import MassBasis
import tempfile
import os
import sys
import random



def higgs_basis_check(MyBasis,param_card,tolerance=1e-4):
    tmpdir = tempfile.mkdtemp(prefix = 'rosetta_temp_', dir=os.getcwd())
    out_card = '{}/output_card.dat'.format(tmpdir)
    myinstance = MyBasis(param_card=param_card, 
                         keep_old=False, 
                         output_basis='mass')
    # print 'cll1221:',myinstance.coeffs.cll1221
    myinstance.write_param_card(out_card)

    HB_instance = HiggsBasis(param_card=out_card, 
                         keep_old=False, 
                         output_basis='mass')
    
    wrong_inputs = []
    for _input,value in HB_instance.input.items():
        my_value = myinstance.input[_input]
        if abs(value - my_value) > tolerance*abs(value):
            wrong_inputs.append([_input,
                                 ('{}:'.format(myinstance.__class__.__name__),
                                  my_value,
                                  'HiggsBasis:',
                                  value)
                                ])
 
    if wrong_inputs:
        print 'Some inconsistent SM inputs detected:'
        for w in wrong_inputs:
            print w[0]
            print w[1]
 
    wrong_coeffs = []
    for coeff,value in myinstance.newpar.items():
        value_HB = HB_instance.newpar[coeff]
        if abs(value - value_HB) > tolerance*abs(value):
            wrong_coeffs.append([coeff,
                                 ('{}:'.format(myinstance.__class__.__name__),
                                  value,
                                  'HiggsBasis:',
                                  value_HB)
                                ])

    if wrong_coeffs:
        print 'Some inconsistent coefficients detected:'
        for w in wrong_coeffs:
            print w[0]
            print w[1]

    os.remove(out_card)
    os.rmdir(tmpdir)

def generate_coeffs(basis_class, rand=False):
    myinstance = basis_class()
    for i,coeff in enumerate(myinstance.independent):
        val = random.uniform(-1.,1.) if rand else 0.1
        print '    {:<2} {:.5e} # {}'.format(i,val,coeff)

if __name__=='__main__':
    # generate_coeffs(WarsawBasis,rand=True)
    # higgs_basis_check(WarsawBasis,'../Cards/param_card_WarsawBasis.dat')
    higgs_basis_check(SILHBasis,'../Cards/param_card_SILHBasis.dat')
