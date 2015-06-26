#!/usr/bin/env python
from HiggsBasis import HiggsBasis
from WarsawBasis import WarsawBasis
from SILHBasis import SILHBasis
from MassBasis import MassBasis
import SLHA
import tempfile
import os
import sys
import random

def compare_inputs(basis1, basis2, tolerance=1e-4):
    wrong_inputs = []
    for _input,value in basis1.input.items():
        other_value = basis2.input[_input]
        if abs(value - other_value) > tolerance*abs(value):
            wrong_inputs.append([_input,
                                 ('{}:'.format(basis1.__class__.__name__),
                                  other_value,
                                  '{}:'.format(basis2.__class__.__name__),
                                  value)
                                ])
 
    if wrong_inputs:
        print 'Some inconsistent SM inputs detected:'
        for w in wrong_inputs:
            print w[0]
            print w[1]
    else:
        print 'Inputs are OK!'
        
def compare_coeffs(basis1, basis2, tolerance=1e-4):
    wrong_coeffs = []
    for coeff,value in basis1.newpar.items():
        other_value = basis2.newpar[coeff]
        if abs(value - other_value) > tolerance*abs(value):
            wrong_coeffs.append([coeff,
                                 ('{}:'.format(basis1.__class__.__name__),
                                  other_value,
                                  '{}:'.format(basis2.__class__.__name__),
                                  value)
                                ])

    if wrong_coeffs:
        print 'Some inconsistent coefficients detected:'
        for w in wrong_coeffs:
            print w[0]
            print w[1]
    else:
        print 'Coeffs are OK!'

def higgs_basis_check(MyBasis,param_card,tolerance=1e-4):
    tmpdir = tempfile.mkdtemp(prefix = 'rosetta_temp_', dir=os.getcwd())
    out_card = '{}/output_card.dat'.format(tmpdir)
    myinstance = MyBasis(param_card=param_card, 
                         keep_old=False, 
                         output_basis='mass')
    myinstance.write_param_card(out_card)

    HB_instance = HiggsBasis(param_card=out_card, 
                         keep_old=False, 
                         output_basis='mass')
                         
    compare_inputs(HB_instance, myinstance, tolerance=tolerance)
    compare_coeffs(HB_instance, myinstance, tolerance=tolerance)

    os.remove(out_card)
    os.rmdir(tmpdir)
    
def SILH_Warsaw_triangle(param_card,tolerance=1e-4):
    tmpdir = tempfile.mkdtemp(prefix = 'rosetta_temp_', dir=os.getcwd())
    out_card = '{}/output_card.dat'.format(tmpdir)
    to_warsaw = SILHBasis(param_card=param_card, 
                         keep_old=False, 
                         output_basis='warsaw')
    to_warsaw.write_param_card(out_card)
    WB_instance = WarsawBasis(param_card=out_card, 
                         keep_old=False, 
                         output_basis='mass')

    to_mass = SILHBasis(param_card=param_card, 
                         keep_old=False, 
                         output_basis='mass')

    compare_inputs(WB_instance, to_mass, tolerance=tolerance)
    compare_coeffs(WB_instance, to_mass, tolerance=tolerance)
    
    os.remove(out_card)
    os.rmdir(tmpdir)
    
    
def generate_coeffs(basis_class, val, rand=False):
    myinstance = basis_class()
    SLHA_card = myinstance.card
    # print myinstance.independent
    for name, blk in SLHA_card.blocks.iteritems():
        for coeff in blk:
            if blk.get_name(coeff) not in myinstance.independent:
                del blk[coeff]
            else:
                blk[coeff]=val if not rand else random.uniform(-1.,1.)
    for name, blk in myinstance.card.blocks.iteritems():
        print blk
        
def generate_frdef(basis_class,filename):
    myinstance = basis_class()
    SLHA_card = myinstance.card
    for name, blk in SLHA_card.blocks.iteritems():
        for coeff in blk:
            if blk.get_name(coeff) not in myinstance.independent:
                del blk[coeff]
        # print repr(blk)
    with open(filename,'w') as out:
        for name, blk in myinstance.card.blocks.iteritems():
            bname = blk.name
            for index, value in blk.iteritems():
                cname = blk.get_name(index)
                out.write('{} == {{ ParameterType -> External,\n'.format(cname))
                out.write('   Value -> {},\n'.format(value))
                out.write('   InteractionOrder -> {QNP, 1},\n')
                out.write('   TeX -> {},\n'.format(cname))
                out.write('   BlockName -> {},\n'.format(bname))
                out.write('   OrderBlock -> {},\n'.format(index))
                out.write('   Description -> "{} coupling parameter"}},\n'.format(cname))
                out.write('\n')

if __name__=='__main__':
    # generate_coeffs(WarsawBasis,0.,rand=False)
    # generate_coeffs(SILHBasis,0.,rand=True)
    generate_coeffs(HiggsBasis,1.,rand=True)
    # generate_coeffs(MassBasis,1.,rand=True)
    # generate_frdef(HiggsBasis,'HB_definitions.fr')
    # generate_frdef(MassBasis,'MB_definitions.fr')
    # higgs_basis_check(SILHBasis,'../Cards/param_card_SILHBasis.dat',tolerance=1e-4)
    # higgs_basis_check(SILHBasis,'../Cards/param_card_SILHBasis.dat')
    # SILH_Warsaw_triangle('../Cards/param_card_SILHBasis_random.dat', tolerance=1e-4)
