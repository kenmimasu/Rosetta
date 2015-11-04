#!/usr/bin/env python
from Rosetta import HiggsBasis as HB
from Rosetta import WarsawBasis as WB
from Rosetta import SILHBasis as SB
from Rosetta import BSMCharacterisation as BC
from Rosetta import TemplateBasis as TB
# from Rosetta import MufBasis as MUF
import SLHA
import tempfile
import os
import sys
import re
import random

__doc__='''
Some tester functions to validate translations, generate sample input cards etc.
'''

def compare_inputs(basis1, basis2, tolerance=1e-4):
    wrong_inputs = []
    for _input,value in basis1.inputs.items():
        other_value = basis2.inputs[_input]
        if abs(value - other_value) > tolerance*abs(value):
            wrong_inputs.append([basis1.inputs._names[_input],
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
    for coeff in basis1.all_coeffs:
        c1, c2 = basis1[coeff], basis2[coeff]
        if abs(c1 - c2) > tolerance*abs(c1) and abs(c1 - c2) > 1e-15:
            wrong_coeffs.append([coeff,
                                 ('{}:'.format(basis1.__class__.__name__), c1,
                                  '{}:'.format(basis2.__class__.__name__), c2)
                                ])
    for matrix in basis1.flavored:
        if matrix not in basis1.card: continue
        b1, b2 = basis1.card.matrices[matrix], basis2.card.matrices[matrix]
        for k,c1 in b1.iteritems():
            c2 = b2[k]
            coeff = b1._names[k]
            if abs(c1 - c2) > tolerance*abs(c1) and abs(c1 - c2) > 1e-15:
                wrong_coeffs.append([coeff,
                                     ('{}:'.format(basis1.__class__.__name__), c1,
                                      '{}:'.format(basis2.__class__.__name__), c2)
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
                         output_basis='bsmc',silent=True)
    myinstance.write_param_card(out_card)
        
    HB_instance = HB.HiggsBasis(param_card=out_card, silent=True, 
                                translate=False)
    BC_instance = BC.BSMCharacterisation(param_card=out_card, silent=True)
    
    compare_inputs(HB_instance, BC_instance, tolerance=tolerance)
    compare_coeffs(HB_instance, BC_instance, tolerance=tolerance)
    
    os.remove(out_card)
    os.rmdir(tmpdir)

def two_way_test(basis, card ,target, tolerance=1e-4, flavor='general'):
    # not working for Higgs -> x two way
    first = basis(param_card=card, silent=True, translate=False, 
                  flavor=flavor)
    first.init_dependent()
    first.calculate_dependent()
    
    intermediate = first.translate(target=target)
    intermediate.set_flavor(intermediate.flavor, 'general')
    intermediate.flavor='general'
    intermediate.init_dependent()
    intermediate.calculate_dependent()
    
    second = intermediate.translate(target=basis.name)
    
    first.set_flavor('general', flavor)
    second.set_flavor('general', flavor)
    
    compare_coeffs(first , second, tolerance=tolerance)
    compare_inputs(first , second, tolerance=tolerance)

def circle_test(basis, card, tolerance=1e-4, reverse=False, flavor='general'):
    others = [x for x in ('silh','warsaw', 'higgs') if x is not basis.name]
    
    one = basis(param_card=card, silent=True, translate=False, flavor=flavor)
    if not reverse:
        two = one.translate(target=others[0])
        two.set_flavor(two.flavor, 'general')
        three = two.translate(target=others[1])
    else:
        two = one.translate(target=others[1])
        two.set_flavor(two.flavor, 'general')
        three = two.translate(target=others[0])
    
    three.set_flavor(three.flavor, 'general')
    four = three.translate(target=basis.name)
    four.set_flavor('general',flavor)
    
    compare_inputs(one , four, tolerance=tolerance)
    compare_coeffs(one , four, tolerance=tolerance)
    
def triangle_test(basis, card, target, tolerance=1e-4):
    other = [x for x in ('silh','warsaw', 'higgs') 
             if x not in (basis.name,target)][0]
    
    one = basis(param_card=card, silent=True, translate=False)
    two = one.translate(target=target)
    three = one.translate(target=other)
    three.set_flavor(three.flavor, 'general')
    four = three.translate(target=target)
    four.set_flavor(four.flavor, 'general')
    
    compare_coeffs(two , four, tolerance=tolerance)

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
    for name, blk in SLHA_card.blocks.iteritems():
        if len(blk)==0:
            del SLHA_card.blocks[name]
    for name, blk in myinstance.card.blocks.iteritems():
        print blk
        
def generate_frdef(basis_class,filename):
    myinstance = basis_class()
    SLHA_card = myinstance.card
    for name, blk in SLHA_card.blocks.iteritems():
        for coeff in blk:
            if blk.get_name(coeff) not in myinstance.independent:
                del blk[coeff]
    with open(filename,'w') as out:
        for name, blk in myinstance.card.blocks.iteritems():
            bname = blk.name
            for index, value in blk.iteritems():
                cname = blk.get_name(index)
                is_sin = re.match(r'S[ude]\d\d',cname)
                out.write('{} == {{ ParameterType -> External,\n'.format(cname))
                out.write('   Value -> {},\n'.format(value))
                if not is_sin:
                    out.write('   InteractionOrder -> {QNP, 1},\n')
                # out.write('   TeX -> {},\n'.format(cname))
                out.write('   BlockName -> {},\n'.format(bname))
                out.write('   OrderBlock -> {},\n'.format(index))
                out.write('   Description -> "{} coupling parameter"}},\n'.format(cname))
                out.write('\n')
                if is_sin:
                    cos = cname.replace('S','C')
                    out.write('{} == {{ ParameterType -> Internal,\n'.format(cos))
                    out.write('   Value -> Sqrt[1-{}^2],\n'.format(cname))
                    out.write('   Description -> "{} coupling parameter"}},\n'.format(cos))
                    out.write('\n')
            for real in [blk.get_name(x) for x in blk.keys() 
                                if blk.get_name(x)[-2:]=='Re' ]:
                name = real.replace('Re','')
                imag = real.replace('Re','Im')
                out.write('{} == {{ ParameterType -> Internal,\n'.format(name))
                out.write('   ComplexParamter -> True,\n')
                out.write('   Value -> {{{}+I*{}}},\n'.format(real,imag))
                out.write('   InteractionOrder -> {QNP, 1},\n')
                out.write('   Description -> "{} coupling parameter"}},\n'.format(name))
                out.write('\n')
                
    print 'wrote ',filename

if __name__=='__main__':
    # pass

################################################################################
# # card writers
#     # for flav in ('general',):
#     # for flav in ('universal',):
#     # for flav in ('diagonal',):
#     for flav in ('general','diagonal','universal'):
#         if flav =='general':
#             tail = ''
#         else:
#             tail = '_%s' % flav
#
#         instance = TB.TemplateBasis(flavor=flav)
#         instance.write_template_card('Cards/TemplateBasis{}.dat'.format(tail))
#
#
#         instance = HB.HiggsBasis(flavor=flav)
#         instance.write_template_card('Cards/HiggsBasis{}.dat'.format(tail))
#         instance.write_template_card('HiggsBasis{}_rand.dat'.format(tail), value='random')
#         instance.write_template_card('HiggsBasis{}_1e-3.dat'.format(tail), value=0.001)
#
#         instance = WB.WarsawBasis(flavor=flav)
#         instance.write_template_card('Cards/WarsawBasis{}.dat'.format(tail))
#         instance.write_template_card('WarsawBasis{}_rand.dat'.format(tail), value='random')
#         instance.write_template_card('WarsawBasis{}_1e-3.dat'.format(tail), value=0.001)
#     #
#         instance = SB.SILHBasis(flavor=flav)
#         instance.write_template_card('Cards/SILHBasis{}.dat'.format(tail))
#         instance.write_template_card('SILHBasis{}_rand.dat'.format(tail), value='random')
#         instance.write_template_card('SILHBasis{}_1e-3.dat'.format(tail), value=0.001)
#     instance = BC.BSMCharacterisation()
#     instance.write_template_card('BSMCharacterisation.dat')
    #
################################################################################
# translation testers   
    # two_way_test(WB.WarsawBasis,'WarsawBasis_rand.dat','higgs',tolerance=1e-14)

    # two_way_test(WB.WarsawBasis,'WarsawBasis_diagonal_rand.dat','higgs',tolerance=1e-10,flavor='diagonal')
    # two_way_test(WB.WarsawBasis,'WarsawBasis_universal_rand.dat','higgs',tolerance=1e-14,flavor='universal')
    #
    # two_way_test(HB.HiggsBasis,'HiggsBasis_rand.dat','warsaw')
    #
    # two_way_test(WB.WarsawBasis,'WarsawBasis_rand.dat','silh')
    # two_way_test(WB.WarsawBasis,'WarsawBasis_diagonal_rand.dat','silh',flavor='diagonal')
    # two_way_test(WB.WarsawBasis,'WarsawBasis_universal_rand.dat','silh',flavor='universal')
    # two_way_test(SB.SILHBasis,'SILHBasis_rand.dat','warsaw')

    # two_way_test(SB.SILHBasis,'SILHBasis_rand.dat','higgs',flavor='general')
    two_way_test(SB.SILHBasis,'SILHBasis_diagonal_rand.dat','higgs',flavor='diagonal')
    two_way_test(SB.SILHBasis,'SILHBasis_universal_rand.dat','higgs',flavor='universal')
    two_way_test(HB.HiggsBasis,'HiggsBasis_rand.dat','silh')
    
    # circle_test(HB.HiggsBasis,'HiggsBasis_rand.dat')
    # circle_test(SB.SILHBasis,'SILHBasis_rand.dat')
    # circle_test(WB.WarsawBasis,'WarsawBasis_rand.dat')
    # circle_test(HB.HiggsBasis,'HiggsBasis_rand.dat',reverse=True)
    # circle_test(SB.SILHBasis,'SILHBasis_rand.dat',reverse=True)
    # circle_test(WB.WarsawBasis,'WarsawBasis_rand.dat',reverse=True)
    # circle_test(HB.HiggsBasis,'Cards/param_card_HiggsBasis.dat',reverse=True)
    # circle_test(SB.SILHBasis,'Cards/param_card_SILHBasis.dat')
    # circle_test(SB.SILHBasis,'Cards/param_card_SILHBasis.dat',reverse=True)
    # circle_test(WB.WarsawBasis,'Cards/param_card_WarsawBasis.dat')
    # circle_test(WB.WarsawBasis,'Cards/param_card_WarsawBasis.dat',reverse=True)
    # triangle_test(HB.HiggsBasis,'HiggsBasis_rand.dat','silh')
    # triangle_test(HB.HiggsBasis,'HiggsBasis_rand.dat','warsaw')
    # triangle_test(SB.SILHBasis,'Cards/param_card_SILHBasis.dat','higgs')
    # triangle_test(SB.SILHBasis,'Cards/param_card_SILHBasis.dat','warsaw')
    # triangle_test(WB.WarsawBasis,'Cards/param_card_WarsawBasis.dat','higgs')
    # triangle_test(WB.WarsawBasis,'Cards/param_card_WarsawBasis.dat','silh')
################################################################################
    
    
    # generate_coeffs(WB.WarsawBasis,0.,rand=True)
    # generate_coeffs(SB.SILHBasis,0.,rand=True)
    # generate_coeffs(HB.HiggsBasis,0.,rand=True)
    # generate_coeffs(TB.TemplateBasis,0.,rand=True)
    # generate_coeffs(MassBasis,1.,rand=True)
    # generate_frdef(HiggsBasis,'HB_definitions.fr')
    # generate_frdef(MB.MassBasis,'MB_definitions.fr')
    # higgs_basis_check(SB.SILHBasis,'Cards/param_card_SILHBasis.dat',tolerance=1e-3)
    # higgs_basis_check(SILHBasis,'../Cards/param_card_SILHBasis.dat')
    # SILH_Warsaw_triangle('Cards/param_card_SILHBasis.dat', tolerance=1e-4)

