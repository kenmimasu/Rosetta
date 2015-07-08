#!/usr/bin/env python
import HiggsBasis as HB
import WarsawBasis as WB
import SILHBasis as SB
import MassBasis as MB
import TemplateBasis as TB
import SLHA
import tempfile
import os
import sys
import re
import random

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
        if abs(c1 - c2) > tolerance*abs(c1):
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
                         output_basis='higgs',silent=True)

    myinstance.write_param_card(out_card)
    try:
        HB_instance = HB.HiggsBasis(param_card=out_card, 
                             output_basis='mass',silent=True)
        MB_instance = MB.MassBasis(param_card=out_card, silent=True)
        compare_inputs(HB_instance, MB_instance, tolerance=tolerance)
        compare_coeffs(HB_instance, MB_instance, tolerance=tolerance)

    except Exception as e:
        print 'Error!!!'
        print e
    finally:
        os.remove(out_card)
        os.rmdir(tmpdir)
    
def SILH_Warsaw_triangle(param_card,tolerance=1e-4, keep=False):
    tmpdir = tempfile.mkdtemp(prefix = 'rosetta_temp_', dir=os.getcwd())
    SB2MB_out = '{}/SB2MB.dat'.format(tmpdir)
    SB2WB_out = '{}/SB2WB.dat'.format(tmpdir)
    WB2MB_out = '{}/WB2MB.dat'.format(tmpdir)
    try:
        # translate directly to Mass Basis
        SB2MB = SB.SILHBasis(param_card=param_card, 
                             output_basis='mass')
        SB2MB.write_param_card(SB2MB_out)
        
        # translate to Warsaw Basis    
        SB2WB = SB.SILHBasis(param_card=param_card, 
                             output_basis='warsaw')
        SB2WB.write_param_card(SB2WB_out)
        # translate to Mass basis  
        WB2MB = WB.WarsawBasis(param_card=SB2WB_out, 
                             output_basis='mass')
        WB2MB.write_param_card(WB2MB_out)

        # read & compare both Mass Basis translations
        from_SB = MB.MassBasis(param_card=SB2MB_out)
        via_WB = MBMassBasis(param_card=WB2MB_out)

        compare_inputs(from_SB, via_WB, tolerance=tolerance)
        compare_coeffs(from_SB, via_WB, tolerance=tolerance)
    except Exception as e:
        print 'Error!!!'
        print e
    finally:
        if not keep:
            os.remove(SB2WB_out)
            os.remove(SB2MB_out)
            os.remove(WB2MB_out)
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
                out.write('   Description -> "{} coupling parameter"}},\n'.format(name))
                out.write('\n')
                
    print 'wrote ',filename

if __name__=='__main__':
    # generate_coeffs(WB.WarsawBasis,0.,rand=True)
    # generate_coeffs(SB.SILHBasis,0.,rand=True)
    # generate_coeffs(HB.HiggsBasis,0.,rand=True)
    # generate_coeffs(TB.TemplateBasis,0.,rand=True)
    # generate_coeffs(MassBasis,1.,rand=True)
    # generate_frdef(HiggsBasis,'HB_definitions.fr')
    generate_frdef(MB.MassBasis,'MB_definitions.fr')
    # higgs_basis_check(SILHBasis,'../Cards/param_card_SILHBasis_blocks.dat',tolerance=1e-4)
    # higgs_basis_check(SILHBasis,'../Cards/param_card_SILHBasis.dat')
    # SILH_Warsaw_triangle('../Cards/param_card_SILHBasis_blocks.dat', tolerance=1e-4, keep=True)

