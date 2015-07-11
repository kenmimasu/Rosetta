import os
__doc__ ='''
    Reads config.txt and sets some variables relevant to running Rosetta.
'''
################################################################################
dirname= os.path.dirname(__file__)
eHDECAY_dir=''
################################################################################
with open('{}/../config.txt'.format(dirname)) as config:
    lines = [x.strip() for x in config.readlines()]
    
for i,l in enumerate(lines):
    if l and not l.startswith('#'):
        try:
            field, value = tuple(l.split())
        except Exception:
            print ('Rosetta ignored line {} of config.txt'.format(i+1))
            field, value = None, None
    else:
        field, value = None, None
    
    if field == 'eHDECAY_dir':
        eHDECAY_dir = value
    else:
        print ('Rosetta did not recognise option "{}" '.format(field) +
               'in line {} of config.txt'.format(i+1))
################################################################################
        
        
    
    