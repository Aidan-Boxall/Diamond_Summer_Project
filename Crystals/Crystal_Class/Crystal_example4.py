
# EXAMPLE 4: LOAD AND SAVE WITH CIF AND JSON FILES

import sys
sys.path.insert(0,'/dls_sw/i16/software/python/crystal')

import Crystal as Cr

# Create an instance

mycrys = Cr.Crystal()

# Load from CIF file
# If no 'data code' such as '53058-ICSD' is given, it loads the first crystal structure in the CIF file

mycrys.load_cif('/dls_sw/i16/software/python/crystal/CoO.cif','53058-ICSD')

# Save a copy of the CIF file from the python object

mycrys.save_cif('/dls_sw/i16/software/python/crystal/test_save.cif')

# Recreate the python object from the copy of the CIF file; test by printing reflection table

test = Cr.Crystal()

test.load_cif('/dls_sw/i16/software/python/crystal/test_save.cif','53058-ICSD')

print test.general.cifstring

mycrys.reflection_list(5,refl='sym',anomalous_flag=False,sort='two_theta',print_list=True)

print '\n'

# Save the python object in a JSON file

mycrys.save('/dls_sw/i16/software/python/crystal/test.json')

# Recreate the python object from the JSON file; test by printing reflection table

test2 = Cr.Crystal()

test2.load('/dls_sw/i16/software/python/crystal/test.json')

test2.reflection_list(5,refl='sym',anomalous_flag=False,sort='two_theta',print_list=True)





