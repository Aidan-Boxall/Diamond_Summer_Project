
# EXAMPLE 5: POINT GROUP (STRUCTURE AND SITE-SPECIFIC)

import sys
sys.path.insert(0,'/dls_sw/i16/software/python/crystal')

import Crystal as Cr

# Create an instance

mycrys = Cr.Crystal()

# Set lattice parameters, space group and structure manually

mycrys.lattice.parameters([5,5,10,90,90,90])

mycrys.structure.space_group(79)

mycrys.structure.site(label='Fe',coordinates=[0,0.5,0.77])

# Extract the point group of the structure

print mycrys.structure.point_group(), '\n'

# Extract the operators of the point group (here they are printed with a loop for nicer printing)

ops = mycrys.structure.point_group_operators()
for op in ops:
    print op, '\n'

# Extract the point group of a specific site

print mycrys.structure.site_point_group('Fe'), '\n'

# Extract the operators of the point group of the site

print mycrys.structure.site_point_group_operators('Fe')

