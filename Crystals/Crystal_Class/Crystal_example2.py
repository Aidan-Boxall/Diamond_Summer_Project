
# EXAMPLE 2: SITE PROPERTIES

import sys
sys.path.insert(0,'/dls_sw/i16/software/python/crystal')

import Crystal as Cr

# Create an instance

mycrys = Cr.Crystal()

# Set lattice parameters and space group manually

mycrys.lattice.parameters([5,5,10,90,90,90])

mycrys.structure.space_group(79)

# Initialise one site

mycrys.structure.site(label='Fe')

print mycrys.structure.site(), '\n'

# Initialise a list of sites

mycrys.structure.site(label=['Fe','O'])

print mycrys.structure.site(), '\n'

# Remove one site

del mycrys.structure.site()['O']

print mycrys.structure.site(), '\n'

# Initialise/modify sites and properties

mycrys.structure.site(label='Fe',coordinates=[0,0.5,0.77])

print mycrys.structure.site(), '\n'

# Add/modify a property of an existing site

mycrys.structure.site()['Fe']['coordinates'] = [0,0,0.35]

print mycrys.structure.site(), '\n'

print '\n'

# Extract info for a specific site

print mycrys.structure.site_multiplicity('Fe'),  mycrys.structure.site_wyckoff('Fe')
print mycrys.structure.site_coordinates_exact('Fe'), '\n'

print mycrys.structure.site_equivalent_coordinates('Fe'), '\n'
print mycrys.structure.site_symmetry_points('Fe')

print '\n'

# This is identified as a special position
mycrys.structure.site(label=['Fe','O'],coordinates=[[0.05,0.05,0.35],[0.27,0.42,0.77]])
print mycrys.structure.site_multiplicity('Fe'),  mycrys.structure.site_wyckoff('Fe')
print mycrys.structure.site_coordinates_exact('Fe'), '\n'

# This is identified as a general position
mycrys.structure.site(label=['Fe','O'],coordinates=[[0.07,0.07,0.35],[0.27,0.42,0.77]])  
print mycrys.structure.site_multiplicity('Fe'),  mycrys.structure.site_wyckoff('Fe')
print mycrys.structure.site_coordinates_exact('Fe'), '\n'



