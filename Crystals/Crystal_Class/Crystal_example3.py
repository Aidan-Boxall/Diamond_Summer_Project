
# EXAMPLE 3: REFLECTION TABLE

import sys
sys.path.insert(0,'/dls_sw/i16/software/python/crystal')

import Crystal as Cr

# Create an instance

mycrys = Cr.Crystal()

# Set lattice parameters, space group and structure manually

mycrys.lattice.parameters([5,5,10,90,90,90])

mycrys.structure.space_group(79)

mycrys.structure.site(label=['Fe','O','Ge'],coordinates=[[0.15,0.25,0.35],[0,0,0.35],[0,0.5,0.77]],atom_type=['Fe','O','Ge'])

# 1) Calculate and print reflection table, only symmetrically equivalent reflections
# First, required argument is the Energy
# Use sort = 'two_theta','int','d_spacing','index' to sort by two theta, intensity, d spacing, hkl

# Considering Friedel pairs as equivalent
print 'Symmetrically equivalent, No Anomalous'
mycrys.reflection_list(2,refl='sym',anomalous_flag=False,sort='two_theta')

# Considering Friedel pairs as NOT equivalent
print '\nSymmetrically equivalent, Anomalous'
mycrys.reflection_list(2,refl='sym',anomalous_flag=True,sort='two_theta')
print '\n'

# 2) Calculate and print reflection table, all allowed reflections

# Considering Friedel pairs as equivalent
print 'Allowed, No Anomalous'
mycrys.reflection_list(2,refl='allowed',anomalous_flag=False,sort='two_theta')

# Considering Friedel pairs as NOT equivalent
print '\nAllowed, Anomalous'
mycrys.reflection_list(2,refl='allowed',anomalous_flag=True,sort='two_theta')
print '\n'

# 3) Calculate and print reflection table, all reflections (allowed and extinct)

# Considering Friedel pairs as equivalent
print 'All, No Anomalous'
mycrys.reflection_list(2,refl='all',anomalous_flag=False,sort='two_theta')

# Considering Friedel pairs as NOT equivalent
print '\nAll, Anomalous'
mycrys.reflection_list(2,refl='all',anomalous_flag=True,sort='two_theta')

