
# EXAMPLE 1: LATTICE PROPERTIES AND SPACE GROUP

import sys
sys.path.insert(0,'/dls_sw/i16/software/python/crystal')

import Crystal as Cr

# Create an instance

mycrys = Cr.Crystal()

# Set lattice parameters manually

mycrys.lattice.parameters([5,5,10,90,90,90])

# Get lattice parameters

print mycrys.lattice.parameters(), '\n'

# Extract lattice info

print mycrys.lattice.reciprocal_parameters(), '\n'

print mycrys.lattice.b_matrix(), '\n'

print mycrys.lattice.g_matrix(), '\n'

print mycrys.lattice.g_star_matrix(), '\n'

# Set space group from its number in the International Tables (you can also use the Hermann-Mauguin symbol, especially for non standard setting)

mycrys.structure.space_group(79)

# Extract space group as symbol (give the optional argument format='number' to extract as IT number)

print mycrys.structure.space_group(), '\n'

# Extract space group operators (here they are printed with a loop for nicer printing)

ops = mycrys.structure.space_group_operators()
for op in ops:
    print op[0]
    print op[1], '\n'

# Extract equivalent points for the general position

eps = mycrys.structure.equivalent_points()
for ep in eps:
    print ep

