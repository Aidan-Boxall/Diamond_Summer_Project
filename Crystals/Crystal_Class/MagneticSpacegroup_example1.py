
import sys
sys.path.insert(0,'C:\\Users\\gpk42638\\Desktop\\PD2\\Project software\\Crystal class v10')

from MagneticSpacegroup import *


# Create an instance

mycrys = MagneticSpacegroup()

# Load a magCIF file

mycrys.load_cif('C:\\Users\\gpk42638\\Desktop\\PD2\\Project software\\Crystal class v10\\Ag2CrO2.mcif')

#########################

print '\n', 'cell parameters: ', mycrys.lattice.parameters(), '\n'
print 'magnetic space group:', mycrys.structure.magnetic_space_group(), '\n'
print 'crystallographic space group (extracted from the magnetic space group):', mycrys.structure.space_group(), '\n'

print 'sites (label, coordinates, magnetic moments): \n'
for site_label in mycrys.structure.site():
    print site_label, mycrys.structure.site()[site_label]['coordinates'], mycrys.structure.site()[site_label]['magnetic_moment']

#########################

print '\n', 'Magnetic Space group operators with cell centering (rotation, translation, time inversion): \n'

m_operators = mycrys.structure.magnetic_space_group_operators()

for op in m_operators:
    print op[0], op[1], op[2], '\n'

#########################

print 'We apply one magnetic space group operator to the coordinates and the magnetic moment of a site.'

operator = m_operators[5]

print 'The magnetic operator (rot, trasl, time rev) is: \n', operator[0], operator[1], operator[2]

site_coord = mycrys.structure.site()['Cr3']['coordinates']
site_moment = mycrys.structure.site()['Cr3']['magnetic_moment']

print 'that is, 2-fold axis along y + centering translation + time reversal'
print '\n', 'we apply it to the coordinates of site Cr3:', site_coord
print 'and the magnetic moment of site Cr3:', site_moment

new_coord = mycrys.structure.apply_to_coordinates(operator,site_coord)
new_moment = mycrys.structure.apply_to_magnetic_moment(operator,site_moment)

print '\n', 'and we get the transformed coordinates:', new_coord
print 'and the transformed magnetic moment:', new_moment, '(reversed twice by the 2-fold axis and the time reversal)'

#########################

print '\n', 'Choose a magnetic site, e.g. Cr2, and find its magnetic point group operators. \n'
print 'Magnetic Point group operators (rotation, time inversion): \n'

point_operators = mycrys.structure.site_magnetic_point_group_operators('Cr2')

for op in point_operators:
    print op[0], op[1], '\n'

print 'Use this information to calculate the allowed magnetic components. \n'
print 'Start with component along crystallographic axis x, and apply the Magnetic point group operators of the site. \n'
print 'The transformed moments are: \n'

moment = [1,0,0]
moment_transformed = []

for op in point_operators:
    moment_transformed.append(mycrys.structure.apply_to_magnetic_moment(op,moment))

print moment_transformed

print '\n', 'therefore a moment along x is allowed. \n'
print 'Magnetic component along y: The transformed moments are \n'


moment = [0,1,0]
moment_transformed = []

for op in point_operators:
    moment_transformed.append(mycrys.structure.apply_to_magnetic_moment(op,moment))

print moment_transformed

print '\n', 'therefore a moment along y is NOT allowed. \n'
print 'Magnetic component along z: The transformed moments are \n'

moment = [0,0,1]
moment_transformed = []

for op in point_operators:
    moment_transformed.append(mycrys.structure.apply_to_magnetic_moment(op,moment))

print moment_transformed

print '\n', 'therefore a moment along z is allowed. \n'
print 'These results are confirmed in the magCIF file, in the tag _atom_site_magnetic_moment_symmetry_constraints_mxmymz,'
print 'but we have found them here independently. The actual magnetic moment in the site is: \n'

print mycrys.structure.site()['Cr2']['magnetic_moment']

print '\n', 'and now it can be seen that the component along y is exactly 0, while the component along x needs not to be.'

#########################

print '\n', 'The only magnetic info that is needed from the magCIF file (and that needs to be saved in the json file) is:'
print '1) which space group operators are primed (= multiplied by time inversion)'
print '2) the magnetic moments in every site'
print "It is easy to construct the same object 'mycrys' by uploading only the crystallographic info, and setting the magnetic info from command line."

mycrys.save('C:\\Users\\gpk42638\\Desktop\\PD2\\Project software\\Crystal class v10\\mycrys.json')

mycrys2 = MagneticSpacegroup()

mycrys2.load('C:\\Users\\gpk42638\\Desktop\\PD2\\Project software\\Crystal class v10\\mycrys.json')

