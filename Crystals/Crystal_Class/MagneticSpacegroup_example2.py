
import sys
sys.path.insert(0,'C:\\Users\\gpk42638\\Desktop\\PD2\\Project software\\Crystal class v10')        

from MagneticSpacegroup import *


mycr = MagneticSpacegroup()

##
mycr.load_cif('C:\\Users\\gpk42638\\Desktop\\PD2\\Project software\\Crystal class v10\\Ag2CrO2.mcif')

ppcr=mycr.parentphase.get_parentphase_crystal()


print '\n', 'Magnetic phase:'
print 'Cell ', mycr.lattice.parameters()
print 'Magnetic space group ', mycr.structure.magnetic_space_group()

print '\n', 'Parent phase:'
print 'Cell ', ppcr.lattice.parameters()
print 'Space group ', ppcr.structure.space_group()


print '\n', 'Reflection list (magnetic structure factor only)', '\n'

mycr.reflection_list(1.2)

print '\n', 'Reflection list (magnetic amplitude, azimuth = [1,0,0])', '\n'

mycr.reflection_list(1.2,azimuth_ref=[1,0,0])

print '\n', 'Reflection list (magnetic amplitude, azimuth = [0,0,1])', '\n'

mycr.reflection_list(1.2,azimuth_ref=[0,0,1])


