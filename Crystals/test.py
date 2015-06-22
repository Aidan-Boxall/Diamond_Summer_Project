import cctbx
import Crystal as c
import math as mt
import numpy as np

import rotation as r

l = 1  # l is the length between the sample and the detector.
w = 0.01  # w is the width of the detector
energy = 8  # photon energy in keV

def change_in_momentum(momentum_in, axis_of_rotation, two_theta):
    """Finds the change in momentum vector given momentum in, the axis of 
    rotation and the two theta angle
    
    Args:
        momentum_in: The momentum vector of the incident beam. Units keV/c.
        
        axis_of_rotion: The axis perpendicular to the difraction plane centred
            on the sample.
        
        two_theta: The two_theta angle of difraction. Angle must be passed in
            degrees.
        
    Returns:
        change_in_momentum: The change in momentum vector.
    """
    _q = r.rotation_to_quaternion(axis_of_rotation, two_theta)
    momentum_out = r.rotate(momentum_in, _q)
    return momentum_out - momentum_in

# Loads the crystal data
mycrys = c.Crystal()
mycrys.load_cif('NiCO3_icsd_61067.cif')

# Produces delta_chi for each reflection
mylist = mycrys.reflection_list(energy)
two_thetas = []

for lists in mylist:
    two_thetas.append(lists[4])

thetas = [two_theta / 2.0 * (mt.pi / 180) for two_theta in two_thetas]

radii = [l * mt.sin(theta) for theta in thetas]

delta_chis = [w / radius for radius in radii]
for i, theta in enumerate(two_thetas):
    print i + 1, theta, delta_chis[i] * (180 / mt.pi)

# Tests the change in momentum code.
# Uses the beamline as the x-axis and the vertical as the z axis.
# For testing purposes assumes the axis of rotation is the y axis and the 
#     diffracted beam is the x-z plane.

y_axis = np.array([0, -1.0, 0])  # The negative y axis is the axis of rotation 
#                                     for the momentum.
momentum_in = np.array([energy, 0, 0])
two_theta = 45

print change_in_momentum(momentum_in, y_axis, two_theta)
