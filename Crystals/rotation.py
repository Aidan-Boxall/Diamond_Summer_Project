import numpy as np

import quaternion

def unit_vector(vector):
    """Returns a unit vector in the same direction as vector"""
    return vector / np.linalg.norm(vector)
        
def convert_to_radians(angle):
    """Converts an angle in degrees to radians"""
    return angle * (np.pi/180.0)

def rotation_to_quaternion(vector, angle):
    """Converts a rotation in terms of an angle and unit vector into a 
    quaternion.
    
    Args:
        axis: The vector in the direction of the axis of rotation.
    
        angle The angle through which the rotation is made. It has to 
            be passed in degrees.
    
    Returns:
        A quaternion that performs the rotation.
    """
    axis = unit_vector(vector)
    angle = convert_to_radians(angle)
    
    _ = axis * np.sin(angle/2)
    
    q = np.quaternion(np.cos(angle/2), _[0], _[1], _[2])
    
    return q



def rotate(vector, quaternion):
    """Rotates a vector using a quaternion.
    
    Args:
        vector: The vector which is to be rotated.
    
        quaternion: The quaternion which describes the rotation.
    
    Returns:
        The rotated vector.
    """
    
    vector = np.quaternion(0, vector[0], vector[1], vector[2])
    new_vector = quaternion * vector * np.conjugate(quaternion)
    new_vector = np.array(new_vector.imag)
    return new_vector
    
axis = np.array([0,0,1])
x = np.array([1,0,0])
angle = 45


q = rotation_to_quaternion(axis,angle)

print rotate(x,q)

