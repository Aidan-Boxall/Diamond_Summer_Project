__doc_title__ = "Quaternion dtype for NumPy"
__doc__ = "Adds a quaternion dtype to NumPy."

__all__ = ['quaternion']

import numpy
from numpy_quaternion import quaternion

if numpy.__dict__.get('quaternion') is not None:
    raise RuntimeError('The NumPy package already has a quaternion type')

numpy.quaternion = quaternion
numpy.typeDict['quaternion'] = numpy.dtype(quaternion)

