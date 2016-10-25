#
import math
import numpy as np

#square magnitude
def sqmag(a):
    return a*a

# Magnitude of a vector

def mag(a):
    """
    return the magnitude of a vector
    :param a: 3 dimesnionsal vector a
    :return: magnitude of a
    """ 
    b = sqmag(a)
    return (b[0] + b[1] + b[2])**0.5

# Scalar Multiple of a vector
def scalar(a, scalar):
    return a*scalar

#The vector difference between two vectors
def sub(a, b):
    """
    subtract one vector from another 
    :param a: 3 dimesnionsal vector a
    :param b: 3 dimesnionsal vector b
    :return: vector diffreence between a and b 
    """ 
    d = a-b
    return d

#The vector sum between two vectors
def add(a,b):
    d = a+b
    return d

#Cross product of three vectors
def cross(a,b)
    """
    return the cross product of two vectors
    :param a: 3 dimesnionsal vector a
    :param b: 3 dimesnionsal vector b
    :return: cross product axb
    """ 
    d = np.array[(a[2]*b[3] - a[3]*b[2]),(a[3]*b[1] - a[1]*b[3]),(a[1]*b[2] - a[2]*b[1])]
    return d

#Dot product of two vectors
def dot(a,b):
    d = a[1]*b[1] + a[2]*b[2] + a[3]*b[3]
    return d
