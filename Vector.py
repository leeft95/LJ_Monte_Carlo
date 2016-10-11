import math
import numpy as np

#square magnitude
def sqmag(a):
    return a*a

# Magnitude of a vector
def mag(a):
    b = sqmag(a)
    return (b[0] + b[1] + b[2])**0.5

# Scalar Multiple of a vector
def scalar(a, scalar):
    return a*scalar

#The vector difference between two vectors
def sub(a, b):
    c=a-b
    return np.array[c]

#Cross product of two vectors
def cross(a,b):
    c= np.cross(a,b)
    return np.array[c]

#Dot product of two vectors
def dot(a,b):
    c=np.inner(a,b)
    return np.array[c]
