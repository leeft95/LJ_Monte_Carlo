"""
 CMod Ex2: methods for complex numbers 
 represented as Numpy 1*2 arrays
"""

import math
import numpy as np

# Complex conjugate
def conj(c):
    return np.array([c[0], -c[1]])

# Square modulus
def normSq(c):
    return c[0]**2 + c[1]**2

# Modulus
def norm(c):
    return math.sqrt(normSq(c))

# Multiplication of complex with scalar
def scale(c, scalar):
    return c*scalar

# Complex addition
def add(a, b):
    realPart = a[0] + b[0]
    imagPart = a[1] + b[1]
    return np.array([realPart, imagPart ])

# Complex subtraction
def sub(a, b):
    realPart = a[0] - b[0]
    imagPart = a[1] - b[1]
    return np.array([realPart, imagPart])

# Complex multiplication
def mul(a, b):
    realPart = a[0]*b[0] - a[1]*b[1]
    imagPart = a[0]*b[1] + a[1]*b[0]
    return np.array([realPart, imagPart])

# Complex division
def div(a, b):
   	z = mul(a, conj(b))
	return scale(z, 1.0/normSq(b))

