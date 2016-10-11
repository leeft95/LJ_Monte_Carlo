"""
CMod Ex2: tester of Complex number operations
"""
import Vector as vec
import numpy as np

"""
Open file with two complex variables,
read them in line-by-line and assign
them to variables x1 and x2.
"""
inFile = open("vectors.in","r")
line = inFile.readline()
tokens = line.split()
v1 = np.array([float(tokens[0]), float(tokens[1]), float(tokens[2])])
line = inFile.readline()
tokens = line.split()
v2 = np.array([float(tokens[0]), float(tokens[1]), float(tokens[2])])
# Print out both numbers
print "x1 = " + str(x1)
print "x2 = " + str(x2)
print ""

# Test complex number methods:
# conjugation, modulus (squared), scaling
print "conj(x2) = " + str(cplx.conj(x2))
print "|x2|^2 = " + str(cplx.normSq(x2))
print "|x2|   = " + str(cplx.norm(x2))
print "3*x2   = " + str(cplx.scale(x2, 3.0))

# Test arithmetic operations

# addition
z2 = cplx.add(x1,x2)
print "z2 = x1+x2 = " + str(z2)
# subtraction
z3 = cplx.sub(x1, x2)
print "z3 = x1-x2 = " + str(z3)
# complex multiplication
z4 = cplx.mul(x1, x2)
print "z4 = x1*x2 = " + str(z4)
# complex division
z5 = cplx.div(x2, x1)
print "z5 = x2/x1 = " + str(z5)
