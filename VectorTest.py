"""
CMod Ex2: Testing of Vector methods 
"""
import Vector as vec
import numpy as np

"""
Open file with two vectors,
read them in line-by-line and assign
them to vectors v1,v2 and v3.
"""
inFile = open("vector.in","r")
line = inFile.readline()
tokens = line.split()
v1 = np.array([float(tokens[0]), float(tokens[1]), float(tokens[2])])
line = inFile.readline()
tokens = line.split()
v2 = np.array([float(tokens[0]), float(tokens[1]), float(tokens[2])])
line = inFile.readline()
tokens = line.split()
v3 = np.array([float(tokens[0]), float(tokens[1]), float(tokens[2])])
# Print out all 3 numbers
print "v1 = " + str(v1)
print "v2 = " + str(v2)
print "v3 = " + str(v3)
print ""

# Test vector methods:
# magnitude, suming , dot product, cross product 
print "magnitude of:\n v1 = " + str(vec.mag(v1))
print "Vector sum of: \n v1, v2 and v3 = " + str(vec.add(vec.add(v1,v2),v3))
print "dot product of:\n v1.v2 = " + str(vec.dot(v1,v2))
print "cross product of: \n v1 x v2 =" + str(vec.cross(v1,v2))

# Test vector operations:
print " \nShow that the cross product of v1 x v2 = -v2 x v1"

v4 = vec.scalar(v2,-1)
print "-v2 = " + str(v4)
print "-v2 x v1 = " + str(vec.cross(v4,v1))
print " v1 x v2 = " + str(vec.cross(v1,v2))

print "\nShow v1 x (v2 + v3) = (v1 x v2) +(v1 x v3)"
print "v1 x (v2 + v3) = " + str(vec.cross(v1,vec.add(v2,v3)))
print "(v1 x v2) = " + str(vec.cross(v1,v2))
print "(v1 x v3) = " + str(vec.cross(v1,v3))
v4 = vec.cross(v1,v2) + vec.cross(v1,v3)
print "(v1 x v2) +(v1 x v3) = " + str(v4)

print "\nShow v1 x (v2 x v3) = (v1.v3)v2 - (v1.v2)v3"
print "v1 x (v2 x v3) = " + str(vec.cross(v1,vec.cross(v2,v3)))
print "(v1.v3)v2 = " + str(vec.scalar(v2,vec.dot(v1,v3)))
print "(v1.v2)v3 = " + str(vec.scalar(v3,vec.dot(v1,v2)))
v4 = vec.scalar(v2,vec.dot(v1,v3))
v5 = vec.scalar(v3,vec.dot(v1,v2))
v6 = vec.sub(v4,v5)
print "(v1.v3)v2 - (v1.v2)v3 = " + str(v6)

