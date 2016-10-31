"""
CMod Ex2: Testing of Vector methods 
"""
import Vector as vec
import numpy as np


inFile = open("vector.in","r")
"""
Open file with two vectors
extracts data from file by reading line-by-line
assign data to 3 dimensional numpy arrays: v1, v2 and v3
"""
line = inFile.readline()
tokens = line.split()
v1 = np.array([float(tokens[0]), float(tokens[1]), float(tokens[2])])
line = inFile.readline()
tokens = line.split()
v2 = np.array([float(tokens[0]), float(tokens[1]), float(tokens[2])])
line = inFile.readline()
tokens = line.split()
v3 = np.array([float(tokens[0]), float(tokens[1]), float(tokens[2])])
# Print out all 3 elements of all 3 numpy arrays
print "v1 = " + str(v1)
print "v2 = " + str(v2)
print "v3 = " + str(v3)
print ""

# Test vector methods:
# magnitude, summing , dot product, cross product 

print "magnitude of:\n v1 = " + str(vec.mag(v1))
print "Vector sum of: \n v1, v2 and v3 = " + str(vec.add(vec.add(v1,v2),v3))
print "dot product of:\n v1.v2 = " + str(vec.dot(v1,v2))
print "cross product of: \n v1 x v2 =" + str(vec.cross(v1,v2))

# Test vector operations:
"""The numpy vector operations are printed to the terminal using strings
  There are three "Show that..." statements below, which are demonstrated 
  step by step and then printed to the terminal"""

"""
first "show that.." statement
Print results of v1 x v2 and -v2 x v1 respectively to console to show that they are equal
"""  
print " \nShow that the cross product of v1 x v2 = -v2 x v1"
v4 = vec.scalar(v2,-1)
print "-v2 = " + str(v4)
print "-v2 x v1 = " + str(vec.cross(v4,v1))
print " v1 x v2 = " + str(vec.cross(v1,v2))

"""
Second "show that.." statement
Print results of v1 x (v2 + v3) and (v1 x v2) +(v1 x v3)  respectively to console to show that they are equal
"""  
print "\nShow v1 x (v2 + v3) = (v1 x v2) +(v1 x v3)"
print "v1 x (v2 + v3) = " + str(vec.cross(v1,vec.add(v2,v3)))
print "(v1 x v2) = " + str(vec.cross(v1,v2))
print "(v1 x v3) = " + str(vec.cross(v1,v3))
v4 = vec.cross(v1,v2) + vec.cross(v1,v3)
print "(v1 x v2) +(v1 x v3) = " + str(v4)



"""
Third "show that.." statement
Print results of v1 x (v2 x v3) and (v1.v3)v2 - (v1.v2)v3  respectively to console to show that they are equal
"""  
#Third "show that..." statement
print "\nShow v1 x (v2 x v3) = (v1.v3)v2 - (v1.v2)v3"
print "v1 x (v2 x v3) = " + str(vec.cross(v1,vec.cross(v2,v3)))
print "(v1.v3)v2 = " + str(vec.scalar(v2,vec.dot(v1,v3)))
print "(v1.v2)v3 = " + str(vec.scalar(v3,vec.dot(v1,v2)))
v4 = vec.scalar(v2,vec.dot(v1,v3))
v5 = vec.scalar(v3,vec.dot(v1,v2))
v6 = vec.sub(v4,v5)
print "(v1.v3)v2 - (v1.v2)v3 = " + str(v6)

