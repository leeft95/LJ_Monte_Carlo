"""
CMod Ex3: symplectic Euler time integration of
particle in anharmonic potential.

Produce both a plot of the amplitude and an
output file.
"""

import sys
import matplotlib.pyplot as pyplot
from Particle3D import Particle3D
import numpy as np

# Read name of output file from command line
if len(sys.argv)!=3:
    print "Wrong number of arguments."
    print "Usage: " + sys.argv[0] + "<input file>" + "<output file>"
    quit()
else:
    filename = sys.argv[1]    
    outfileName = sys.argv[2]


# Open output file for writing
outfile = open(outfileName, "w")
infile = open(filename,"r")
# Set up particle
p1 = Particle3D.from_file(infile)
p2 = Particle3D.from_file(infile)

# Set up simulation parameters
numstep = 100
time = 0.0
dt = 0.1

# Set up force constants
fc2 = p1.mass
fc4 = p2.mass

# Set up data lists
tValue = [time]
posValue_x = [p1.position[0]]
posValue_y = [p1.position[1]]
outfile.write("{0:f} {1:f} {2:f}\n".format(time, p1.position[0], p1.position[1]))

# Start the time integration loop

for i in range(numstep):
    # Update particle position
    p1.leapPos1st(dt)
    # Update force
    force = -(fc2*fc4)/(np.linalg.norm(p1.position-p2.position)**3)*Particle3D.vec_sep(p1.position,p2.position)
    # Update particle velocity
    p1.leapVelocity(dt, force)

    # Increase time
    time = time + dt
    
    # Output particle information
    posValue_y.append(p1.position[1])
    posValue_x.append(p1.position[0])
	
    outfile.write("{0:f} {1:f} {2:f}\n".format(time, p1.position[0], p1.position[1]))

# Close output file
outfile.close()

# Plot graph of amplitude vs time
pyplot.plot(posValue_x,posValue_y)
pyplot.show()

