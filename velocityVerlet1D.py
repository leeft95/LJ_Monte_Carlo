"""
CMod Ex3: velocity Verlet time integration of
particle in anharmonic potential.

Produce both a plot of the amplitude and an
output file.
"""

import sys
import matplotlib.pyplot as pyplot
from Particle3D import Particle3D
import sys



# Read name of output file from command line
if len(sys.argv)!=3:
    print "Wrong number of arguments."
    print "Usage: " + sys.argv[0] + " <input file>" + "<output file>"
    quit()
else:
    filename = sys.argv[1]    
    outfileName = sys.argv[2]

# Open output file for writing
outfile = open(outfileName, "w")
infile = open(filename,"r")
# Set up particle
p1 = Particle3D.from_file(infile)
particle2 = Particle3D.from_file(infile)
print p1
print particle2
print Particle3D.kineticEnergy(p1)
"""
# Set up simulation parameters
numstep = 100
time = 0.0
dt = 0.1

# Set up anharmonic potential constants
fc2 = 1.0
fc4 = 1.0
force = -fc2*p1.position -fc4*p1.position**3

# Set up data lists
tValue = [time]
posValue = [p1.position]

outfile.write("{0:} {1:}\n".format(time, p1.position))

# Start the time integration loop

for i in range(numstep):
    # Update particle position
    p1.leapPos2nd(dt, force)
    # Update force
    force_new = -fc2*p1.position -fc4*p1.position**3
    # Update particle velocity, based on average
    # of current and new forces
    p1.leapVelocity(dt, 0.5*(force+force_new))

    # Reset force variable
    force = force_new

    # Increase time
    time = time + dt
    
    # Output particle information
    tValue.append(time)
    posValue.append(p1.position)
    outfile.write("{0:} {1:}\n".format(time, p1.position))

# Close output file
outfile.close()

# Plot graph of amplitude vs time
pyplot.plot(tValue,posValue)
pyplot.show()
"""
