"""
CMod Ex3: velocity Verlet time integration of
particle 
"""

import sys
import matplotlib.pyplot as pyplot
import numpy as np
import math as math
from Particle3D import Particle3D
from copy import copy


""" Read name of input and output file from command line
"""

if len(sys.argv)!=3:
    print "Wrong number of arguments."
    print "Usage: " + sys.argv[0] + "<input file>" + "<output file>"
    quit()
else:
    filename = sys.argv[1]    
    outfileName = sys.argv[2]

""" Open output and input file for input and output
"""

outfile = open(outfileName, "w")
infile = open(filename,"r")

""" Set up particle
"""

p1 = Particle3D.from_file(infile)
p2 = Particle3D.from_file(infile)



""" Set up simulation parameters
"""
numstep = 100
time = 0.0
dt = 0.1
pe = 0.0
e = 0.0


""" Set up data lists
"""

tValue = []
posValue_x = [p1.position[0]]
posValue_y = [p1.position[1]]
kE = []
outfile.write("Time" + "     X" + "        Y" + "        Total Energy" "\n")
outfile.write("{0:f} {1:f} {2:f} {3}\n".format(time, p1.position[0], p1.position[1],e))

""" Set up force constants
"""

fc2 = p1.mass
fc4 = p2.mass
force = -(fc2*fc4)/(np.linalg.norm(p1.position-p2.position)**3)*Particle3D.vec_sep(p1.position,p2.position)

""" Start the time integration loop
"""

for i in range(numstep):
    """ Update particle position """
    p1.leapPos2nd(dt, force)
    """ Update force """
    force_new = -(fc2*fc4)/(np.linalg.norm(p1.position-p2.position)**3)*Particle3D.vec_sep(p1.position,p2.position)
    """ Update particle velocity, based on average
    of current and new forces """
    v = p1.leapVelocity(dt, 0.5*(force+force_new))

    """ Reset force variable """
    force = copy(force_new)
    """ update particle potential energy """
    pe = float(-(fc2*fc4)/(np.linalg.norm(p1.position-p2.position)))
    e = pe + p1.kineticEnergy(v)
    """ Increase time """
    time = time + dt
    
    """ Output particle information """
    	
    posValue_y.append(p1.position[1])
    posValue_x.append(p1.position[0])
    kE.append(e)
    tValue.append(time/math.pi)	

    outfile.write("{0:f} {1:f} {2:f} {3}\n".format(time, p1.position[0], p1.position[1],e))

""" Close output file """
outfile.close()
infile.close()

""" Plot graph of x and y position """
pyplot.figure()
pyplot.subplot(111)
pyplot.plot(posValue_x,posValue_y)
pyplot.title('Trajectory using Velocity Verlet')
pyplot.xlabel('x position')
pyplot.ylabel('y position')
pyplot.savefig('TrajectoryVV.png')
pyplot.figure()
pyplot.subplot(111)
pyplot.plot(tValue,kE)
pyplot.title('Total Energy of the particle against time')
pyplot.xlabel('Time(s)')
pyplot.ylabel('Energy(J)')
pyplot.savefig('EnergyVV.png')

