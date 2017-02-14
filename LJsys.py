import sys
import matplotlib.pyplot as pyplot
import numpy as np
import math as math
from Particle3D import Particle3D
import MDUtilities as MDUtilities
from copy import copy


if len(sys.argv)!=3:
    print "Wrong number of arguments."
    print "Usage: " + sys.argv[0] + "<input file>" + "<output file>"
    quit()
else:
    filename = sys.argv[1]    
    outfileName = sys.argv[2]

outfile = open(outfileName, "w")
infile = open(filename,"r")


pos = np.array([0.0,0.0,0.0],float)
vel = np.array([0.0,0.0,0.0],float)
mass = float(1)
particles = [Particle3D.base(pos,vel,mass) for i in range(2)]

rho = 1
temp = 273
MDUtilities.setInitialPositions(rho, particles)
MDUtilities.setInitialVelocities(temp, particles)

for i in range(2):
	particles[i].name = i 

for i in range(2):
	print particles[i]
	


numstep = 1000
time = 0.0
dt = 0.001
pe = 0.0
e = 0.0

force = 48.0*((1/(np.linalg.norm(particles[0].position-particles[1].position)**14))-(1/(2*(np.linalg.norm(particles[0].position-particles[1].position)**8))))*Particle3D.vec_sep(particles[0].position,particles[1].position)


tValue = []	
posValue_x = [particles[0].position[0]]
posValue_y = [particles[0].position[1]]
kE = []
outfile.write("Time" + "     X" + "        Y" + "        Total Energy" "\n")
outfile.write("{0:f} {1:f} {2:f} {3}\n".format(time, particles[0].position[0], particles[0].position[1],e))

for i in range(numstep):


    particles[0].leapPos2nd(dt, force)
 
    force_new = 48.0*((1/(np.linalg.norm(particles[0].position-particles[1].position)**14))-(1/(2*(np.linalg.norm(particles[0].position-particles[1].position)**8))))*Particle3D.vec_sep(particles[0].position,particles[1].position)
    
    v = particles[0].leapVelocity(dt, 0.5*(force+force_new))

  
    force = copy(force_new)
  
    pe = 4.0*((1/(np.linalg.norm(particles[0].position-particles[1].position)**12))-(1/(np.linalg.norm(particles[0].position-particles[1].position)**6)))
    e = pe + particles[0].kineticEnergy(v)

    time = time + dt
    
    posValue_y.append(particles[0].position[1])
    posValue_x.append(particles[0].position[0])
    kE.append(e)
    tValue.append(time)	

    outfile.write("{0:f} {1:f} {2:f} {3}\n".format(time, particles[0].position[0], particles[0].position[1],e))


outfile.close()
infile.close()

pyplot.figure()
pyplot.subplot(111)
pyplot.plot(posValue_y,posValue_x)
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

