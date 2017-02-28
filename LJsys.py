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
particles = [Particle3D.base(pos,vel,mass) for i in range(4)]

rho = 1
temp = 0.1
MDUtilities.setInitialPositions(rho, particles)
MDUtilities.setInitialVelocities(temp, particles)

for i in range(4):
	particles[i].name = "s" + str(i + 1)
	print particles[i]

numstep = 100
time = 0.0
dt = 0.01
pe = 0.0
e = 0.0

tValue = []	
posValue_x = [particles[0].position[0]]
posValue_y = [particles[0].position[1]]
kE = []
force = []
force_new = []
for j in range(0,4):
	print j
	f = 0.0
	print f
	for i in range(0,4):
		if i != j:
			f = f + 48.0*((1/(np.linalg.norm(particles[j].position-particles[i].position)**14))-(1/(2*(np.linalg.norm(particles[j].position-particles[i].position)**8))))*Particle3D.vec_sep(particles[j].position,particles[i].position)
	force.append(f)
	print force[j]
"""
for i in range(numstep):
	point = i+1
	outfile.write(str(len(particles)) + "\nPoint = " + str(point))

	for j in range(len(particles)):
	    particles[j].leapPos2nd(dt, force)
	   
	    for y in range(len(particles) - 1):
			z = y + 1
			x = 48.0*((1/(np.linalg.norm(particles[y].position-particles[z].position)**14))-(1/(2*(np.linalg.norm(particles[y].position-particles[z].position)**8))))*Particle3D.vec_sep(particles[y].position,particles[z].position)	
			force_new += x
	    
	    		v = particles[y].leapVelocity(dt, 0.5*(force+force_new))

			force = force_new
	  
	    		pe = 4.0*((1/(np.linalg.norm(particles[y].position-particles[z].position)**12))-(1/(np.linalg.norm(particles[y].position-particles[z].position)**6)))
	    		e = pe + particles[y].kineticEnergy(v)
	
	time = time + dt
	posValue_y.append(particles[0].position[1])
	posValue_x.append(particles[0].position[0])
	kE.append(e)
	tValue.append(time)   		
	for i in range(len(particles)):
	    	outfile.write("\n" + str(particles[i].name) + str(point) + " " + str(particles[i].position[0]) + str(point) + " " + str(particles[i].position[1]) + str(point) + " " + str(particles[i].position[2]) + str(point)) 
	outfile.write("\n")
			

    


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
"""
