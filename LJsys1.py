"""Importing necessary packages and classes
"""
import sys
import matplotlib.pyplot as pyplot
import numpy as np
import math as math
from Particle3D import Particle3D
import MDUtilities as MDUtilities
import copy

""" Ensure user inputs correct number of arguments
"""
if len(sys.argv)!=4:
    print "Wrong number of arguments."
    print "Usage: " + sys.argv[0] + " <input file> " + " <output file 1> " + " <output file 2> "
    quit()
else:
    filename = sys.argv[1]    
    outfile1 = sys.argv[2]
    outfile2 = sys.argv[3]
"""Open files to write on and read from
"""
outfile1_0 = open(outfile1, "w")
outfile2_0 = open(outfile2,"w")
infile = open(filename,"r")
line = infile.readline()
tokens = line.split(",")
"""defining values for variables:
Where 	rho= number density of particles
	temp= temperature
	n= number of particles in box
	rcutoff= Cut off radius beyond which force between particles is not calculated (negligible)
	nAt= number of particles converted into float form
	boxsize= length of each side of the box
"""
n = int(tokens[0])
temp = float(tokens[1])
rho = float(tokens[2])
nAt = float(n)
rcutoff = float(tokens[3])
boxSize = (nAt/rho)**(1.0/3.0)
dr = 0.1
maxr = 3
"""set 3D arrays for position and velocity of particles
where 	pos= numpy array for position
	vel= numpy array for velocity
	mass= mass of particle
	particles"""
pos = np.array([0.0,0.0,0.0],float)
vel = np.array([0.0,0.0,0.0],float)
mass = float(1)
particles = [Particle3D.base(pos,vel,mass) for i in range(n)]

"""Use MDUtilities class to set initial positions and velocities for particles"""
MDUtilities.setInitialPositions(rho, particles)
MDUtilities.setInitialVelocities(temp, particles)
#attempt at function to histogram Radial Distribution Function
"""
def RDF(particles,n,maxr,dr,count,rho):
	i = 0	
	j = 0	
	rmin = 0
	rmax = dr
	ci = 0
	while i < n:
		for j in range (n):
			rmin = 0
			rmax = dr
			ci = 0	
			if i!= j:
				vec_sep = Particle3D.vec_sep(particles[i],particles[j],boxSize)
				sep_mag = np.linalg.norm(vec_sep)
				while rmax < maxr:
					if sep_mag > rmin and sep_mag < rmax:
						#print len(count)
						count[ci] = count[ci] + (1.0)/(n*4*rho*math.pi*dr*(rmin+(dr/2))**2)
						rmax = maxr + 1
					else:
						ci = ci +1
						rmin = rmin + dr
						rmax = rmax + dr
		i = i+1
	return count

"""
"""
	Define the timescale for the plot
	numstep= number of steps in the program
	time= time (initial value set to zero)
	dt= increment by which time is increased for each step
	tlog= number of time steps before position vector is updated
	pe= potential energy (initial value set to zero)
	e= total energy (initial value set to zero)
"""



numstep = 10000
time = 0.0
dt = 0.001
tlog = 100
pe = 0.0
e = 0.0
seppy = 0

particles1 = copy.deepcopy(particles)
""" 
	Creating lists for each numstep
	tvalue= list of time values of each step
	
	kE = list of kinetic energy values for each step
	force = force between particles before positional update
	force_new = between particles after positional update
	tE = 
"""
tvalue = []
MSSD = []
rds_list = []
count = np.array([30])
posValue_x = [particles[0].position[0]]
posValue_y = [particles[0].position[1]]
tE = []
kE = []
pel = []
force = 0.0
force_new = 0.0
potential_new = 0.0
gg = 0
for i in range(n):
	print particles[i]

#main time integration loop 
for i in range(numstep):
	point = i+1

	#to limit output file size readings will be out put only every 100 timesteps, embedded for loops to make sure that the same particle is not counted twice, and the interaction between each particle is taken into account

	if i%100 == 0:
		outfile1_0.write(str(n) + "\nPoint = " + str(point) + "\n")
	for j in range(len(particles)):
		for k in range(len(particles)):
			if j!= k:
				#minimum image convention is obeyed everytime a vector calculation is done
				vec_sep = Particle3D.vec_sep(particles[j],particles[k],boxSize)
				
				vec_sepm = np.linalg.norm(vec_sep)
				#force and potential before position update is computed 				
				if vec_sepm<=rcutoff:
					pe_mod = Particle3D.Lj_pot(particles[j],vec_sepm)
					f_mod = Particle3D.Lj_force(particles[j],vec_sep,vec_sepm)	
					force = force + f_mod
					
					pe = pe + pe_mod
		#position update using the velocity verlet algorithim is done		
		p1 = j+1
	    	particles[j].leapPos2nd(dt,force)
		particles[j].position = Particle3D.PBC(particles[j],boxSize)	
		#As the particle position is calcualted the mean square displacement is calculated and put into a array to be plotted, data calcualted every 100 timesteps		
		if i%100 == 0:
			bb = 0
			for aa in range(len(particles)):
				seppy = Particle3D.MsD(particles[aa],particles1[bb])			
				sum_seppy = seppy + seppy
				bb = bb + 1			
			#data is written and outputted for the vmd program and is formated and outputted			
	            	outfile1_0.write(str(p1) + " " + str(particles[j].position[0]) + " " + str(particles[j].position[1]) + " " + str(particles[j].position[2]) + "\n")
	
	if i%100 == 0:
		MsD = (1/nAt)*sum_seppy
		tvalue.append(gg)
		gg = gg + 1
		MSSD.append(MsD)
		#count = copy.deepcopy(RDF(particles,n,maxr,dr,count,rho))
		#print count
		#The forces and potential are calcualted for all the particles after the position update, this is then used to calucalte the new velocity of the particles, which has then be used to determine the Ke
	for l in range(len(particles)): 
        	for m in range(len(particles)): 
		    	if l!=m: 
				vec_sep1 = Particle3D.vec_sep(particles[l],particles[m],boxSize)
				vec_sepm1 = np.linalg.norm(vec_sep1)

				if vec_sepm1<=rcutoff:
					pot_mod1 = Particle3D.Lj_pot(particles[j],vec_sepm1)
					f_mod1 = Particle3D.Lj_force(particles[j],vec_sep1,vec_sepm1)	
					force_new = force_new + f_mod1
					potential_new = pe + pot_mod1 
					peo = potential_new/point  
		v = particles[l].leapVelocity(dt, 0.5*(force+force_new))
		ke = Particle3D.kineticEnergy(particles[l])
		total_energy = potential_new/point + ke	
		#the data is dumped to an outputfile for all the data calucalted 
		if l%100 == 0:		
			outfile2_0.write(str(total_energy) + " te \n" + str(ke) +  " ke \n" + str(peo) + " pe\n" + str(MsD) + " MSD\n\n" + str(force) + " force \n\n")
		force_new = 0
		force = 0
		potential = 0
		time = time + dt
	#lists used is drawing the energy graphs is compiled every 100 timesteps
	if i%100 == 0:		
		tE.append(total_energy)
		pel.append(peo)
		kE.append(ke)

#output and input files closed
outfile1_0.close()
outfile2_0.close()
infile.close()
tValue = np.array(tvalue)
#All the graphs are drawn and outputed with axis labled
MSD = np.array(MSSD)
pyplot.figure()
pyplot.subplot(111)
pyplot.plot(tValue,MSD)
pyplot.title('MSD')
pyplot.xlabel('Timestep')
pyplot.ylabel('MSD')
pyplot.savefig('MSD.png')
pyplot.figure()
pyplot.subplot(111)
pyplot.plot(tValue,tE)
pyplot.title('Total Energy of the particle against timestep')
pyplot.xlabel('Timestep')
pyplot.ylabel('Energy(J)')
pyplot.savefig('EnergyT.png')
pyplot.figure()
pyplot.subplot(111)
pyplot.plot(tValue,pel)
pyplot.title('Potential Energy of the particle against timestep')
pyplot.xlabel('Timestep')
pyplot.ylabel('Energy(J)')
pyplot.savefig('EnergyP.png')
pyplot.figure()
pyplot.subplot(111)
pyplot.plot(tValue,kE)
pyplot.title('Kinetic Energy of the particle against timestep')
pyplot.xlabel('Timestep')
pyplot.ylabel('Energy(J)')
pyplot.savefig('EnergykE.png')

