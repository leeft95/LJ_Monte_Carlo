import sys
import matplotlib.pyplot as pyplot
import numpy as np
import math as math
from Particle3D import Particle3D
import MDUtilities as MDUtilities
import copy


if len(sys.argv)!=4:
    print "Wrong number of arguments."
    print "Usage: " + sys.argv[0] + " <input file> " + " <output file 1> " + " <output file 2> "
    quit()
else:
    filename = sys.argv[1]    
    outfile1 = sys.argv[2]
    outfile2 = sys.argv[3]

outfile1_0 = open(outfile1, "w")
outfile2_0 = open(outfile2,"w")
infile = open(filename,"r")
line = infile.readline()
tokens = line.split(",")

n = int(tokens[0])
temp = float(tokens[1])
rho = float(tokens[2])
nAt = float(n)
rcutoff = float(tokens[3])
boxSize = (nAt/rho)**(1.0/3.0)
x_rsize = 1.0 /boxSize
y_rsize = 1.0 /boxSize
z_rsize = 1.0 /boxSize

pos = np.array([0.0,0.0,0.0],float)
vel = np.array([0.0,0.0,0.0],float)
mass = float(1)
particles = [Particle3D.base(pos,vel,mass) for i in range(n)]


MDUtilities.setInitialPositions(rho, particles)
MDUtilities.setInitialVelocities(temp, particles)





numstep = 10000
time = 0.0
dt = 0.00001
tlog = 100
pe = 0.0
e = 0.0
seppy = 0

particles1 = copy.deepcopy(particles)

tvalue = []
MSSD = []
posValue_x = [particles[0].position[0]]
posValue_y = [particles[0].position[1]]
kE = []
force = 0.0
force_new = 0.0
potential_new = 0.0
gg = 0
print particles[1]
for i in range(numstep):
	point = i+1
	if i%100 == 0:
		outfile1_0.write(str(n) + "\nPoint = " + str(point) + "\n")
	for j in range(len(particles)):
		for k in range(len(particles)):
			if j!= k:
				vecsep_x = abs(particles[j].position[0]-particles[k].position[0])
				vecsep_x = vecsep_x - int(vecsep_x*x_rsize + 0.5)*boxSize
				vecsep_y = abs(particles[j].position[1]-particles[k].position[1])
				vecsep_y = vecsep_y - int(vecsep_y*y_rsize + 0.5)*boxSize
				vecsep_z = abs(particles[j].position[2]-particles[k].position[2])
				vecsep_z = vecsep_z - int(vecsep_z*z_rsize + 0.5)*boxSize
				#print vecsep_x
				#print vecsep_y	
				#print vecsep_z	
				"""
				if (vecsep_x <=-0.5*boxSize): 
					vecsep_x = vecsep_x+boxSize
				if (vecsep_x > 0.5*boxSize):
					vecsep_x = vecsep_x-boxSize
				if (vecsep_y<=-0.5*boxSize): 
					vecsep_y = vecsep_y+boxSize
				if (vecsep_y> 0.5*boxSize):
					vecsep_y = vecsep_y-boxSize
				if (vecsep_z <=-0.5*boxSize): 
					vecsep_z = vecsep_z+boxSize
				if (vecsep_z > 0.5*boxSize):
					vecsep_z = vecsep_z-boxSize
				"""
				#print vecsep_x
				#print vecsep_y	
				#print vecsep_z	
				"""				
				pos_img = np.array((xpos_im,ypos_im,zpos_im),float)
				img = Particle3D(pos_img, particles[k].velocity, particles[k].mass)
				"""
				vec_sep = Particle3D.vec_sep(particles[j],particles[k])
				#print vec_sep
				vec_sepm = math.sqrt(sum(vec_sep*vec_sep))
				#print vec_sepm
				if vec_sepm<=rcutoff:
					pe_mod = Particle3D.Lj_pot(particles[j],vec_sepm)
					f_mod = Particle3D.Lj_force(particles[j],vec_sep,vec_sepm)	
					force = force + f_mod
					#print force
					pe = pe + pe_mod
		print '1'
		p1 = j+1
	    	particles[j].leapPos2nd(dt,force)
		while particles[j].position[0]>=boxSize*0.5:
			print 'a'
           		x_pos = particles[j].position[0]-boxSize
            		particles[j].position[0] = x_pos
			print 'a1'
        	while particles[j].position[0]<-boxSize*0.5:
			print 'b'
           		x_pos = particles[j].position[0]+boxSize
            		particles[j].position[0] = x_pos
			print 'b1'
		while particles[j].position[1]>=boxSize*0.5:
			print 'c'
           		y_pos = particles[j].position[1]-boxSize
            		particles[j].position[1] = y_pos
			print 'c1'
        	while particles[j].position[1]<-boxSize*0.5:
			print 'd'
           		y_pos = particles[j].position[1]+boxSize
            		particles[j].position[1] = y_pos
			print 'd1'
		while particles[j].position[2]>=boxSize*0.5:
			print 'e'
           		z_pos = particles[j].position[2]-boxSize
            		particles[j].position[2] = z_pos
			print 'e1'
        	while particles[j].position[2]<-boxSize*0.5:
			print 'f'
           		z_pos = particles[j].position[2]+boxSize
            		particles[j].position[2] = z_pos
			print 'f1'
		print '2'
		if i%100 == 0:
			bb = 0
			for aa in range(len(particles)):
				seppy = Particle3D.MsD(particles[aa],particles1[bb])
				sum_seppy = seppy + seppy
				bb = bb + 1			
	            	outfile1_0.write(str(p1) + " " + str(particles[j].position[0]) + " " + str(particles[j].position[1]) + " " + str(particles[j].position[2]) + "\n")
	
	if i%100 == 0:
		MsD = (1/nAt)*sum_seppy
		tvalue.append(gg)
		gg = gg + 1
		MSSD.append(MsD)
	for l in range(len(particles)): 
        	for m in range(len(particles)): 
		    	if l!=m: 
				x_sep= abs(particles[l].position[0]-particles[m].position[0])
				x_sep=  x_sep - int(x_sep*x_rsize + 0.5)*boxSize
				y_sep= abs(particles[l].position[1]-particles[m].position[1])
				y_sep=  y_sep - int(y_sep*y_rsize + 0.5)*boxSize
				z_sep= abs(particles[l].position[2]-particles[m].position[2])
				z_sep=  z_sep - int(z_sep*z_rsize + 0.5)*boxSize
				if (x_sep<=0.5*-boxSize): 
					x_sep = x_sep+boxSize
				if (x_sep> 0.5*boxSize):
					x_sep = x_sep-boxSize
				if (y_sep<=0.5*-boxSize): 
					y_sep = y_sep+boxSize
				if (y_sep> 0.5*boxSize):
					y_sep = y_sep-boxSize
				if (z_sep<=0.5*-boxSize): 
					z_sep = z_sep+boxSize
				if (z_sep> 0.5*boxSize):
					z_sep = z_sep-boxSize
				"""
				img_pos= np.array((xposim, yposim, zposim), float) 
				img = Particle3D(img_pos, particles[m].velocity, particles[m].mass)
				img_sep = Particle3D.vec_sep(particles[l], img)
				img_sqmag = math.sqrt(sum(img_sep*img_sep))
				"""
				vec_sep1 = Particle3D.vec_sep(particles[l],particles[m])
				vec_sepm1 = math.sqrt(sum(vec_sep*vec_sep))
				if vec_sepm<=rcutoff:
					pot_mod1 = Particle3D.Lj_pot(particles[j],vec_sepm1)
					f_mod1 = Particle3D.Lj_force(particles[j],vec_sep1,vec_sepm1)	
					force_new = force_new + f_mod1
					potential_new = pe + pot_mod1   
		v = particles[l].leapVelocity(dt, 0.5*(force+force_new))
		ke = Particle3D.kineticEnergy(particles[l])
		total_energy = potential_new + ke		
		if l%32 == 0:		
			outfile2_0.write(str(total_energy) + " te \n" + str(ke) +  " ke \n" + str(potential_new) + " pe\n" + str(MsD) + " MSD\n\n" + str(force) + " force \n\n")
		force_new = 0
		force = 0
		potential = 0
		time = time + dt
	


outfile1_0.close()
outfile2_0.close()
infile.close()
tValue = np.array(tvalue)
print tValue
print len(tValue)
MSD = np.array(MSSD)
print MSD
print len(MSD)
pyplot.figure()
pyplot.subplot(111)
pyplot.plot(tValue,MSD)
pyplot.title('MSD')
pyplot.xlabel('Time')
pyplot.ylabel('MSD')
pyplot.savefig('MSD.png')
"""
pyplot.figure()
pyplot.subplot(111)
pyplot.plot(tValue,kE)
pyplot.title('Total Energy of the particle against time')
pyplot.xlabel('Time(s)')
pyplot.ylabel('Energy(J)')
pyplot.savefig('EnergyVV.png')
"""
