"""
 CMod Ex3: Particle3D, a class to describe 3D particles
"""
import math
import numpy as np
class Particle3D(object):
    """
from_file:
This static method inputs data from an input file called particle.in where each value of the file corresponds to an attribute of the particle

vec_sep:      
This static method returns the vector seperation of two particles

__init__:
From the static method the attributes of the particles are assigned and initialized

__str__:
Formated output of the particle attributes

kineticEnergy:
The kinetic energy of the particle is updated

leapPos1st:
The first order positon calculation method

leapPos2nd:
The second order position calcualtion method

leapVelocity: 
The velocity of the particle is updated using this method
    
"""
    @staticmethod
    def base(pos,vel,mass):
	position = pos
	velocity = vel
	mass = mass
	return Particle3D(pos,vel,mass)

    
    def Lj_pot(particles,vec_sqmag):
	pot = 4*(1/((vec_sqmag)**12)-1/((vec_sqmag)**6))
	return pot

    def Lj_force(particles,vec_sep,vec_sqmag):
	force = 48.0*(1/((vec_sqmag)**14)-1/(2*(vec_sqmag)**8))*vec_sep
  	return force
	
	

    def MsD(particles,particles1):
	vec_sep = particles.position-particles1.position
	sqmag = math.sqrt(sum(vec_sep*vec_sep))
	return sqmag
  
    @staticmethod
    def vec_sep(p1,p2,boxSize):
	vecsep_x = p1.position[0]-p2.position[0]
				
	vecsep_y = p1.position[1]-p2.position[1]
				
	vecsep_z = p1.position[2]-p2.position[2]
				
				
	if (vecsep_x <= -0.5*boxSize): 
		vecsep_x = vecsep_x+boxSize
	if (vecsep_x > 0.5*boxSize):
		vecsep_x = vecsep_x-boxSize
	if (vecsep_y<= -0.5*boxSize): 
		vecsep_y = vecsep_y+boxSize
	if (vecsep_y> 0.5*boxSize):
		vecsep_y = vecsep_y-boxSize
	if (vecsep_z <= -0.5*boxSize): 
		vecsep_z = vecsep_z+boxSize
	if (vecsep_z > 0.5*boxSize):
		vecsep_z = vecsep_z-boxSize
								
	vec_sep = np.array([vecsep_x,vecsep_y,vecsep_z])
	return vec_sep


    @staticmethod
    def PBC(particles,boxSize):
    	if particles.position[0]>=boxSize:
          	x_pos = particles.position[0]-boxSize
            	particles.position[0] = x_pos
			
        if particles.position[0]<0.0:
           	x_pos = particles.position[0]+boxSize
            	particles.position[0] = x_pos
			
	if particles.position[1]>=boxSize:
           	y_pos = particles.position[1]-boxSize
            	particles.position[1] = y_pos
			
       	if particles.position[1]<0.0:
           	y_pos = particles.position[1]+boxSize
            	particles.position[1] = y_pos
			
	if particles.position[2]>=boxSize:	
           	z_pos = particles.position[2]-boxSize
            	particles.position[2] = z_pos
			
        if particles.position[2]<0.0:
           	z_pos = particles.position[2]+boxSize
            	particles.position[2] = z_pos
	return particles.position
    



    def __init__(self,pos, vel, mass):
        self.position = pos
	self.velocity = vel
	self.mass = mass
	


    def __str__(self):
        return  " postion = " + str(self.position) + "\n velocity = " + str(self.velocity) + "\n mass = " + str(self.mass)
    
   
    def kineticEnergy(particles):
	e = 0.0
	for i in range(0,3):
		b = float(particles.velocity[i])
		c = 0.5*particles.mass*b**2
		e = e + c	
	return e


    def leapVelocity(self, dt, force):
	self.velocity = self.velocity + dt*force/self.mass
	return self.velocity

    def leapPos1st(self, dt):
        self.position = self.position + dt*self.velocity

 
    def leapPos2nd(self, dt, force):
        self.position = self.position + dt*self.velocity + 0.5*dt**2*force/self.mass
	return self.position
