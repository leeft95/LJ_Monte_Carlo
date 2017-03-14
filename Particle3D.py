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

    
    def Lj_pot(particles,img):
	vec_sep = Particle3D.vec_sep(particles,img)
	
	vec_sqmag = math.sqrt(sum(vec_sep*vec_sep))
	print vec_sqmag
	pot = 4*((1/(vec_sqmag)**12)-(1/(vec_sqmag)**6))
	return pot

    def Lj_force(particles,img):
	vec_sep = Particle3D.vec_sep(particles,img)
	vec_sqmag = math.sqrt(sum(vec_sep*vec_sep))
	force = 48.0*((1/(vec_sqmag)**14))-(1/(2*(vec_sqmag)**8))*Particle3D.vec_sep(particles,img)
  	return force
  
    @staticmethod
    def vec_sep(p1,p2):
	return p1.position - p2.position

    def __init__(self,pos, vel, mass):
        self.position = pos
	self.velocity = vel
	self.mass = mass
	


    def __str__(self):
        return  " postion = " + str(self.position) + "\n velocity = " + str(self.velocity) + "\n mass = " + str(self.mass)
    
   
    def kineticEnergy(self, velocity):
	e = 0.0
	for i in range(0,3):
		b = float(velocity[i])
		c = 0.5*self.mass*b**2
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
