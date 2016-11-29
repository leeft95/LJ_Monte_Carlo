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
The second order positoin calcualtion method

leapVelocity: 
The velocity of the particle is updated using this method
    """

    @staticmethod
    def from_file(infile):
	line = infile.readline()
	tokens = line.split(",")
	name = str(tokens[0])
	pos = np.array([float(tokens[1]),float(tokens[2]),float(tokens[3])],float)
	vel = np.array([float(tokens[4]),float(tokens[5]),float(tokens[6])],float)
	mass = float(tokens[7])
	return Particle3D(pos,vel,mass,name)

    @staticmethod
    def vec_sep(p1,p2):
	return p1 - p2

    def __init__(self,pos, vel, mass,name):
        self.position = pos
	self.velocity = vel
	self.mass = mass
	self.name = name


    def __str__(self):
        return " Name = " + str(self.name) + "\n postion = " + str(self.position) + "\n velocity = " + str(self.velocity) + "\n mass = " + str(self.mass)
    
   
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
