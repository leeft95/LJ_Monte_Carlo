"""
 CMod Ex3: Particle3D, a class to describe 3D particles
"""
import math
import numpy as np
class Particle3D(object):

    @staticmethod
    def from_file(infile):
	line = infile.readline()
	tokens = line.split(",")
	name = str(tokens[0])
	pos = np.array([float(tokens[1]),float(tokens[2]),float(tokens[3])])
	
	vel = np.array([float(tokens[4]),float(tokens[5]),float(tokens[6])])
	
	mass = float(tokens[7])
	return Particle3D(pos,vel,mass,name)

    # Initialise a Particle3D instance
    def __init__(self,pos, vel, mass,name):
        self.position = pos
	self.velocity = vel
	self.mass = mass
	self.name = name

    # Formatted output as String
    def __str__(self):
        return "Name = " + str(self.name) + "\n postion = " + str(self.position) + "\n velocity = " + str(self.velocity) + "\n mass = " + str(self.mass)
    
    # Kinetic energy, which is the each components KE summed up
    def kineticEnergy(self):
	e = 0.0
	for i in range(0,3):
		b = float(self.velocity[i])
		print i
		print b
		e = 0.5*self.mass*(b**2)
		print e
		e += e
		
	return e

    # Time integration methods
    # First-order velocity update
    def leapVelocity(self, dt, force):
        self.velocity = self.velocity + dt*force/self.mass

    # First-order position update
    def leapPos1st(self, dt):
        self.position = self.position + dt*self.velocity

    # Second-order position update
    def leapPos2nd(self, dt, force):
        self.position = self.position + dt*self.velocity + 0.5*dt**2*force/self.mass
