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
	name = [str(tokens[0])]
	x_pos = [float(tokens[1])]
	y_pos = [float(tokens[2])]
	z_pos = [float(tokens[3])]
	x_vel = [float(tokens[4])]
	y_vel = [float(tokens[5])]
	z_vel = [float(tokens[6])]
	mass = [float(tokens[7])]
	return Particle3D(x_pos,x_vel,y_pos,y_vel,z_pos,z_vel,mass,name)

    # Initialise a Particle3D instance
    def __init__(self, x_pos, x_vel, y_pos, y_vel, z_pos, z_vel, mass,name):
        self.position = np.array([x_pos,y_pos,z_pos],float)
	self.velocity = np.array([x_vel,y_vel,z_vel],float)
	self.mass = mass
	self.name = name

    # Formatted output as String
    def __str__(self):
        return "Name = " + str(self.name) + "\n postion = " + str(self.position) + "\n velocity = " + str(self.velocity) + "\n mass = " + str(self.mass)
    
    # Kinetic energy, which is the each components KE summed up
    def kineticEnergy(self):
	e = 0.0
	step = 2
	for i in range (step):
		b = float(self.velocity[i])
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
