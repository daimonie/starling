from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import matplotlib.animation as animation
from fbehaviour import ethology
class Flock:
    number = 10
    tau = None
    positions = None
    velocities = None
    sensitivities = None
    eta = None
    
    speed = None
    figure = None
    axis = None
    ani = None
    
    mode = 0
    def __init__(self):
        """Mostly empty, except for a joke."""
        print "We are boid."
        
        self.figure = plt.figure()
        self.axis = self.figure.gca(projection='3d')
    def initBoids(self):
        """Function that initialises the necessary matrices and such. """
        self.positions = np.random.rand(self.number, 3)* self.boxsize
        
        randomness = 0.95;
        
        self.velocities = np.random.rand(self.number,3) * randomness
        self.velocities[:, 0] += 1 - randomness
        
        self.velocities /= np.sqrt( np.sum( np.square(self.velocities) ) )
        
        self.sensitivities = np.ones((self.number,1), dtype=np.float)
    def show(self):
        """Function that sets the animation running""" 
        self.ani = animation.FuncAnimation(self.figure, self.evolve,interval=100) 
        plt.show()
    def behaviour(self):
        if self.mode == 0:
            return ethology.simplebehaviour( positions = self.positions, velocities = self.velocities, tau = self.tau, eta=self.eta, sensitivities=self.sensitivities, number=self.number )
        
    def evolve(self, r):
        """Function that calculates the next frame"""
        newpositions = self.positions
        diff = self.velocities
        
        newpositions, diff = self.behaviour()
                
        self.axis.cla ()
        self.axis.quiver(newpositions[:,0], newpositions[:,1], newpositions[:,2], diff[:,0], diff[:,1], diff[:,2], length=self.boxsize*2.0/self.number);
        self.axis.set_xlim(np.min( newpositions[:,0]), np.max( newpositions[:,0]))
        self.axis.set_ylim(np.min( newpositions[:,1]), np.max( newpositions[:,1]))
        self.axis.set_zlim(np.min( newpositions[:,2]), np.max( newpositions[:,2]))
        
        self.positions = newpositions
        self.velocities = diff * self.speed
        
        print "Evolve, generation %d ." % (r)