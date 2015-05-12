from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import matplotlib.animation as animation
from fbehaviour import ethology
class Flock:
    number = 10
    positions = None
    velocities = None
    sensitivities = None
    
    figure = None
    axis = None
    ani = None
    def __init__(self):
        """Mostly empty, except for a joke."""
        print "We are boid."
        
        self.figure = plt.figure()
        self.axis = self.figure.gca(projection='3d')
    def initBoids(self):
        """Function that initialises the necessary matrices and such. """
        self.positions = np.random.rand(self.number, 3)
        self.velocities = np.random.rand(self.number,3) * 0.05
        self.velocities[:, 0] += 0.95
        
        self.velocities /= np.sqrt( np.sum( np.square(self.velocities) ) )
        
        self.sensitivities = np.ones((self.number,1), dtype=np.float)
    def show(self):
        """Function that sets the animation running"""
        self.ani = animation.FuncAnimation(self.figure, self.evolve,interval=100) 
        plt.show()
    def evolve(self, r):
        """Function that calculates the next frame"""
        newpositions = self.positions
        diff = self.velocities
        
        newpositions, diff = ethology.simplebehaviour( positions = self.positions, velocities = self.velocities, sensitivities=self.sensitivities, number=self.number )
        
        self.axis.cla ()
        self.axis.quiver(newpositions[:,0], newpositions[:,1], newpositions[:,2], diff[:,0], diff[:,1], diff[:,2], length=0.1);
        self.axis.set_xlim([-4.0, 4.0])
        self.axis.set_ylim([-4.0, 4.0])
        self.axis.set_zlim([-4.0, 4.0])
        
        self.positions = newpositions
        self.velocities = diff
        
        print "Evolve, generation %d ." % (r)