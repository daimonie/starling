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
    
    boxSize = None
    habitatSize = None
    habitatStrength = None
    
    length = None
    
    display = 1
    def __init__(self):
        """Mostly empty, except for a joke."""
        print "We are boid."
        
        self.figure = plt.figure()
        self.axis = self.figure.gca(projection='3d')
    def initBoids(self):
        """Function that initialises the necessary matrices and such. """
        self.positions = np.random.rand(self.number, 3)* self.boxSize
        
        randomness = 0.75;
        
        self.length = self.boxSize*2.0/self.number
        self.velocities = np.random.rand(self.number,3) * randomness
        self.velocities[:, :] += 1 - randomness
        
        self.velocities /= np.sqrt( np.sum( np.square(self.velocities) ) )
        
        self.sensitivities = np.ones((self.number,1), dtype=np.float)
    def show(self):
        """Function that sets the animation running""" 
        self.ani = animation.FuncAnimation(self.figure, self.evolve, interval=self.display) 
        plt.show()
    def behaviour(self):
        """This returns the proper calculation given self.mode"""
        if self.mode == 0: #Simplest behaviour
            return ethology.simplebehaviour( positions = self.positions, velocities = self.velocities, tau = self.tau, eta=self.eta, sensitivities=self.sensitivities, number=self.number )
        elif self.mode == 1: #Simple + habitat
            return ethology.simplehabitat( positions = self.positions, velocities = self.velocities, tau = self.tau, eta=self.eta, sensitivities=self.sensitivities, number=self.number, habitatsize=self.habitatSize, habitatstrength=self.habitatStrength)
    def setAxes(self):
        """Proper axes given self.mode"""
        if self.mode == 0: 
            self.axis.set_xlim(np.min( newpositions[:,0]), np.max( newpositions[:,0]))
            self.axis.set_ylim(np.min( newpositions[:,1]), np.max( newpositions[:,1]))
            self.axis.set_zlim(np.min( newpositions[:,2]), np.max( newpositions[:,2]))
        elif self.mode == 1: 
            self.axis.set_xlim( -self.habitatSize, self.habitatSize )
            self.axis.set_ylim( -self.habitatSize, self.habitatSize )
            self.axis.set_zlim( -self.habitatSize, self.habitatSize )
    def evolve(self, r):
        """Function that calculates the next frame"""
        newpositions = self.positions
        diff = self.velocities
        
        newpositions, diff = self.behaviour() 
                
        self.axis.cla ()
        self.axis.quiver(newpositions[:,0], newpositions[:,1], newpositions[:,2], diff[:,0], diff[:,1], diff[:,2], length=self.length);
        
        self.setAxes()
        
        self.positions = newpositions
        self.velocities = diff * self.speed
        
        print "Evolve, generation %d ." % (r)