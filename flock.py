from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import matplotlib.animation as animation
from fbehaviour import ethology
from matplotlib.ticker import NullFormatter
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
    axisOrder = None
    axisOrderPrime = None
    ani = None
    
    mode = 0
    
    boxSize = None
    habitatSize = None
    habitatStrength = None
    noSphere = True
    length = None
    
    display = 1
    order = 0.0
    orderPrime = 0.0
     
    i0 = None
    i1 = None
    i2 = None 
    i3 = None   
    i4 = None   
    i5 = None   
    def __init__(self):
        """Mostly empty, except for a joke."""
        print "We are boid."
        
        self.figure = plt.figure(figsize=(8,8), )
        
        left, width = 0.1, 0.65
        bottom, height = 0.1, 0.65
        bottom_h = left_h = left+width+0.02
        rect_sim = [left, bottom, width, height]
        rect_ord = [left, bottom_h, width, 0.2]
        rect_pri = [left_h, bottom, 0.2, height]
        
        
        #self.axis = self.figure.gca(projection='3d')
        self.axis = axes3d.Axes3D(self.figure, rect=rect_sim)
        self.axisOrder = plt.axes(rect_ord)
        self.axisOrderPrime = plt.axes(rect_pri)
        
    def initBoids(self):
        """Function that initialises the necessary matrices and such. """
        self.positions = 1.0*self.boxSize*np.ones((self.number,3)) - 2.0*np.random.rand(self.number, 3)* self.boxSize
                
        self.length = self.boxSize*2.0/self.number
        self.velocities = np.random.rand(self.number,3)*2 - 1 * np.ones((self.number,3), dtype=np.float)
        self.velocities /= np.sqrt( np.sum( np.square(self.velocities) ) )
        self.velocities *= self.speed
        
        self.sensitivities = np.ones((self.number,1), dtype=np.float)
        
        self.orderCalculation()
    def show(self):
        """Function that sets the animation running""" 
        self.ani = animation.FuncAnimation(self.figure, self.evolve, interval=self.display) 
        plt.show()
    def behaviour(self):
        """This returns the proper calculation given self.mode"""
        if self.mode == 0: #Simplest behaviour
            return ethology.simplebehaviour( positions = self.positions, velocities = self.velocities, tau = self.tau,
                eta=self.eta, sensitivities=self.sensitivities, number=self.number )
        elif self.mode == 1: #Simple + habitat
            return ethology.simplehabitat( positions = self.positions, velocities = self.velocities, tau = self.tau,
                eta=self.eta, sensitivities=self.sensitivities, number=self.number, habitatsize=self.habitatSize,
                habitatstrength=self.habitatStrength)
        elif self.mode == 2: #Simple + habitat + nteraction
            return ethology.interactionhabitat( positions = self.positions, velocities = self.velocities, tau = self.tau,
                eta=self.eta, sensitivities=self.sensitivities, number=self.number, habitatsize=self.habitatSize,
                habitatstrength=self.habitatStrength, i0=self.i0, i1=self.i1,i2=self.i2, i3=self.i3, i4=self.i4, i5=self.i5)
    def setAxes(self):
        """Proper axes given self.mode"""
        if self.mode == 0: 
            self.axis.set_xlim(np.min( self.positions[:,0]), np.max( self.positions[:,0]))
            self.axis.set_ylim(np.min( self.positions[:,1]), np.max( self.positions[:,1]))
            self.axis.set_zlim(np.min( self.positions[:,2]), np.max( self.positions[:,2])) 
        elif self.mode == 1 or self.mode == 2: 
            self.axis.set_xlim( -self.habitatSize, self.habitatSize )
            self.axis.set_ylim( -self.habitatSize, self.habitatSize )
            self.axis.set_zlim( -self.habitatSize, self.habitatSize )
            
            if self.noSphere == True:
                u = np.linspace(0, 2 * np.pi, 100)
                v = np.linspace(0, np.pi, 100)

                self.sphere_x = self.habitatSize * np.outer(np.cos(u), np.sin(v))
                self.sphere_y = self.habitatSize * np.outer(np.sin(u), np.sin(v))
                self.sphere_z = self.habitatSize * np.outer(np.ones(np.size(u)), np.cos(v))
                
                self.noSphere = False;
            self.axis.plot_wireframe(self.sphere_x, self.sphere_y, self.sphere_z,  rstride=13, cstride=13, color='r', alpha=0.3)
                 
    def orderCalculation(self):
        """Calculate the order parameter and its derivative"""
        vsum = np.sum( self.velocities, axis=0) 
        newOrder = 1.0/self.number/self.speed * np.sqrt( np.dot(vsum, vsum))
        
        self.orderPrime = ( newOrder - self.order) / self.tau
        self.order = newOrder 
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
        
        self.orderCalculation()
        
        self.axisOrder.plot(r*self.tau, self.order)
        
        
        
        print "%d\t%2.3e\t%2.3e." % (r, self.order, self.orderPrime)