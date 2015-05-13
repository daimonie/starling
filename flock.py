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
    lineOrder = None
    lineOrderPrime = None
    timeData = None
    orderData = None
    primeData = None
    
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
    predatorSense = None
    predatorStrength = None
    predatorLocation = None
    sharkPath    = None
    predatorNumber = 10
    
    scatterPlot = False
    def __init__(self):
        """Mostly empty, except for a joke."""
        print "We are boid."
        
        self.figure = plt.figure(figsize=(8,8), )
        
        left, width = 0.1, 0.65
        bottom, height = 0.1, 0.65
        bottom_h = left_h = left+width+0.05
        rect_sim = [left, bottom, width, height]
        rect_ord = [left, bottom_h, width, 0.2]
        rect_pri = [left_h, bottom, 0.2, height]
        
        
        #self.axis = self.figure.gca(projection='3d')
        self.axis = axes3d.Axes3D(self.figure, rect=rect_sim)

        self.axisOrder = plt.axes(rect_ord)
        self.axisOrderPrime = plt.axes(rect_pri)

        #line, = ax.plot([], [], lw=2)        
        self.lineOrder, = self.axisOrder.plot([],[],lw=2)
        self.lineOrderPrime, = self.axisOrderPrime.plot([],[],lw=2)
        
        self.timeData = []
        self.orderData = []
        self.primeData = []
        
        self.axisOrder.set_xlim( 0, 10.0)
        self.axisOrder.set_ylim( 0, 1.0)
        
        self.axisOrder.grid()
        
        self.axisOrderPrime.set_xlim( 0, 10.0)
        self.axisOrderPrime.set_ylim( 0, 1.0)
         
        self.axisOrderPrime.grid()
        self.axisOrder.set_title("Order parameter")
        self.axisOrder.set_xlabel("Time")
        self.axisOrder.set_ylabel("psi")
        
        
        self.axisOrderPrime.set_title("Order parameter derivative")
        self.axisOrderPrime.set_xlabel("Time")
        self.axisOrderPrime.set_ylabel("psi'")
        
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
        self.axisOrderPrime.set_ylim( 0, 1.0/self.tau)
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
        elif self.mode == 3: #Simple + habitat + nteraction + SHARK
            return ethology.shark( positions = self.positions, velocities = self.velocities, tau = self.tau,
                eta=self.eta, sensitivities=self.sensitivities, number=self.number, habitatsize=self.habitatSize,
                habitatstrength=self.habitatStrength, i0=self.i0, i1=self.i1,i2=self.i2, i3=self.i3, i4=self.i4, i5=self.i5,
                predatorsense=self.predatorSense, predatorstrength=self.predatorStrength, predatorlocation = self.predatorLocation, predatornumber=self.predatorNumber)
        
        #predatorSense, predatorStrength, predatorLocation)
    def setAxes(self):
        """Proper axes given self.mode"""
        if self.mode == 0: 
            self.axis.set_xlim(np.min( self.positions[:,0]), np.max( self.positions[:,0]))
            self.axis.set_ylim(np.min( self.positions[:,1]), np.max( self.positions[:,1]))
            self.axis.set_zlim(np.min( self.positions[:,2]), np.max( self.positions[:,2])) 
        elif self.mode == 1 or self.mode == 2 or self.mode == 3: 
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
        if self.mode == 3:
            #Draw shark
            self.axis.scatter( self.predatorLocation[:,0], self.predatorLocation[:,1], self.predatorLocation[:,2], color='r')
    def orderCalculation(self):
        """Calculate the order parameter and its derivative"""
        vsum = np.sum( self.velocities, axis=0) 
        newOrder = 1.0/self.number/self.speed * np.sqrt( np.dot(vsum, vsum))
        
        self.orderPrime = ( newOrder - self.order) / self.tau
        self.order = newOrder 
    def initShark(self):
        """Initialises THE SHARK"""
        self.predatorNumber = 100
        self.predatorLocation = np.ones((self.predatorNumber, 3)) * 3.0 * self.habitatSize; 
        self.predatorSense = np.ones((self.predatorNumber,1))
        self.predatorStrength = np.ones((self.predatorNumber,1))
        
        self.predatorLocation[0, :] = [-2.50, 0.00, 0.00]
        self.predatorLocation[1, :] = [-2.00, 0.00, 0.00]
        self.predatorLocation[2, :] = [-1.50, 0.00, 0.00]
        self.predatorLocation[3, :] = [-1.00, 0.00, 0.00]
        self.predatorLocation[4, :] = [-0.50, 0.00, 0.00]
        self.predatorLocation[5, :] = [0.00, 0.00, 0.00]
        self.predatorLocation[6, :] = [0.50, 0.00, 0.00]
        self.predatorLocation[7, :] = [1.00, 0.00, 0.00]
        self.predatorLocation[8, :] = [1.50, 0.00, 0.00]
        self.predatorLocation[9, :] = [2.00, 0.00, 0.00]
        self.predatorLocation[10, :] = [2.50, 0.00, 0.00]
        
        
        self.predatorLocation[11, :] = [0.50, 0.30, 0.00]
        self.predatorLocation[12, :] = [1.00, 0.30, 0.00]
        self.predatorLocation[13, :] = [1.50, 0.30, 0.00]
        self.predatorLocation[14, :] = [2.00, 0.30, 0.00]
        self.predatorLocation[15, :] = [0.50, -0.30, 0.00]
        self.predatorLocation[16, :] = [1.00, -0.30, 0.00]
        self.predatorLocation[17, :] = [1.50, -0.30, 0.00]
        self.predatorLocation[18, :] = [2.00, -0.30, 0.00]
        
        
        self.predatorLocation[19, :] = [-2.50, 0.00, 0.50]
        self.predatorLocation[20, :] = [-2.00, 0.00, 0.50]
        self.predatorLocation[21, :] = [-1.50, 0.00, 0.50]
        self.predatorLocation[22, :] = [-1.00, 0.00, 0.50]
        self.predatorLocation[23, :] = [-0.50, 0.00, 0.50]
        self.predatorLocation[24, :] = [0.00, 0.00, 0.50]
        self.predatorLocation[25, :] = [0.50, 0.00, 0.50]
        self.predatorLocation[26, :] = [1.00, 0.00, 0.50]
        self.predatorLocation[27, :] = [1.50, 0.00, 0.50]
        self.predatorLocation[28, :] = [2.00, 0.00, 0.50]
        self.predatorLocation[29, :] = [2.50, 0.00, 0.50]
        
        
        self.predatorLocation[30, :] = [0.50, 0.30, 0.30]
        self.predatorLocation[31, :] = [1.00, 0.30, 0.30]
        self.predatorLocation[32, :] = [1.50, 0.30, 0.30]
        self.predatorLocation[33, :] = [2.00, 0.30, 0.30]
        self.predatorLocation[34, :] = [0.50, -0.30, 0.30]
        self.predatorLocation[35, :] = [1.00, -0.30, 0.30]
        self.predatorLocation[36, :] = [1.50, -0.30, 0.30]
        self.predatorLocation[37, :] = [2.00, -0.30, 0.30]
        
        
        self.predatorLocation[38, :] = [-2.50, -0.00, 0.70]
        self.predatorLocation[39, :] = [-2.50, -0.00, 0.90]
         
        self.predatorLocation[40, :] = [-2.70, -0.00, 0.70]
        self.predatorLocation[41, :] = [-2.70, -0.00, 0.90]
        
        self.predatorLocation[42, :] = [-2.50, -0.00, 0.50] 
        self.predatorLocation[43, :] = [-2.70, -0.00, 0.50] 
        self.predatorLocation[46, :] = [-2.50, -0.00, 0.10] 
        self.predatorLocation[47, :] = [-2.70, -0.00, 0.10]
        self.predatorLocation[48, :] = [-2.50, -0.00, -0.10] 
        self.predatorLocation[49, :] = [-2.70, -0.00, -0.10] 
        self.predatorLocation[52, :] = [-2.50, -0.00, -0.50] 
        self.predatorLocation[53, :] = [-2.70, -0.00, -0.50] 
        self.predatorLocation[54, :] = [-2.70, -0.00, -0.90]
        self.predatorLocation[55, :] = [-2.90, -0.00, -0.90]
        
        self.predatorLocation[56, :] = [0.70, -0.00, 1.00]
        self.predatorLocation[57, :] = [0.70, -0.00, 1.50]
        self.predatorLocation[57, :] = [0.50, -0.00, 1.50]

        self.predatorLocation[58, :] = [-2.90, 0.00, 1.10]
        
        self.predatorLocation *= 10.00;
        
    def evolve(self, r):
        """Function that calculates the next frame""" 
        newpositions, diff = self.behaviour() 
                
        self.axis.cla ()
        if self.scatterPlot == True:
            self.axis.scatter(newpositions[:,0], newpositions[:,1], newpositions[:,2], color='b');
        else:
            self.axis.quiver(newpositions[:,0], newpositions[:,1], newpositions[:,2], diff[:,0], diff[:,1], diff[:,2], length=self.length);
        
        self.setAxes()
        
        self.positions = newpositions
        self.velocities = diff * self.speed
        
        self.orderCalculation()
        #Plot order data
        
        self.timeData.append(r*self.tau)
        self.orderData.append(self.order)
        self.primeData.append(self.orderPrime)
        
        self.axisOrder.set_xlim(0, int(1+r*self.tau/10.00) * 10.00)
        self.axisOrderPrime.set_xlim(0, int(1+r*self.tau/10.00) * 10.00)
        
        self.lineOrder.set_data( self.timeData, self.orderData )
        self.lineOrderPrime.set_data( self.timeData, self.primeData )
  
        print "%d\t%2.3e\t%2.3e." % (r, self.order, self.orderPrime)