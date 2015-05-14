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
    #Torque parameters:  
    orientations = None   
    viscosity = None       
    torqueCutOff = None     
    alignWithNeighbors = None  
    seeTheFlock1 = None    
    seeTheFlock2 = None     
    seeTheFlock3 = None     
    
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
    predatorNumber = 10
    predatorLocationPrevious = None
    
    sharkSpeed = None
    sharkOmega = None
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
        
        self.orientations = self.velocities / self.speed 
        
    def show(self):
        """Function that sets the animation running""" 
        self.axisOrderPrime.set_ylim( 0, 1.0/self.tau)
        self.ani = animation.FuncAnimation(self.figure, self.evolve, interval=self.display) 
        plt.show()
    def behaviour(self):
        """This returns the proper calculation given self.mode"""
        if self.mode == 0: #Simplest behaviour
            newpositions, diff = ethology.simplebehaviour( positions = self.positions, velocities = self.velocities, tau = self.tau,
                eta=self.eta, sensitivities=self.sensitivities, number=self.number )
            newpositions -= np.sum( newpositions, axis=0)/self.number
            
            return newpositions, diff
        elif self.mode == 1: #Simple + habitat
            return ethology.bowlhabitat( positions = self.positions, velocities = self.velocities, tau = self.tau,
                eta=self.eta, sensitivities=self.sensitivities, number=self.number, habitatsize=self.habitatSize,
                habitatstrength=self.habitatStrength)
        elif self.mode == 2: #Simple + habitat + nteraction
            return ethology.interactionbowlhabitat( positions = self.positions, velocities = self.velocities, tau = self.tau,
                eta=self.eta, sensitivities=self.sensitivities, number=self.number, habitatsize=self.habitatSize,
                habitatstrength=self.habitatStrength, i0=self.i0, i1=self.i1,i2=self.i2, i3=self.i3, i4=self.i4, i5=self.i5)
        elif self.mode == 3: #Simple + habitat + nteraction + SHARK 
            return ethology.sharkbowl( positions = self.positions, velocities = self.velocities, tau = self.tau,
                eta=self.eta, sensitivities=self.sensitivities, number=self.number, habitatsize=self.habitatSize,
                habitatstrength=self.habitatStrength, i0=self.i0, i1=self.i1,i2=self.i2, i3=self.i3, i4=self.i4, i5=self.i5,
                predatorsense=self.predatorSense, predatorstrength=self.predatorStrength, predatorlocation = self.predatorLocation, predatornumber=self.predatorNumber)
        
        #predatorSense, predatorStrength, predatorLocation)
    def calculateRotation(self):        
        '''Calculates the new orientations of the individuals'''
        neworientations = ethology.simplerotation(number = self.number, viscosity = self.viscosity, cutoff = self.torqueCutOff, 
                                                  tau = self.tau, tboundary1 = self.seeTheFlock1, tboundary2 = self.seeTheFlock2, 
                                                  tboundary3 = self.seeTheFlock3, talign = self.alignWithNeighbors, 
                                                  positions = self.positions, orientations = self.orientations)
        return neworientations
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
            
    def orderCalculation(self):
        """Calculate the order parameter and its derivative"""
        vsum = np.sum( self.velocities, axis=0) 
        newOrder = 1.0/self.number/self.speed * np.sqrt( np.dot(vsum, vsum))
        
        self.orderPrime = ( newOrder - self.order) / self.tau
        self.order = newOrder 
    def updateShark(self, r):
        """The idea is that the shark can move around and orient itself.""" 
        if self.predatorLocationPrevious == None:
            self.predatorLocationPrevious = np.sum(self.predatorLocation)/self.predatorNumber
        if r < 2:
            return;
        radius = 0.00
        
        if r > self.sharkSpeed: 
            radius = self.habitatSize * 0.75  
            #Recall velocities go at speed omega * r, so total shark speed should be (omega theta^2 + omega phi^2
            
            
            theta   = np.pi * 2 / 100.00 * (r-self.sharkSpeed) * self.tau * self.sharkOmega
            phi     = np.pi/2.0 + np.pi / 200.00 * (r-self.sharkSpeed) * self.tau * self.sharkOmega
            #print "ANGLE\t%2.3f\t%2.3f" % (theta, phi)
            sharkCentre = np.array([ radius * np.cos(theta) * np.sin(phi), radius * np.sin(theta) * np.sin(phi), radius * np.cos(phi) ])
        else:
            sharkCentre = np.array([r / self.sharkSpeed * self.habitatSize*0.75, 0, 0])
        
        sharkDiff = sharkCentre - self.predatorLocationPrevious   
        
        #Arrow from shark belly to nose
        sharkDirection = ( self.predatorLocation[60,:] - self.predatorLocation[5,:])
        sharkDirection /= np.sqrt( np.sum( np.square( sharkDirection )))
        
        self.predatorLocation += sharkDiff
        
        velocityCentre = sharkCentre - self.predatorLocationPrevious
        velocityCentreLength = np.sqrt( np.sum( np.square( velocityCentre )))
        if velocityCentreLength > 0:
            velocityCentre /= velocityCentreLength
            #Velocity direction is now normalised. Note that if it is zero, there is no change in direction
            #and thus no need for rotation. 
            axis = np.cross(sharkDirection, velocityCentre)
            axisLength =  np.sqrt(np.sum( np.square( axis)))
            if axisLength > 0:
                axis /= axisLength
                theta = np.arccos( np.inner(sharkDirection, velocityCentre))
                print "Rotating over [%2.3e]" % theta
                self.predatorLocation = ethology.rotatepoints(number = self.predatorNumber, points = self.predatorLocation, axis=axis, theta=theta)
        self.predatorLocationPrevious = sharkCentre
         
    def evolveDraw(self, r):
        """ This can contain triggers for things to be drawn, e.g. the shark."""
        if self.mode == 3:
            #Draw shark
            self.updateShark(r)
             
            self.axis.scatter( self.predatorLocation[:,0], self.predatorLocation[:,1], self.predatorLocation[:,2], color='r', s=4*self.length)
        
            #self.axis.set_xlim(np.min( self.predatorLocation[:,0]), np.max( self.predatorLocation[:,0]))
            #self.axis.set_ylim(np.min( self.predatorLocation[:,1]), np.max( self.predatorLocation[:,1]))
            #self.axis.set_zlim(np.min( self.predatorLocation[:,2]), np.max( self.predatorLocation[:,2])) 
        
        if self.mode == 1 or self.mode == 2 or self.mode == 3:  
            
            if self.noSphere == True:
                u = np.linspace(0, 2 * np.pi, 100)
                v = np.linspace(0, np.pi, 100)

                self.sphere_x = self.habitatSize * np.outer(np.cos(u), np.sin(v))
                self.sphere_y = self.habitatSize * np.outer(np.sin(u), np.sin(v))
                self.sphere_z = self.habitatSize * np.outer(np.ones(np.size(u)), np.cos(v))
                
                self.noSphere = False;
            self.axis.plot_wireframe(self.sphere_x, self.sphere_y, self.sphere_z,  rstride=13, cstride=13, color='r', alpha=0.3)
    def evolve(self, r):
        """Function that calculates the next frame""" 
        newpositions, diff = self.behaviour() 
        neworientations = self.calculateRotation()
                
        self.axis.cla ()
        if self.scatterPlot == True:
            self.axis.scatter(newpositions[:,0], newpositions[:,1], newpositions[:,2], color='b', s=self.length);
        else:
            self.axis.quiver(newpositions[:,0], newpositions[:,1], newpositions[:,2], neworientations[:,0], neworientations[:,1], neworientations[:,2], length=self.length);
        
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
  
        self.evolveDraw(r)
        print "%d\t%2.3e\t%2.3e." % (r, self.order, self.orderPrime)
    
    def initShark(self):
        """Initialises THE SHARK"""
        self.predatorNumber = 61
        self.predatorLocation = np.ones((self.predatorNumber, 3)); 
        self.predatorSense = np.ones((self.predatorNumber,1))
        self.predatorStrength = np.ones((self.predatorNumber,1)) 
        
        pixelNumber = 0
        #One of those moments where you miss ++pixelNumber :(
        self.predatorLocation[pixelNumber, :] = [-2.50, 0.00, 0.00] 
        pixelNumber += 1 
        
        self.predatorLocation[pixelNumber, :] = [-2.00, 0.00, 0.00]
        pixelNumber += 1 
        
        self.predatorLocation[pixelNumber, :] = [-1.50, 0.00, 0.00]
        pixelNumber += 1 
        
        self.predatorLocation[pixelNumber, :] = [-1.00, 0.00, 0.00]
        pixelNumber += 1 
        
        self.predatorLocation[pixelNumber, :] = [-0.50, 0.00, 0.00]
        pixelNumber += 1 
        
        self.predatorLocation[pixelNumber, :] = [0.00, 0.00, 0.00] 
        pixelNumber += 1 
        
        self.predatorLocation[pixelNumber, :] = [0.50, 0.00, 0.00]
        pixelNumber += 1 
        
        self.predatorLocation[pixelNumber, :] = [1.00, 0.00, 0.00]
        pixelNumber += 1 
        
        self.predatorLocation[pixelNumber, :] = [1.50, 0.00, 0.00]
        pixelNumber += 1 
        
        self.predatorLocation[pixelNumber, :] = [2.00, 0.00, 0.00]
        pixelNumber += 1 
        
        self.predatorLocation[pixelNumber, :] = [2.50, 0.00, 0.00]
        pixelNumber += 1 
                     
        self.predatorLocation[pixelNumber, :] = [0.50, 0.10, 0.00]
        pixelNumber += 1 
        
        self.predatorLocation[pixelNumber, :] = [1.00, 0.10, 0.00]
        pixelNumber += 1 
        
        self.predatorLocation[pixelNumber, :] = [1.50, 0.10, 0.00]
        pixelNumber += 1 
        
        self.predatorLocation[pixelNumber, :] = [2.00, 0.10, 0.00]
        pixelNumber += 1 
        
        self.predatorLocation[pixelNumber, :] = [0.50, -0.10, 0.00]
        pixelNumber += 1 
        
        self.predatorLocation[pixelNumber, :] = [1.00, -0.10, 0.00]
        pixelNumber += 1 
        
        self.predatorLocation[pixelNumber, :] = [1.50, -0.10, 0.00]
        pixelNumber += 1 
        
        self.predatorLocation[pixelNumber, :] = [2.00, -0.10, 0.00] 
        pixelNumber += 1 
        
        self.predatorLocation[pixelNumber, :] = [-2.50, 0.00, 0.50]
        pixelNumber += 1 
        
        self.predatorLocation[pixelNumber, :] = [-2.00, 0.00, 0.50]
        pixelNumber += 1 
        
        self.predatorLocation[pixelNumber, :] = [-1.50, 0.00, 0.50]
        pixelNumber += 1 
        
        self.predatorLocation[pixelNumber, :] = [-1.00, 0.00, 0.50]
        pixelNumber += 1 
        
        self.predatorLocation[pixelNumber, :] = [-0.50, 0.00, 0.50]
        pixelNumber += 1 
        
        self.predatorLocation[pixelNumber, :] = [0.00, 0.00, 0.50]
        pixelNumber += 1 
        
        self.predatorLocation[pixelNumber, :] = [0.50, 0.00, 0.50]
        pixelNumber += 1 
        
        self.predatorLocation[pixelNumber, :] = [1.00, 0.00, 0.50]
        pixelNumber += 1 
        
        self.predatorLocation[pixelNumber, :] = [1.50, 0.00, 0.50]
        pixelNumber += 1 
        
        self.predatorLocation[pixelNumber, :] = [2.00, 0.00, 0.50]
        pixelNumber += 1 
        
        self.predatorLocation[pixelNumber, :] = [2.50, 0.00, 0.50] 
        pixelNumber += 1 
        
        self.predatorLocation[pixelNumber, :] = [0.50, 0.10, 0.30]
        pixelNumber += 1 
        
        self.predatorLocation[pixelNumber, :] = [1.00, 0.10, 0.30]
        pixelNumber += 1 
        
        self.predatorLocation[pixelNumber, :] = [1.50, 0.10, 0.30]
        pixelNumber += 1 
        
        self.predatorLocation[pixelNumber, :] = [2.00, 0.10, 0.30]
        pixelNumber += 1 
        
        self.predatorLocation[pixelNumber, :] = [0.50, -0.10, 0.30]
        pixelNumber += 1 
        
        self.predatorLocation[pixelNumber, :] = [1.00, -0.10, 0.30]
        pixelNumber += 1 
        
        self.predatorLocation[pixelNumber, :] = [1.50, -0.10, 0.30]
        pixelNumber += 1 
        
        self.predatorLocation[pixelNumber, :] = [2.00, -0.10, 0.30] 
        pixelNumber += 1 
        
        self.predatorLocation[pixelNumber, :] = [-2.50, -0.00, 0.70]
        pixelNumber += 1 
        
        self.predatorLocation[pixelNumber, :] = [-2.50, -0.00, 0.90] 
        pixelNumber += 1 
        
        self.predatorLocation[pixelNumber, :] = [-2.70, -0.00, 0.70]
        pixelNumber += 1 
        
        self.predatorLocation[pixelNumber, :] = [-2.70, -0.00, 0.90] 
        pixelNumber += 1 
        
        
        self.predatorLocation[pixelNumber, :] = [-3.00, -0.00, -1.10] 
        pixelNumber += 1 
        
        
        self.predatorLocation[pixelNumber, :] = [-3.10, -0.00, -1.10] 
        pixelNumber += 1 
        
        self.predatorLocation[pixelNumber, :] = [-2.50, -0.00, 0.50] 
        pixelNumber += 1 
        
        self.predatorLocation[pixelNumber, :] = [-2.70, -0.00, 0.50]   
        pixelNumber += 1 
        
        self.predatorLocation[pixelNumber, :] = [-2.50, -0.00, 0.10] 
        pixelNumber += 1 
        
        self.predatorLocation[pixelNumber, :] = [-2.70, -0.00, 0.10]
        pixelNumber += 1 
        
        self.predatorLocation[pixelNumber, :] = [-2.50, -0.00, -0.10] 
        pixelNumber += 1 
        
        self.predatorLocation[pixelNumber, :] = [-2.70, -0.00, -0.10] 
        pixelNumber += 1 
        
        self.predatorLocation[pixelNumber, :] = [-2.50, -0.00, -0.50] 
        pixelNumber += 1 
        
        self.predatorLocation[pixelNumber, :] = [-2.70, -0.00, -0.50] 
        pixelNumber += 1 
        
        self.predatorLocation[pixelNumber, :] = [-2.70, -0.00, -0.90]
        pixelNumber += 1 
        
        self.predatorLocation[pixelNumber, :] = [-2.90, -0.00, -0.90] 
        pixelNumber += 1 
        
        self.predatorLocation[pixelNumber, :] = [0.70, -0.00, 1.00]
        pixelNumber += 1 
        
        self.predatorLocation[pixelNumber, :] = [0.10, -0.00, 1.50]
        pixelNumber += 1 
        
        self.predatorLocation[pixelNumber, :] = [0.30, -0.00, 1.30]
        pixelNumber += 1 
        
        
        self.predatorLocation[pixelNumber, :] = [0.55, -0.00, 0.80]
        pixelNumber += 1 
        
        self.predatorLocation[pixelNumber, :] = [0.50, -0.00, 1.10] 
        pixelNumber += 1 
        
        self.predatorLocation[pixelNumber, :] = [-2.90, 0.00, 1.10] 
        pixelNumber += 1 
        
        
        self.predatorLocation[pixelNumber, :] = [2.70, 0.00, 0.00]
        pixelNumber += 1 
        print "The shark consists of %d dots." % pixelNumber
        self.predatorLocation *= 10.00;  