#So, given the latest matplotlib version, quiver3D and quiver should support a full 3D vector field, like we want.
#
#Let's make a minimal example.
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import matplotlib.animation as animation
 

fig = plt.figure()
ax = fig.gca(projection='3d')
print "Minimal example of a quiver animation."

#Pretty much an altered version of the standard demo.
def newFrame(r): 
    global fig, ax
    x, y, z = np.meshgrid(np.arange(-0.8, 1.2, 0.2),
                        np.arange(-0.8, 1.2, 0.2),
                        np.arange(-0.8, 1.2, 0.2)) 
    omega = r * 0.1
    u = np.sin(np.pi * x* (0.01  + omega)) * np.cos(np.pi * y * (0.01  + omega)) * np.cos(np.pi * z* (0.01  + omega))
    v = np.cos(np.pi * x* (0.01  + omega)) * np.sin(np.pi * y* (0.01  + omega)) * np.cos(np.pi * z* (0.01  + omega))
    w = np.cos(np.pi * x* (0.01  + omega)) * np.cos(np.pi * y* (0.01  + omega)) *np.sin(np.pi * z* (0.01  + omega))
    
    ax.cla()
    ax.quiver(x, y, z, u, v, w, length=0.1)
    
iMax = 150
lastpositions = np.zeros((iMax,3));
for i in range(0,iMax):
    r = 0
    
    
    theta = 2.0 * np.pi / iMax * (i + r * 0.5 )
    phi = np.pi / 2.0 * (1 + r * 0.01 * i/ iMax )
    lastpositions[i,0] = 3.0 * np.cos( theta ) * np.sin( phi )
    lastpositions[i,1] = 3.0 * np.sin( theta ) * np.sin( phi )
    lastpositions[i,2] = 3.0 * np.cos( phi )
def chase(r):
    global fig, ax, lastpositions, iMax
    
    newpositions = np.zeros((iMax,3));
    for i in range(0,iMax):
        theta = 2.0 * np.pi / iMax * (i + r * 0.5 )
        phi = np.pi / 2.0 * (1 + r * 0.01 * i/ iMax )
        newpositions[i,0] = 3.0 * np.cos( theta ) * np.sin( phi )
        newpositions[i,1] = 3.0 * np.sin( theta ) * np.sin( phi )
        newpositions[i,2] = 3.0 * np.cos( phi )
        
    diff = newpositions - lastpositions 
    
    length = np.sqrt( np.sum( np.sum( np.square( diff ) ) ) )
    if length > 0.05:
        diff /= length
    
    ax.cla ()
    ax.quiver(newpositions[:,0], newpositions[:,1], newpositions[:,2], diff[:,0], diff[:,1], diff[:,2], length=0.1);
    ax.set_xlim([-4.0, 4.0])
    ax.set_ylim([-4.0, 4.0])
    ax.set_zlim([-4.0, 4.0])
    
    
    lastpositions = newpositions
    
#A small function of a few little boids flying in circles
ani = animation.FuncAnimation(fig, chase,interval=100) 
plt.show()

