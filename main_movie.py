from flock import Flock
import numpy as np
import argparse as argparse;

parser = argparse.ArgumentParser(prog="python main.py",
  description = "Main.py handles everything.  Use python main.py -h to see your options.");

parser.add_argument('-m', '--mode', help='Set the simulation mode.', default = 3, action='store', type = int); 

parser.add_argument('-s', '--snapshots', help='Set the number of frames for the video', default = 1200, action='store', type = int); 

parser.add_argument('-f', '--filename', help='Set the simulation mode.', default = "default", action='store'); 

starlings = Flock()

args = parser.parse_args();
 
flocksim = args.mode
filename = args.filename + ".mp4"
snapshots = args.snapshots #number of frames
if flocksim == 0:
    starlings.mode = 1
    starlings.number = 40
    starlings.tau = 0.1
    starlings.eta = 0.25
    starlings.speed = 25.0
    starlings.boxSize = 100.0
    starlings.habitatSize = 200.00
    starlings.habitatStrength = 20.0
    starlings.initBoids()
    
    starlings.velocities = np.ones((starlings.number,3), dtype=np.float)
    starlings.positions[ int(starlings.number/2):,:] += 0.8 * starlings.habitatSize / (2.0**0.5)

    starlings.sensitivities *= 30.00
    starlings.display = 1
    starlings.length = 10.0
    starlings.record(filename, snapshots)
elif flocksim == 1:
    starlings.mode = 2
     
    starlings.number = 40
    starlings.tau = 0.1
    starlings.eta = 0.25
    starlings.speed = 25.0
    starlings.boxSize = 100.0
    starlings.habitatSize = 200.00
    starlings.habitatStrength = 20.0
    starlings.initBoids()
        
    starlings.i0 = 0.0
    starlings.i1 = 0.9
    starlings.i2 = 0.35
    starlings.i3 = 3.5
    starlings.i4 = -0.20
    starlings.i5 = 0.65 
    
    starlings.velocities = np.ones((starlings.number,3), dtype=np.float)
    starlings.positions[ int(starlings.number/2):,:] += 0.8 * starlings.habitatSize / (2.0**0.5)

    starlings.sensitivities *= 30.00
    starlings.display = 1
    starlings.length = 10.0
    starlings.record(filename, snapshots)
elif flocksim == 2:
    print "SHARK TIME";
    starlings.mode = 3
    
    
    starlings.scatterPlot = True
    starlings.number = 800
    starlings.tau = 0.1
    starlings.eta = 0.45
    starlings.habitatSize = 100.00
    starlings.boxSize = starlings.habitatSize
    starlings.speed = 10.0
    starlings.habitatStrength = 3.0
    starlings.initBoids()
    
    starlings.eta = 0.25
    starlings.boxSize = 125.0
    starlings.habitatSize = 200.00
    starlings.habitatStrength = 20.0
    starlings.initBoids()
        
    starlings.i0 = 0.0
    starlings.i1 = 0.9
    starlings.i2 = 0.35
    starlings.i3 = 3.5
    starlings.i4 = -0.20
    starlings.i5 = 1.65 
     
    #starlings.velocities[:,0] = starlings.velocities[:,1] = 0.0
    #starlings.positions -= 0.2* starlings.habitatSize

    starlings.sensitivities *= 25.00
    starlings.display = 1
    starlings.length = 15.0
    
    starlings.initShark()
    starlings.predatorLocation *= 3.00;  
    starlings.predatorSense *= starlings.speed
    starlings.predatorStrength *= 5.00
    starlings.sharkSpeed = 25.0
    starlings.sharkOmega = 25.0
    
    starlings.record(filename, snapshots)
    
    
elif flocksim == 4: 
    starlings.mode = 2
    
    
    starlings.scatterPlot = True
    starlings.number = 800
    starlings.tau = 0.1
    starlings.eta = 0.45
    starlings.habitatSize = 100.00
    starlings.boxSize = starlings.habitatSize
    starlings.speed = 10.0
    starlings.habitatStrength = 3.0
    starlings.initBoids()
    
    starlings.eta = 0.25
    starlings.boxSize = 125.0
    starlings.habitatSize = 200.00
    starlings.habitatStrength = 20.0
    starlings.initBoids()
        
    starlings.i0 = 0.0
    starlings.i1 = 0.9
    starlings.i2 = 0.35
    starlings.i3 = 3.5
    starlings.i4 = -0.20
    starlings.i5 = 1.65 
     
    #starlings.velocities[:,0] = starlings.velocities[:,1] = 0.0
    #starlings.positions -= 0.2* starlings.habitatSize

    starlings.sensitivities *= 25.00
    starlings.display = 1
    starlings.length = 15.0 
    
    starlings.record(filename, snapshots)
elif flocksim == 3:
    print "SHARK TIME 2";
    starlings.mode = 4
    starlings.scatterPlot = True
    starlings.number = 800
    starlings.tau = 0.1
    starlings.eta = 0.45
    starlings.habitatSize = 100.00
    starlings.boxSize = starlings.habitatSize
    starlings.speed = 10.0
    starlings.habitatStrength = 3.0
    starlings.initBoids()
        
    starlings.i0 = 0.0
    starlings.i1 = 0.8
    starlings.i2 = 0.25
    starlings.i3 = 3.5
    starlings.i4 = -0.20
    starlings.i5 = 1.65 
     
    #starlings.velocities[:,0] = starlings.velocities[:,1] = 0.0
    #starlings.positions -= 0.2* starlings.habitatSize

    starlings.sensitivities *= 3.00
    starlings.display = 1
    starlings.length = 15.0
    
    starlings.initShark()
    starlings.predatorLocation *= 3.00;  
    starlings.predatorSense *= starlings.speed
    starlings.predatorStrength *= 5.00
    starlings.sharkSpeed = 25.0
    starlings.sharkOmega = 25.0
    
    starlings.record(filename, snapshots)
elif flocksim == 5: 
    starlings.mode = 5
    starlings.scatterPlot = True
    starlings.number = 800
    starlings.tau = 0.1
    starlings.eta = 0.45
    starlings.habitatSize = 100.00
    starlings.boxSize = starlings.habitatSize
    starlings.speed = 10.0
    starlings.habitatStrength = 3.0
    starlings.initBoids()
        
    starlings.i0 = 0.0
    starlings.i1 = 0.8
    starlings.i2 = 0.25
    starlings.i3 = 3.5
    starlings.i4 = -0.20
    starlings.i5 = 1.65 
     
    #starlings.velocities[:,0] = starlings.velocities[:,1] = 0.0
    #starlings.positions -= 0.2* starlings.habitatSize

    starlings.sensitivities *= 3.00
    starlings.display = 1
    starlings.length = 15.0
    
    starlings.initShark()
    starlings.predatorLocation *= 3.00;  
    starlings.predatorSense *= starlings.speed
    starlings.predatorStrength *= 5.00
    starlings.sharkSpeed = 25.0
    starlings.sharkOmega = 25.0
    
    starlings.record(filename, snapshots)
else: 
    starlings.mode = 0
    starlings.number = 60
    starlings.tau = 0.1
    starlings.eta = 0.15
    starlings.speed = 30.0
    starlings.boxSize = 10.0
    starlings.habitatSize = 200.00
    starlings.habitatStrength = starlings.number / 0.3
    starlings.initBoids()
    
    starlings.velocities = np.ones((starlings.number,3), dtype=np.float)
    starlings.positions[ int(starlings.number/2):,:] += 0.8 * starlings.habitatSize / (2.0**0.5)

    starlings.sensitivities *= 30.00
    starlings.display = 1
    starlings.length = 10.0
    starlings.record(filename, snapshots)