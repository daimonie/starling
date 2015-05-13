from flock import Flock
import numpy as np

starlings = Flock()
flocksim = 1
if flocksim == 0:
    starlings.mode = 1
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
    starlings.show()
elif flocksim == 1:
    starlings.mode = 2
     
    starlings.number = 60
    starlings.tau = 0.1
    starlings.eta = 0.15
    starlings.speed = 30.0
    starlings.boxSize = 10.0
    starlings.habitatSize = 200.00
    starlings.habitatStrength = starlings.number / 0.3
    starlings.initBoids()
        
    starlings.i0 = 2.0
    starlings.i1 = 0.4
    starlings.i2 = 0.35
    starlings.i3 = 3.5
    starlings.i4 = -0.20
    starlings.i5 = 0.65 
    
    starlings.velocities = np.ones((starlings.number,3), dtype=np.float)
    starlings.positions[ int(starlings.number/2):,:] += 0.8 * starlings.habitatSize / (2.0**0.5)

    starlings.sensitivities *= 30.00
    starlings.display = 1
    starlings.length = 10.0
    starlings.show()
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
    starlings.show()