from flock import Flock
import numpy as np

starlings = Flock()
flocksim = 1
if flocksim == 0:
    starlings.mode = 1
    starlings.number = 600
    starlings.tau = 0.1
    starlings.eta = 0.15
    starlings.speed = 30.0
    starlings.boxSize = 50.0
    starlings.habitatSize = 200.00
    starlings.habitatStrength = starlings.number / 0.3
    starlings.initBoids()

    starlings.positions[ int(starlings.number/2):,:] += 0.8 * starlings.habitatSize / (2.0**0.5)

    starlings.sensitivities *= 0.5 * starlings.boxSize
    starlings.display = 1
    starlings.length = 10.0
    starlings.show()
elif flocksim == 1:
    starlings.mode = 2
    starlings.number = 100
    starlings.tau = 0.1
    starlings.eta = 0.5 
    starlings.boxSize = 100.0
    starlings.habitatSize = 200.00
    starlings.habitatStrength = 1.05 * starlings.number
    starlings.speed = 25.00
     
    starlings.i0 = 2.0
    starlings.i1 = 2.0
    starlings.i2 = 5.5
    starlings.i3 = -0.50
    
    
    starlings.initBoids()
    starlings.velocities = np.ones((starlings.number,3), dtype=np.float)
    starlings.sensitivities *= 5.0
    starlings.display = 1
    starlings.length = 5.00
    starlings.show()    