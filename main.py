from flock import Flock

starlings = Flock()
flocksim = 0
if flocksim == 0:
    starlings.mode = 1
    starlings.number = 200
    starlings.tau = 0.1
    starlings.eta = 0.15
    starlings.speed = 20.0
    starlings.boxSize = 15.0
    starlings.habitatSize = 190.00
    starlings.habitatStrength = starlings.number / 0.05
    starlings.initBoids()

    starlings.positions[ int(starlings.number/2):,:] += 0.8 * starlings.habitatSize / (2.0**0.5)

    starlings.sensitivities *= 0.5 * starlings.boxSize
    starlings.display = 1
    starlings.length = 0.5
    starlings.show()