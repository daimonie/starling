from flock import Flock

starlings = Flock()
starlings.number = 50
starlings.tau = 0.1
starlings.eta = 0.5
starlings.speed = 5.0
starlings.boxsize = 5.0
starlings.initBoids()
starlings.sensitivities *= 5.0

starlings.show()