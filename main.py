from flock import Flock

starlings = Flock()
starlings.mode = 1
starlings.number = 50
starlings.tau = 0.1
starlings.eta = 0.5
starlings.speed = 10.0
starlings.boxSize = 15.0
starlings.habitatSize = 3.0 * starlings.boxSize
starlings.habitatStrength = 5.0
starlings.initBoids()
starlings.sensitivities *= starlings.habitatSize
starlings.display = 1
starlings.length = 0.5
starlings.show()