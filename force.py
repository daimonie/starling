import matplotlib.pyplot as plt
import numpy as np


x = np.linspace(0.0,25,1000);

y = x*0;

s = 5.0


i0 = 0.0
i1 = 1.25
i2 = 0.35
i3 = 4.5
i4 = -0.50
i5 = 0.65

for i in range(0,1000):
    thisx = x[i]; 
    y[i] = i0 + i1/( 1.0 + np.exp((thisx-i1)/(2*i2**2))) 
    if thisx > i3:
        y[i] += i4*(thisx-i3)*1.0 / ( 1.0 + np.exp( (thisx-s)/(2*i5**2)))
    #The fermi-dirac term takes care of the 'heaviside' behaviour in a way that is continous

plt.plot(x,y) 
plt.show()