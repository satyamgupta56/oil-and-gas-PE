import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

"""Darcy law"""

re=3500   #ft
rwf=0.6   #ft
pwf=1400  #psi
h=50      #ft
k=145     #md
q=200     #STB/day
mu=15     #cp
b=1

r=np.linspace(rwf,re,700)   #700 discreate points
# print(r)

pressure=[]
for i in range(len(r)):
  p=pwf + (141.2*q*mu*b*(np.log(r[i]/rwf))/k/h)
  pressure.append(p)

# print(pressure)

plt.style.use("classic")

plt.figure(figsize=[8,6])
plt.plot(r,pressure)
plt.xlabel('r (ft)',size=20)
plt.ylabel('P (r), Psi', size=20)
plt.title('Reservoir Pressure Profile',size=20)
plt.grid(True)
plt.show()