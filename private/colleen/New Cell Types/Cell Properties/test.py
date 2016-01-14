import matplotlib.pyplot as plt
import numpy as np
#import pandas #comment this out or move one line down - and it works!
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()    
ax = fig.gca(projection='3d')
theta = np.linspace(-4 * np.pi, 4 * np.pi, 100)
z = np.linspace(-2, 2, 100)
r = z**2 + 1
x = r * np.sin(theta)
y = r * np.cos(theta)
ax.plot(x, y, z, label='parametric curve')
ax.legend()
plt.show()