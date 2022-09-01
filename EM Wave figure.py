import scipy as sp
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import FancyArrowPatch
import matplotlib.pyplot as plt
fig=plt.figure()
ax=fig.gca(projection='3d')
y=sp.linspace(0, 8*sp.pi, 2000)
x=sp.sin(y)
z=sp.sin(y)
r=sp.linspace(0, 0, 2000)
ax.plot(x, y, r, color='blue', linewidth=5, label="Magnetic component")
ax.plot(r, y, z, color='red', linewidth=5, label="Electric component")
ax.plot(r, y, r, color="black", linewidth=5, label="Direction of propagation")
ax.set_title("Electromagentic wave", fontsize=50)
ax.set_xlabel('$x$', fontsize=20)
ax.set_ylabel('$z$', fontsize=20)
ax.set_zlabel('$y$', fontsize=20)
ax.legend(fontsize=20)