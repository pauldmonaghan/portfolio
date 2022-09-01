#%% Preamble
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
pi = np.pi
cos = np.cos
sin = np.sin
exp = np.exp
re = np.real
im = np.imag
sqrt = np.sqrt
j = 1j
#%% Constants
freq0 = 340                                     # Frequency [Hz]
L = 1.0                                         # Length of pipe [m]
Up = 1.0                                        # Piston velocity [ms^-1]
z = 1.0 + 0.0*j                                 # Non-dimensional impedance at x=L (z=1 is a rho_c impedance)
N0 = 20                                         # Number of equally spaced elements (n is the number of nodes)                                      
c0 = 340                                        # Speed of sound [ms^-1]
rho0 = 1.2                                      # Fluid density [kgm^-3]
#%% Calculations and Functions
def calc_p(N,freq):
    n = int(N) + 1    
    omega = 2*pi*freq                           # Angular Frequency [rads^-1]
    k = omega/c0                                # Wavenumber [radm^-1]
    y = 1/z                                     # Non-dimensional admittance
    delta = L/N                                 # Element length [m]
    x = np.linspace(0,L,n)                      # Position vector [m]
    xex = np.linspace(0,L,1000)                 # Position vector for exact solution [m]
    repex = rho0*c0*Up*cos(k*xex)               # Real part of exact pressure at plotting points [Pa]
    pexnode = rho0*c0*Up*exp(-j*k*x)            # Complex exact pressures at nodes [Pa]
    pexnode = np.transpose(pexnode)             # Converting into a collumn vector
    M = np.zeros((n,n))                         # Mass matrix [kg]
    K = np.zeros_like(M)                        # Stiffness matrix [Nm^-1]
    C = np.zeros((n,n), dtype=np.csingle)       # Damping matrix [Nsm^-1]
    f = np.zeros((n,1), dtype=np.csingle)       # Force vector
    facM = delta/((c0**2)*6)                    # Coefficient for mass Matrix
    facK = 1/delta                              # Coefficient for stiffness matrix
    M[0,0], K[0,0] = facM*2, facK
    M[0,1], K[0,1] = facM, -facK
    M[n-1,n-1], K[n-1,n-1] = facM*2, facK
    M[n-1,n-2], K[n-1,n-2] = facM, -facK 
    for i in range(2,n):
        M[i-1,i-1], K[i-1,i-1] = facM*4, 2*facK
        M[i-1,i], K[i-1,i] = facM, -facK
        M[i-1,i-2], K[i-1,i-2] = facM, -facK
    facC = y/c0                                 # Coefficient for damping matrix
    C[n-1,n-1] = facC
    f[0,0] = j*omega*rho0*Up
    A = K + (j*omega*C) - (omega**2)*M          # Coefficient Matrix
    p = np.matmul(np.linalg.inv(A), f)          # Nodal Pressures [Pa]
    rep = re(p)
    return x, rep, xex, repex, pexnode, p , n
def update(val):
    N = sN.val
    freq = sfreq.val
    l1.set_xdata(calc_p(N,freq)[0])
    l1.set_ydata(calc_p(N,freq)[1])
    l2.set_ydata(calc_p(N,freq)[3])
    fig.canvas.draw_idle()
#%% Plotting
x, rep, xex, repex, pexnode, p, n = calc_p(N0,freq0)
fig, ax = plt.subplots()
plt.subplots_adjust(left=0.03, right=0.955, top=0.975, bottom=0.15)
l1, = plt.plot(x,rep,'ko')
l2, = plt.plot(xex, repex, 'b-', linewidth=3)
plt.xlabel('Distance Along Duct [m]', fontsize=20)
plt.ylabel('Real Acoustic Pressure [Pa]', fontsize=20)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
axN = plt.axes([0.03, 0.075, 0.925, 0.03])
axfreq = plt.axes([0.03, 0.025, 0.925, 0.03])
sN = Slider(axN, '$N$', 10, 400, valinit=N0, valfmt='%0.0f')
sfreq = Slider(axfreq, '$f$', 100, 6000, valinit=freq0)
sN.on_changed(update)
sfreq.on_changed(update)
ax.legend(('Finite Element Approximation','Exact Pressure Field'), fontsize=20)
ax.grid(which='both')
#%% Error analysis
nodalerror = 0.0
pexnorm = 0.0
perror = p-pexnode
for i in range(0,n):
    nodalerror += abs(perror[i])**2
    pexnorm += abs(pexnode[i])**2
nodalerror = sqrt(nodalerror)
pexnorm = sqrt(pexnorm)
percent_nodal_error = (1-(nodalerror/pexnorm))*100
print(percent_nodal_error)