#%% Preamble
import numpy as np
import matplotlib.pyplot as plt
pi = np.pi
cos = np.cos
sin = np.sin
exp = np.exp
re = np.real
im = np.imag
sqrt = np.sqrt
J = 1j
#%% Constants
L = 1.0                                                     # Length of domain [m]
W = 0.4                                                     # Width of domain [m]
Z = 1.0 + 0.0*J                                             # Non-dimensional impedance at x=L (z=1 is a rho_c impedance)                                     
c0 = 340                                                    # Speed of sound [ms^-1]
rho0 = 1.2                                                  # Fluid density [kgm^-3]
U0 = 1.0                                                    # Vibrating surface velocity [ms^-1]
#%% Calculations and Functions
def e_to_g(Ge, n):                                          # Ge is element matrix to be used, n is side 
    G = np.zeros([n,n], dtype=np.csingle)                   # dimension of global matrix to be formed
    m = len(Ge[:,0])                                        # Side dimemsion of element matrix
    for i in range(4,n+1):
        for j in range(1,m+1):
            for k in range(j,m+1):
                G[i-k+j-1,i-k] += Ge[j-1,0]
                G[i-k,i-k+j-1] = G[i-k+j-1,i-k]
    return G
def calc_p(col, row, freq):
    n = col*row
    omega = 2*pi*freq                                       # Angular Frequency [rads^-1]
    k = omega/c0                                            # Wavenumber [radm^-1]
    ky = 0
    kx = np.sqrt(k**2 - ky**2)
    a = L/(2*(col-1))                                       # 1/2 Nodal spacing/Element Length in x-direction [m]
    b = W/(2*(row-1))                                       # 1/2 Nodal spacing/Element Length in y-direction [m]
    Ae = 2*a*2*b                                            # Element Area [m^2]
    Y = 1/Z                                                 # Non-dimensional admittance
    x = np.linspace(0,L,n)
    y = np.linspace(0,W,n)
    Me = np.zeros([4,4], dtype=np.csingle)
    Ke = np.zeros([4,4], dtype=np.csingle)
    facM = (a*b)/(9*(c0**2))
    facK = 1/(6*a*b)
    facC = (Ae*b)/(3*c0)
    for i in range(2,4):                                    # Element mass matrix construction
        for j in range(0,i+1):
            Me[j,j] = 4*facM
            Me[j-1,j] = 2*facM
            Me[j-2,j] = facM
            Me[j-3,j] = 2*facM
    for i in range(0,4):                                    # Element stiffness matrix construction
        for j in [3,2,1,0]:
            Ke[i,j] = (2*(b**2) - (a**2))*facK
            Ke[i,i] = (2*((a**2) + (b**2)))*facK
    for i in range(0,2):
        Ke[i+2,i] = (-((a**2) + (b**2)))*facK
        Ke[i,i+2] = Ke[i+2,i]
    Ke[0,3] = (2*(b**2) - (a**2))*facK
    Ke[3,0] = Ke[0,3]
    M = e_to_g(Me, n)                                       # Global mass matrix construction
    K = e_to_g(Ke, n)                                       # Global stiffness matrix construction
    C = np.zeros([n,n], dtype=np.csingle)
    # Global damping matrix construction
    # Acoustically hard walls, so expect zero damping
    # Only expect finite damping with impedance
    # boundary condition
    fQ = np.zeros([n,1], dtype=np.csingle)                  # Force vector (Distributed source). No sources in domain, therefore zero    
    fU = np.zeros([n,1], dtype=np.csingle)                  # Force vector (Structural load)
    facfU = -J*k*rho0*c0*b*U0                               # Global Force vector construction
    fU[0] = facfU
    fU[n-row] = fU[0]
    for i in range(1,col-1):
        fU[(row*i)] = 2*facfU
    A = K + J*omega*C - (omega**2)*M                        # Coefficient Matrix (currently without damping)
    p = np.zeros([n,1], dtype=np.csingle)
    p = np.matmul(np.linalg.inv(A), (fU+fQ))                            # Nodal Pressures [Pa]
    xex = np.linspace(0,L,1000)
    pex = -(J*cos(kx*xex))/(sin(kx*L))*U0*rho0*c0
    pexnode = -(J*cos(kx*x))/(sin(kx*L))*U0*rho0*c0
    return x, y, kx, p, xex, pex, pexnode, n

x, y, kx, p, xex, pex, pexnode, n = calc_p(4, 3, 400)
print(p)
#%% Plotting
plt.plot(xex, re(pex), 'b-', linewidth=5)
plt.plot(x, re(p), 'k.', markersize=20)
#plt.plot(x, re(pexnode), 'r.', linewidth=3)
plt.xlabel('Distance along length of tube, $x/L$')
plt.ylabel('Pressure, $|P|$ [Pa]')
plt.grid(which='both')
plt.legend(('Exact analytic solution','Finite element approximation'))