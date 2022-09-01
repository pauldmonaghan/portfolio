#%% Preamble
import numpy as np
#import scipy as sp
import matplotlib.pyplot as plt
import math
pi = np.pi
cos = np.cos
sin = np.sin
exp = np.exp
re = np.real
im = np.imag
sqrt = np.sqrt
J = 1j
#
# Input data from user
L = 1.0                                                     # Length of domain [m]
W = 0.4                                                     # Width of domain [m]
c0 = 340                                                    # Speed of sound [ms^-1]
rho0 = 1.2                                                  # Fluid density [kgm^-3]
# =============================================================================
freq = 320
U0 = 1.0 + 0*J                                              # Vibrating surface velocity [ms^-1]
#Z = 1.0 - 0.0*J                                            # Non-dimensional impedance at x=L (z=1 is a rho_c impedance)                                     
Z = 1.0e+10 - 0.0*J                                             # Non-dimensional impedance at x=L (z<10^6~inf - Hard wall)                                     
piston_end = 1                                               # 1: Left, 2: Bottom, 3: Right, 4: Top
absorbing_end = 3                                           # 1: Left, 2: Bottom, 3: Right, 4: Top, 0: None
elem_type = 1                                            # 1: linear rectangle (NOT YET)2: linear triangle 
# =============================================================================
nex = 20                                                    # number of elements in x direction
ney = 8                                                    # number of elements in y direction
# =============================================================================
# =============================================================================
omega = 2*pi*freq                                       # Angular Frequency [rads^-1]
k = omega/c0                                            # Wavenumber [radm^-1]
ky = 0
kx = np.sqrt(k**2 - ky**2)
# =============================================================================
#from Helmholtz_FE_functions import e_to_g2, calc_p 
from Helmholtz_FE_functions import *                        # import functions from a separate file
# =============================================================================
Y = calc_nondim_admittance(Z)                                 # get non-dim admittance Y from Z
# =============================================================================
#
if elem_type == 1:
    nn_elem = 4                                                # Number of nodes per element (linear rectangular)
elif elem_type == 2:
    nn_elem = 3                                                # Number of nodes per element (linear triangular)
else:
    print("ERROR!  Element type must be 1 or 2!!")    
    quit()
n_elem_edges = nn_elem                                         # Number of edges of each element          
#
#==============================================================================
# mesh creation, assign BC to element edges, and set up necessary arrays
#==============================================================================
xn,yn,elem_nodes,a,b,Ae,nnx,nny,nn,nel = set_mesh(L,W,nex,ney,nn_elem)
# =============================================================================
# Set array for element edges with excitation boundary (piston)
elem_edges_piston_num, elem_edges_piston_nodes = set_bc_excite(piston_end,nnx, nny, nex, ney)
#print("elem_edges_piston_num")
#print(elem_edges_piston_num)    
#print(elem_edges_piston_nodes)    
# =============================================================================
# Set array for element edges with absorbing boundary
elem_edges_absorb_num, elem_edges_absorb_nodes = set_bc_absorb(absorbing_end,nnx, nny, nex, ney)
#print(elem_edges_absorb_num)    
#print(elem_edges_absorb_nodes)
# =============================================================================
# Initialise arrays for element edge BCs (NOT used at the moment)
#elem_edges_bc_U = np.zeros([nel,n_elem_edges], dtype=np.csingle)  # Complex array of element edge U (normal velocity)
#elem_edges_bc_Y = np.zeros([nel,n_elem_edges], dtype=np.csingle)  # Complex array of element edge U (normal velocity)
#print("elem_nodes")
#print(elem_nodes)
#print("")
#
p_FE = calc_p(nex,ney,k,omega,elem_nodes, a,b,c0,rho0,U0,Y,elem_edges_piston_num,elem_edges_piston_nodes,elem_edges_absorb_num,elem_edges_absorb_nodes)
# =============================================================================
# Call analytical solution for 1D propagation at x-coordinates of nodes 
xex, pex = analytic_soln_1D(xn, L, k, c0, rho0, U0, Y)
# =============================================================================
#
plot_p_vs_x(nnx,nny,xn,p_FE,xex,pex) 
# =============================================================================
