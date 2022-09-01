import numpy as np
col = 5           # Number of collumns in mesh (or number of nodes in x-direction)
row = 4           # Number of rows in mesh (of number of nodes in y-direction)
n = col*row       # Total number of nodes when arranged in rectangular mesh
C = np.zeros([n,n])
a = 1
b = 1
Ae = 2*a*2*b
A = 1e9
rho0 = 1.2
c0 = 343
J = 1j
freq = 340
omega = 2*np.pi*freq
k = omega/c0
U0 = 1
# Suppose we want the nodes arranged in a i x j arrangement
facC = (A*b)/(3*c0)
for i in range(1, row):
    C[(i*col)-1, (i*col)-1] = 4*facC
    C[(i*col)-1, ((i+1)*col)-1] = facC
    C[((i+1)*col)-1, (i*col)-1] = C[(i*col)-1, ((i+1)*col)-1]
C[col-1, col-1] = 2*facC
C[n-1, n-1] = C[col-1, col-1]
print(C/facC)
fU = np.zeros([n,1], dtype=np.csingle)                  # Force vector (Structural load)
facfU = -J*k*rho0*c0*b*U0
for i in range(1,row):
    fU[(col*i)] = 2*facfU
fU[0] = facfU
fU[n-col] = fU[0]
print(fU/facfU)