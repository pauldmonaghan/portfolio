import scipy as sp
import matplotlib.pyplot as plt
def electron_motion (px, py, pz, x, y, z, t, dt):
    gamma=1.0
    q=-1.60e-19
    E=2.74e12
    c=3.0e8
    l=8e-7
    m=9.11e-31
    w=(2*(sp.pi)*c)/l
    k=(-1)*(w/c)
    envelope=list()
    time=list()
    ydisplacement=list()
    zdisplacement=list()
    tau=30e-15
    sigma=c*tau
    z0=4*sigma
    for i in range(0,t):
        r=sp.sqrt((z**2)+(y**2))
        env=sp.exp(-(((r/sigma)**2)+((((z-z0)+(c*i*(dt)))/sigma)**2)))
        gamma=(sp.sqrt(1+(((sp.sqrt(py**2+pz**2))/(m*c))**2)))**(-1)
        dpy=env*(q*E*(sp.sin((k*z)-(w*i*(dt))))*dt*((pz/(m*c*gamma))+1))
        dpz=env*(q*E*(sp.sin((k*z)-(w*i*(dt))))*dt*(py/(m*c*gamma))*(-1))
        py+=dpy
        pz+=dpz
        dy=((py)/(m*gamma))*dt
        dz=((pz)/(m*gamma))*dt
        y+=dy
        z+=dz
        envelope.append(env)
        time.append((i*(dt)))
        ydisplacement.append(y)
        zdisplacement.append(z)
        i+=1
    params = {
    'axes.labelsize': 18,# Set to size 8 to start with
    'font.size': 18,# Set to size 8 to start with
    'legend.fontsize': 18,# Set to size 8 to start with
    'xtick.labelsize': 16,# Set to size 8 to start with
    'ytick.labelsize': 16, # Set to size 8 to start with
    'figure.figsize': [8.8, 8.8/1.618] # Set to [6,4] to start with
    }
    plt.rcParams.update(params)
    plt.plot(zdisplacement, ydisplacement, '-', color="black")
    plt.title("Motion of an electron under the Lorentz Force Law")
    plt.xlabel("Z-Displacement, $z$ (m)")
    plt.ylabel("Y-Displacement, $y$ (m)")
    plt.grid()
    return plt.savefig("Momentum_Draft_10.png")
electron_motion(0, 0, 0, 0, 0, 0, 99999999, 1e-19)
