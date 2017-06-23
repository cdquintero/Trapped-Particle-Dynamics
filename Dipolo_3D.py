from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np


x = np.linspace(-4,4,9)
y = np.linspace(-4,4,9)
z = np.linspace(-4,4,9)

x,y,z = np.meshgrid(x,y,z)

# 3d figure
fig = plt.figure()
ax = fig.gca(projection='3d')

def B(x,y,z):
    alpha = -1/((x**2 + y**2 + z**2)**(5/2.0))
    Bx = (3*x*z*alpha)
    By = 3*y*z*alpha
    Bz = (2*(z**2) - x**2 - y**2)*alpha
    return Bx,By,Bz

def cylinder(u,v):
    x = np.outer(np.cos(u), np.sin(v))
    y = np.outer(np.sin(u), np.sin(v))
    z = np.outer(np.ones(np.size(u)), np.cos(v))
    return x,y,z

u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)

# Plot of the fields
Bx,By,Bz = B(x,y,z)                                   #Magnetic field
cx,cy,cz = cylinder(u,v)                               #Wire

ax.set_xlim3d(-5, 5)
ax.set_ylim3d(-5, 5)
ax.set_zlim3d(-5, 5)  

# Plot of the 3d vector field
ax.quiver(x,y,z,Bx,By,Bz,color='#4682b4', length=0.3,arrow_length_ratio=0.4,linestyles='solid',linewidths=0.5,zorder=1)        #Plot the magnetic field
ax.plot_surface(cx,cy, cz,label='Tierra', rstride=6, cstride=6, linewidth=0.5, color='cornflowerblue') 
ax.set_aspect('equal')

ax.set_xlabel('Eje X')
ax.set_ylabel('Eje Y')
ax.set_zlabel('Eje Z')  
plt.show()