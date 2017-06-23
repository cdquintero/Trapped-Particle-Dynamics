'''
===========================
More triangular 3D surfaces
===========================

Two additional examples of plotting surfaces with triangular mesh.

The first demonstrates use of plot_trisurf's triangles argument, and the
second sets a Triangulation object's mask and passes the object directly
to plot_trisurf.
'''

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.tri as mtri
import math  as mth
from copy import copy, deepcopy

#============
# First plot
#============

def B(x,y,z):
    alpha = -1/((x**2 + y**2 + z**2)**(5/2.0))
    Bx = (3*x*z*alpha)
    By = 3*y*z*alpha
    Bz = (2*(z**2) - x**2 - y**2)*alpha
    return Bx,By,Bz
    
    
def CAPAS(L):

    Distancia = []
    Latitud = []
    BE = 3.11*10E-5 
    B= []
    
    # Make a mesh in the space of parameterisation variables u and v
    Phi = np.linspace(0, 2*np.pi, endpoint=True, num=50)
    Theta = np.linspace(0, np.pi/2, endpoint=True, num=50)
    ii=np.arange(len(Theta))

    flag=1
    while flag:
        for li in zip(Theta, ii):
            Latitud.append(li[0]) 
            Alt_rad = li[0]
#            L=4
            B0=BE/(L**3)
            Bm=B0/((mth.sin(pitch_angle))**2)                                                                                                                  
            B_Dipolo = B0*mth.sqrt(3.0*(mth.sin(Alt_rad)**2)+1.0)/(mth.cos(Alt_rad)**6)
            B.append(B_Dipolo)  
            if B_Dipolo >= Bm:
                flag=0
                break

    Lat=deepcopy(Latitud)
    jj=np.arange(len(Latitud)-1,0-1,-1)

    for i in jj:
        Lat[i]=-Lat[i]

    Lat.reverse()
    Lat.pop(-1)
    Latitud=Lat+Latitud

    Theta=np.array(Latitud)

    Phi, Theta = np.meshgrid(Phi, Theta)
    Phi, Theta = Phi.flatten(), Theta.flatten()
    lm=max(Theta)
    lmp=min(Theta)

    r= L*((np.cos(Theta))**2)
    rm=L*((np.cos(lm))**2)
    rmp=L*((np.cos(lmp))**2)

    # This is the Mobius mapping, taking a u, v pair and returning an x, y, z
    # triple
    x = r*np.cos(Theta)* np.cos(Phi)
    y = r*np.cos(Theta)* np.sin(Phi)
    z = r*np.sin(Theta)

    xm = rm*np.cos(lm)* np.cos(Phi)
    ym = rm*np.cos(lm)* np.sin(Phi)
    zm = rm*np.sin(lm)

    xmp = rmp*np.cos(lmp)* np.cos(Phi)
    ymp = rmp*np.cos(lmp)* np.sin(Phi)
    zmp = rmp*np.sin(lmp)

    return x,y,z,xm,ym,zm,xmp,ymp,zmp,Phi, Theta

    
  
xxx = np.linspace(-4,4,9)
yyy = np.linspace(-4,4,9)
zzz = np.linspace(-4,4,9)

xxx,yyy,zzz = np.meshgrid(xxx,yyy,zzz)

# Plot of the fields  
Bx,By,Bz = B(xxx,yyy,zzz)     


pitch_angle = 45.0                                                             #  initial angle between velocity and mag.field (degrees)
pitch_angle = pitch_angle*np.pi/180

L=4
x,y,z,xm,ym,zm,xmp,ymp,zmp, Phi, Theta= CAPAS(L)

# Triangulate parameter space to determine the triangles
tri = mtri.Triangulation(Phi, Theta)

# Plot the surface.  The triangles in parameter space determine which x, y, z
# points are connected by an edge.

fig = plt.figure()
ax = fig.gca(projection='3d')      

ri=-5
rf=5
ax.set_xlim3d(ri,rf)
ax.set_ylim3d(ri,rf)
ax.set_zlim3d(ri,rf) 

w = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)

xx = np.outer(np.cos(w), np.sin(v))
yy = np.outer(np.sin(w), np.sin(v))
zz = np.outer(np.ones(np.size(w)), np.cos(v))

ax.plot_surface(xx, yy, zz, rstride=6, cstride=6, linewidth=0.3, color='cornflowerblue') 
ax.plot_trisurf(x, y, z,linewidth=0.2, triangles=tri.triangles, color='#e6e6fa',alpha=0.4)
ax.scatter(xm,ym,zm, c='#2f4f4f',linewidth=0.1, s=12)
ax.scatter(xmp,ymp,zmp, c='#2f4f4f',linewidth=0.1,s=12)
ax.quiver(xxx,yyy,zzz,Bx,By,Bz,color='#4682b4', length=0.5,arrow_length_ratio=0.4,linestyles='dotted',linewidths=0.5)        #Plot the magnetic field


ax.set_aspect('equal')

ax.set_xlabel('Eje X')
ax.set_ylabel('Eje Y')
ax.set_zlabel('Eje Z')  
plt.show()
