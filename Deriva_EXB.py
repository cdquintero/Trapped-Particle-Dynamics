# -*- coding: utf-8 -*-
"""
Created on Sat Jan 07 17:52:26 2017

@author: cd67
"""

import math  as mth
import numpy as np                                                             #Paquete Arrays
import matplotlib.pyplot as plt 
import sympy as sym      
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
from copy import copy, deepcopy
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d

class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0, 0), (0, 0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        FancyArrowPatch.draw(self, renderer)

class V_Static(object) :                                                       #Se define una clase que contiene todas las variables estáticas del programa
    e = 1.602176565e-19                                                        #Carga de las partículas
    m_pr = 1.672621777e-27                                                     #Masa del protón [kg]
    m_el = 9.10938291e-31                                                      #Masa del electrón [kg]
    
    B0=3.07e-5                                                                 #Campo magnético ecuatorial
    qom=0                                                                  
    n=6 
    ng=4     
    qomp = 1
    qome=-3    
    Bz = 1.0 
    Bx = 0.0 
    By = 0.0 
    Ex = 0.2 
    Ey = 0.0 
    Ez = 0.0 
    
    B=mth.sqrt(Bx**2+By**2+Bz**2)
    Omega=B
    vperp0=0
    vpar0=0

#______________________________________________________________________________

def RK4(u0,time,ht,St):

    t=time[0]
    for jj in j:
        sr[0,jj]=u0[jj]

    for ii in i:
        for jj in j:
            K1[jj]=ht*DERIVA_EXB(t,sr[ii,0],sr[ii,1],sr[ii,2],sr[ii,3],sr[ii,4],sr[ii,5],jj,St)
        for jj in j:   
            K2[jj]=ht*DERIVA_EXB((t+ht/2),(sr[ii,0]+K1[0]/2),(sr[ii,1]+K1[1]/2),(sr[ii,2]+K1[2]/2),(sr[ii,3]+K1[3]/2),(sr[ii,4]+K1[4]/2),(sr[ii,5]+K1[5]/2),jj,St)
        for jj in j:   
            K3[jj]=ht*DERIVA_EXB((t+ht/2),(sr[ii,0]+K2[0]/2),(sr[ii,1]+K2[1]/2),(sr[ii,2]+K2[2]/2),(sr[ii,3]+K2[3]/2),(sr[ii,4]+K2[4]/2),(sr[ii,5]+K2[5]/2),jj,St)
        for jj in j:
            K4[jj]=ht*DERIVA_EXB(t+ht,(sr[ii,0]+K3[0]),(sr[ii,1]+K3[1]),(sr[ii,2]+K3[2]),(sr[ii,3]+K3[3]),(sr[ii,4]+K3[4]),(sr[ii,5]+K3[5]),jj,St)
        for jj in j:    
            sr[ii+1,jj]=sr[ii,jj]+(K1[jj]+2*K2[jj]+2*K3[jj]+K4[jj])/6.0
    
        t=time[ii+1]
    return(sr)
#______________________________________________________________________________
def DERIVA_EXB(t,u0,u1,u2,u3,u4,u5,j,St):
  
    if j==0 :
        y1 = u3                                                                # dx/dt = u
    if j==1 :                                                                 
        y1 = u4                                                                # dy/dt = v
    if j==2 :
        y1 = u5                                                                # dz/dt = w
    if j==3 :    
        y1 = St.qom*(St.Ex + u4*St.Bz - u5*St.By)                                    # du/dt = q/m*(vel x B)_x
    if j==4 :     
        y1 = St.qom*(St.Ey + u5*St.Bx - u3*St.Bz)                                         # dv/dt = q/m*(vel x B)_y
    if j==5 :    
        y1 = St.qom*(St.Ez+ u3*St.By - u4*St.Bx)                                         # dw/dt = q/m*(vel x B)_z

    return y1

#______________________________________________________________________________

def RK4gc(u0,time,ht,St):

    t=time[0]
    for jj in jg:
        srgc[0,jj]=u0[jj]

    for ii in i:
        for jj in jg:
            K1[jj]=ht*EXBGC(t,srgc[ii,0],srgc[ii,1],srgc[ii,2],srgc[ii,3],jj,St)
        for jj in jg:   
            K2[jj]=ht*EXBGC((t+ht/2),(srgc[ii,0]+K1[0]/2),(srgc[ii,1]+K1[1]/2),(srgc[ii,2]+K1[2]/2),(srgc[ii,3]+K1[3]/2),jj,St)
        for jj in jg:   
            K3[jj]=ht*EXBGC((t+ht/2),(srgc[ii,0]+K2[0]/2),(srgc[ii,1]+K2[1]/2),(srgc[ii,2]+K2[2]/2),(srgc[ii,3]+K2[3]/2),jj,St)
        for jj in jg:
            K4[jj]=ht*EXBGC(t+ht,(srgc[ii,0]+K3[0]),(srgc[ii,1]+K3[1]),(srgc[ii,2]+K3[2]),(srgc[ii,3]+K3[3]),jj,St)
        for jj in jg:    
            srgc[ii+1,jj]=srgc[ii,jj]+(K1[jj]+2*K2[jj]+2*K3[jj]+K4[jj])/6.0
    
        t=time[ii+1]
    return(srgc)
#______________________________________________________________________________
def EXBGC(t,u0,u1,u2,u3,j,St):
    
    # b unit vector
    b_unit_x = St.Bx/St.B
    b_unit_y = St.By/St.B
    b_unit_z = St.Bz/St.B  

    #U
    Ux = St.Ey*b_unit_z - St.Ez*b_unit_y;
    Uy = St.Ez*b_unit_x - St.Ex*b_unit_z;
    Uz = St.Ex*b_unit_y - St.Ey*b_unit_x;
    
    if j==0 :
        y2 = Ux/(St.B) + u3*b_unit_x                                                                # dx/dt = u
    if j==1 :                                                                 
        y2 = Uy/(St.B) +u3*b_unit_y                                                                # dy/dt = v
    if j==2 :
        y2 = Uz/(St.B) +u3*b_unit_z                                                                # dz/dt = w
    if j==3 :    
        y2 = St.e*St.Ez                                    # du/dt = q/m*(vel x B)_x

    return y2
   
#______________________________________________________________________________


def B(x,y,z):
    Bx = 0
    By = 0
    Bz = np.abs(z)
    
    return Bx,By,Bz,    

def E(x,y,z):    
    Ex = -np.abs(x)
    Ey = 0
    Ez = 0
    return Ex,Ey,Ez
    
Bxx = np.linspace(-4,4,3)
Byy = np.linspace(-4,4,3)
Bzz = np.linspace(-4,4,3)/1.7

Exx = np.linspace(-4,4,3)-0.5
Eyy = np.linspace(-4,4,3)
Ezz = np.linspace(-4,4,2)/1.7-0.5

Bxx,Byy,Bzz = np.meshgrid(Bxx,Byy,Bzz)
Exx,Eyy,Ezz = np.meshgrid(Exx,Eyy,Ezz)


# Plot of the fields
Bx,By,Bz = B(Bxx,Byy,Bzz)        
Ex,Ey,Ez = E(Exx,Eyy,Ezz)

St = V_Static()   

#Parámetros de entrada del método RK4
ti=0
tf=28 
ht=0.01                                                
N=int((tf-ti)/ht) 
sr=np.zeros([N,St.n])                                                
srgc=np.zeros([N,St.ng])                                                
vpr=np.zeros(N)                                                
vprgc=np.zeros(N)                                                         
srpgc=np.zeros([N,St.ng]) 
sregc=np.zeros([N,St.ng]) 
K1=np.zeros(St.n)
K2=np.zeros(St.n)
K3=np.zeros(St.n)
K4=np.zeros(St.n)
u0p=np.zeros(St.n)
u0e=np.zeros(St.n)
u0pgc=np.zeros(St.ng)
u0egc=np.zeros(St.ng)
j=np.arange(St.n)
jg=np.arange(St.ng)
i=np.arange(N-1)
k=np.arange(N)
r=3
time=np.arange(ti,tf,ht)

# Condiciones iniciales para la posición y la velocidad
u0p[0] = -2.0                                                                    #x
u0p[1] = 2.0                                                                    #y
u0p[2] = -2.0                                                                    #z
u0p[3] = 0.5                                                                    #vx
u0p[4] = 0.5                                                                    #vy
u0p[5] = 0.0    

xlp=u0p[3]/St.B
ylp=u0p[4]/St.B 
St.vperp0=mth.sqrt(u0p[3]**2+u0p[4]**2)
St.vpar0=u0p[5]

u0pgc[0] = -2.0 +xlp                                                                    #x
u0pgc[1] = 2.0 -ylp                                                                     #y
u0pgc[2] = -2.0                                                                    #z
u0pgc[3] = St.vpar0                                                                    #vx

St.qom=St.qomp
srp=deepcopy(RK4(u0p,time,ht,St))
srpgc=deepcopy(RK4gc(u0pgc,time,ht,St))

# Condiciones iniciales para la posición y la velocidad
rl=St.vperp0/St.B

u0e[0] = 2.0                                                                    #x
u0e[1] = 2.0                                                                    #y
u0e[2] = -2.0                                                                    #z
u0e[3] = 0.5                                                                    #vx
u0e[4] = 0.5                                                                    #vy
u0e[5] = 0.0                                                                     #vz

xle=u0e[3]/(np.abs(St.qome)*St.B)
yle=u0e[4]/(np.abs(St.qome)*St.B)  
St.vperp0=mth.sqrt(u0e[3]**2+u0e[4]**2)
St.vpar0=u0e[5]

u0egc[0] = 2.0-xle                                                                   #x
u0egc[1] = 2.0+yle                                                                   #y
u0egc[2] = -2.0                                                                    #z
u0egc[3] = St.vpar0                                                                     #vx

St.qom=St.qome
sre=deepcopy(RK4(u0e,time,ht,St))
sregc=deepcopy(RK4gc(u0egc,time,ht,St))

if  (u0p[5]==0 and u0e[5]==0):
    for indice in zip(srp[:,3],srp[:,4],k):
        vpr[indice[2]]=mth.sqrt(indice[0]**2+indice[1]**2)
    
    fig, ax = plt.subplots()
    fig2, ax2 = plt.subplots()      
    ax2.plot(time,vpr)
    ax2.plot(time,srpgc[:,3], color='g')
    ax2.set_ylim(-0.1, 1.2)
    ax.plot(srp[:,0], srp[:,1], color='r')
    ax.plot(sre[:,0], sre[:,1], color='b')
    ax.plot(srpgc[:,0],srpgc[:,1], color='orange')
    ax.plot(sregc[:,0],sregc[:,1], color='g') 

    ax.scatter(0, 0)   
    ax.arrow(0.2, 0, 0.5, 0, color='r',head_width=0.05, head_length=0.1) 
    ax.arrow(0, -0.5, 0 , -1.5,head_width=0.05, head_length=0.1)  

    ax.scatter(0, 0)   
    ax.text(x = -0.1, y = 0.1, s = u'B', fontsize = 12)
    ax.text(x = 0.2, y = 0.1, s = u'E', fontsize = 12)
    ax.text(x = 0.1, y = -0.9, s = u'$V_{E}$', fontsize = 15)
    ax.scatter(u0pgc[0], u0pgc[1], color='black')   
    ax.text(u0pgc[0], u0pgc[1]+0.1, s = u'$r_{cg}$', fontsize = 15)
    ax.scatter(u0egc[0], u0egc[1], color='black')    
    ax.text(u0egc[0],u0egc[1]+0.1, s = u'$r_{cg}$', fontsize = 15)

    ax.scatter(u0p[0], u0p[1], s = 60, color='r')   
    ax.text(u0p[0]-0.3, u0p[1]+0.1, s = u'+', fontsize = 15)
    ax.scatter(u0e[0], u0e[1], color='b',s = 30)    
    ax.text(u0e[0]+0.2, u0e[1]-0.1, s = u'-', fontsize = 20)
    
    ax.set_xlabel('Eje X')
    ax.set_ylabel('Eje Y')

else:
    for indice in zip(srp[:,3],srp[:,4],srp[:,5],k):
        vpr[indice[3]]=mth.sqrt(indice[0]**2+indice[1]**2+indice[2]**2)

    fig = plt.figure()
    fig2, ax2 = plt.subplots()      
    ax = fig.gca(projection='3d')      
    ax.set_xlim3d(-r, r)
    ax.set_ylim3d(-r, r)
    ax.set_zlim3d(-r, r)  
    ax.set_aspect('equal')

    ax2.plot(time,vpr)
    ax2.plot(time,srpgc[:,3], color='g')
    ax.plot(srp[:,0], srp[:,1], srp[:,2], color='r')
    ax.plot(sre[:,0], sre[:,1], sre[:,2], color='b')
    ax.quiver(Bxx,Byy,Bzz,Bx,By,Bz,color='#4682b4', length=0.5,arrow_length_ratio=0.3,linestyles='solid',linewidths=0.5)        #Plot the magnetic field
    ax.quiver(Exx,Eyy,Ezz,Ex,Ey,Ez,color='r', length=0.5,arrow_length_ratio=0.3,linestyles='solid',linewidths=0.5)        #Plot the magnetic field
    VD = Arrow3D([0, 0], [2, 1], [-2, -2], mutation_scale=20, lw=1, arrowstyle="->", color="black")
    ax.text(0+0.2, 2+0.5, -2, s = u'$V_E$', fontsize = 15)

    ax.add_artist(VD)
    ax.plot(srpgc[:,0],srpgc[:,1],srpgc[:,2], color='orange')
    ax.plot(sregc[:,0],sregc[:,1],sregc[:,2], color='g') 

    ax.scatter(u0p[0], u0p[1], u0p[2], s = 60, c='r',color='r')   
    ax.text(u0p[0]-0.3, u0p[1]+0.1, u0p[2], s = u'+', fontsize = 15)
    ax.scatter(u0e[0], u0e[1], u0e[2], c='b',color='b',s = 40)    
    ax.text(u0e[0]+0.2, u0e[1]-0.1,u0e[2], s = u'-', fontsize = 20)

    ax.set_xlabel('Eje X')
    ax.set_ylabel('Eje Y')
    ax.set_zlabel('Eje Z')       
