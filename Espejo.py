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
    Bz = 2 
    Bx = 0.0 
    By = 0.0 
    Ex = 0.0 
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
    
    St.Bz = u2
  
    if j==0 :
        y1 = u3                                                                # dx/dt = u
    if j==1 :                                                                 
        y1 = u4                                                                # dy/dt = v
    if j==2 :
        y1 = u5                                                                # dz/dt = w
    if j==3 :    
        y1 = St.qom*(u4*St.Bz - u5*St.By)                                    # du/dt = q/m*(vel x B)_x
    if j==4 :     
        y1 = St.qom*(u5*St.Bx - u3*St.Bz)                                         # dv/dt = q/m*(vel x B)_y
    if j==5 :    
        y1 = St.qom*(u3*St.By - u4*St.Bx)                                         # dw/dt = q/m*(vel x B)_z

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
    
    St.Bz = u2
    St.B=mth.sqrt(St.Bx**2+St.By**2+St.Bz**2)       
       
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
    Bz = z
  
    return Bx,By,Bz, 
    
def B_mod(x,y,z):
    Bx = 0
    By = 0
    Bz = z
    B=mth.sqrt(Bx**2+By**2+Bz**2)
    
    return B
  
ndata=2
x=np.linspace(1,4,4)
y= np.linspace(0,3,3)
z=np.linspace(1.5,4.5,4)
Bxx,Byy,Bzz=  np.meshgrid(x,y,z)
                                                                                                
Bxxx, Byyy = np.meshgrid(np.linspace(0,3,4), np.linspace(0.5,2.5,3))
Bzzz = 1/Bxxx 
u, v = np.gradient(Bxxx,abs(Bzzz))

Bx=np.zeros(len(x))
By=np.zeros(len(y)) 
Bz=v[1,:]

Bx,By,Bz= np.meshgrid(Bx,By,Bz)

St = V_Static()   

#Parámetros de entrada del método RK4
ti=0
tf=15 
ht=0.01                                                
N=int((tf-ti)/ht) 
sr=np.zeros([N,St.n])                                                
srgc=np.zeros([N,St.ng])                                                
vpr=np.zeros(N)                                                
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
r=5
time=np.arange(ti,tf,ht)

# Condiciones iniciales para la posición y la velocidad
u0p[0] = 1.5                                                                    #x
u0p[1] = 1.5                                                                    #y
u0p[2] = 1.5                                                                    #z
u0p[3] = 1.0                                                                    #vx
u0p[4] = 1.0                                                                    #vy
u0p[5] = 0.2    

Bmod =B_mod(u0p[0],u0p[1],u0p[2])

xlp=u0p[3]/Bmod
ylp=u0p[4]/Bmod 
St.vperp0=mth.sqrt(u0p[3]**2+u0p[4]**2)
St.vpar0=u0p[5]

u0pgc[0] = u0p[0] + xlp                                                                     #x
u0pgc[1] = u0p[1] - ylp                                                                  #y
u0pgc[2] = u0p[2]                                                                   #z
u0pgc[3] = St.vpar0                                                                    #vx

St.qom=St.qomp
srp=deepcopy(RK4(u0p,time,ht,St))
srpgc=deepcopy(RK4gc(u0pgc,time,ht,St))

# Condiciones iniciales para la posición y la velocidad

u0e[0] = 1.5                                                                  #x
u0e[1] = 1.5                                                                #y
u0e[2] = 1.5                                                                    #z
u0e[3] = 1.0                                                                    #vx
u0e[4] = 1.0                                                                    #vy
u0e[5] = 0.2                                                                     #vz

xle=u0e[3]/(np.abs(St.qome)*Bmod)
yle=u0e[4]/(np.abs(St.qome)*Bmod)  
St.vperp0=mth.sqrt(u0e[3]**2+u0e[4]**2)
St.vpar0=u0e[5]

u0egc[0] = u0e[0] - xle                                                                   #x
u0egc[1] = u0e[1] + yle                                                                   #y
u0egc[2] = u0e[2]                                                                    #z
u0egc[3] = St.vpar0                                                                     #vx

St.qom=St.qome
sre=deepcopy(RK4(u0e,time,ht,St))
sregc=deepcopy(RK4gc(u0egc,time,ht,St))

for indice in zip(srp[:,3],srp[:,4],srp[:,5],k):
    vpr[indice[3]]=mth.sqrt(indice[0]**2+indice[1]**2+indice[2]**2)

fig = plt.figure()
ax = fig.gca(projection='3d')      
ax.set_xlim3d(1, 3.5)
#ax.set_ylim3d(-r, r)
ax.set_zlim3d(1, 5)  
#ax.set_aspect('equal')

ax.plot(srp[:,0], srp[:,1], srp[:,2], color='r')
ax.plot(sre[:,0], sre[:,1], sre[:,2], color='b')
ax.plot(srpgc[:,0],srpgc[:,1],srpgc[:,2], color='orange')
ax.plot(sregc[:,0],sregc[:,1],sregc[:,2], color='g') 
ax.quiver(Bxx,Byy,Bzz,Bx,By,Bz,color='#4682b4', length=0.2,arrow_length_ratio=0.3,linestyles='solid',linewidths=0.5)        #Plot the magnetic field

 
ax.set_xlabel('Eje X')
ax.set_ylabel('Eje Y')
ax.set_zlabel('Eje Z')       
