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
    n=6 
    ng=4     
    qomp = 1
    qome=-3  
    Bz = 0.7 
    Bx = 0.0 
    By = 0.0 
    Ex = 0.0 
    Ey = 0.0 
    Ez = 0.0 
    
    v_mod=0
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
            K1[jj]=ht*MOV_B(t,sr[ii,0],sr[ii,1],sr[ii,2],sr[ii,3],sr[ii,4],sr[ii,5],jj,St)
        for jj in j:   
            K2[jj]=ht*MOV_B((t+ht/2),(sr[ii,0]+K1[0]/2),(sr[ii,1]+K1[1]/2),(sr[ii,2]+K1[2]/2),(sr[ii,3]+K1[3]/2),(sr[ii,4]+K1[4]/2),(sr[ii,5]+K1[5]/2),jj,St)
        for jj in j:   
            K3[jj]=ht*MOV_B((t+ht/2),(sr[ii,0]+K2[0]/2),(sr[ii,1]+K2[1]/2),(sr[ii,2]+K2[2]/2),(sr[ii,3]+K2[3]/2),(sr[ii,4]+K2[4]/2),(sr[ii,5]+K2[5]/2),jj,St)
        for jj in j:
            K4[jj]=ht*MOV_B(t+ht,(sr[ii,0]+K3[0]),(sr[ii,1]+K3[1]),(sr[ii,2]+K3[2]),(sr[ii,3]+K3[3]),(sr[ii,4]+K3[4]),(sr[ii,5]+K3[5]),jj,St)
        for jj in j:    
            sr[ii+1,jj]=sr[ii,jj]+(K1[jj]+2*K2[jj]+2*K3[jj]+K4[jj])/6.0
    
        t=time[ii+1]
    return(sr)
#______________________________________________________________________________
def MOV_B(t,u0,u1,u2,u3,u4,u5,j,St):
    v_mod=mth.sqrt(u3**2+u4**2+u5**2)
    pitch_angle=90
    v_perp= v_mod*mth.sin(pitch_angle*np.pi/180) 
    rl=v_perp/St.B  
  
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

def RGC(srp):
    # b unit vector
    b_unit_x = St.Bx/St.B
    b_unit_y = St.By/St.B
    b_unit_z = St.Bz/St.B    

    # b unit vector cross gradB
    bxgB_x = srp[4]*b_unit_z - srp[5]*b_unit_y;
    bxgB_y = srp[5]*b_unit_x - srp[3]*b_unit_z;
    bxgB_z = srp[3]*b_unit_y - srp[4]*b_unit_x;
    
    xgc=srp[0]+bxgB_x/St.B
    ygc=srp[1]+bxgB_y/St.B
    zgc=srp[2]+bxgB_z/St.B
    return xgc,ygc,zgc

def RK4gc(u0,time,ht,St):

    t=time[0]
    for jj in jg:
        sr[0,jj]=u0[jj]

    for ii in i:
        for jj in jg:
            K1[jj]=ht*MOV_B_GC(t,sr[ii,0],sr[ii,1],sr[ii,2],sr[ii,3],jj,St)
        for jj in jg:   
            K2[jj]=ht*MOV_B_GC((t+ht/2),(sr[ii,0]+K1[0]/2),(sr[ii,1]+K1[1]/2),(sr[ii,2]+K1[2]/2),(sr[ii,3]+K1[3]/2),jj,St)
        for jj in jg:   
            K3[jj]=ht*MOV_B_GC((t+ht/2),(sr[ii,0]+K2[0]/2),(sr[ii,1]+K2[1]/2),(sr[ii,2]+K2[2]/2),(sr[ii,3]+K2[3]/2),jj,St)
        for jj in jg:
            K4[jj]=ht*MOV_B_GC(t+ht,(sr[ii,0]+K3[0]),(sr[ii,1]+K3[1]),(sr[ii,2]+K3[2]),(sr[ii,3]+K3[3]),jj,St)
        for jj in jg:    
            sr[ii+1,jj]=sr[ii,jj]+(K1[jj]+2*K2[jj]+2*K3[jj]+K4[jj])/6.0
    
        t=time[ii+1]
    return(sr)

#______________________________________________________________________________
def MOV_B_GC(t,u0,u1,u2,u3,j,St):

    # b unit vector
    b_unit_x = St.Bx/St.B
    b_unit_y = St.By/St.B
    b_unit_z = St.Bz/St.B  
    
    if j==0 :
        y2 = u3*b_unit_x                                                                # dx/dt = u
    if j==1 :                                                                 
        y2 = u3*b_unit_y                                                                # dy/dt = v
    if j==2 :
        y2 = u3*b_unit_z                                                                # dz/dt = w
    if j==3 :    
        y2 = 0                                    # du/dt = q/m*(vel x B)_x

    return y2
    
#______________________________________________________________________________

def B(x,y,z):
    Bx = 0
    By = 0
    Bz = np.abs(z)
    return Bx,By,Bz

xx = np.linspace(-4,4,3)
yy = np.linspace(-4,4,3)
zz = np.linspace(-4,4,3)

xx,yy,zz = np.meshgrid(xx,yy,zz)

# Plot of the fields
Bx,By,Bz = B(xx,yy,zz)              

St = V_Static()   

#Parámetros de entrada del método RK4
ti=0
tf=8  
ht=0.01                                                
N=int((tf-ti)/ht) 
sr=np.zeros([N,St.n]) 
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
time=np.arange(ti,tf,ht)

# Condiciones iniciales para la posición y la velocidad
u0p[0] = -2.0                                                                    #x
u0p[1] = 0.0                                                                    #y
u0p[2] = -2.0                                                                    #z
u0p[3] = 0.5                                                                    #vx
u0p[4] = 0.5                                                                    #vy
u0p[5] = 0.0    

St.v_mod=mth.sqrt(u0p[3]**2+u0p[4]**2+u0p[5]**2)
pitch_angle=90
v_perp= mth.sqrt(u0p[3]**2+u0p[4]**2)
v_par= u0p[5]

xlp=u0p[3]/St.B
ylp=u0p[4]/St.B 
St.vperp0=mth.sqrt(u0p[3]**2+u0p[4]**2)
St.vpar0=u0p[5]

u0pgc[0]= u0p[0]+xlp                                                                    #x
u0pgc[1]= u0p[1]-ylp                                                                    #y
u0pgc[2]= u0p[2]                                                                  #z
u0pgc[3]= St.vpar0 

St.qom=St.qomp
srp=deepcopy(RK4(u0p,time,ht,St))
srpgc=deepcopy(RK4gc(u0pgc,time,ht,St))

# Condiciones iniciales para la posición y la velocidad

u0e[0] = 2.0                                                                    #x
u0e[1] = 0.0                                                                    #y
u0e[2] = -2.0                                                                    #z
u0e[3] = 0.5                                                                    #vx
u0e[4] = 0.5                                                                    #vy
u0e[5] = 0.0      

xle=u0e[3]/(np.abs(St.qome)*St.B)
yle=u0e[4]/(np.abs(St.qome)*St.B)  
St.vperp0=mth.sqrt(u0e[3]**2+u0e[4]**2)
St.vpar0=u0e[5]

u0egc[0]= u0e[0]-xle                                                                    #x
u0egc[1]= u0e[1]+yle                                                                    #y
u0egc[2]= u0e[2]                                                                   #z
u0egc[3]= St.vpar0 

St.qom=St.qome
sre=deepcopy(RK4(u0e,time,ht,St))
sregc=deepcopy(RK4gc(u0egc,time,ht,St))

for indice in zip(srp[:,3],srp[:,4],srp[:,5],k):
     vpr[indice[3]]=mth.sqrt(indice[0]**2+indice[1]**2+indice[2]**2)


if  (u0p[5]==0 and u0e[5]==0):
    
    fig, ax = plt.subplots()
    fig2, ax2 = plt.subplots()      
    ax.plot(srp[:,0], srp[:,1], color='r')
    ax.plot(sre[0:250,0], sre[0:250,1], color='b')
    ax.set_xlim(-r, r)
    ax.set_ylim(-r, r)
    ax.set_aspect('equal')
  
    ax2.plot(time,vpr, color='r')
    ax2.plot(time,srpgc[:,3], color='g')
    ax2.set_ylim(-0.1, 0.8)

    ax.scatter(0, 0)   
    ax.text(x = 0.0, y = 0.1, s = u'B')
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
    fig = plt.figure()
    ax = fig.gca(projection='3d')      
    fig2, ax2 = plt.subplots()          
    
    ax.set_xlim3d(-r, r)
    ax.set_ylim3d(-r, r)
    ax.set_zlim3d(-r, r)  
    ax.set_aspect('equal')
    ax.plot(srp[:,0], srp[:,1], srp[:,2], color='r')
    ax.plot(sre[:,0], sre[:,1], sre[:,2], color='b')
   
    ax2.plot(time,vpr, color='r')
    ax2.plot(time,srpgc[:,3], color='g')
    ax2.set_ylim(-0.1, 0.8)
    ax.quiver(xx,yy,zz,Bx,By,Bz,color='#4682b4', length=0.7,arrow_length_ratio=0.4,linestyles='solid',linewidths=0.5)        #Plot the magnetic field
    ax.plot(srpgc[:,0],srpgc[:,1],srpgc[:,2], color='orange')
    ax.plot(sregc[:,0],sregc[:,1],sregc[:,2], color='g') 

    ax.scatter(u0p[0], u0p[1], u0p[2], s = 60, c='r',color='r')   
    ax.text(u0p[0]-0.3, u0p[1]+0.1, u0p[2], s = u'+', fontsize = 15)
    ax.scatter(u0e[0], u0e[1], u0e[2], c='b',color='b',s = 40)    
    ax.text(u0e[0]+0.2, u0e[1]-0.1,u0e[2], s = u'-', fontsize = 20)

    ax.set_xlabel('Eje X')
    ax.set_ylabel('Eje Y')
    ax.set_zlabel('Eje Z')       
