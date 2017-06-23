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
import mpl_toolkits.mplot3d.axes3d

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
    ng=6     
    qomp = 1
    qome=-3    
    Bz = 1.0 
    Bx = 0.0 
    By = 0.0 
    Ex = 0.0 
    Ey = 0.0 
    Ez = 0.0 
    
    B=mth.sqrt(Bx**2+By**2+Bz**2)
    Omega=B
    vperp0=0
    vpar0=0
    vpar=0
    v_mod=0
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
    
    St.Bz = 1*u0
  
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
            K1[jj]=ht*EXBGC(t,srgc[ii,0],srgc[ii,1],srgc[ii,2],srgc[ii,3],srgc[ii,4],srgc[ii,5],jj,St)
        for jj in jg:   
            K2[jj]=ht*EXBGC((t+ht/2),(srgc[ii,0]+K1[0]/2),(srgc[ii,1]+K1[1]/2),(srgc[ii,2]+K1[2]/2),(srgc[ii,3]+K1[3]/2),(srgc[ii,4]+K1[4]/2),(srgc[ii,5]+K1[5]/2),jj,St)
        for jj in jg:   
            K3[jj]=ht*EXBGC((t+ht/2),(srgc[ii,0]+K2[0]/2),(srgc[ii,1]+K2[1]/2),(srgc[ii,2]+K2[2]/2),(srgc[ii,3]+K2[3]/2),(srgc[ii,4]+K2[4]/2),(srgc[ii,5]+K2[5]/2),jj,St)
        for jj in jg:
            K4[jj]=ht*EXBGC(t+ht,(srgc[ii,0]+K3[0]),(srgc[ii,1]+K3[1]),(srgc[ii,2]+K3[2]),(srgc[ii,3]+K3[3]),(srgc[ii,4]+K3[4]),(srgc[ii,5]+K3[5]),jj,St)
        for jj in jg:    
            srgc[ii+1,jj]=srgc[ii,jj]+(K1[jj]+2*K2[jj]+2*K3[jj]+K4[jj])/6.0
    
        t=time[ii+1]
    return(srgc)    
#______________________________________________________________________________
def EXBGC(t,u0,u1,u2,u3,u4,u5,j,St):

    St.Bz = 1*u0
    delta = 0.001
    B=mth.sqrt(St.Bx**2+St.By**2+St.Bz**2)

    # calculate gradBx
    gradB_x =  (-getBmod((u0+2*delta),u1,u2) + 8*getBmod((u0+delta),u1,u2) - 8*getBmod((u0-delta),u1,u2) + getBmod((u0-2*delta),u1,u2))/(12.0*delta)
    gradB_y =  (-getBmod(u0,(u1+2*delta),u2) + 8*getBmod(u0,(u1+delta),u2) - 8*getBmod(u0,(u1-delta),u2) + getBmod(u0,(u1-2*delta),u2))/(12.0*delta)
    gradB_z =  (-getBmod(u0,u1,(u2+2*delta)) + 8*getBmod(u0,u1,(u2+delta)) - 8*getBmod(u0,u1,(u2-delta)) + getBmod(u0,u1,(u2-2*delta)))/(12.0*delta)

    # b unit vector
    b_unit_x = St.Bx/B
    b_unit_y = St.By/B
    b_unit_z = St.Bz/B    

    # b unit vector cross gradB
    bxgB_x = b_unit_y*gradB_z - b_unit_z*gradB_y;
    bxgB_y = b_unit_z*gradB_x - b_unit_x*gradB_z;
    bxgB_z = b_unit_x*gradB_y - b_unit_y*gradB_x;
    
    beta = (St.v_mod**2-u5**2)/(2*St.qom*B**2)
    Vx = beta*bxgB_x
    Vy = beta*bxgB_y
    Vz = beta*bxgB_z

    if j==0 :
        y1 = Vx + u5*b_unit_x  
    if j==1 :                                                      # dx/dt = u
        y1 = Vy + u5*b_unit_y                                                  # dy/dt = v
    if j==2 :
        y1 = Vz + u5*b_unit_z                                                     # dz/dt = w
    if j==3 :    
        y1 = St.qom*(St.By*Vz- St.Bz*Vy)                                  # du/dt = qom*(vel x B)_x
    if j==4 :    
        y1 = St.qom*(St.Bz*Vx- St.Bx*Vz)                                 # du/dt = qom*(vel x B)_x
    if j==5 :    
        y1 = St.qom*(St.Bx*Vy- St.By*Vx)                                   # du/dt = qom*(vel x B)_x
   
    return y1
   
#______________________________________________________________________________
   
def getBmod(x,y,z) :
    
    Bx = St.Bx
    By = St.By
    Bz = 1*x
    B = mth.sqrt(Bx**2 + By**2 + Bz**2);

    return B
   
#______________________________________________________________________________   
   
def B(x,y,z):
    Bx = 0
    By = 0
    Bz = x
  
    return Bx,By,Bz,    

ndata=2
x=np.linspace(0,5,5)
y= np.linspace(0.5,2.5,ndata)
z=np.linspace(0,3,5)
Bxx,Byy,Bzz=  np.meshgrid(x,y,z)
                                                                                                
Bxxx, Byyy = np.meshgrid(np.linspace(0,3,5), np.linspace(0.5,2.5,ndata))
Bzzz = 1/Bxxx 
u, v = np.gradient(Bxxx,abs(Bzzz))

Bx=np.zeros(len(x))
By=np.zeros(len(y)) 
Bz=v[1,:]

Bx,By,Bz= np.meshgrid(Bx,By,Bz)

St = V_Static()   

#Parámetros de entrada del método RK4
ti=0
tf=12
ht=0.1                                                
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
time=np.arange(ti,tf,ht)

# Condiciones iniciales para la posición y la velocidad
u0p[0] = 1.0                                                                    #x
u0p[1] = 2.0                                                                    #y
u0p[2] = 0.0                                                                    #z
u0p[3] = 0.4                                                                    #vx
u0p[4] = 0.4                                                                    #vy
u0p[5] = 0.4    
St.vperp0=mth.sqrt(u0p[3]**2+u0p[4]**2)
St.vpar0=u0p[5]
St.vpar=St.vpar0
St.v_mod=mth.sqrt(u0p[3]**2+u0p[4]**2+u0p[5]**2)

St.qom=St.qomp
srp=deepcopy(RK4(u0p,time,ht,St))

xlp=u0p[3]/St.B
ylp=u0p[4]/St.B 

u0pgc[0] = u0p[0]+xlp                                                                    #x
u0pgc[1] = u0p[1]-ylp                                                                    #y
u0pgc[2] = u0p[2]                                                                   #z
u0pgc[3] = u0p[3]                                                                     #vx
u0pgc[4] = u0p[4]                                                                     #vx
u0pgc[5] = u0p[5]                                                                   #vx

St.vperp0=mth.sqrt(u0pgc[3]**2+u0pgc[4]**2)
St.vpar0=u0pgc[5]

St.Bz=1
srpgc=deepcopy(RK4gc(u0pgc,time,ht,St))

# Condiciones iniciales para la posición y la velocidad
rl=St.vperp0/St.B

u0e[0] = 1.0                                                                    #x
u0e[1] = 1.0                                                                    #y
u0e[2] = 0.0                                                                    #z
u0e[3] = 0.4                                                                    #vx
u0e[4] = 0.4                                                                    #vy
u0e[5] = 0.4                                                                     #vz

xle=u0e[3]/(np.abs(St.qome)*St.B)
yle=u0e[4]/(np.abs(St.qome)*St.B)  
St.vperp0=mth.sqrt(u0e[3]**2+u0e[4]**2)
St.vpar0=u0e[5]
St.vpar=St.vpar0

u0egc[0] = u0e[0]-xle                                                                   #x
u0egc[1] = u0e[1]+yle                                                                   #y
u0egc[2] = u0e[2]                                                                    #z
u0egc[3] = u0e[3]                                                                    #vx
u0egc[4] = u0e[4]                                                                     #vx
u0egc[5] = u0e[5]                                                                     #vx

St.qom=St.qome
sre=deepcopy(RK4(u0e,time,ht,St))
sregc=deepcopy(RK4gc(u0egc,time,ht,St))

for indice in zip(srp[:,3],srp[:,4],srp[:,5],k):
    vpr[indice[3]]=mth.sqrt(indice[0]**2+indice[1]**2+indice[2]**2)

for indice in zip(srpgc[:,3],srpgc[:,4],srpgc[:,5],k):
    vprgc[indice[3]]=mth.sqrt(indice[0]**2+indice[1]**2+indice[2]**2)


r=3

if  (u0p[5]==0 and u0e[5]==0):
    for indice in zip(srp[:,3],srp[:,4],k):
        vpr[indice[2]]=mth.sqrt(indice[0]**2+indice[1]**2)
    
    fig, ax = plt.subplots()
    #fig2, ax2 = plt.subplots()      
    #ax2.plot(time,vpr)
    #ax2.plot(time,srpgc[:,3], color='g')
    #ax2.set_ylim(-0.1, 1.2)
    ax.plot(srp[:,0], srp[:,1], color='r')
    ax.plot(sre[:,0], sre[:,1], color='b')
    ax.plot(srpgc[:,0],srpgc[:,1], color='orange')
    ax.plot(sregc[:,0],sregc[:,1], color='g') 
  
    ax.arrow(u0e[0]+0.3, u0e[1]-0.3, 0, -0.5,lw=0.6, head_width=0.03, head_length=0.06)  
    ax.arrow(u0p[0]-0.4, u0p[1]+0.3, 0,  0.5,lw=0.6,head_width=0.03, head_length=0.06)  
    ax.arrow(-0.1, 0, 0.2, 0.0,lw=0.6,head_width=0.05, head_length=0.03)  
    
    ax.scatter([-0.3,-0.3,-0.3,-0.3]+0.25*np.ones(4), [1.0, 1.2, 1.4, 1.6]-0.25*np.ones(4))   
    ax.scatter([-0.2,-0.2,-0.2,-0.2]+0.25*np.ones(4), [0.85, 1.15, 1.45, 1.75]-0.25*np.ones(4))   
    ax.scatter([-0.1,-0.1,-0.1,-0.1]+0.25*np.ones(4), [0.625, 1.05, 1.425, 1.85]-0.25*np.ones(4))   


    ax.text(x = -0.1, y = 0.1, s = r'$\nabla B$', fontsize = 12)
    ax.text(u0e[0]+0.35, u0e[1]-0.3, s = u'$V_{G_e}$', fontsize = 15)
    ax.text(u0p[0]-0.6, u0p[1]+0.3, s = u'$V_{G_p}$', fontsize = 15)
    ax.text(0.06+0.25, 1.2-0.25, s = u'$B$', fontsize = 15)

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
    #fig2, ax2 = plt.subplots()      
    ax = fig.gca(projection='3d')      
    ax.set_xlim3d(-0, r)
    ax.set_ylim3d(-0, r)
    ax.set_zlim3d(0, 5)  
    ax.set_aspect('equal')

    #ax2.plot(time,vpr)
    #ax2.plot(time,vprgc, color='g')
    ax.plot(srp[:,0], srp[:,1],srp[:,2],  color='r')
    ax.plot(sre[:,0], sre[:,1],sre[:,2],  color='b')
    ax.quiver(Bzz,Byy,Bxx,Bx,By,Bz,color='#4682b4', length=0.2,arrow_length_ratio=0.3,linestyles='solid',linewidths=0.5)        #Plot the magnetic field
    VDp = Arrow3D([2, 2], [2, 2.5], [0, 0], mutation_scale=20, lw=1, arrowstyle="->", color="black")
    VDe = Arrow3D([1.5, 1.5], [1, 0.5], [0, 0], mutation_scale=20, lw=1, arrowstyle="->", color="black")
    VG = Arrow3D([2.25, 2.75], [0.75, 0.75], [0, 0], mutation_scale=20, lw=1, arrowstyle="->", color="black")

    ax.add_artist(VDp)
    ax.add_artist(VDe)
    ax.add_artist(VG)
       
    ax.text(2+0.2, 1.7, 0, s = r'$V_{G_p}$', fontsize = 15)
    ax.text(1.5+0.2, 1, 0, s = r'$V_{G_e}$', fontsize = 15)
    ax.text(2.5, 0.8, 0, s = r'$\nabla B$', fontsize = 15)


    ax.plot(srpgc[:,0],srpgc[:,1],srpgc[:,2], color='orange')
    ax.plot(sregc[:,0],sregc[:,1],sregc[:,2], color='g') 

    ax.scatter(u0p[0], u0p[1], u0p[2], s = 60, c='r',color='r')   
    ax.text(u0p[0]-0.3, u0p[1]+0.1, u0p[2], s = u'+', fontsize = 15)
    ax.scatter(u0e[0], u0e[1], u0e[2], c='b',color='b',s = 40)    
    ax.text(u0e[0]+0.2, u0e[1]-0.1,u0e[2], s = u'-', fontsize = 20)
 
    ax.set_xlabel('Eje X')
    ax.set_ylabel('Eje Y')
    ax.set_zlabel('Eje Z')       
