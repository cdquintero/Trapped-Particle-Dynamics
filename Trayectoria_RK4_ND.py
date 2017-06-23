# -*- coding: utf-8 -*-
"""
Created on Mon Dec 19 11:34:01 2016

@author: cd67
"""

#%matplotlib inline
import math  as mth
import numpy as np                                                             #Paquete Arrays
import matplotlib.pyplot as plt 
import sympy as sym      
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D 

class V_Static(object) :                                                       #Se define una clase que contiene todas las variables estáticas del programa
    e = 1.602176565e-19                                                        #Carga de las partículas
    m_pr = 1.672621777e-27                                                     #Masa del protón [kg]
    m_el = 9.10938291e-31                                                      #Masa del electrón [kg]
    
    B0=3.07e-5                                                                 #Campo magnético ecuatorial
    q=e                                                                             
    m=m_pr                                                                  
    Re=6378137                                                                 #Radio de la tierra [m]
    n=6                                                                        #Número de Ecuaciones

    c = 299792458                                                              #Velocidad de la Luz
    K = 1e7                                                                    #Energía cinética [eV]
    K = K*e;
    K=K/(m*c**2)

    pitch_angle = 30.0                                                         #Ángulo de paso inicial (grados)    
    v_mod = mth.sqrt(1-(1/(1+K))**2)                                          #Velocidad de la partícula     
    v_par0 = v_mod*mth.cos(pitch_angle*np.pi/180)                            #Componente paralela de la velocidad

 #    mu=(m*v_perp0**2)/(2*B0)                                                   #Momento magnético
    

def EC_MOV(t,u0,u1,u2,u3,u4,u5,j) :    

    St = V_Static()    
    
    # Componentes cartesianas del campo magnético
    alpha = -((St.q*St.Re*St.B0)/(St.m*St.c))/((u0**2 + u1**2 + u2**2)**(5/2.0))
    Bx = 3*u0*u2*alpha
    By = 3*u1*u2*alpha
    Bz = (2*(u2**2) - u0**2 - u1**2)*alpha
    B = mth.sqrt(Bx**2 + By**2 + Bz**2)
    mu = (St.v_mod**2-u3**2)/(2.0*B)  
    v = mth.sqrt(u3**2 + u4**2 + u5**2)
    E = (1/(1-(v**2/St.c**2)))

    Tau_c = 1/B   
    
    if j==0 :
        y1 = u3                                                                # dx/dt = u
    if j==1 :                                                                 
        y1 = u4                                                                # dy/dt = v
    if j==2 :
        y1 = u5                                                                # dz/dt = w
    if j==3 :    
        y1 = (u4*Bz - u5*By)                                         # du/dt = q/m*(vel x B)_x
    if j==4 :     
        y1 = (u5*Bx - u3*Bz)                                         # dv/dt = q/m*(vel x B)_y
    if j==5 :    
        y1 = (u3*By - u4*Bx)                                         # dw/dt = q/m*(vel x B)_z

    return y1, Tau_c, mu, E

#______________________________________________________________________________

def B(x,y,z,St):
    alpha = -((St.q*St.Re*St.B0)/(St.m*St.c))/((x**2 + y**2 + z**2)**(5/2.0))
    Bx = (3*x*z*alpha)
    By = 3*y*z*alpha
    Bz = (2*(z**2) - x**2 - y**2)*alpha
    return Bx,By,Bz


St = V_Static()   

xx = np.linspace(-4,4,9)
yy = np.linspace(-4,4,9)
zz = np.linspace(-4,4,9)

xx,yy,zz = np.meshgrid(xx,yy,zz)

# Plot of the fields
Bx,By,Bz = B(xx,yy,zz,St)              


#Parámetros de entrada del método RK4
ti=0
tf=(St.c/St.Re)*20
ht=(St.c/St.Re)*0.001                                                                      #Paso temporal                                        
N=int((tf-ti)/ht)                                                              #Número de pasos temporales 
sr=np.zeros([N,St.n])                                                      
K1=np.zeros(St.n)
K2=np.zeros(St.n)
K3=np.zeros(St.n)
K4=np.zeros(St.n)
u0=np.zeros(St.n)
t=ti
j=np.arange(St.n)
i=np.arange(N-1)
Tau=np.zeros(N) 
Mu=np.zeros(N) 
En=np.zeros(N) 

# initial velocity
# initial position: equatorial plane 4Re from Earth
u0[0] = 4
u0[1] = 0 
u0[2] = 0                                                                      #Para que el ángulo de paso sea ecuatorial z0=0
u0[3] = 0.0
u0[5] = St.v_par0
u0[4] = mth.sqrt(

St.v_mod**2-u0[5]**2)
Tau[0]= St.q*St.B0/St.m

time=np.arange(ti,tf,ht)

for jj in j:
    sr[0,j]=u0[j]

for ii in i:
    for jj in j:
        y,Tau_c, mu, E= EC_MOV(t,sr[ii,0],sr[ii,1],sr[ii,2],sr[ii,3],sr[ii,4],sr[ii,5],jj)
        K1[jj]=ht*y
    for jj in j:
        y,Tau_c, mu, E=EC_MOV((t+ht/2),(sr[ii,0]+K1[0]/2),(sr[ii,1]+K1[1]/2),(sr[ii,2]+K1[2]/2),(sr[ii,3]+K1[3]/2),(sr[ii,4]+K1[4]/2),(sr[ii,5]+K1[5]/2),jj)
        K2[jj]=ht*y
    for jj in j:
        y,Tau_c, mu, E=EC_MOV((t+ht/2),(sr[ii,0]+K2[0]/2),(sr[ii,1]+K2[1]/2),(sr[ii,2]+K2[2]/2),(sr[ii,3]+K2[3]/2),(sr[ii,4]+K2[4]/2),(sr[ii,5]+K2[5]/2),jj)
        K3[jj]=ht*y
    for jj in j:
        y,Tau_c, mu, E=EC_MOV(t+ht,(sr[ii,0]+K3[0]),(sr[ii,1]+K3[1]),(sr[ii,2]+K3[2]),(sr[ii,3]+K3[3]),(sr[ii,4]+K3[4]),(sr[ii,5]+K3[5]),jj)
        K4[jj]=ht*y
    for jj in j:    
        sr[ii+1,jj]=sr[ii,jj]+(K1[jj]+2*K2[jj]+2*K3[jj]+K4[jj])/6.0
    Tau[ii+1]=Tau_c
    Mu[ii+1]=mu
    En[ii+1]=E    
    
    t=time[ii+1]

fig = plt.figure()
ax = fig.gca(projection='3d')      
Mu=Mu*St.q*St.Re*St.c       
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)

x = np.outer(np.cos(u), np.sin(v))
y = np.outer(np.sin(u), np.sin(v))
z = np.outer(np.ones(np.size(u)), np.cos(v))


ax.set_xlim3d(-4.2, 4.2)
ax.set_ylim3d(-4.2, 4.2)
ax.set_zlim3d(-4.2, 4.2)  
ax.plot_surface(x, y, z, rstride=6, cstride=6, linewidth=0.5, color='cornflowerblue') 
ax.set_aspect('equal')
ax.plot(sr[:,0], sr[:,1], sr[:,2], color='r')
ax.quiver(xx,yy,zz,Bx,By,Bz,color='#4682b4', length=0.5,arrow_length_ratio=0.4,linestyles='dotted',linewidths=0.5,alpha=0.8)        #Plot the magnetic field

ax.set_xlabel('X axis')
ax.set_ylabel('Y axis')
ax.set_zlabel('Z axis')       

fig2, ax2 = plt.subplots()
#ax2.plot(time,Mu*St.q*St.Re*St.c)
ax2.plot(time,Mu)
ax2.set_ylim(-.1, .1)
ax2.set_xlabel('$Tiempo$')
ax2.set_ylabel('$Momento\quad magnetico$')

fig3, ax3 = plt.subplots()
#ax3.plot(time,En*St.m*St.c**2/St.e)
ax3.plot(time,En)
ax3.set_ylim(0, 2)
ax3.set_xlabel('$Tiempo$')
ax3.set_ylabel('$Energia\quad cinetica$')
#ax3.semilogy()
   
plt.show()       



