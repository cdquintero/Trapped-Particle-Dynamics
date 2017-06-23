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
    m=m_el                                                                   
    Re=6378137                                                                 #Radio de la tierra [m]
    n=4                                                                        #Número de Ecuaciones

    c = 299792458                                                              #Velocidad de la Luz
    K = 1e7                                                                    #Energía cinética [eV]
    K = K*e   

    pitch_angle = 6.0                                                         #Ángulo de paso inicial (grados)    
    v_mod = mth.sqrt(1-(1/(1+(K/(m*(c**2))))**2))                                          #Velocidad de la partícula     
    v_par0 = v_mod*mth.cos(pitch_angle*np.pi/180)                              #Componente paralela de la velocidad

def EC_MOV(t,u0,u1,u2,u3,j) :    
    
    St = V_Static()   
     
    # calculate 3 Cartesian components of the magnetic field
    alpha = -((St.q*St.Re*St.B0)/(St.m*St.c))/((u0**2 + u1**2 + u2**2)**(5/2.0))
    Bx = 3*u0*u2*alpha
    By = 3*u1*u2*alpha
    Bz = (2*(u2**2) - u0**2 - u1**2)*alpha
    B = mth.sqrt(Bx**2 + By**2 + Bz**2)
    mu = (St.v_mod**2-u3**2)/(2.0*B)  
    delta = 0.001

    # calculate gradBx
    gradB_x =  (-getBmod((u0+2*delta),u1,u2) + 8*getBmod((u0+delta),u1,u2) - 8*getBmod((u0-delta),u1,u2) + getBmod((u0-2*delta),u1,u2))/(12.0*delta)
    gradB_y =  (-getBmod(u0,(u1+2*delta),u2) + 8*getBmod(u0,(u1+delta),u2) - 8*getBmod(u0,(u1-delta),u2) + getBmod(u0,(u1-2*delta),u2))/(12.0*delta)
    gradB_z =  (-getBmod(u0,u1,(u2+2*delta)) + 8*getBmod(u0,u1,(u2+delta)) - 8*getBmod(u0,u1,(u2-delta)) + getBmod(u0,u1,(u2-2*delta)))/(12.0*delta)
 
    # b unit vector
    b_unit_x = Bx/B
    b_unit_y = By/B
    b_unit_z = Bz/B    

    # b unit vector cross gradB
    bxgB_x = b_unit_y*gradB_z - b_unit_z*gradB_y;
    bxgB_y = b_unit_z*gradB_x - b_unit_x*gradB_z;
    bxgB_z = b_unit_x*gradB_y - b_unit_y*gradB_x;
    
    # b unit vector inner product gradB
    dotpr = b_unit_x*gradB_x +  b_unit_y*gradB_y + b_unit_z*gradB_z;

    beta = (St.v_mod**2 + u3**2)/B**2
 
    if j==0 :
        y1 = beta*bxgB_x + u3*b_unit_x
    if j==1 :                                                      # dx/dt = u
        y1 = beta*bxgB_y + u3*b_unit_y                                                          # dy/dt = v
    if j==2 :
        y1 = beta*bxgB_z + u3*b_unit_z                                                        # dz/dt = w
    if j==3 :    
        y1 = -mu*dotpr                                   # du/dt = qom*(vel x B)_x
    return y1

#______________________________________________________________________________
   
def getBmod(x,y,z) :
    
    alpha = -((St.q*St.Re*St.B0)/(St.m*St.c))/((x**2 + y**2 + z**2)**(5/2.0))
    Bx = 3*x*z*alpha
    By = 3*y*z*alpha
    Bz = (2*(z**2) - x**2 - y**2)*alpha
    B = mth.sqrt(Bx**2 + By**2 + Bz**2)

    return B
   
#______________________________________________________________________________   
   
# Parámetros
St = V_Static()   

#Parámetros de entrada del método RK4
ti=0
tf=(St.c/St.Re)*20
ht=(St.c/St.Re)*0.001                                                    
N=int((tf-ti)/ht) 
sr=np.zeros([N,St.n])
K1=np.zeros(St.n)
K2=np.zeros(St.n)
K3=np.zeros(St.n)
K4=np.zeros(St.n)
u0=np.zeros(St.n)
t=ti
j=np.arange(St.n)
i=np.arange(N-1)

# initial velocity
# initial position: equatorial plane 4Re from Earth
u0[0] = 4
u0[1] = 0 
u0[2] = 0                                                                      #Para que el ángulo de paso sea ecuatorial z0=0
u0[3] = St.v_par0

time=np.arange(ti,tf,ht)

for jj in j:
    sr[0,j]=u0[j]

for ii in i:
    for jj in j:
        K1[jj]=ht*EC_MOV(t,sr[ii,0],sr[ii,1],sr[ii,2],sr[ii,3],jj)
    for jj in j:   
        K2[jj]=ht*EC_MOV((t+ht/2),(sr[ii,0]+K1[0]/2),(sr[ii,1]+K1[1]/2),(sr[ii,2]+K1[2]/2),(sr[ii,3]+K1[3]/2),jj)
    for jj in j:   
        K3[jj]=ht*EC_MOV((t+ht/2),(sr[ii,0]+K2[0]/2),(sr[ii,1]+K2[1]/2),(sr[ii,2]+K2[2]/2),(sr[ii,3]+K2[3]/2),jj)
    for jj in j:
        K4[jj]=ht*EC_MOV(t+ht,(sr[ii,0]+K3[0]),(sr[ii,1]+K3[1]),(sr[ii,2]+K3[2]),(sr[ii,3]+K3[3]),jj)
    for jj in j:    
        sr[ii+1,jj]=sr[ii,jj]+(K1[jj]+2*K2[jj]+2*K3[jj]+K4[jj])/6.0
    
    t=time[ii+1]

fig = plt.figure()
ax = fig.gca(projection='3d')      
       
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
ax.plot(sr[:,0], sr[:,1], sr[:,2], color='b')
  
ax.set_xlabel('X axis')
ax.set_ylabel('Y axis')
ax.set_zlabel('Z axis')       
   
plt.show()       
