# -*- coding: utf-8 -*-
"""
Created on Mon Dec 19 11:34:01 2016

@author: cd67
"""

#%matplotlib inline
from scipy import integrate

import math  as mth
import numpy as np                                                             #Paquete Arrays
import sympy as sym      
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class V_Static(object) :                                                       #Se define una clase que contiene todas las variables estáticas del programa
    e = 1.602176565e-19                                                        #Carga de las partículas
    m_pr = 1.672621777e-27                                                     #Masa del protón [kg]
    m_el = 9.10938291e-31                                                      #Masa del electrón [kg]
#    m_el = m_pr/1800 
    
    B0=3.07e-5                                                                 #Campo magnético ecuatorial
    q=e                                                                             
    m=m_pr                                                                  
    Re=6378137                                                                 #Radio de la tierra [m]
    n=4                                                                        #Número de Ecuaciones

    c = 299792458                                                              #Velocidad de la Luz
    K = 1e7                                                                  #Energía cinética [eV]
    K = K*e   

    pitch_angle = 45.0                                                         #Ángulo de paso inicial (grados)    
    v_mod = c/mth.sqrt(1-(1/(1+(K/(m*(c**2))))**2))                                           #Velocidad de la partícula     
    v_par0 = v_mod*mth.cos(pitch_angle*np.pi/180)                              #Componente paralela de la velocidad
#______________________________________________________________________________   

def EC_MOV(t,u) :    
    
    St = V_Static()   
       
    # calculate 3 Cartesian components of the magnetic field
    alpha = -(St.B0*St.Re**3)/((u[0]**2 + u[1]**2 + u[2]**2)**(2.5))
    Bx = 3*u[0]*u[2]*alpha
    By = 3*u[1]*u[2]*alpha
    Bz = (2*(u[2]**2) - u[0]**2 - u[1]**2)*alpha
    B = mth.sqrt(Bx**2 + By**2 + Bz**2)
    mu = St.m*(St.v_mod**2-u[3]**2)/(2.0*B)                                      # magnetic moment is an adiabatic invariant

    delta = 0.01*St.Re
    # calculate gradBx
    gradB_x =  (-getBmod((u[0]+2*delta),u[1],u[2]) + 8*getBmod((u[0]+delta),u[1],u[2]) - 8*getBmod((u[0]-delta),u[1],u[2]) + getBmod((u[0]-2*delta),u[1],u[2]))/(12.0*delta)
    gradB_y =  (-getBmod(u[0],(u[1]+2*delta),u[2]) + 8*getBmod(u[0],(u[1]+delta),u[2]) - 8*getBmod(u[0],(u[1]-delta),u[2]) + getBmod(u[0],(u[1]-2*delta),u[2]))/(12.0*delta)
    gradB_z =  (-getBmod(u[0],u[1],(u[2]+2*delta)) + 8*getBmod(u[0],u[1],(u[2]+delta)) - 8*getBmod(u[0],u[1],(u[2]-delta)) + getBmod(u[0],u[1],(u[2]-2*delta)))/(12.0*delta)

    # b unit vector
    b_unit_x = Bx/B
    b_unit_y = By/B
    b_unit_z = Bz/B    

    # b unit vector cross gradB
    bxgB_x = b_unit_y*gradB_z - b_unit_z*gradB_y
    bxgB_y = b_unit_z*gradB_x - b_unit_x*gradB_z
    bxgB_z = b_unit_x*gradB_y - b_unit_y*gradB_x
    
    # b unit vector inner product gradB
    dotpr = b_unit_x*gradB_x +  b_unit_y*gradB_y + b_unit_z*gradB_z;

    beta = (St.m/(2*St.q*B**2))*(St.v_mod**2 + u[3]**2)
    
    dudt=np.zeros((St.n,1))

    dudt[0] = beta*bxgB_x + u[3]*b_unit_x 
    dudt[1] = beta*bxgB_y + u[3]*b_unit_y 
    dudt[2] = beta*bxgB_z + u[3]*b_unit_z                                                            # dz/dt = w 
    dudt[3] = -mu/St.m*dotpr
    
    return dudt

def EC_MOVD(t,u) :    
    
    St = V_Static()   
       
    # calculate 3 Cartesian components of the magnetic field
    alpha = -(St.B0*St.Re**3)/((u[0]**2 + u[1]**2 + u[2]**2)**(2.5))
    Bx = 3*u[0]*u[2]*alpha
    By = 3*u[1]*u[2]*alpha
    Bz = (2*(u[2]**2) - u[0]**2 - u[1]**2)*alpha
    B = mth.sqrt(Bx**2 + By**2 + Bz**2)

    delta = 0.01*St.Re
    # calculate gradBx
    gradB_x =  (-getBmod((u[0]+2*delta),u[1],u[2]) + 8*getBmod((u[0]+delta),u[1],u[2]) - 8*getBmod((u[0]-delta),u[1],u[2]) + getBmod((u[0]-2*delta),u[1],u[2]))/(12.0*delta)
    gradB_y =  (-getBmod(u[0],(u[1]+2*delta),u[2]) + 8*getBmod(u[0],(u[1]+delta),u[2]) - 8*getBmod(u[0],(u[1]-delta),u[2]) + getBmod(u[0],(u[1]-2*delta),u[2]))/(12.0*delta)
    gradB_z =  (-getBmod(u[0],u[1],(u[2]+2*delta)) + 8*getBmod(u[0],u[1],(u[2]+delta)) - 8*getBmod(u[0],u[1],(u[2]-delta)) + getBmod(u[0],u[1],(u[2]-2*delta)))/(12.0*delta)

    # b unit vector
    b_unit_x = Bx/B
    b_unit_y = By/B
    b_unit_z = Bz/B    

    # b unit vector cross gradB
    bxgB_x = b_unit_y*gradB_z - b_unit_z*gradB_y
    bxgB_y = b_unit_z*gradB_x - b_unit_x*gradB_z
    bxgB_z = b_unit_x*gradB_y - b_unit_y*gradB_x
    
    # b unit vector inner product gradB
    beta = (St.m/(2*St.q*B**2))*(St.v_mod**2 + u[3]**2)
    
    dvdt=np.zeros((St.n,1))
    
    dvdt[0] = beta*bxgB_x 
    dvdt[1] = beta*bxgB_y 
    dvdt[2] = beta*bxgB_z                                                            # dz/dt = w 
    dvdt[3] = 0
    return dvdt


    
#______________________________________________________________________________

def getBmod(x,y,z) :
    
    alpha = -(St.B0*St.Re**3)/((x**2 + y**2 + z**2)**(2.5))
    Bx = 3*x*z*alpha
    By = 3*y*z*alpha
    Bz = (2*(z**2) - x**2 - y**2)*alpha
    B = mth.sqrt(Bx**2 + By**2 + Bz**2)

    return B
#______________________________________________________________________________

def B(x,y,z):
    alpha = -1/((x**2 + y**2 + z**2)**(5/2.0))
    Bx = (3*x*z*alpha)
    By = 3*y*z*alpha
    Bz = (2*(z**2) - x**2 - y**2)*alpha
    return Bx,By,Bz
  
def BDF(L,St,num_steps,ti,ht):
    
    r = integrate.ode(EC_MOV).set_integrator('vode', method='bdf')
    rd= integrate.ode(EC_MOVD).set_integrator('vode', method='bdf')   
    

# initial velocity
# initial position: equatorial plane 4Re from Earth
    u0 = L*St.Re
    u1 = 0 
    u2 = 0*St.Re
    u3 = St.v_par0
    
    r.set_initial_value([u0,u1,u2,u3], ti)
    rd.set_initial_value([u0,u1,u2,u3], ti)

    t = np.zeros((num_steps, 1))
    CA1 = np.zeros((num_steps, 1))
    CA2 = np.zeros((num_steps, 1))
    CA3 = np.zeros((num_steps, 1))
    CA4 = np.zeros((num_steps, 1))
    CA1d = np.zeros((num_steps, 1))
    CA2d = np.zeros((num_steps, 1))
    CA3d = np.zeros((num_steps, 1))
    CA4d = np.zeros((num_steps, 1))

    t[0] = ti
    CA1[0] = u0
    CA2[0] = u1
    CA3[0] = u2
    CA4[0] = u3
    CA1d[0] = u0
    CA2d[0] = u1
    CA3d[0] = u2
    CA4d[0] = u3


    k = 1
    while r.successful() and k < num_steps:
        r.integrate(r.t + ht)
 
        # Store the results to plot later
        t[k] = r.t
        CA1[k] = r.y[0]
        CA2[k] = r.y[1]
        CA3[k] = r.y[2]        
        CA4[k] = r.y[3] 
        k += 1
        
    k = 1
    while rd.successful() and k < num_steps:
        rd.integrate(rd.t + ht)
 
        # Store the results to plot later
        t[k] = rd.t
        CA1d[k] = rd.y[0]
        CA2d[k] = rd.y[1]
        CA3d[k] = rd.y[2]        
        CA4d[k] = rd.y[3] 
        k += 1
        
    return CA1,CA2,CA3,CA4,CA1d,CA2d,CA3d,CA4d,t


xx = np.linspace(-4,4,9)
yy = np.linspace(-4,4,9)
zz = np.linspace(-4,4,9)

xx,yy,zz = np.meshgrid(xx,yy,zz)

# Plot of the fields  
Bx,By,Bz = B(xx,yy,zz)     

St = V_Static()   

ti = 0.0
tf = 50
ht = 0.01 
L=4

u0 = L*St.Re
u1 = 0 
u2 = 0*St.Re
u3 = St.v_par0


# Number of time steps: 1 extra for initial condition
num_steps = int(np.floor((tf - ti)/ht)) + 1  

CA1,CA2,CA3,CA4,CA1d,CA2d,CA3d,CA4d,t= BDF(L,St,num_steps,ti,ht)


fig = plt.figure()
ax = fig.gca(projection='3d')      
       
w = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)

x = np.outer(np.cos(w), np.sin(v))
y = np.outer(np.sin(w), np.sin(v))
z = np.outer(np.ones(np.size(w)), np.cos(v))

xsol = CA1.reshape(num_steps)
ysol = CA2.reshape(num_steps)
zsol = CA3.reshape(num_steps)   
vsol = CA4.reshape(num_steps)     

xdsol = CA1d.reshape(num_steps)
ydsol = CA2d.reshape(num_steps)
zdsol = CA3d.reshape(num_steps)   
vdsol = CA4d.reshape(num_steps)     

ri=-5
rf=5
ax.set_xlim3d(ri,rf)
ax.set_ylim3d(ri,rf)
ax.set_zlim3d(ri,rf)  
ax.plot_surface(x, y, z, rstride=6, cstride=6, linewidth=0.5, color='cornflowerblue') 
ax.set_aspect('equal')
ax.plot(xsol/St.Re, ysol/St.Re, zsol/St.Re, color='#a0522d')
#ax.plot(xdsol[0:30000]/St.Re, ydsol[0:30000]/St.Re, zdsol[0:30000]/St.Re, color='#ffa500')
ax.quiver(xx,yy,zz,Bx,By,Bz,color='#4682b4', length=0.5,arrow_length_ratio=0.4,linestyles='dotted',linewidths=0.5)        #Plot the magnetic field

ax.scatter(u0/St.Re, u1/St.Re, u2/St.Re, s = 20, c='r',color='r')   
ax.text(u0/St.Re-0.3, u1/St.Re+0.1, u2/St.Re, s = u'+', fontsize = 15)
 
ax.set_xlabel('Eje X')
ax.set_ylabel('Eje Y')
ax.set_zlabel('Eje Z')  
  
plt.show()  
     
