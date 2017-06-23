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
from functools import partial
import scipy.optimize
import matplotlib.pyplot as plt

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
    K = 1e6                                                                  #Energía cinética [eV]
    K = K*e   
    N=15
    pitch_angle = np.linspace(0.1, 90.0, num=N)                                                 #Ángulo de paso inicial (grados)    
    v_mod = c*mth.sqrt(1-(1/(1+(K/(m*(c**2))))**2))
    v_par0 = v_mod*np.cos(pitch_angle*np.pi/180)                              #Componente paralela de la velocidad
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

    F=-50000000
    dudt[0] = beta*bxgB_x + u[3]*b_unit_x 
    dudt[1] = beta*bxgB_y + u[3]*b_unit_y 
    dudt[2] = beta*bxgB_z + u[3]*b_unit_z                                                            # dz/dt = w 
    dudt[3] = -mu/St.m*dotpr
    
    return dudt
    
#______________________________________________________________________________

def getBmod(x,y,z) :
    
    alpha = -(St.B0*St.Re**3)/((x**2 + y**2 + z**2)**(2.5))
    Bx = 3*x*z*alpha
    By = 3*y*z*alpha
    Bz = (2*(z**2) - x**2 - y**2)*alpha
    B = mth.sqrt(Bx**2 + By**2 + Bz**2)

    return B
#______________________________________________________________________________

  
def BDF(L,St,num_steps,ti,ht):
    
    r = integrate.ode(EC_MOV).set_integrator('vode', method='bdf')
    lx = np.zeros(St.N)
    ly = np.zeros(St.N)
    lz = np.zeros(St.N)
    
    
    
# initial velocity
# initial position: equatorial plane 4Re from Earth
   
    for i in zip(St.v_par0, range(0, St.N)):
        
        u0 = L*St.Re
        u1 = 0 
        u2 = 0*St.Re
        u3 = i[0]
    
        r.set_initial_value([u0,u1,u2,u3], ti)
        
        fl1=0
        fl2=1
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
            if (r.y[3]>0):
                fl1=1
            if (r.y[3]<0 and fl1 and fl2):
                latx=r.y[0]
                laty=r.y[1]
                latz=r.y[2]
                fl1=0
                fl2=0
            k += 1
     
        lx[i[1]] = latx
        ly[i[1]] = laty
        lz[i[1]] = latz

        
    return lx,ly,lz

def z(x, y):
    return (np.sin(x))**2-(np.cos(y))**6/(1+3*np.sin(y)**2)**0.5


St = V_Static()   

ti = 0.0
tf = 5
ht = 0.01 
L=4

u0 = L*St.Re
u1 = 0 
u2 = 0*St.Re
u3 = St.v_par0


# Number of time steps: 1 extra for initial condition
num_steps = int(np.floor((tf - ti)/ht)) + 1  


lx,ly,lz= BDF(L,St,num_steps,ti,ht)

r=np.zeros(St.N)
theta=np.zeros(St.N)
lat=np.zeros(St.N)
for i in range(0,St.N):
    r[i]=mth.sqrt(lx[i]**2 + ly[i]**2)
    theta[i]=mth.atan(r[i]/lz[i])
    theta[i]=theta[i]*180/np.pi
    lat[i]=90-theta[i]

fig, ax = plt.subplots()

x_window = 0.1*np.pi/180, 90*np.pi/180
y_window = 0.1*np.pi/180, 90*np.pi/180

xs = []
ys = []
for x in np.linspace(*x_window, num=15):
    try:
        # A more efficient technique would use the last-found-y-value as a 
        # starting point
        y = scipy.optimize.brentq(partial(z, x), *y_window)
    except ValueError:
        # Should we not be able to find a solution in this window.
        pass
    else:
        xs.append(x*180/np.pi)
        ys.append(y*180/np.pi)

ax.plot(xs, ys)
ax.scatter(St.pitch_angle,lat, c='g', s=40, linewidths='0.7')
   
plt.rc('font', family='serif')
plt.ylabel(r'Latitud del punto espejo en grados', size='large',  family='monospace', weight= 'light')
plt.xlabel(r'Angulo de paso ecuatorial en grados',  size='large',  family='monospace', weight= 'light')
plt.legend()
ax.set_ylim(0.1,90)
ax.set_xlim(0.1,90)

plt.show()  
    
import os

outputname = 'P_esp_teo.txt'

archivo=open(outputname,'w')

for i in ys:
    archivo.write(str(i) + '\n')

archivo.close()
          
