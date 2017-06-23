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
    K = 3e7                                                              #Energía cinética [eV]
    K = K*e   

    pitch_angle = 89.0                                                         #Ángulo de paso inicial (grados)    
    v_mod = c*mth.sqrt(1-(1/(1+(K/(m*(c**2))))**2))                                                      #Ángulo de paso inicial (grados)    
    v_par0 = v_mod*mth.cos(pitch_angle*np.pi/180)                          #Componente paralela de la velocidad
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
def T_BOUNCE(domL,St):
    T_Bounce=np.empty(len(domL))    
    Li=np.arange(len(domL))     
    
    for ii in zip(Li,domL):
        CA1,CA2,CA3,CA4,t,Tau_exp = BDF(ii[1],St,num_steps,ti,ht)

        Tau=np.copy(np.array(Tau_exp))
        ll=np.arange(1,len(Tau_exp))

        for j in zip(ll,Tau_exp):
            Tau[j[0]]= Tau_exp[j[0]]-Tau_exp[j[0]-1]
     
        Tau1=np.mean(Tau)
        T_Bounce[ii[0]]=Tau1

    return T_Bounce  

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
    
# initial velocity
# initial position: equatorial plane 4Re from Earth
    u0 = L*St.Re
    u1 = 0 
    u2 = 0
    u3 = St.v_par0
    
    r.set_initial_value([u0,u1,u2,u3], ti)

    t = np.zeros((num_steps, 1))
    T_Bounce=[]
    CA1 = np.zeros((num_steps, 1))
    CA2 = np.zeros((num_steps, 1))
    CA3 = np.zeros((num_steps, 1))
    CA4 = np.zeros((num_steps, 1))

    t[0] = ti
    CA1[0] = u0
    CA2[0] = u1
    CA3[0] = u2
    CA4[0] = u3

    fl3=1
    fl4=1
    fl1=0
    hr=0
    fl5=1

    k = 1
    while r.successful() and k < num_steps:
        r.integrate(r.t + ht)
 
        # Store the results to plot later
        t[k] = r.t
        CA1[k] = r.y[0]
        CA2[k] = r.y[1]
        CA3[k] = r.y[2] 
        
        if (r.y[2]<=0 and fl3):
            hr=mth.sqrt((u0-r.y[0])**2+(u1-r.y[1])**2)             
            fl3=0
            fl1=1
  
        if (r.y[2]>=0 and fl1 and fl4 and fl5):
            sr=mth.sqrt((u0-r.y[0])**2+(u1-r.y[1])**2)           
            if sr<=hr:
                T_Bounce.append(r.t)
                fl5=0
            fl4=0
            
        if (r.y[2]<=0 and fl1 and not(fl4) and fl5):
            sr=mth.sqrt((u0-r.y[0])**2+(u1-r.y[1])**2)           
            if sr<=hr:
                T_Bounce.append(r.t)
                fl5=0
            fl4=1
            
        CA4[k] = r.y[3] 
        k += 1
        
        
    return CA1,CA2,CA3,CA4,t,T_Bounce


St = V_Static()   

domL=np.linspace(2,6,20)
T=St.K
Tau_bounce=((np.pi*St.q*St.B0*St.Re**2)/(3*domL*T))*((0.35 + 0.15*mth.sin(St.pitch_angle*np.pi/180))**(-1))

ti = 0.0
tf = 140
ht = 0.01 
L=4

u0 = L*St.Re
u1 = 0 
u2 = 0*St.Re
u3 = St.v_par0


# Number of time steps: 1 extra for initial condition
num_steps = int(np.floor((tf - ti)/ht)) + 1  
Tau_Bounce= T_BOUNCE(domL,St)


fig, ax = plt.subplots()

ax.plot(domL,Tau_bounce, color = 'b',  label= '%1.1f MeV' % 10.0)
ax.scatter(domL,Tau_Bounce, c='g', s=40, linewidths='0.7')

#ax.plot(domL,Tau_bounce2, color = 'r', label= '%1.1f MeV' % 20.0)
#ax.scatter(domL,Tau_Bounce2, c='cyan', s=40, linewidths='0.7')          

#ax.plot(domL,Tau_bounce3, color = '#ff69b4', label= '%1.1f MeV' % 30.)
#ax.scatter(domL,Tau_Bounce3, c='y', s=40, linewidths='0.7')
      
plt.rc('font', family='serif')
plt.ylabel(r'Periodo de deriva [s]', size='large',  family='monospace', weight= 'light')
plt.xlabel(r'Capa L',  size='large',  family='monospace', weight= 'light')
plt.legend()
#ax.set_ylim(1,100)
ax.set_xlim(2,6)
#ax.semilogy()    

plt.show()  

#import os

#outputname = 'Periodo_deriva_simGC3.txt'

#archivo=open(outputname,'w')

#for i in Tau_Bounce3:
#    archivo.write(str(i) + '\n')

#archivo.close()
          
