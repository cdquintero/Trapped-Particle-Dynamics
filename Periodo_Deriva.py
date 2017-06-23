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
import matplotlib.lines as mlines
import matplotlib.patches as mpatches


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
    K1 = 1e7                                                            #Energía cinética [eV]
    K1 = K1*e   
    K2 = 2e7                                                              #Energía cinética [eV]
    K2 = K2*e   
    K3 = 3e7                                                              #Energía cinética [eV]
    K3 = K3*e   

    pitch_angle = 90.0  
    v_mod = 0                                                       #Ángulo de paso inicial (grados)    
    v_mod1 = c*mth.sqrt(1-(1/(1+(K1/(m*(c**2))))**2))                                           #Velocidad de la partícula     
    v_mod2 = c*mth.sqrt(1-(1/(1+(K2/(m*(c**2))))**2))                                           #Velocidad de la partícula     
    v_mod3 = c*mth.sqrt(1-(1/(1+(K3/(m*(c**2))))**2))                                           #Velocidad de la partícula     
#______________________________________________________________________________   

def EC_MOV(t,u) :    

    St = V_Static()    
    
    # calculate 3 Cartesian components of the magnetic field
    alpha = -St.B0*St.Re**3/((u[0]**2 + u[1]**2 + u[2]**2)**(2.5))
    Bx = 3*u[0]*u[2]*alpha
    By = 3*u[1]*u[2]*alpha
    Bz = (2*(u[2]**2) - u[0]**2 - u[1]**2)*alpha

    dudt=np.zeros((St.n,1))

    F=-50000000
    dudt[0] = u[3] 
    dudt[1] = u[4]
    dudt[2] = u[5]                                                         # dz/dt = w 
    dudt[3] = St.q/St.m*(u[4]*Bz - u[5]*By)
    dudt[4] = St.q/St.m*(u[5]*Bx - u[3]*Bz) 
    dudt[5] = St.q/St.m*(u[3]*By - u[4]*Bx)
   
    return dudt
    
#______________________________________________________________________________

def T_BOUNCE(domL,vmod,St):
    T_Bounce=np.empty([len(domL),len(vmod)])    
    Li=np.arange(len(domL))     
    Ti=np.arange(len(vmod))    
    
    for kk in zip(Ti,vmod): 
        for ii in zip(Li,domL):
            CA1,CA2,CA3,CA4,CA5,CA6,t,Tau_exp = BDF(ii[1],kk[1],St,num_steps,ti,ht)

            Tau=np.copy(np.array(Tau_exp))
            ll=np.arange(1,len(Tau_exp))

            for j in zip(ll,Tau_exp):
                Tau[j[0]]= Tau_exp[j[0]]-Tau_exp[j[0]-1]
     
            Tau1=np.mean(Tau)
            T_Bounce[ii[0],kk[0]]=Tau1

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
#______________________________________________________________________________
    
def BDF(L,vmod,St,num_steps,ti,ht):    
    
    r = integrate.ode(EC_MOV).set_integrator('vode', method='bdf')  

    # initial velocity
    # initial position: equatorial plane 4Re from Earth
    u0 = L*St.Re
    u1 = 0*St.Re
    u2 = 0.0                                                                      #Para que el ángulo de paso sea ecuatorial z0=0
    u3 = 0.0
    u4 = vmod*mth.sin(St.pitch_angle*np.pi/180)
    u5 = vmod*mth.cos(St.pitch_angle*np.pi/180)
    
    r.set_initial_value([u0,u1,u2,u3,u4,u5], ti)

    t = np.zeros((num_steps, 1))
    T_Bounce=[]
    CA1 = np.zeros((num_steps, 1))
    CA2 = np.zeros((num_steps, 1))
    CA3 = np.zeros((num_steps, 1))
    CA4 = np.zeros((num_steps, 1))
    CA5 = np.zeros((num_steps, 1))
    CA6 = np.zeros((num_steps, 1))
    t[0] = ti
    CA1[0] = u0
    CA2[0] = u1
    CA3[0] = u2
    CA4[0] = u3
    CA5[0] = u4
    CA6[0] = u5
    fl3=1
    fl4=1
    fl1=0
    hr=0
    fl5=1
    
    
    fPl1=0
    fPl2=0
    fPl3=1

    k = 1
    while r.successful() and k < num_steps:
        r.integrate(r.t + ht)
         
    # Store the results to plot later
        t[k] = r.t
        CA1[k] = r.y[0]
        CA2[k] = r.y[1]
        CA3[k] = r.y[2]        

        if (r.y[2]<=0 and fl3):
            xh=r.y[0]
            yh=r.y[1]             
            fl3=0
            fl1=1
  
        if (r.y[2]>=0 and fl1 and fl4 and fl5): 
            if (r.y[1]>=0):
                fPl1=1
            if (r.y[1]<=0 and fPl1):
                fPl2=1     
            if (r.y[1]>=0 and fPl2):
                fPl1=0
                fPl2=0            
                T_Bounce.append(r.t)
            
        if (r.y[2]<=0 and fl1 and not(fl4) and fl5):
            if (r.y[1]>=0):
                fPl1=1
            if (r.y[1]<=0 and fPl1):
                fPl2=1     
            if (r.y[1]>=0 and fPl2):
                fPl1=0
                fPl2=0            
                T_Bounce.append(r.t)

        CA4[k] = r.y[3] 
        CA5[k] = r.y[4]         
        CA6[k] = r.y[5] 
        k += 1

    return CA1,CA2,CA3,CA4,CA5,CA6,t,T_Bounce
#______________________________________________________________________________    


#______________________________________________________________________________ 
 
St = V_Static()   
domL=np.linspace(2,6,20)
Tau_bounce=np.empty([len(domL),3])
T1=St.K1
T2=St.K2
T3=St.K3

Tau_bounce1=((np.pi*St.q*St.B0*St.Re**2)/(3*domL*T1))*((0.35 + 0.15*mth.sin(St.pitch_angle*np.pi/180))**(-1))
Tau_bounce2=((np.pi*St.q*St.B0*St.Re**2)/(3*domL*T2))*((0.35 + 0.15*mth.sin(St.pitch_angle*np.pi/180))**(-1))
Tau_bounce3=((np.pi*St.q*St.B0*St.Re**2)/(3*domL*T3))*((0.35 + 0.15*mth.sin(St.pitch_angle*np.pi/180))**(-1))

Tau_bounce[:,0]=Tau_bounce1
Tau_bounce[:,1]=Tau_bounce2
Tau_bounce[:,2]=Tau_bounce3

T=np.zeros(3)
K=np.zeros(3)
vmod=np.zeros(3)

T[0] = T1
T[1] = T2
T[2] = T3

K[0] = 10
K[1] = 20
K[2] = 30

vmod[0] = St.v_mod1
vmod[1] = St.v_mod2
vmod[2] = St.v_mod3


ti = 0.0
tf = 140
ht = 0.1   

# Number of time steps: 1 extra for initial condition
num_steps = int(np.floor((tf - ti)/ht)) + 1
Tau_Bounce= T_BOUNCE(domL,vmod,St)

iT=np.arange(3)

fig, ax = plt.subplots()

for i in iT :    
    ax.plot(domL,Tau_bounce[:,i], label= '%1.1f MeV' % (K[i]))
    ax.plot(domL,Tau_Bounce[:,i], marker='o')
    
#ax.semilogy()
    
#ax.scatter(domL,Tau_Bounce[:,i], color='b')     
#ax.scatter(domL,Tau_Bounce[:,i], color='g')     
#ax.scatter(domL,Tau_Bounce[:,i], color='r')     
#ax.scatter(domL,Tau_Bounce[:,i], color='Cyan')     
     
plt.rc('font', family='serif')
plt.ylabel(r'Periodo de deriva [s]', size='large',  family='monospace', weight= 'light')
plt.xlabel(r'Capa L',  size='large',  family='monospace', weight= 'light')
plt.legend()
#ax.set_ylim(1,100)
ax.set_xlim(2,6)

#import os

#outputname = 'Periodo_Deriva_teo2.txt'

#archivo=open(outputname,'w')

#for i in Tau_bounce:
#   archivo.write(str(i) + '\n')

#archivo.close()

