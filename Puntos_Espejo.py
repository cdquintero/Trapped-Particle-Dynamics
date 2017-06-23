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
from copy import copy, deepcopy

#______________________________________________________________________________

def TRAYECTORIA(L,pitch_angle): 
    Distancia2 = []
    Latitud2 = []

    flag=1

    while flag:
        for li in zip(l, ii):
            Latitud2.append(li[0]) 
            Alt_rad = li[0]
            Distancia2.append(L*((mth.cos(Alt_rad))**2))    
            B0=BE/(L**3)
            Bm=B0/((mth.sin(pitch_angle))**2)                                                                                                                
            B_Dipolo = B0*mth.sqrt(3.0*(mth.sin(Alt_rad)**2)+1.0)/(mth.cos(Alt_rad)**6)  
            if B_Dipolo >= Bm:
                flag=0
                break

    Dist=deepcopy(Distancia2)
    jj=np.arange(len(Distancia2)-1,0-1,-1)

    Dist.reverse()
    Dist.pop(-1)
    Distancia2=Dist+Distancia2

    Lat=deepcopy(Latitud2)
    jj=np.arange(len(Latitud2)-1,0-1,-1)

    for i in jj:
        Lat[i]=-Lat[i]

    Lat.reverse()
    Lat.pop(-1)
    Latitud2=Lat+Latitud2

    return Latitud2,Distancia2 
#______________________________________________________________________________

def ESPEJO(pitch_angle):
    lm=[]
    rm=[]

    for jj in ri:
        flag=1
        L=jj
        B0=BE/(L**3)
        Bm=B0/((mth.sin(pitch_angle))**2)
    
        while flag:
            for li in zip(l, ii):
                Alt_rad = li[0]
                Dist=L*((mth.cos(Alt_rad))**2)
                B_Dipolo = BE*mth.sqrt(3.0*(mth.sin(Alt_rad)**2)+1.0)/((L**3)*(mth.cos(Alt_rad)**6))                                                                                                                           
                if B_Dipolo >= Bm:
                    lm.append(Alt_rad)               
                    rm.append(Dist)         
                    flag=0
                    break

    return lm,rm
#______________________________________________________________________________

Range=100
RE=20
r_start= 1                                                                     
r_end= 10
q_start= 0
q_end= 2*np.pi
l = np.linspace(q_start, q_end, num=Range)  

Distancia = np.empty([Range,6])
Latitud =  np.linspace(q_start, q_end, num=Range)  
ii=np.arange(Range)
ri=np.arange(2,7)
BE = 3.11*10E-5 
B=np.empty([Range,6])

pitch_angle = 60.0                                                             #  initial angle between velocity and mag.field (degrees)
pitch_angle = pitch_angle*np.pi/180

fig, ax = plt.subplots(figsize=(10, 10), subplot_kw=dict(projection='polar'))  
       
Rs, Qs =  np.ogrid[r_start:r_end:100j, q_start:q_end:100j]                                     #Malla
ax.plot(Latitud, Distancia, color='cornflowerblue')
circle = plt.Circle((0, 0), 1, transform=ax.transData._b, color='cornflowerblue', linewidth=0)   
ax.add_artist(circle)   

for jj in ri:
    ax.text(0, jj-0.8, s = u'$L=%d$' %jj, fontsize = 10)

for jj in ri:
    ax.text(np.pi, jj-0.1, s = u'$L=%d$' %jj, fontsize = 10)

lm,rm= ESPEJO(pitch_angle)
lm=np.array(lm) 

pitch_angle2 = 30.0                                                             #  initial angle between velocity and mag.field (degrees)
pitch_angle2 = pitch_angle2*np.pi/180

lm2,rm2= ESPEJO(pitch_angle2)
lm2=np.array(lm2) 

   
ax.plot(lm, rm, color='r')
ax.plot(-lm, rm, color='r')
ax.plot(lm2+np.pi, rm2, color='r')
ax.plot(-lm2+np.pi, rm2, color='r')

for i in ri:
    Latitud2, Distancia2=TRAYECTORIA(i+0.15,pitch_angle)
    ax.plot(Latitud2, Distancia2, color='black')
    ax.scatter(min(Latitud2), min(Distancia2), color='black',s=8)
    ax.scatter(max(Latitud2), min(Distancia2), color='black',s=8)

ax.text(max(Latitud2), min(Distancia2)+0.1, s = r'$\alpha_0=%10.0f ^\circ$' %(pitch_angle*180/np.pi), fontsize = 12)

for i in ri:   
    Latitud2, Distancia2=np.array(TRAYECTORIA(i+0.15,pitch_angle2))
    Latitud2= Latitud2+np.pi 
    ax.plot(Latitud2, Distancia2, color='black')
    ax.scatter(min(Latitud2), min(Distancia2), color='black',s=8)
    ax.scatter(max(Latitud2), min(Distancia2), color='black',s=8)

ax.text(min(Latitud2), min(Distancia2)+0.1, s = r'$\alpha_0=%10.0f ^\circ$' %(pitch_angle2*180/np.pi), fontsize = 12)


plt.show()  


