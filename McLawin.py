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

pitch_angle = 30.0                                                             #  initial angle between velocity and mag.field (degrees)
pitch_angle = pitch_angle*np.pi/180

Range=100
RE=20
r_start= 1                                                                     
r_end= 10
q_start= 0
q_end= 2*np.pi

Distancia = np.empty([Range,6])
Latitud =  np.linspace(q_start, q_end, num=Range)  
ii=np.arange(Range)
ri=np.arange(2,7)
BE = 3.11*10E-5 
B=np.empty([Range,6])

for jj in ri:
    for li in zip(Latitud, ii):
        Alt_rad = li[0]
        L=jj
        Distancia[li[1],jj-2]=2    
        B_Dipolo = BE*mth.sqrt(3.0*(mth.sin(Alt_rad)**2)+1.0)/((L**3)*(mth.cos(Alt_rad)**6))
        B[li[1],jj-2] = B_Dipolo  
 
fig, ax = plt.subplots(figsize=(10, 10), subplot_kw=dict(projection='polar'))  
       
Rs, Qs =  np.ogrid[r_start:r_end:100j, q_start:q_end:100j]                                     #Malla
ax.plot(Latitud, Distancia, color='cornflowerblue')
circle = plt.Circle((0, 0), 1, transform=ax.transData._b, color='cornflowerblue', linewidth=0)   
ax.add_artist(circle)   
for jj in ri:
    ax.text(0, jj-1, s = u'$L=%d$' %jj, fontsize = 15)
plt.show()  


