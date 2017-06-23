# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 20:54:29 2016

@author: cd67
"""
#%matplotlib inline
import numpy as np                                                             #Paquete Arrays
import matplotlib.pyplot as plt                                                #Paquete Gráficas
import sympy as sym                                                            #Paquete Simbólicas
from sympy.abc import r, q

def PotencialB(BE = 3.11**(-5), RE = 6378.137):                                #Función Potencial Magnético
    return  BE*(RE**3)*sym.sin(q)/r**2                                             #BE(Campo magnético ecuatorial en la superficie de la tierra)
                                                                                   #RE(Distancia normalizada a radios Terrestres)
def CampoB(psi):                                                               #Función Campo Magnético
    u = sym.lambdify((r, q),  psi.diff(r), 'numpy')                                #Componente radial
    v = sym.lambdify((r, q),  (1/r)*psi.diff(q), 'numpy')                          #Componente polar
    return u, v
    
def plot_streamlines(ax, u, v, rlim, qlim):                                    #Función lineas de Campo vectorial
    r0, r1 = rlim                                                                  #Fronteras r                    
    q0, q1 = qlim                                                                  #Fronteras theta
    R, Q =  np.ogrid[r0:r1:100j, q0:q1:100j]                                       #Malla
    ax.streamplot(Q, R, v(R, Q), u(R, Q))
    
def format_axes(ax, r_end):                                                    #Función formato de los ejes
    ax.set_aspect('equal')                         
    ax.figure.subplots_adjust(bottom=0, top=1, left=0, right=1)
    ax.xaxis.set_ticks(np.arange(0, 2*np.pi, np.pi/4))                             #Marcas theta
    ax.yaxis.set_ticks(np.arange(0, r_end, 1))                                   #Marcas r
    for spine in ax.spines.itervalues():
        spine.set_visible(False)                                                
   
r_start= 0.1
r_end= 6.6
q_start= 0
q_end= 2*np.pi
rlim = (r_start, r_end)                                                        #Fronteras r
qlim = (q_start, q_end)                                                        #Fronteras theta
    
PotB = PotencialB()                                                            
u, v = CampoB(PotB)
fig, ax = plt.subplots(figsize=(10, 10), subplot_kw=dict(projection='polar'))  
plot_streamlines(ax, u, v, rlim, qlim)
circle = plt.Circle((0, 0), 1, transform=ax.transData._b, color='cornflowerblue', linewidth=0)   
ax.add_artist(circle)  
format_axes(ax, r_end)     

                            
