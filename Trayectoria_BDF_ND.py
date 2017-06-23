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
    K = 1e5                                                               #Energía cinética [eV]
    K = K*e;
    K=K/(m*c**2)

    pitch_angle = 30.0                                                         #Ángulo de paso inicial (grados)    
    v_mod = mth.sqrt(1-(1/(1+K))**2)                                          #Velocidad de la partícula     
    v_par0 = v_mod*mth.cos(pitch_angle*np.pi/180)                              #Componente paralela de la velocidad
#______________________________________________________________________________   

def EC_MOV(t,u) :    

    St = V_Static()    
    
    # calculate 3 Cartesian components of the magnetic field
    alpha = -((St.q*St.Re*St.B0)/(St.m*St.c))/((u[0]**2 + u[1]**2 + u[2]**2)**(5/2.0))
    Bx = 3*u[0]*u[2]*alpha
    By = 3*u[1]*u[2]*alpha
    Bz = (2*(u[2]**2) - u[0]**2 - u[1]**2)*alpha
    B = mth.sqrt(Bx**2 + By**2 + Bz**2)
    mu = (St.v_mod**2-u[3]**2)/(2.0*B)  
    Tau_c = 1/B   

    dudt=np.zeros((St.n,1))

    dudt[0] = u[3] 
    dudt[1] = u[4]
    dudt[2] = u[5]                                                         # dz/dt = w 
    dudt[3] = (u[4]*Bz - u[5]*By)
    dudt[4] = (u[5]*Bx - u[3]*Bz) 
    dudt[5] = (u[3]*By - u[4]*Bx)
   
    return dudt
    
#______________________________________________________________________________
def B(x,y,z,St):
    alpha = -((St.q*St.Re*St.B0)/(St.m*St.c))/((x**2 + y**2 + z**2)**(5/2.0))
    Bx = (3*x*z*alpha)
    By = 3*y*z*alpha
    Bz = (2*(z**2) - x**2 - y**2)*alpha
    return Bx,By,Bz
#______________________________________________________________________________
    
def BDF(L,St,num_steps,ti,ht):    
    
    r = integrate.ode(EC_MOV).set_integrator('vode', method='bdf') 
    Tau=np.zeros(num_steps) 
    Mu=np.zeros(num_steps) 

    # initial velocity
    # initial position: equatorial plane 4Re from Earth
    u0 = L
    u1 = 0.0
    u2 = 0.0                                                                      #Para que el ángulo de paso sea ecuatorial z0=0
    u3 = 0.0
    u5 = St.v_par0      
    u4 = mth.sqrt(St.v_mod**2-u5**2)
    
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

    k = 1
    while r.successful() and k < num_steps:
        r.integrate(r.t + ht)
         
    # Store the results to plot later
        t[k] = r.t
        CA1[k] = r.y[0]
        CA2[k] = r.y[1]
        CA3[k] = r.y[2]        
        CA4[k] = r.y[3] 
        CA5[k] = r.y[4]         
        CA6[k] = r.y[5]
        
        k += 1

    return CA1,CA2,CA3,CA4,CA5,CA6,t
#______________________________________________________________________________    

def add_arrow_to_line2D(
    axes, line, arrow_locs=[0.2, 0.4, 0.6, 0.8],
    arrowstyle='-|>', arrowsize=1.0, transform=None):
    """
    Add arrows to a matplotlib.lines.Line2D at selected locations.

    Parameters:
    -----------
    axes: 
    line: Line2D object as returned by plot command
    arrow_locs: list of locations where to insert arrows, % of total length
    arrowstyle: style of the arrow
    arrowsize: size of the arrow
    transform: a matplotlib transform instance, default to data coordinates

    Returns:
    --------
    arrows: list of arrows
    """
    if not isinstance(line, mlines.Line2D):
        raise ValueError("expected a matplotlib.lines.Line2D object")
    x, y = line.get_xdata(), line.get_ydata()

    arrow_kw = {
        "arrowstyle": arrowstyle,
        "mutation_scale": 10 * arrowsize,
    }

    color = line.get_color()
    use_multicolor_lines = isinstance(color, np.ndarray)
    if use_multicolor_lines:
        raise NotImplementedError("multicolor lines not supported")
    else:
        arrow_kw['color'] = color

    linewidth = line.get_linewidth()
    if isinstance(linewidth, np.ndarray):
        raise NotImplementedError("multiwidth lines not supported")
    else:
        arrow_kw['linewidth'] = linewidth

    if transform is None:
        transform = axes.transData

    arrows = []
    for loc in arrow_locs:
        s = np.cumsum(np.sqrt(np.diff(x) ** 2 + np.diff(y) ** 2))
        n = np.searchsorted(s, s[-1] * loc)
        arrow_tail = (x[n], y[n])
        arrow_head = (np.mean(x[n:n + 2]), np.mean(y[n:n + 2]))
        p = mpatches.FancyArrowPatch(
            arrow_tail, arrow_head, transform=transform,
            **arrow_kw)
        axes.add_patch(p)
        arrows.append(p)
    return arrows

#______________________________________________________________________________ 
 

St= V_Static()  
 
xx = np.linspace(-4,4,9)
yy = np.linspace(-4,4,9)
zz = np.linspace(-4,4,9)

xx,yy,zz = np.meshgrid(xx,yy,zz)

# Plot of the fields  
Bx,By,Bz = B(xx,yy,zz,St)     
 
ti = 0.0
tf = (St.c/St.Re)*2000
ht = (St.c/St.Re)*0.01   
# Number of time steps: 1 extra for initial condition
num_steps = int(np.floor((tf - ti)/ht)) + 1

L=4.0
CA1,CA2,CA3,CA4,CA5,CA6,t = BDF(L,St,num_steps,ti,ht)

#Tau_Bounce= T_BOUNCE(domL)

#fig3, ax3 = plt.subplots()
#ax3.set_ylim(1,100)
#ax3.plot(domL,Tau_bounce)
#ax3.plot(domL,Tau_Bounce,marker='o')
#ax3.semilogy()

u0 = L
u1 = 0.0
u2 = 0.0                                                                      #Para que el ángulo de paso sea ecuatorial z0=0
u3 = 0.0
u5 = St.v_par0
u4 = mth.sqrt(St.v_mod**2-u5**2)
               
w = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)

x = np.outer(np.cos(w), np.sin(v))
y = np.outer(np.sin(w), np.sin(v))
z = np.outer(np.ones(np.size(w)), np.cos(v))

xsol = CA1.reshape(num_steps)
ysol = CA2.reshape(num_steps)
zsol = CA3.reshape(num_steps)   
vxsol = CA4.reshape(num_steps)     
vysol = CA5.reshape(num_steps)
vzsol = CA6.reshape(num_steps)

flag=1

ri=-5
rf=5

if flag:
    
    fig = plt.figure()
    ax = fig.gca(projection='3d')      

    ax.set_xlim3d(ri,rf)
    ax.set_ylim3d(ri,rf)
    ax.set_zlim3d(ri,rf)  
    ax.plot_surface(x, y, z, rstride=6, cstride=6, linewidth=0.5, color='cornflowerblue') 
    ax.set_aspect('equal')
    ax.plot(xsol, ysol, zsol, color='b', alpha=0.8)

    ax.quiver(xx,yy,zz,Bx,By,Bz,color='#4682b4', length=0.5,arrow_length_ratio=0.4,linestyles='dotted',linewidths=0.5,alpha=0.8)        #Plot the magnetic field

    ax.scatter(u0, u1, u2, s = 20, c='black',color='black', alpha=1)   
    ax.text(u0-0.3, u1+0.1, u2, s = u'+', fontsize = 15)
 
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')  

vpr=np.zeros(num_steps )
k=np.arange(num_steps )

if (not flag):

    L=3.0
    CA11,CA21,CA31,CA41,CA51,CA61,t1,Tau_exp1 = BDF(L,St,num_steps,ti,ht)

    L=2.0
    CA12,CA22,CA32,CA42,CA52,CA62,t2,Tau_exp2 = BDF(L,St,num_steps,ti,ht)

    x1sol = CA11.reshape(num_steps)
    y1sol = CA21.reshape(num_steps)
    z1sol = CA31.reshape(num_steps)   
    v1xsol = CA41.reshape(num_steps)     
    v1ysol = CA51.reshape(num_steps)
    v1zsol = CA61.reshape(num_steps)

    x2sol = CA12.reshape(num_steps)
    y2sol = CA22.reshape(num_steps)
    z2sol = CA32.reshape(num_steps)   
    v2xsol = CA42.reshape(num_steps)     
    v2ysol = CA52.reshape(num_steps)
    v2zsol = CA62.reshape(num_steps)

    rl=(St.m*St.v_mod)/(St.q*St.B0)    
    
    for indice in zip(vxsol,vysol,vzsol,k):
        vpr[indice[3]]=mth.sqrt(indice[0]**2+indice[1]**2+indice[2]**2)

    fig, ax = plt.subplots()
    ax.plot(xsol/St.Re, ysol/St.Re, color='brown')
#    line, = ax.plot(xsol/St.Re, ysol/St.Re, 'k-')
#    add_arrow_to_line2D(ax, line, arrow_locs=np.arange(0., 10, 0.565*np.pi*rl/St.Re),arrowstyle='-|>')
    ax.plot(x1sol/St.Re, y1sol/St.Re, color='r')
    ax.plot(x2sol/St.Re, y2sol/St.Re, color='tomato')
    ax.set_xlim(ri, rf)
    ax.set_ylim(ri, rf)
    ax.set_aspect('equal')
    circle = plt.Circle((0, 0), 1, transform=ax.transData._b, color='cornflowerblue', linewidth=0)   
    circleout = plt.Circle((0, 0), 1+0.05, transform=ax.transData._b, color='#4682b4', linewidth=0)   
    ax.add_artist(circleout)   
    ax.add_artist(circle)
#    ax.scatter(u0/St.Re, u1/St.Re, s = 60, c='r',color='black')   
#    ax.text(u0/St.Re-0.5, u1/St.Re, s = u'+', fontsize = 15)
  
    n=6
    a=4
    b=np.arange(np.pi/(2*n), np.pi/2,np.pi/(2*n))
    c=a*np.cos(b)    
    d=np.sqrt(16-c**2)     
#    e=n+1    
#    while e<len(d):    
#        d[e]=-1*d[e]
#        e=e+1

    ax.scatter(c, d, s = 15, c='b',color='b',)     
    ax.scatter(4, -4, marker='x', s = 60, c='b',color='b',)   
    ax.scatter(4, 4, s = 60, c='b',color='b',)   
    ax.text(4, 4-0.5,s = r'$\nabla B$', fontsize = 12)
    ax.text(4, -4+0.3, s = r'$\nabla B$', fontsize = 12)
    ax.text(-0.2, 3.5, s = u'$B$', fontsize = 12)

    fig2, ax2 = plt.subplots(3)   
    ax2[0].plot(t, zsol/St.Re, color='brown', marker='o')
    ax2[1].plot(t1, z1sol/St.Re, color='r')
    ax2[2].plot(t2, z2sol/St.Re, color='tomato')
    ax.text(2-0.5, 0+0.5, s = u'$L=2$', fontsize = 10)
    ax.text(3-0.5, 0+0.5, s = u'$L=3$', fontsize = 10)
    ax.text(4-0.5, 0+0.5, s = u'$L=4$', fontsize = 10)
    fig2.text(0.4,0.01,'Tiempo',ha='center')
    fig2.text(0.04,0.5,'Eje Z',va='center',rotation='vertical')    
    plt.show()  
     
