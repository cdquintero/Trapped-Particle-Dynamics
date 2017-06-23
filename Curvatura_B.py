# -*- coding: utf-8 -*-
"""
Created on Mon Dec 19 11:34:01 2016

@author: cd67
"""

#%matplotlib inline
from scipy import integrate

import math  as mth
import numpy as np                                                             #Paquete Arrays
import matplotlib.pyplot as plt 
import sympy as sym      
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
from copy import copy, deepcopy
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
import mpl_toolkits.mplot3d.axes3d
class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0, 0), (0, 0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        FancyArrowPatch.draw(self, renderer)

class V_Static(object) :                                                       #Se define una clase que contiene todas las variables estáticas del programa
    e = 1.602176565e-19                                                        #Carga de las partículas
    m_pr = 1.672621777e-27                                                     #Masa del protón [kg]
    m_el = 9.10938291e-31                                                      #Masa del electrón [kg]
    
    B0=3.07e-5                                                                 #Campo magnético ecuatorial
    q=e                                                                            
    m=3*m_pr                                                                  
    Re=6378137                                                                 #Radio de la tierra [m]
    n=6  
    ngc=4                                                                       #Número de Ecuaciones

    c = 299792458                                                              #Velocidad de la Luz
    K = 1e5                                                               #Energía cinética [eV]
    K = K*e;   

    pitch_angle = 40.0                                                         #Ángulo de paso inicial (grados)    
    v_mod = c*mth.sqrt(1-(1/(1+(K/(m*(c**2))))**2))                                          #Velocidad de la partícula     

    v_perp0 = v_mod*mth.sin(pitch_angle*np.pi/180)                             #Componente perpendicular de la velocidad
    v_par0 = v_mod*mth.cos(pitch_angle*np.pi/180)                              #Componente paralela de la velocidad
#_____________________________________________________________________________   

def EC_MOV(t,u) :    

    St = V_Static()    
    
    # calculate 3 Cartesian components of the magnetic field
    alpha = -St.B0*St.Re**3/((u[0]**2 + u[1]**2 + u[2]**2)**(2.5))
    Bx = 3*u[0]*u[2]*alpha
    By = 3*u[1]*u[2]*alpha
    Bz = (2*(u[2]**2) - u[0]**2 - u[1]**2)*alpha

    dudt=np.zeros((St.n,1))

    dudt[0] = u[3] 
    dudt[1] = u[4]
    dudt[2] = u[5]                                                         # dz/dt = w 
    dudt[3] = St.q/St.m*(u[4]*Bz - u[5]*By)
    dudt[4] = St.q/St.m*(u[5]*Bx - u[3]*Bz) 
    dudt[5] = St.q/St.m*(u[3]*By - u[4]*Bx) 
   
    return dudt
    
#______________________________________________________________________________

def EC_MOVGC(t,u) :    
    
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
    
    dudt=np.zeros((St.ngc,1))

    dudt[0] = beta*bxgB_x + u[3]*b_unit_x
    dudt[1] = beta*bxgB_y + u[3]*b_unit_y
    dudt[2] = beta*bxgB_z + u[3]*b_unit_z                                                            # dz/dt = w 
    dudt[3] = -mu/St.m*dotpr 
    
    return dudt   
    
#______________________________________________________________________________

def T_BOUNCE(domL):
    T_Bounce=np.zeros(len(domL))    
    Ti=np.arange(len(domL))     
     
    for ii in zip(Ti,domL):
        CA1,CA2,CA3,CA4,CA5,CA6,t,Tau_exp = BDF(ii[1],St,num_steps,ti,ht)
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
    By = 0
    Bz = (2*(z**2) - x**2 - y**2)*alpha
    B = mth.sqrt(Bx**2 + By**2 + Bz**2)

    return B
#______________________________________________________________________________
def B(x,y,z):
    alpha = -1/((x**2 + y**2 + z**2)**(5/2.0))
    Bx = (3*x*z*alpha)
    By = 0
    Bz = (2*(z**2) - x**2 - y**2)*alpha
    return Bx,By,Bz
#______________________________________________________________________________
    
def BDF(L,St,num_steps,ti,ht):    
    
    r = integrate.ode(EC_MOV).set_integrator('vode', method='bdf')  

    # initial velocity
    # initial position: equatorial plane 4Re from Earth
    u0 = L*St.Re
    u1 = 0*St.Re
    u2 = 0.0                                                                      #Para que el ángulo de paso sea ecuatorial z0=0
    u3 = 0.0
    u4 = St.v_mod*mth.sin(St.pitch_angle*np.pi/180)
    u5 = St.v_mod*mth.cos(St.pitch_angle*np.pi/180)
    
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
    fl1=0
    fl2=0
    fl3=1

    k = 1
    while r.successful() and k < num_steps:
        r.integrate(r.t + ht)
         
    # Store the results to plot later
        t[k] = r.t
        CA1[k] = r.y[0]
        CA2[k] = r.y[1]
        CA3[k] = r.y[2]        

        if (r.y[2]>=0):
            fl1=1
        if (r.y[2]<=0 and fl1):
            fl2=1     
        if (r.y[2]>=0 and fl2):
            fl1=0
            fl2=0            

            T_Bounce.append(r.t)

        CA4[k] = r.y[3] 
        CA5[k] = r.y[4]         
        CA6[k] = r.y[5] 
        k += 1

    return CA1,CA2,CA3,CA4,CA5,CA6,t
#______________________________________________________________________________    

def BDFGC(L,St,num_steps,ti,ht):
    
    r = integrate.ode(EC_MOVGC).set_integrator('vode', method='bdf')
   

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

    t[0] = ti
    CA1[0] = u0
    CA2[0] = u1
    CA3[0] = u2
    CA4[0] = u3

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
        
    return CA1,CA2,CA3,CA4,t


St = V_Static()   

ti = 0.0
tf = 50
ht = 0.01 
L=4

u0 = L*St.Re
u1 = 0 
u2 = 0*St.Re
u3 = St.v_par0

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
#______________________________________________________________________________

Range=100
RE=20
r_start= 1                                                                     
r_end= 10
q_start= 0
q_end= np.pi/2
p_start= 0
p_end= 2*np.pi
l = np.linspace(q_start, q_end, num=Range)  

pitch_angle = 50.0                                                             #  initial angle between velocity and mag.field (degrees)
pitch_angle = pitch_angle*np.pi/180

Distancia = []
Latitud = []

ii=np.arange(Range)
ri=np.arange(2,7)
BE = 3.11*10E-5 
Bb= []

flag=1

while flag:
    for li in zip(l, ii):
        Latitud.append(li[0]) 
        Alt_rad = li[0]
        L=5
        Distancia.append(L*((mth.cos(Alt_rad))**2))    
        B0=BE/(L**3)
        Bm=B0/((mth.sin(pitch_angle))**2)
        Bmp=Bm                                                                                                                  
        B_Dipolo = B0*mth.sqrt(3.0*(mth.sin(Alt_rad)**2)+1.0)/(mth.cos(Alt_rad)**6)
        Bb.append(B_Dipolo)  
        if B_Dipolo >= Bm:
            flag=0
            break

Dist=deepcopy(Distancia)
jj=np.arange(len(Distancia)-1,0-1,-1)

Dist.reverse()
Dist.pop(-1)
Distancia=Dist

Lat=deepcopy(Latitud)
jj=np.arange(len(Latitud)-1,0-1,-1)

for i in jj:
    Lat[i]=Lat[i]

Lat.reverse()
Lat.pop(-1)
Latitud=Lat

#______________________________________________________________________________   
 
 
xx = np.linspace(-4,4,9)
yy = np.linspace(-4,4,9)
zz = np.linspace(-4,4,9)

xx,yy,zz = np.meshgrid(xx,yy,zz)

# Plot of the fields  
Bx,By,Bz = B(xx,yy,zz)     
 
St = V_Static()   
domL=np.linspace(2,8,20)
T=1/2.0*St.m*St.v_mod**2
Tau_bounce=(domL*St.Re)/(mth.sqrt(T/St.m))*(3.7-1.6*mth.sin(St.pitch_angle*np.pi/180))
Tau=(np.pi*St.e*St.B0*St.Re**2)/((3*4.0*(1/2.0*St.m_pr*St.v_mod**2))*(0.35+0.15*mth.sin(St.pitch_angle*np.pi/180)))

ti = 0.0
tf = 7
ht = 0.01   
# Number of time steps: 1 extra for initial condition
num_steps = int(np.floor((tf - ti)/ht)) + 1

L=5.0
CA1,CA2,CA3,CA4,CA5,CA6,t = BDF(L,St,num_steps,ti,ht)
CA1cg,CA2cg,CA3cg,CA4cg,t= BDFGC(L,St,num_steps,ti,ht)


#Tau_Bounce= T_BOUNCE(domL)

#fig3, ax3 = plt.subplots()
#ax3.set_ylim(1,100)
#ax3.plot(domL,Tau_bounce)
#ax3.plot(domL,Tau_Bounce,marker='o')
#ax3.semilogy()

u0 = L*St.Re
u1 = 0*St.Re
u2 = 0.0                                                                      #Para que el ángulo de paso sea ecuatorial z0=0
u3 = 0.0
u4 = St.v_mod*mth.sin(St.pitch_angle*np.pi/180)
u5 = St.v_mod*mth.cos(St.pitch_angle*np.pi/180)
               
xsol = CA1.reshape(num_steps)
ysol = CA2.reshape(num_steps)
zsol = CA3.reshape(num_steps)   
vxsol = CA4.reshape(num_steps)     
vysol = CA5.reshape(num_steps)
vzsol = CA6.reshape(num_steps)

xsolcg = CA1cg.reshape(num_steps)
ysolcg = CA2cg.reshape(num_steps)
zsolcg = CA3cg.reshape(num_steps)   
vsolcg = CA4cg.reshape(num_steps)     

flag=1

ri=-3
rf=3

if flag:
    
    fig = plt.figure()
    ax = fig.gca(projection='3d')      

#    ax.set_xlim3d(2,5.5)
    ax.set_ylim3d(-.4,.2)
 #  ax.set_zlim3d(ri,rf)  
    ax.plot(xsol/St.Re, ysol/St.Re, zsol/St.Re, color='g', alpha=0.8)
    ax.plot(xsolcg/St.Re, ysolcg/St.Re, zsolcg/St.Re, color='brown', alpha=0.8)


    ax.scatter(u0/St.Re, u1/St.Re, u2/St.Re, s = 20, c='black',color='black', alpha=1)   
 
    ax.set_xlabel('Eje X')
    ax.set_ylabel('Eje Y')
    ax.set_zlabel('Eje Z')  

Longitud = np.linspace(p_start, p_end, len(Distancia))

Phi, Theta = np.mgrid[0:2*np.pi:200j, min(Latitud):max(Latitud):75j]

Long = 0

Rc = Arrow3D([4, 5], [0, 0], [0, 0], mutation_scale=20, lw=1, arrowstyle="->", color="black")
VR = Arrow3D([5, 5], [0, -0.2], [0, 0], mutation_scale=20, lw=1, arrowstyle="->", color="black")

ax.add_artist(Rc)
ax.add_artist(VR)
       
ax.text(5.2, -0.17, 0, s = r'$V_{R}$', fontsize = 15)
ax.text(4.6, 0.05, 0, s = r'$Rc$', fontsize = 15)
ax.text(4.2, 0.05, 1.4, s = r'$B$', fontsize = 15)
   
xx = Distancia*np.sin(Latitud)*np.cos(Long)
yy = Distancia*np.sin(Latitud)*np.sin(Long)  
zz = Distancia*np.cos(Latitud)   
ax.plot(zz, yy, xx, color='cornflowerblue')



