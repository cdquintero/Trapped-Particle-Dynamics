# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 15:28:02 2016
@author: cd67

Archivos externos

AE8MIN.txt        Coeficientes del modelo de electrones  (una línea más corta que AE8MIN.ASC)
AP8MIN.txt        Coeficientes del modelo de protones  (una línea más corta que AP8MIN.ASC)


"""

#%matplotlib inline
import math  as mth
import numpy as np                                                             #Paquete Arrays
import pylab as pl                                                             
import matplotlib.pyplot as plt                                                #Paquete Gráficas
import sympy as sym                                                            #Paquete Simbólicas
from sympy.abc import r, q
from matplotlib import ticker                                                    
#_________________________________________________________________________________________________________

class V_Static(object) :                                                       #Se define una clase que contiene todas las variables estáticas del programa
    FI_Step = 0.0
    I1=0
    X_W1=0.0; X_W2=0.0; X_INCR2=0.0; X_INCR1=0.0
    X_PCM=0.0; Z_PCM=0.0; SL2=0.0; FNB=0.0; DFL=0.0
    J1=0; J2=0; I_Time=0; L1=0; L2=0
    I2=0; FLOG1=0; FLOG2=0
    XB_W1=0.0; XB_W1=0.0; SL1=0.0
    X_PCP=0.0; FLOG=0.0
    FLAG=0


def TRARA1 (FL, BB0, E, Descr, Map) :                                          #TRARA1 encuentra el logaritmo del flujo de partículas para una energía E, 
                                                                               #intensidad de campo normalizado BB0 y un valor de L dados
    St = V_Static()
    S0 = 0
    I0 = F0 = E0 = 0.0                                                         #Se instancia la clase
    F1 = 1.001; F2 = 1.002 
    St.FI_Step = float(Descr[6])/Descr[1]                                      #Se rescatan el los factores de escala de E,L,B y F del vector de cabecera
    E_Scale = float(Descr[3])
    F_Scale = float(Descr[6])
    XNL = min(15.6,abs(FL))                                                    #Se establece un límite superior para L
    NL = XNL*Descr[4]                                                          #Valor L con escala                 
    if (BB0 < 1.0 ) :                                                          #Se establece un límite inferior par BB0
        BB0 = 1.0                                                             
    NB = ( BB0-1.0)*Descr[5]                                                   #BB0 con escala  
  
    St.I1 = 0                                                         
    I2 = Map[0]                                                                #Número de elementos en el primer submapa de energías                                          
    I3 = I2+Map[I2]                                                            #Índice del último elemento del submapa 2  
    L3 = Map[I3]                                                               #Número de elementos en el tercer submapa de energías                                                  
    E1 = Map[St.I1+1]/E_Scale                                                  #Energía del primer submapa de energías (Sin escala)  
    E2 = Map[I2+1]/E_Scale                                                     #Energía del segundo submapa de energías (Sin escala) 

    S1 = 1                                                                     #Variables lógicas que indican si un flujo f(E,B,L) ya ha sido calculado 
    S2 = 1                                                                     #antes de interpolar (TRARA2). Si no Si=True
               
                                                                           
    while (not((E <= E2) or (L3==0))) :                                        #Para cada valor de energía E(I) recorre los sucesivos submapas de energía 
                                                                               #para encontrar los valores E0, E1 y E2 que obedecen  E0 < E1 < E(I) < E2
        I0 = St.I1                                                             
        St.I1 = I2                                                            
        I2 = I3                                                               
        I3 = I3+L3                                                             
        L3 = Map[I3]                                                                                      
                                                    
        E0 = E1                                                                
        E1 = E2                                                               
        E2 = Map[I2+1]/E_Scale                                                
     
        S0 = S1                                                            
        S1 = S2                                                            
        S2 = 1                                                        
      
        F0 = F1                                                             
        F1 = F2

    if S1 :                                                                    #Llama la función TRARA2 para interpolar los mapas de flujo para E1,E2 en el espacio L-B/B0
        F1 = TRARA2((St.I1+2), NL, NB, Map,St)/F_Scale                         #para los flujos F1, F2

    if S2 :
        F2 = TRARA2((I2+2), NL, NB, Map,St)/F_Scale
    
    S1 = 0                                                                     #Como ya se calcularon  F1 y F2: S1 = S2 = False 
    S2 = 0                                                

    F=F1+(F2-F1)/(E2-E1)*(E-E1)                                                #Finalmente se interpola linealmente 
       
    if (F2 <= 0.0 or St.I1 == 0 ) :
        if (S0):
            F0=TRARA2((I0+3-1), NL, NB, Map, St)/F_Scale
        S0=0;                                       
        F = min(F , F0+(F1-F0)*(E-E0)/(E1-E0))                   
    
    F = max(F, 0.0)                                                            #El flujo logarítmico debe ser siempre mayor o igual a 0   
           
    return F

#_________________________________________________________________________________________________________

def TRARA2 (I, IL, IB, Sub_Map, St) :                                          #* Parametros arbitrarios
    
    St.I2 = St.Z_W1 = St.Z_W2 = St.J1 = St.J2 = St.I_Time = St.L1 = St.L2 = 0
    St.X_W1 = St.X_W2 = St.X_INCR2 = St.X_INCR1 = St.X_PCM = St.Z_PCM = St.SL2 = St.FNB = St.DFL = 0

    St.J2 = 4
    St.I2 = I
 
    St.FNB = IB
    St.FNL = IL 
    St.L2 = Sub_Map[St.I2]
        
        
    while Sub_Map[St.I2+1] <= IL: 
        St.I1 = St.I2                                                            
        St.L1 = St.L2
        St.I2 = St.I2+St.L2
        St.L2 = Sub_Map[St.I2]
                 
    if St.L1 < 4  and St.L2 < 4:
        return(TRARA5(St))
   
    if Sub_Map[St.I2+2] <= Sub_Map[St.I1+2] :
        reversal(St)
         
    while St.I_Time <= 1 :
        FLL1 = float(Sub_Map[St.I1+1])
        FLL2 = float(Sub_Map[St.I2+1])
        St.DFL = (St.FNL-FLL1)/(FLL2-FLL1)
        St.Z_W1 = Sub_Map[St.I1+2]
        St.Z_W2 = Sub_Map[St.I2+2]
        St.X_W1 = 0.0
        St.X_W2 = 0.0 

        if St.L1 < 4 : 
           return(Case2(I, Sub_Map, St))
           
        while St.J2 <= St.L2 : #for                
            St.X_INCR2 = Sub_Map[St.I2+St.J2-1]
                 
            if ((St.X_W2+St.X_INCR2) > St.FNB):
                return (TRARA3(I, Sub_Map, St))
            St.X_W2 = St.X_W2+St.X_INCR2
            St.Z_W2 = St.Z_W2-St.FI_Step
            St.J2 = St.J2+1
             
        St.I_Time = St.I_Time+1
        reversal(St)
     
    return 0.0

#_________________________________________________________________________________________________________
    
def TRARA3 (I, Sub_Map, St) :                                                #* Parametros arbitrarios
    
    St.XB_W1 = St.XB_W1 = St.SL1 = St.FLAG= 0
    St.J1 = 4

    #Caso (Goto 30): Caso en que las curvas se han invertido y el punto AT 
    #se ha encontrado, se trata como el caso 3.
    if St.I_Time == 1 :
        St.X_W2=0.0 
        St.J2=4 
        St.X_INCR2 = Sub_Map[St.I2+St.J2-1]
        St.Z_W2 = Sub_Map[St.I2+2]
        St.Z_W1 = Sub_Map[St.I1+2]
        return(Case3(I, Sub_Map, St))
    
    if St.J2 == 4 : 
        return(Case3(I, Sub_Map, St))


    return(Case4(I, Sub_Map, St))
    #Tercer caso (Goto 28): En el caso de que XAT se encuentre en el primer paso (FNX está cerca del ecuador).   
   

#_________________________________________________________________________________________________________

def TRARA4 (I, Sub_Map, St) :                                               #* Parametros arbitrarios

    St.X_PCP = St.Z_PCP = 0.0
                                          
    while 1 :
        if  St.SL1 < St.SL2 : 
            St.XB_W1 = ((St.Z_W1/St.FI_Step)*St.X_INCR1+St.X_W1)/((St.X_INCR1/St.FI_Step)*St.SL2+1.0)
            St.X_PCP = St.XB_W1+(St.X_W2-St.XB_W1)*St.DFL
            St.Z_PCP = St.X_PCP*St.SL2
            if(St.X_PCP >= St.FNB) :
                return(TRARA5(St))
            St.X_PCM = St.X_PCP
            St.Z_PCM = St.Z_PCP
            if(St.J2 >= St.L2) :
                return 0
            St.J2 = St.J2+1;
            St.X_INCR2 = Sub_Map[St.I2+St.J2-1]
            St.Z_W2 = St.Z_W2-St.FI_Step
            St.X_W2 = St.X_W2+St.X_INCR2
            St.SL2 = St.Z_W2/St.X_W2
   
        else :
            St.XB_W1 = ((St.Z_W2/St.FI_Step)*St.X_INCR2+St.X_W2)/((St.X_INCR2/St.FI_Step)*St.SL1+1.0)
            St.X_PCP = St.X_W1+(St.XB_W1-St.X_W1)*St.DFL
            St.Z_PCP = St.X_PCP*St.SL1
            if St.X_PCP >= St.FNB :
               return(TRARA5(St))
            St.X_PCM = St.X_PCP
            St.Z_PCM = St.Z_PCP
            if(St.J1 >= St.L1) :
               return 0
            St.J1 = St.J1+1
            St.X_INCR1 = Sub_Map[St.I1+St.J1-1]
            St.Z_W1 = St.Z_W1-St.FI_Step
            St.X_W1 = St.X_W1+St.X_INCR1
            St.SL1 = St.Z_W1/St.X_W1    

#_________________________________________________________________________________________________________
    
def TRARA5 (St) :                                                                #* Parametros arbitrarios
    if St.X_PCP < (St.X_PCM+1.0e-10): 
        return 0.0
    retval = St.Z_PCM+(St.Z_PCP-St.Z_PCM)*((St.FNB-St.X_PCM)/(St.X_PCP-St.X_PCM))
    retval = max(retval,0.0)
   
    return(retval)


#_________________________________________________________________________________________________________
    
def reversal(St) :
    
    KT = St.I1
    St.I1 = St.I2
    St.I2 = KT
    KT = St.L1
    St.L1 = St.L2
    St.L2 = KT

            
def Case2(I, Sub_Map, St) :
    
    St.J2=4 
    St.X_INCR2 = Sub_Map[St.I2+St.J2-1]
    St.Z_W2 = Sub_Map[St.I2+2]
    St.Z_W1 = Sub_Map[St.I1+2]
    
    #Valor de Z en X=0
    St.Z_PCM = St.Z_W1+(St.Z_W2-St.Z_W1)*St.DFL
    St.X_PCM = 0.0
    #Primer paso en la curva 2 (BT)  
    St.X_W2 = St.X_W2+St.X_INCR2 
    St.Z_W2 = St.Z_W2-St.FI_Step 
    #pendiente OBT
    St.SL2 = St.Z_W2/St.X_W2
    St.X_INCR1 = 0
    St.SL1 = -900000
    return (TRARA4 (I, Sub_Map, St))
           
    
def Case3(I, Sub_Map, St) :

    #Valor de Z en X=0
    St.Z_PCM = St.Z_W1+(St.Z_W2-St.Z_W1)*St.DFL
    St.X_PCM = 0.0
    #Primer paso en la curva 2 (BT)  
    St.X_W2 = St.X_W2+St.X_INCR2 
    St.Z_W2 = St.Z_W2-St.FI_Step 
    #pendiente OBT
    St.SL2 = St.Z_W2/St.X_W2
    St.J1=4
    St.X_INCR1 = Sub_Map[St.I1+St.J1 +1]                                                
    #Primer paso en la curva 1 (CB)   
    St.X_W1 = St.X_W1+St.X_INCR1
    St.Z_W1 = St.Z_W1-St.FI_Step
    St.SL1 = St.Z_W1/St.X_W1
    return(TRARA4(I, Sub_Map, St))


def Case4(I, Sub_Map, St) :

    St.SL2 = St.Z_W2/St.X_W2
    

    while St.J1 <= St.L1 : #for        
        St.X_INCR1 = Sub_Map[St.I1+St.J1-1]
        St.X_W1 = St.X_W1+St.X_INCR1
        St.Z_W1 = St.Z_W1-St.FI_Step
        St.XB_W1 = ((St.Z_W1/St.FI_Step)*St.X_INCR1+St.X_W1)/((St.X_INCR1/St.FI_Step)*St.SL2+1.0)
        
        if St.XB_W1 <= St.X_W1 : #(Goto 31)

            St.X_PCM = St.XB_W1+(St.X_W2-St.XB_W1)*St.DFL
            St.Z_PCM = St.X_PCM*St.SL2
            St.Z_W2 = St.Z_W2-St.FI_Step
            St.X_W2 = St.X_W2+St.X_INCR2
            St.SL1 = St.Z_W1/St.X_W1
            St.SL2 = St.Z_W2/St.X_W2 
            return(TRARA4(I, Sub_Map, St))

            if St.XB_W1 >= St.X_W2 : #(Goto 31)
           
                St.X_W1=0.0
                St.X_W2=0.0 
                St.J2=4 
                St.X_INCR2 = Sub_Map[St.I2+St.J2-1]
                St.Z_W2 = Sub_Map[St.I2+2]
                St.Z_W1 = Sub_Map[St.I1+2] 
                return(Case3(Sub_Map, St))   
           
        St.J1=St.J1+1
    
   
    if St.XB_W1  <= St.X_W2 : 
        return(0.0) 

    
def POPULATE_ARRAYS(n,m) :
    
    import os
    
    Map_prtns= 0
    Prtns_name = os.path.normpath(r'C:\Users\cd67\Desktop\David\1 Documentos\Proyecto de Grado\Program\3 Python\ap8min.txt')
    Map_prtns= np.array(np.loadtxt(Prtns_name, int))                           #Se cargan los coeficientes del mapa 
    Map_prtns= Map_prtns.reshape(m)                                            #Se redimenciona el vector para trabajarlo como una lista de datos
    
    if not(Map_prtns.all):
        print('No se pudo abir el archivo AP8MIN') 

    Map_eltns = 0
    Eltns_name = os.path.normpath(r'C:\Users\cd67\Desktop\David\1 Documentos\Proyecto de Grado\Program\3 Python\ae8min.txt')
    Map_eltns = np.loadtxt(Eltns_name, int) 
    
#    Map_eltns= Map_eltns.reshape(n)                                            #Se redimenciona el vector para trabajarlo como una lista de datos
                                 #Se cargan los coeficientes del mapa   
    Map_eltns= Map_eltns.reshape(n+8)                                          #Se redimenciona el vector para trabajarlo como una lista de datos
    Map_eltns= Map_eltns[0:-8]                                                 #Se eliminan lo 8 últimos elementos del mapa los cuales no existen en el mapa original 
    
    if not(Map_eltns.all):
        print('No su pudo abir el archivo AE8MIN') 
        
    return Map_eltns, Map_prtns
    
#_________________________________________________________________________________________________________    

def PLOT_FLUX(ax, Om_Flux, D, L ) :                                             #Función Gráficas                
    Q, R =  np.meshgrid(L, D)                                                      #Malla
    plt.contourf(Q, R, Om_Flux, locator=ticker.LogLocator())
    plt.colorbar()

def CONTOUR_FLUX(ax, Om_Flux, D, L, r_start, r_end ) :
    Q, R =  np.meshgrid(L, D)                                                   #Curvas de nivel del flujo omnidireccional
    cs = plt.contour(Q, R, Om_Flux, locator=ticker.LogLocator(), Polar= True)
#    manual_locations = [(0.05*np.pi, 6.5), (0*np.pi, 5.2), (0.05*np.pi,4.8),(0*np.pi, 4.2), (0.05*np.pi, 3.5), (-0.05*np.pi+np.pi, 6.5), (0*np.pi+np.pi, 5.2), (-0.05*np.pi+np.pi,4.8), (0*np.pi+np.pi, 4.2), (-0.05*np.pi+np.pi, 3.5)]
#    manual_locations = [(0.05*np.pi, 3.5), (0*np.pi, 3.3), (-0.05*np.pi,3.1),(0.01*np.pi, 2.8), (0.05*np.pi, 2.1), (-0.05*np.pi+np.pi, 3.5), (0*np.pi+np.pi, 3.3), (0.05*np.pi+np.pi,3.1), (-0.01*np.pi+np.pi, 2.8), (-0.05*np.pi+np.pi, 2.1)]
#    manual_locations = [(0*np.pi, 2.8), (0.04*np.pi, 2.5), (0*np.pi,2.1),(0.04*np.pi, 1.6), (0*np.pi+np.pi, 2.8), (-0.04*np.pi+np.pi, 2.5), (0*np.pi+np.pi,2.1), (-0.04*np.pi+np.pi, 1.6)]
    fmt=ticker.LogFormatterMathtext()
    labels= plt.clabel(cs, inline=1, fmt=fmt)
    # manual=manual_locations
    for l in labels:
        l.set_rotation(0)
    
def FORMAT_AXES(ax, r_end):                                                    #Función formato de los ejes
    ax.set_aspect('equal')                         
    ax.figure.subplots_adjust(bottom=0, top=1, left=0, right=1)
    ax.xaxis.set_ticks(np.arange(0, 2*np.pi, np.pi/4))                             #Marcas theta
    ax.yaxis.set_ticks(np.arange(0, r_end+1, 1))                                   #Marcas r
    for spine in ax.spines.itervalues():
        spine.set_visible(False)                                                
#_________________________________________________________________________________________________________    

def B_POTENTIAL(BE = 3.11**(-5), RE = 6378.137):                                #Función Potencial Magnético
    return  BE*(RE**3)*sym.sin(q)/r**2                                             #BE(Campo magnético ecuatorial en la superficie de la tierra)
                                                                                   #RE(Distancia normalizada a radios Terrestres)
def B_FIELD(psi):                                                               #Función Campo Magnético
    u = sym.lambdify((r, q),  psi.diff(r), 'numpy')                                #Componente radial
    v = sym.lambdify((r, q),  (1/r)*psi.diff(q), 'numpy')                          #Componente polar
    return u, v
    
def PLOT_STREAMLINES(ax, u, v, rlim, qlim, Scale):                             #Función lineas de Campo vectorial
    r0, r1 = rlim                                                                  #Fronteras r                    
    q0, q1 = qlim                                                                  #Fronteras theta
    Rs, Qs =  np.ogrid[r0:r1:100j, q0:q1:100j]                                     #Malla
    ax.streamplot(Qs, Rs, v(Rs, Qs), u(Rs, Qs), density=0.5, color='cornflowerblue')

#_________________________________________________________________________________________________________    
    
Range=100
Om_Flux= np.empty([Range,Range])                                                   #Matriz de Flujo Omnidireccional en función de L y BB0                                        
mElect= []                                                                     
mProt= []                                                    

Descr_eltns = np.array([8, 4, 1964, 6400, 2100, 1024, 1024, 13168])            #Vector de cabecera del mapa AEMIN 
#Descr_eltns = np.array([7, 4, 1970, 6400, 2100, 1024, 1024, 13548])
Descr_prtns = np.array([2, 4, 1964, 100, 2048, 2048, 1024, 16584])             #Vector de cabecera del mapa APMIN
#Descr_prtns = np.array([2, 4, 1970, 100, 2048, 2048, 1024, 16296])

r_start= 1.0                                                                     
r_end= 3.0
q_start= 0
q_end= 2*np.pi
rlim = (r_start, r_end)                                                        #Fronteras r
qlim = (q_start, q_end)                                                        # #Fronteras theta


Distancia = np.linspace(r_start,r_end, num=Range)
Latitud =  np.linspace(q_start, q_end, num=Range)  
#Distancia = input('Ingrese el rango principal en "Km": ')
#Latitud = input('Ingrese la latitud magnetica en "Grados": ')
Map_eltns, Map_prtns = POPULATE_ARRAYS(Descr_eltns[7],Descr_prtns[7])

BE = 31000.0                                                                   #Campo magnético ecuatorial en la superficie de la tierra
                             
Map = Map_prtns

En = 10.0                                                                 #Se define el valor de la energía
ii=np.arange(Range)

for ri in zip(Distancia, ii) :                                                 
    for li in zip(Latitud, ii) :
        RE = ri[0]                                                             #Distancia normalizada a radios Terrestres
        Alt_rad = li[0]                                                        #Latitud Magnética
        B_Dipolo = BE*mth.sqrt(3.0*(mth.sin(Alt_rad)**2)+1.0)/(RE**3)          #Aproximación del campo magnético terrestre a un dipolo [nT]
        L_value = RE/(mth.cos(Alt_rad)**2)                                     #Cálculo del valor de McIlwain
        BB0 = (B_Dipolo*(L_value**3))/BE                                       #Campo magnético normalizado al valor ecuatorial
        Flux=TRARA1(L_value, BB0, En, Descr_prtns, Map)                        #Flujo logarítmico
        mProt=10.0**(Flux)                                                     #Flujo Omndirecconal
        Om_Flux[ri[1],li[1]] = mProt              
        
#        print('=============================================================')
 #       print('RE es: ', ri, 'Theta: ', li)
  #      print('BB0 es: ', BB0, 'Lvalue es: ', L_value)
   #     print('El flujo de protones con energia', En, 'es', mProt)


PotB = B_POTENTIAL()                                                            
u, v = B_FIELD(PotB)
Scale=8         
fig, ax = plt.subplots(figsize=(Scale, Scale), subplot_kw=dict(projection='polar'))  
PLOT_FLUX(ax, Om_Flux, Distancia, Latitud)
#CONTOUR_FLUX(ax, Om_Flux, Distancia, Latitud, r_start, r_end)
#PLOT_STREAMLINES(ax, u, v, rlim, qlim, Scale)
circle = pl.Circle((0, 0), 1, transform=ax.transData._b, color='cornflowerblue', linewidth=0)   
ax.add_artist(circle)                                                          #Representación de la tierra
FORMAT_AXES(ax, r_end)    
fig.savefig(r'AP8MIN20.0.png', dpi=720) 

