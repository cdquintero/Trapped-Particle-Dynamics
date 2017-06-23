# -*- coding: utf-8 -*-
"""Created on Tue Sep 27 17:48:58 2016

Blablablablablablalbalblalblalblablalabla

@author: cd67"""

#%matplotlib inline
import math  as mth
import numpy as np                                                             #Paquete Arrays
import matplotlib.pyplot as plt
                                               #Paquete Gráficas
#_________________________________________________________________________________________________________

class V_Static(object) :
    FI_Step = 0.0
    I1=0
    FKB1=0.0; FKB2=0.0; F_Incr2=0.0; F_Incr1=0.0
    FKBM=0.0; F_LogM=0.0; SL2=0.0; FNB=0.0; DFL=0.0
    J1=0; J2=0; I_Time=0; L1=0; L2=0
    I2=0; FLOG1=0; FLOG2=0
    FKBJ1=0.0; FKBJ2=0.0; SL1=0.0
    FKB=0.0; FLOG=0.0

def TRARA1 (FL, BB0, E, F, N, Descr, Map) :                                    #* Parametros arbitrarios
      
    St = V_Static()
    S0 = 0
    I0 = F0 = E0 = 0.0
    IE=1    
    F1 = 1.001; F2 = 1.002 
    St.FI_Step = float(Descr[6])/Descr[1]
    E_Scale = float(Descr[3])
    F_Scale = float(Descr[6])
    XNL = min(15.6,abs(FL))
    NL = XNL*Descr[4]                             
    if (BB0 < 1.0 ): 
        BB0 = 1.0                                             
    NB = ( BB0-1.0)*Descr[5]  
 
    St.I1 = 0                                                         
    I2 = Map[0]                                                                   #Número de elementos en el primer submapa de energías                                          
    I3 = I2+Map[I2]                                                               #Índice del último elemento del submapa 2  
    L3 = Map[I3]                                                                  #Número de elementos en el tercer submapa de energías                                                  
    E1 = Map[St.I1+1]/E_Scale                                                        #Energía del primer submapa de energías (Sin escala)  
    E2 = Map[I2+1]/E_Scale                                                        #Energía del segundo submapa de energías (Sin escala) 

    S1 = 1                                                                      #Variables lógicas que indican si un flujo f(E,B,L) ya ha sido calculado 
    S2 = 1                                                                      #antes de interpolar (TRARA2). Si no Si=True
               
    while IE<=N : #for
                                                            
        while (not((E[IE-1] <= E2) or (L3==0))) :
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
            
        if S1 :
            F1 = TRARA2((St.I1+2), NL, NB, Map,St)/F_Scale

        if S2 :
            F2 = TRARA2((I2+2), NL, NB, Map,St)/F_Scale
    
        S1 = 0
        S2 = 0                                                

        F[IE-1]=F1+(F2-F1)*(E[IE-1]-E1)/(E2-E1)
           
        if (F2 <= 0.0 or St.I1 == 0 ):
            if (S0):
                F0=TRARA2((I0+3-1), NL, NB, Map, St)/F_Scale
            S0=0;                                       
            F[IE-1] = min(F[IE-1] , F0+(F1-F0)*(E[IE-1]-E0)/(E1-E0))                   
    
        F[IE-1] = max(F[IE-1], 0.0)
        IE=IE+1
           
    return 

#_________________________________________________________________________________________________________

def TRARA2 (I, IL, IB, Sub_Map, St) :                                                  #* Parametros arbitrarios
    
    St.I2 = St.F_Log1 = St.F_Log2 = St.J1 = St.J2 = St.I_Time = St.L1 = St.L2 = 0
    St.FKB1 = St.FKB2 = St.F_Incr2 = St.F_Incr1 = St.FKBM = St.F_LogM = St.SL2 = St.FNB = St.DFL = 0

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
        KT = St.I1
        St.I1 = St.I2
        St.I2 = KT
        KT = St.L1
        St.L1 = St.L2
        St.L2 = KT
         
    while St.I_Time == 0 :
        FLL1 = float(Sub_Map[St.I1+1])
        FLL2 = float(Sub_Map[St.I2+1])
        St.DFL = (St.FNL-FLL1)/(FLL2-FLL1)
        St.F_Log1 = Sub_Map[St.I1+2]
        St.F_Log2 = Sub_Map[St.I2+2]
        St.FKB1 = 0.0
        St.FKB2 = 0.0       

        if St.L1 < 4 :  
           return(TRARA3(I, 32, Sub_Map, St))
        
        while St.J2 <= St.L2 : #for                
            St.F_Incr2 = Sub_Map[St.I2+St.J2-1]
                 
            if ((St.FKB2+St.F_Incr2) > St.FNB):
                return(TRARA3(I, 23, Sub_Map, St))
            St.FKB2 = St.FKB2+St.F_Incr2
            St.F_Log2 = St.F_Log2-St.FI_Step
            St.J2 = St.J2+1
            
        St.I_Time = St.I_Time+1
        KT = St.I1
        St.I1 = St.I2
        St.I2 = KT
        KT = St.L1
        St.L1 = St.L2
        St.L2 = KT
     
    return 0.0

#_________________________________________________________________________________________________________
    
def TRARA3 (I, Position, Sub_Map, St) :                                                #* Parametros arbitrarios
    
    St.FKBJ1 = St.FKBJ1 = St.SL1 = 0
    St.J1 = 4

    if Position == 23 and St.I_Time != 1 and St.J2 != 4 : 
        St.SL2 = St.F_Log2/St.FKB2
    
        while St.J1 <= St.L1 : #for        
            St.F_Incr1 = Sub_Map[St.I1+St.J1-1]
            St.FKB1 = St.FKB1+St.F_Incr1
            St.F_Log1 = St.F_Log1-St.FI_Step
            St.FKBJ1 = ((St.F_Log1/St.FI_Step)*St.F_Incr1+St.FKB1)/((St.F_Incr1/St.FI_Step)*St.SL2+1.0)
            if St.FKBJ1 <= St.FKB1 :
                return(TRARA4(I, 0, Sub_Map, St))
            St.J1=St.J1+1
        
        if St.FKBJ1 <= St.FKB2 :
            return(TRARA5(St))
        if St.FKBJ1 <= St.FKB2 : 
            return(TRARA4(I, 0, Sub_Map, St))
        St.FKB1 = 0.0
    
    St.FKB2 = 0.0
    
    if Position == 32 :
        St.J2=4 
        St.F_Incr2 = Sub_Map[St.I2+St.J2-1]
        St.F_Log2 = Sub_Map[St.I2+2]
        St.F_Log1 = Sub_Map[St.I1+2]
    
    St.F_LogM = St.F_Log1+(St.F_Log2-St.F_Log1)*St.DFL
    St.FKBM = 0.0
    St.FKB2 = St.FKB2+St.F_Incr2 
    St.F_Log2 = St.F_Log2-St.FI_Step 
    St.SL2 = St.F_Log2/St.FKB2
    if St.L1 < 4 :
        return(TRARA4(I, 35, Sub_Map, St))
    St.J1=4
    St.F_Incr1 = Sub_Map[St.I1+St.J1 +1]                                                
    St.FKB1 = St.FKB1+St.F_Incr1
    St.F_Log1 = St.F_Log1-St.FI_Step
    St.SL1 = St.F_Log1/St.FKB1
    
    return(TRARA4(I, 15, Sub_Map, St))

#_________________________________________________________________________________________________________

def TRARA4 (I, Start_Psn, Sub_Map, St) :                                               #* Parametros arbitrarios

    Pass_To20 = 0.0
    St.FKB = St.F_Log = 0.0

    if Start_Psn == 15 :
        Pass_To20 = 1 
    elif Start_Psn == 35 :
        St.F_Incr1 = 0
        St.SL1 = -900000
        Pass_To20 = 1 
    else :
        St.FKBM = St.FKBJ1+(St.FKB2-St.FKBJ1)*St.DFL
        St.F_LogM = St.FKBM*St.SL2
        St.F_Log2 = St.F_Log2-St.FI_Step
        St.FKB2 = St.FKB2+St.F_Incr2
        St.SL1 = St.F_Log1/St.FKB1
        St.SL2 = St.F_Log2/St.FKB2
                                          
    while 1 :
        if  St.SL1 >= St.SL2 and (not(Pass_To20)) : 
            St.FKBJ2 = ((St.F_Log2/St.FI_Step)*St.F_Incr2+St.FKB2)/((St.F_Incr2/St.FI_Step)*St.SL1+1.0)
            St.FKB = St.FKB1+(St.FKBJ2-St.FKB1)*St.DFL
            St.F_Log = St.FKB*St.SL1
            if St.FKB >= St.FNB :
               return(TRARA5(St))
            St.FKBM = St.FKB
            St.F_LogM = St.F_Log
            if(St.J1 >= St.L1) :
               return 0
            St.J1 = St.J1+1
            St.F_Incr1 = Sub_Map[St.I1+St.J1-1]
            St.F_Log1 = St.F_Log1-St.FI_Step
            St.FKB1 = St.FKB1+St.F_Incr1
            St.SL1 = St.F_Log1/St.FKB1    
           
        Pass_To20 = 0
        St.FKBJ1 = ((St.F_Log1/St.FI_Step)*St.F_Incr1+St.FKB1)/((St.F_Incr1/St.FI_Step)*St.SL2+1.0)
        St.FKB = St.FKBJ1+(St.FKB2-St.FKBJ1)*St.DFL
        St.F_Log = St.FKB*St.SL2
        if(St.FKB >= St.FNB) :
            return(TRARA5(St))
        St.FKBM = St.FKB
        St.F_LogM = St.F_Log
        if(St.J2 >= St.L2) :
            return 0
        St.J2 = St.J2+1;
        St.F_Incr2 = Sub_Map[St.I2+St.J2-1]
        St.F_Log2 = St.F_Log2-St.FI_Step
        St.FKB2 = St.FKB2+St.F_Incr2
        St.SL2 = St.F_Log2/St.FKB2
    
#_________________________________________________________________________________________________________
    
def TRARA5 (St) :                                                                #* Parametros arbitrarios

    if St.FKB < (St.FKBM+1.0e-10): 
        return 0.0
    retval = St.F_LogM+(St.F_Log-St.F_LogM)*((St.FNB-St.FKBM)/(St.FKB-St.FKBM))
    retval = max(retval,0.0)
   
    return(retval)


#_________________________________________________________________________________________________________
    
def POPULATE_ARRAYS(n,m,j,k) :
    
    import os
    
    Map_prtnsmax= 0
    Prtnsmax_name = os.path.normpath(r'C:\Users\cd67\Desktop\David\1 Documentos\Proyecto de Grado\Program\3 Python\ap8max.txt')
    Map_prtnsmax= np.array(np.loadtxt(Prtnsmax_name, int))
    Map_prtnsmax= Map_prtnsmax.reshape(m)
    
    if not(Map_prtnsmax.all):
        print('No se pudo abir el archivo AP8MAX') 

    Map_prtnsmin= 0
    Prtnsmin_name = os.path.normpath(r'C:\Users\cd67\Desktop\David\1 Documentos\Proyecto de Grado\Program\3 Python\ap8min.txt')
    Map_prtnsmin= np.array(np.loadtxt(Prtnsmin_name, int))
    Map_prtnsmin= Map_prtnsmin.reshape(n)
    
    if not(Map_prtnsmin.all):
        print('No se pudo abir el archivo AP8MIN') 
        
        
    Map_eltnsmin = 0
    Eltnsmin_name = os.path.normpath(r'C:\Users\cd67\Desktop\David\1 Documentos\Proyecto de Grado\Program\3 Python\ae8min.txt')
    Map_eltnsmin = np.loadtxt(Eltnsmin_name, int) 
    
    Map_eltnsmin= Map_eltnsmin.reshape(j+8)                                          #Se redimenciona el vector para trabajarlo como una lista de datos
    Map_eltnsmin= Map_eltnsmin[0:-8]                                                 #Se eliminan lo 8 últimos elementos del mapa los cuales no existen en el mapa original 
    
    if not(Map_eltnsmin.all):
        print('No su pudo abir el archivo AE8MIN')     
        
    Map_eltnsmax = 0
    Eltnsmax_name = os.path.normpath(r'C:\Users\cd67\Desktop\David\1 Documentos\Proyecto de Grado\Program\3 Python\ae8max.txt')
    Map_eltnsmax = np.loadtxt(Eltnsmax_name, int) 
    
    Map_eltnsmax= Map_eltnsmax.reshape(k)                                            #Se redimenciona el vector para trabajarlo como una lista de datos
   
    if not(Map_eltnsmax.all):
        print('No su pudo abir el archivo AE8MAX')            
        
    return Map_prtnsmin, Map_prtnsmax, Map_eltnsmin, Map_eltnsmax
    
#_________________________________________________________________________________________________________    
                                                                              
Num_energ = 5  
Num_energe = 5
Range=1000                                                                  #Número de energías
Enp = np.empty(Num_energ)
Ene = np.empty(Num_energe)
Flux= np.empty(Num_energ)
Fluxe= np.empty(Num_energe)
Omp_Fluxmin = np.empty([Range, Num_energ]) 
Omp_Fluxmax = np.empty([Range, Num_energ])     
Ome_Fluxmin = np.empty([Range, Num_energe]) 
Ome_Fluxmax = np.empty([Range, Num_energe])                                                                                                              
mProt= np.empty(Num_energ)                                                   
mElct= np.empty(Num_energe)                                          

Descr_prtnsmax = np.array([2, 4, 1970, 100, 2048, 2048, 1024, 16296])
Descr_prtnsmin = np.array([2, 4, 1964, 100, 2048, 2048, 1024, 16584])

Descr_eltnsmin = np.array([8, 4, 1964, 6400, 2100, 1024, 1024, 13168])            #Vector de cabecera del mapa AEMIN 
Descr_eltnsmax = np.array([7, 4, 1970, 6400, 2100, 1024, 1024, 13548])

Map_prtnsmin, Map_prtnsmax, Map_eltnsmin, Map_eltnsmax = POPULATE_ARRAYS(Descr_prtnsmin[7],Descr_prtnsmax[7], Descr_eltnsmin[7],Descr_eltnsmax[7])

BB0_start= 1                                                                     
BB0_end= 10.0

L_value = 1.5                                         #Cálculo del valor de McIlwain
BB0 = np.linspace(BB0_start, BB0_end, num=Range)                                                 #Campo magnético normalizado al valor ecuatorial

Mapmin = Map_prtnsmin
Mapmax = Map_prtnsmax

Enp[0] = 0.10
Enp[1] = 10.50
Enp[2] = 30.70
Enp[3] = 50.00
Enp[Num_energ-1] = 1.50

Ene[0] = 0.10
Ene[1] = 0.50
Ene[2] = 1.00
Ene[3] = 1.50
Ene[4] = 3.00
#Ene[5] = 6.00
#Ene[6] = 7.00
#Ene[Num_energe-1] = 7.00

ii=np.arange(Range)
iE=np.arange(Num_energ)
    
for jj in zip(BB0, ii) :                                                 
    B = jj[0]                                      #Campo magnético normalizado al valor ecuatorial
    TRARA1(L_value, B, Enp, Flux, Num_energ, Descr_prtnsmin, Mapmin)
    for kk in zip(Flux, iE) :
        mProt[kk[1]]=(10.0**(kk[0]))                                                    #Flujo Omndirecconal
    Omp_Fluxmin[jj[1]] = mProt 
    
for jj in zip(BB0, ii) :                                                 
    B = jj[0]                                      #Campo magnético normalizado al valor ecuatorial
    TRARA1(L_value, B, Enp, Flux, Num_energ, Descr_prtnsmax, Mapmax)
    for kk in zip(Flux, iE) :
        mProt[kk[1]]=(10.0**(kk[0]))                                                    #Flujo Omndirecconal
    Omp_Fluxmax[jj[1]] = mProt        
    
Mapmin = Map_eltnsmin
Mapmax = Map_eltnsmax   

iEe = np.arange(Num_energe)   

for jj in zip(BB0, ii) :                                                 
    B = jj[0]                                       #Campo magnético normalizado al valor ecuatorial
    TRARA1(L_value, B, Ene, Fluxe, Num_energe, Descr_eltnsmin, Mapmin)
    for kk in zip(Fluxe, iEe) :
        mElct[kk[1]]=(10.0**(kk[0]))                                                    #Flujo Omndirecconal
    Ome_Fluxmin[jj[1]] = mElct 
    
for jj in zip(BB0, ii) :                                                 
    B = jj[0]                                        #Campo magnético normalizado al valor ecuatorial
    TRARA1(L_value, B, Ene, Fluxe, Num_energe, Descr_eltnsmax, Mapmax)
    for kk in zip(Fluxe, iEe) :
        mElct[kk[1]]=(10.0**(kk[0]))                                                    #Flujo Omndirecconal
    Ome_Fluxmax[jj[1]] = mElct       
    
fig, ax = plt.subplots()

for i in iE :    
    ax.plot(BB0, Omp_Fluxmin[:,i], label= 'min %10.2f MeV' % (Enp[i]))   
    ax.plot(BB0, Omp_Fluxmax[:,i], linestyle='--', color='black', label= 'max%10.2f MeV' % (Enp[i]))

plt.rc('font', family='serif')
plt.ylabel(r'Flujo omnidireccional $[\frac{Protones}{cm^2 s}]$', size='large',  family='monospace', weight= 'light')
plt.xlabel(r'B/B0',  size='large',  family='monospace', weight= 'light')
plt.legend()
plt.title(r'Flujo omnidireccional Vs B/B0', size='large',  family='monospace', weight= 'light')
plt.legend()
plt.xlim(1,5) 
ax.semilogx()
ax.semilogy()

fig2, ax2 = plt.subplots()

for i in iEe :    
    ax2.plot(BB0, Ome_Fluxmin[:,i], label= 'min %10.2f MeV' % (Ene[i]))
    ax2.plot(BB0, Ome_Fluxmax[:,i], linestyle='--', color='black', label= 'max %10.2f MeV' % (Ene[i]))

plt.rc('font', family='serif')
plt.ylabel(r'Flujo omnidireccional $[\frac{Protones}{cm^2 s}]$', size='large',  family='monospace', weight= 'light')
plt.xlabel(r'B/B0',  size='large',  family='monospace', weight= 'light')
plt.legend()
plt.title(r'Flujo omnidireccional Vs B/B0', size='large',  family='monospace', weight= 'light')
plt.legend()
#plt.xlim(1,15) 
ax2.semilogx()
ax2.semilogy()

