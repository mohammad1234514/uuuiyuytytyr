import pandas as pd
import numpy as np
from math import exp , pi
import matplotlib.pyplot as plt
from scipy.integrate import   solve_ivp , odeint

R=8314*1000


MWCH4= 16.04      #g/mole  or kg/kmole
MWH2O = 18.01     #g/mole  or kg/kmole
MWCO  = 28.01     #g/mole  or kg/kmole
MWH2 = 2.016      #g/mole  or kg/kmole
MWCO2 = 44.01     #g/mole  or kg/kmole
MWN2 = 28.01      #g/mole  or kg/kmole

#catalyst properties
# Ni/MgAl2O4 
dp = 0.0173  #particle size                       #m
rocata = 2355.2  #density of catalyst             #kg/m3
eps = 0.5   #catalyst bed porosity                #dimensionless

# Tube property
din = 0.1  # inner diameter of tube               #m
dout= 0.13 #outer diameter of tube                #m
L = 12     # length of reactor                    #m
Ac = pi *( din**2)/4  # Area of tube              #m2

#import table 4 (heat capacity data)
table4_data = pd.read_excel (r'copoi.excel.xlsx')

A = table4_data.A
B = table4_data.B
C = table4_data.C
D = table4_data.D

# Cp = 8.314 * (A + B*T + C*T**2 + D/T**2)        #kj/kmole

#import table 5 (viscosity data)
table5_data = pd.read_excel (r'yorits.xlsx') 
a = table5_data.A
b = table5_data.B
c = table5_data.C
d = table5_data.D
# Mo = (a*T**b )/(1+c/T + d/T**2)      #kg / m.s

#initial condition
Fin = 21                 #kmole/hr
Fco20 = 4.846154         #kmole/hr
Fch40= 3.230769          #kmole/hr
Fh2o0= 12.923076         #kmole/hr

Nco20 = Fco20/Ac         #kmole/hr.m2
Nch40 = Fch40/Ac         #kmole/hr.m2
Nh2o0 = Fh2o0/Ac         #kmole/hr.m2
Nco0 = 0.001             #kmole/hr.m2
Nh20 = 0.001             #kmole/hr.m2
Nn20 = 0.001             #kmole/hr.m2
T0 = 973                 #kelvin
P0 = 25                  #bar

# radial heat flux
qr= 65                  #kw/m2

def modelofreactor(H,z):

    NCH4 = H[0]
    NH2O = H[1]
    NCO = H[2]
    NH2 = H[3]
    NCO2 = H[4]
    NN2 = H[5]
    T = H[6]
    P =H[7]

    Ntotal = NCH4 + NH2O + NCO + NH2 + NCO2 + NN2 
    
    XCH4 = NCH4 / Ntotal
    XH2O = NH2O / Ntotal
    XCO = NCO / Ntotal
    XH2 = NH2 / Ntotal
    XCO2 = NCO2 / Ntotal
    XN2  = NN2 / Ntotal
    
    PCH4 = XCH4 * P
    PH2O = XH2O * P
    PCO = XCO * P
    PH2 = XH2 * P
    PCO2 = XCO2 * P
    PN2  = XN2 * P
    
    #calculating heat capacity varying Temp
    Cp = 8.314 * (A + B*T + C*T**2 + D/T**2)  #kj/kmole
    Cpch4 = Cp[0]
    Cph2o = Cp[1]
    Cpco= Cp[2]
    Cph2= Cp[3]
    Cpco2= Cp[4]
    Cpn2= Cp[5]

    #calculating delta Cp of reactions
    deltaCpR1 = 3 * Cph2 + Cpco - Cph2o - Cpch4
    deltaCpR2 = Cph2 + Cpco2 - Cph2o - Cpco
    deltaCpR3 = 4 * Cph2 + Cpco2 - Cph2o - Cpch4

    deltaH01 = 206.11*1000   #converting kj/mole.k to Kj/kmole.K
    deltaH02 = -41.17*1000   #converting kj/mole.k to Kj/kmole.K
    deltaH03 = 164.94*1000   #converting kj/mole.k to Kj/kmole.K

    deltaH1 = deltaH01 + deltaCpR1*(T-298)    #kj/kmole
    deltaH2 = deltaH02 + deltaCpR2*(T-298)    #kj/kmole
    deltaH3 = deltaH03 + deltaCpR3*(T-298)    #kj/kmole
    
    #calculating viscosity varying Temp
    Mo = (a*T**b )/(1+c/T + d/T**2)   #kg / (m.s)
    
    M = [MWCH4 , MWH2O ,MWCO, MWH2 , MWCO2, MWN2 ]  #molecular weight g/mole
    X = [XCH4 , XH2O , XCO , XH2 , XCO2 , XN2]
    sigma = 0
    Mogmix= 0
    
    for i in range(6):
        for j in range(6):
            sigma = sigma + X[i]*(M[j]/M[i])**0.5 
        Mogmix = Mogmix +Mo[i] * X[i]/sigma

    
    #calculating mass flux
    # MWave = XCH4*MWCH4 + XH2O*MWH2O + XCO*MWCO + XH2*MWH2 + XCO2*MWCO2 + XN2*MWN2
    
    #calculating mass flux
    Gflux = (Nch40*MWCH4 + Nh2o0*MWH2O + Nco20*MWCO2)
    
    #friction factor calculating
    f= 150 + 1.75 *(1/(1-eps))*(Gflux * dp/(Mogmix*3600))
    
    #calculating density of gas
    rog = 0.30    #data hysys
    
    #Table 6 data
    #Expressions for kinetic, equilibrium and adsorption constants.
    K1 = 9.49 * 10**15 *exp(-240100/(R*T)) 
    K2 = 4.39* 10**16 *exp(-67130/(R*T))
    K3 = 2.29 * 10**15*exp(-243900/(R*T))
    Keq1 = exp((-26830/T)+30.114)
    Keq2 = exp((4400/T)-4.063)
    Keq3 = Keq1 * Keq2
    Kch4 = 6.65 * 10**-4 * exp(38280/(R*T))
    Kh2o = 1.77 * 10 **5 * exp(-88680/(R*T))
    Kco = 8.23 * 10**-5 * exp(70650/(R*T))
    Kh2 = 6.12 * 10**-9 * exp(82900/(R*T))
    
    #Rate expression for the reactions
    DEN = 1 + Kch4*PCH4 + Kco*PCO + Kh2*PH2 +(Kh2o*PH2O)/(PH2)
    R1 =((K1/PH2**2.5) *(PCH4*PH2O-(PH2**3*PCO/Keq1)))/DEN**2
    R2 =((K2/PH2) *(PCO*PH2O-(PH2*PCO2/Keq2)))/DEN**2
    R3 =((K3/PH2**3.5) *(PCH4*PH2O**2-(PH2**4*PCO2/Keq3)))/DEN**2
    

    #effectivness factor for hemogenous reaction
    eta1 = 1
    eta2 = 1
    eta3 = 1
    
    # diffrential equations
    dNCH4dz = rocata*(1-eps)*(-1*R1*eta1 -1*R3*eta3)
    dNH2Odz = rocata*(1-eps)*(-1*R1*eta1 -1*R2*eta2 - 2*R3*eta3)
    dNCOdz = rocata*(1-eps)*(+1*R1*eta1 -1*R2*eta2 )
    dNH2dz = rocata*(1-eps)*(+3*R1*eta1 + 1*R2*eta2 + 4*R3*eta3)
    dNCO2dz = rocata*(1-eps)*(+1*R2*eta2 + 1 *R3*eta3)
    dNN2dz = 0 
    dTdz = (rocata*(1-eps)*(-deltaH1*R1*eta1-deltaH2*R2*eta2-deltaH3*R3*eta3)+4*qr*3600/din)/(NCH4*Cpch4+NH2O*Cph2o+NCO*Cpco+NH2*Cph2+NCO2*Cpco2+NN2*Cpn2)
    dPdz= (f*Gflux*3600*Mogmix*(1-eps)**2)/((rog*dp**2)*eps**3)

    return [dNCH4dz , dNH2Odz , dNCOdz , dNH2dz , dNCO2dz , dNN2dz , dTdz , dPdz]
    return func (*args,  **kwargs)

H0 = [Nch40, Nh2o0 , Nco0*.005, Nh2o0*0.005 , Nco20*.005 , Nn20*0.005 , T0 , P0 ]
z= np . linspace (0,12,50)
H = odeint(modelofreactor , H0 , z )

NCH4=H[:,0]  
NH2O=H[:,1]  
NCO=H[:,2]  
NH2=H[:,3]  
NCO2=H[:,4]
NN2=H[:,5]  
T=H[:,6]  
P=H[:,7] 

plt.ion()

plt.plot(z,T)
plt.show(block=True)
