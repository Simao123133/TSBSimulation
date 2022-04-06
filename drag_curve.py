import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from Boat_Model import boat_update
import time 
from matplotlib.ticker import (AutoMinorLocator, MultipleLocator)
from numpy.polynomial import Polynomial
from scipy.optimize import curve_fit

rho= 1025.00 # Water density(kg/m**3)   <-I
mu= 1.14E-3 # Coef. Viscosidade Dinâmica para calular Reynolds(Pa*s)   <-I
c_FF=0.08 # Corda foil frente(m)

def cdpf(x,y):
    #determining coefficient of pressure drag of the foils
    p00 =     0.01465
    p10 =  -0.0003439
    p01 =  -5.436e-05
    p20 =  -0.0002518
    p11 =   1.629e-06
    p02 =   8.999e-08
    p30 =   1.962e-05
    p21 =   5.111e-07
    p12 =  -2.115e-09
    p03 =  -7.446e-11
    p40 =    6.85e-06
    p31 =  -4.153e-08
    p22 =  -2.604e-10
    p13 =   1.213e-12
    p04 =   3.002e-14
    p50 =  -3.572e-07
    p41 =  -2.379e-09
    p32 =   1.914e-11
    p23 =   3.918e-14
    p14 =  -2.599e-16
    p05 =  -4.682e-18
    
    return (p00 + p10*x + p01*y + p20*x**2 + p11*x*y + p02*y**2 + p30*x**3 + p21*x**2*y
    + p12*x*y**2 + p03*y**3 + p40*x**4 + p31*x**3*y + p22*x**2*y**2 
    + p13*x*y**3 + p04*y**4 + p50*x**5 + p41*x**4*y + p32*x**3*y**2 
    + p23*x**2*y**3 + p14*x*y**4 + p05*y**5)

def f(x, p00,p10,p01,p20,p11,p02,p30,p21,p12,p03,p31,p22,p13,p04,p41,p32,p23,p14,p05):
    aa = x[0]
    Re = x[1]
    x=aa
    y=Re

    return (p00 + p10*x + p01*y + p20*x**2 + p11*x*y + p02*y**2 + p30*x**3 + p21*x**2*y 
        + p12*x*y**2 + p03*y**3 + p31*x**3*y + p22*x**2*y**2 
        + p13*x*y**3 + p04*y**4 + p32*x**3*y**2 
        + p23*x**2*y**3 + p14*x*y**4 + p05*y**5)

ang = 220
g = 0
d_v = np.zeros(ang*170)
aa_v = np.zeros(ang*170)
Re_v = np.zeros(ang*170)

for i in range(0, ang):
    k = i/10 - 5
    for j in range(0, 170):
        l = j/10
        Re = (rho)*l*(c_FF)/(mu) #Num Reynolds foil frente
        d_v[g] = cdpf(k, Re/1000)
        aa_v[g] = k
        Re_v[g] = Re/1000
        g = g + 1

popt, pcov = curve_fit(f, (aa_v, Re_v), d_v)

def f_aa(aa, speed, p00,p10,p01,p20,p11,p02,p30,p21,p12,p03,p31,p22,p13,p04,p41,p32,p23,p14,p05):
    rho= 1025.00 # Water density(kg/m**3)   <-I
    mu= 1.14E-3 # Coef. Viscosidade Dinâmica para calular Reynolds(Pa*s)   <-I
    c_FF=0.08 # Corda foil frente(m)
    Re = (rho)*speed*(c_FF)/(mu)/1000 #Num Reynolds foil frente
    x=aa
    y=Re

    return (p00 + p10*x + p01*y + p20*x**2 + p11*x*y + p02*y**2 + p30*x**3 + p21*x**2*y 
        + p12*x*y**2 + p03*y**3 + p31*x**3*y + p22*x**2*y**2 
        + p13*x*y**3 + p04*y**4 + p32*x**3*y**2 
        + p23*x**2*y**3 + p14*x*y**4 + p05*y**5)


aa = np.zeros(ang)
cd_v_aa = np.zeros(ang)
speed = 7
Re = (rho)*speed*(c_FF)/(mu) #Num Reynolds foil frente

for i in range(0, ang):
    k = i/10 - 5
    aa[i] = k
    cd_v_aa[i] = cdpf(k, Re/1000)

print(*popt)
plt.plot(aa,  f_aa(aa, speed, *popt), aa, cd_v_aa)
plt.show()
