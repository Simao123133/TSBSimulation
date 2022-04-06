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

def cl(x,y): 
        #determining coefficient of lift of the foils for a set angle of attack
        p00 =     0.04647
        p10 =      0.1369
        p01 =    0.001502
        p20 =   0.0005872
        p11 =  -0.0001138
        p02 =  -2.794e-06
        p30 =  -0.0007229
        p21 =  -2.662e-07
        p12 =   1.545e-07
        p03 =   2.428e-09
        p40 =   1.647e-05
        p31 =   6.574e-07
        p22 =  -3.892e-09
        p13 =  -7.221e-11
        p04 =  -9.966e-13
        p50 =  -6.807e-08
        p41 =  -2.534e-09
        p32 =  -2.657e-10
        p23 =   2.267e-12
        p14 =   9.751e-15
        p05 =   1.555e-16
        
        return (p00 + p10*x + p01*y + p20*x**2 + p11*x*y + p02*y**2 + p30*x**3 + p21*x**2*y 
        + p12*x*y**2 + p03*y**3 + p40*x**4 + p31*x**3*y + p22*x**2*y**2 
        + p13*x*y**3 + p04*y**4 + p50*x**5 + p41*x**4*y + p32*x**3*y**2 
        + p23*x**2*y**3 + p14*x*y**4 + p05*y**5)

def f(x, a, b, c, d, e, f, g, h, i, j, k):
    aa = x[0]
    Re = x[1]
    
    return (a + b*aa + c*Re + d*aa*Re + e*Re**2 + f*aa*Re**2 + g*Re**3 + h*aa*Re**3 + i*Re**4 + j*aa*Re**4 + k*Re**5)

g = 0
cl_v = np.zeros(210*170)
aa_v = np.zeros(210*170)
Re_v = np.zeros(210*170)

for i in range(0, 210):
    k = i/10 - 5
    for j in range(0, 170):
        l = j/10
        Re = (rho)*l*(c_FF)/(mu) #Num Reynolds foil frente
        cl_v[g] = cl(k, Re/1000)
        aa_v[g] = k
        Re_v[g] = Re/1000
        g = g + 1

c = 130*170

popt, pcov = curve_fit(f, (aa_v[0:c], Re_v[0:c]), cl_v[0:c])
popt2, pcov2 = curve_fit(f, (aa_v[c:-1], Re_v[c:-1]), cl_v[c:-1])
#plt.plot(aa, f(aa, *popt), aa, cl_v, aa, f(aa, *popt2))
#plt.plot(aa, f(aa, *popt), aa, cl_v)

def f_aa(aa, speed, a, b, c, d, e, f, g, h, i, j, k):
    rho= 1025.00 # Water density(kg/m**3)   <-I
    mu= 1.14E-3 # Coef. Viscosidade Dinâmica para calular Reynolds(Pa*s)   <-I
    c_FF=0.08 # Corda foil frente(m)
    Re = (rho)*speed*(c_FF)/(mu)/1000 #Num Reynolds foil frente
    
    return (a + b*aa + c*Re + d*aa*Re + e*Re**2 + f*aa*Re**2 + g*Re**3 + h*aa*Re**3 + i*Re**4 + j*aa*Re**4 + k*Re**5)


aa = np.zeros(210)
cl_v_aa = np.zeros(210)
speed = 3
point = 8
Re = (rho)*speed*(c_FF)/(mu) #Num Reynolds foil frente

for i in range(0, 210):
    k = i/10 - 5
    aa[i] = k
    cl_v_aa[i] = cl(k, Re/1000)

#plt.plot(aa, f_aa(aa, speed, *popt), aa, cl_v_aa)
#plt.plot(aa, f_aa(aa, 7, *popt2), aa, cl_v_aa)
#plt.show()

## together middle point 6.5

cl_tog = np.zeros(210)

for i in range(0, 210):
    k = i/10 - 5
    if k <= point:
        cl_tog[i] = f_aa(aa[i], speed, *popt)
    else:
        cl_tog[i] = f_aa(aa[i], speed, *popt2)

plt.plot(aa, cl_tog, aa, cl_v_aa)
plt.show()