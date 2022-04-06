from sympy import *
from sympy.utilities.lambdify import lambdify, implemented_function

a37, a29, a30, a38, a39, a40, a41, a42, a10, a11, a43, a45, a44, a46 = symbols('a37, a29, a30, a38, a39, a40, a41, a42, a10, a11, a43, a45, a44, a46')
eff_aa_LF, eff_aa_RF, Re_LF, Re_RF = symbols('eff_aa_LF, eff_aa_RF, Re_LF, Re_RF')

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

cl_LF = cl(eff_aa_LF,Re_LF/1000)
cl_RF = cl(eff_aa_RF,Re_RF/1000)

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

cdp_LF = cdpf(eff_aa_LF,Re_LF/1000)
cdp_RF = cdpf(eff_aa_RF,Re_RF/1000)

f = ((a37 - (a29*cl_LF + a30*cl_RF + a38*cl_LF**2 + a39*cl_RF**2 + a40*cdp_LF + a41*cdp_RF))*180/9027.5/3)**2 + (a42-(a10*cl_LF + a11*cl_RF + a43*cl_LF**2 + a45*cl_RF**2 + a44*cdp_LF + a46*cdp_RF))**2 

gradient_LF = f.diff(eff_aa_LF)
gradient_RF = f.diff(eff_aa_RF)
print(simplify(gradient_RF))


hessian11 = gradient_LF.diff(eff_aa_LF)
hessian12 = gradient_LF.diff(eff_aa_RF)
hessian21 = gradient_RF.diff(eff_aa_LF)
hessian22 = gradient_RF.diff(eff_aa_RF)




