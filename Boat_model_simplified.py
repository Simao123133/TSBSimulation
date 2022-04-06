import numpy as np
from numpy.core.numerictypes import obj2sctype
from sympy import *
from sympy.polys.partfrac import apart_undetermined_coeffs
from sympy.vector import CoordSys3D
import math as mt

Ixx = symbols('Ixx')
Ixz = symbols('Ixz')
Iyy = symbols('Iyy')
Izz = symbols('Izz')

L = symbols('L')
N = symbols('N')

wx = symbols('wx')
wy = symbols('wy')
wz = symbols('wz')

w = Matrix([wx, wy, wz])
J = Matrix([ [Ixx, 0, Ixz] , [0, Iyy , 0], [Ixz, 0, Izz] ])
Tau = Matrix([L, 0, N])
JinvTAU = J.inv()*Tau

print(JinvTAU[0])
#(u_phi + GyrTerm)*(Ixx*Izz - Ixz**2) = -Ixz*N + Izz*L

L = -Y_LS*LS_CoP - Z_LF*dist_FF_cm_y-Y_RS*RS_CoP + Z_RF*dist_FF_cm_y-Y_C*C_CoP - Y_Prop*(dist_ReF_cm_z-0.2) 
Z_LF = L_b_LF + D_b_LF

L_b_LF = L_LF*(-cosd(aa_LF))
D_b_LF = D_LF*(-sind(aa_LF))

Z_RF = L_b_RF + D_b_RF

L_b_RF = L_RF*(-cosd(aa_RF))
D_b_RF = D_RF*(-sind(aa_RF))

(u_phi + GyrTerm)*(Ixx*Izz - Ixz**2) = Izz*(-Y_LS*LS_CoP-Y_RS*RS_CoP - Y_C*C_CoP - Y_Prop*(dist_ReF_cm_z-0.2) ) - Izz*Z_LF*dist_FF_cm_y + Izz*Z_RF*dist_FF_cm_y


(u_phi + GyrTerm)*(Ixx*Izz - Ixz**2) - Izz*(-Y_LS*LS_CoP-Y_RS*RS_CoP - Y_C*C_CoP - Y_Prop*(dist_ReF_cm_z-0.2) ) = - Izz*Z_LF*dist_FF_cm_y + Izz*Z_RF*dist_FF_cm_y

(u_phi + GyrTerm)*(Ixx*Izz - Ixz**2) - Izz*(-Y_LS*LS_CoP-Y_RS*RS_CoP - Y_C*C_CoP - Y_Prop*(dist_ReF_cm_z-0.2) ) + Izz*D_b_LF*dist_FF_cm_y - Izz*D_b_RF*dist_FF_cm_y
= - Izz*L_LF*(-cosd(aa_LF))*dist_FF_cm_y + Izz*L_RF*(-cosd(aa_RF))*dist_FF_cm_y

stuff =  Izz*(CL_3D_LF)*q_LF*S_FF*cosd(aa_LF)*dist_FF_cm_y - Izz*(CL_3D_RF)*q_RF*S_FF*cosd(aa_RF)*dist_FF_cm_y

stuff = Izz*(cl_LF*(AR_FF/(AR_FF+2))*k2_LF)*q_LF*S_FF*cosd(aa_LF)*dist_FF_cm_y - Izz*(cl_RF*(AR_FF/(AR_FF+2))*k2_RF)*q_RF*S_FF*cosd(aa_RF)*dist_FF_cm_y

stuff1 = cl_LF*stuff2 + cl_RF*stuff3

#####################################

u_z = (Z*cosd(Pitch)*cosd(Roll)-X*sind(Pitch)+Y*cosd(Pitch)*sind(Roll))/m
u_z*m - Y*cosd(Pitch)*sind(Roll) = (L_b_LF[2] + D_b_LF[2] + L_b_RF[2] + D_b_RF[2] + Z_ReF + Z_LS + Z_RS + Z_C + W_b[2] + BUOY_b[2]+ D_b_H[2])*cosd(Pitch)*cosd(Roll) - (L_b_LF[0] + D_b_LF[0] + L_b_RF[0] + D_b_RF[0] + X_ReF + X_LS + X_RS + X_C + X_Prop + W_b[0] + BUOY_b[0] + D_b_H[0])*sind(Pitch)
u_z*m - Y*cosd(Pitch)*sind(Roll) - (D_b_LF[2] + D_b_RF[2] + Z_ReF + Z_LS + Z_RS + Z_C + W_b[2] + BUOY_b[2]+ D_b_H[2])*cosd(Pitch)*cosd(Roll) + (D_b_LF[0] + D_b_RF[0] + X_ReF + X_LS + X_RS + X_C + X_Prop + W_b[0] + BUOY_b[0] + D_b_H[0])*sind(Pitch) = (L_b_LF[2] + L_b_RF[2])*cosd(Pitch)*cosd(Roll) - (L_b_LF[0] + L_b_RF[0])*sind(Pitch)
stuff = -(L_LF*cosd(aa_LF) + L_RF*cosd(aa_RF))*cosd(Pitch)*cosd(Roll) - (L_LF*sind(aa_LF) + L_RF*sind(aa_LF))*sind(Pitch)
stuff = L_LF*(-cosd(aa_LF)*cosd(Pitch)*cosd(Roll) -sind(aa_LF)*sind(Pitch)) + L_RF*(-cosd(aa_RF)*cosd(Pitch)*cosd(Roll)- sind(aa_RF)*sind(Pitch))
stuff4 = stuff5*cl_LF + stuff6*cl_RF



## Real stuff:

cl_LF = (stuff4-stuff6*cl_RF)/stuff5
stuff1 = (stuff4-stuff6*cl_RF)/stuff5*stuff2 + cl_RF*stuff3
(=)
cl_RF(stuff3 - stuff6/stuff5*stuff2) = stuff1 - stuff4/stuff5*stuff2
cl_RF = (stuff1 - stuff4/stuff5*stuff2)/(stuff3 - stuff6/stuff5*stuff2)
cl_LF = (stuff4-stuff6*cl_RF)/stuff5

#No drag approximation:

u_z = (Z*cosd(Pitch)*cosd(Roll)-X*sind(Pitch)+Y*cosd(Pitch)*sind(Roll))/m
u_z*m - Y*cosd(Pitch)*sind(Roll) = (L_b_LF[2] + D_b_LF[2] + L_b_RF[2] + D_b_RF[2] + Z_ReF + Z_LS + Z_RS + Z_C + W_b[2] + BUOY_b[2]+ D_b_H[2])*cosd(Pitch)*cosd(Roll) - (L_b_LF[0] + D_b_LF[0] + L_b_RF[0] + D_b_RF[0] + X_ReF + X_LS + X_RS + X_C + X_Prop + W_b[0] + BUOY_b[0] + D_b_H[0])*sind(Pitch)
u_z*m - Y*cosd(Pitch)*sind(Roll) - (Z_ReF + Z_LS + Z_RS + Z_C + W_b[2] + BUOY_b[2]+ D_b_H[2])*cosd(Pitch)*cosd(Roll) + (X_ReF + X_LS + X_RS + X_C + X_Prop + W_b[0] + BUOY_b[0] + D_b_H[0])*sind(Pitch) = (L_b_LF[2] + L_b_RF[2] + D_b_LF[2] + D_b_RF[2])*cosd(Pitch)*cosd(Roll) - (L_b_LF[0] + L_b_RF[0] + D_b_LF[0] + D_b_RF[0])*sind(Pitch)
a1 = (L_b_LF[2] + L_b_RF[2] + D_b_LF[2] + D_b_RF[2])*cosd(Pitch)*cosd(Roll) - (L_b_LF[0] + L_b_RF[0] + D_b_LF[0] + D_b_RF[0])*sind(Pitch)
a1 = -cl_LF*AR_FF/(AR_FF+2)*k2_LF*q_LF*S_FF*cosd(aa_LF)*cosd(Pitch)*cosd(Roll) - cl_RF*AR_FF/(AR_FF+2)*k2_RF*q_RF*S_FF*cosd(aa_RF)*cosd(Pitch)*cosd(Roll) - cl_LF*AR_FF/(AR_FF+2)*k2_LF*q_LF*S_FF*sind(aa_LF)*sind(Pitch) - cl_RF*AR_FF/(AR_FF+2)*k2_LF*q_LF*S_FF*sind(aa_LF)*sind(Pitch) - D_LF*sind(aa_LF)*cosd(Pitch)*cosd(Roll) + D_LF*cosd(aa_LF)*sind(Pitch) - D_RF*sind(aa_RF)*cosd(Pitch)*cosd(Roll) + D_RF*cosd(aa_RF)*sind(Pitch)
a1 = (a2+a3)*cl_LF + (a4+a5)*cl_RF + (a6+a7)*D_LF + (a8+a9)*D_RF
a1 = a10*cl_LF + a11*cl_RF + a12*D_LF + a13*D_RF

##
D_LF = D_induced_LF+D_skf_LF+D_form_LF_left+D_interference_LF+Dw_LF
     = (CL_3D_LF)**2/(pi*e*(AR_FF))*q_LF*(S_FF) + D_skf_LF + (cdp_LF)*q_LF*(S_FF) + D_interference_LF + q_LF*0.5*((CL_3D_LF/(F_h_LS*2))**2)*h_LS*b_FF
     = a14 + ((AR_FF/(AR_FF+2)*k2_LF)**2)/(pi*e*AR_FF)*q_LF*S_FF*cl_LF**2 + cdp_LF*q_LF*S_FF + q_LF*0.5/((F_h_LS*2)**2)*h_LS*b_FF*((AR_FF/(AR_FF+2)*k2_LF)**2)*cl_LF**2
     = a14 + (a15+a16)*cl_LF**2 + a17*cdp_LF

D_RF = a18 + a19*cl_RF**2 + a20*cdp_RF
##

a1 = a10*cl_LF + a11*cl_RF + a12*(a14 + (a15+a16)*cl_LF**2 + a17*cdp_LF) + a13*(a18 + a19*cl_RF**2 + a20*cdp_RF)
a1 = a10*cl_LF + a11*cl_RF + a12*a14 + a12*(a15+a16)*cl_LF**2 + a12*a17*cdp_LF + a13*a18 + a13*a19*cl_RF**2 + a13*a20*cdp_RF
a1 - a12*a14 - a13*a18 = a10*cl_LF + a11*cl_RF + a12*(a15+a16)*cl_LF**2 + a12*a17*cdp_LF + a13*a19*cl_RF**2 + a13*a20*cdp_RF (1º equação)


###### constants to compute
a1 = u_z*m - Y*cosd(Pitch)*sind(Roll) - (Z_ReF + Z_LS + Z_RS + Z_C + W_b[2] + BUOY_b[2]+ D_b_H[2])*cosd(Pitch)*cosd(Roll) + (X_ReF + X_LS + X_RS + X_C + X_Prop + W_b[0] + BUOY_b[0] + D_b_H[0])*sind(Pitch)
a2 = -AR_FF/(AR_FF+2)*k2_LF*q_LF*S_FF*cosd(aa_LF)*cosd(Pitch)*cosd(Roll)
a3 = -AR_FF/(AR_FF+2)*k2_LF*q_LF*S_FF*sind(aa_LF)*sind(Pitch)
a4 = -AR_FF/(AR_FF+2)*k2_RF*q_RF*S_FF*cosd(aa_RF)*cosd(Pitch)*cosd(Roll)
a5 = -AR_FF/(AR_FF+2)*k2_RF*q_RF*S_FF*sind(aa_RF)*sind(Pitch)
a6 = -sind(aa_LF)*cosd(Pitch)*cosd(Roll)
a7 = cosd(aa_LF)*sind(Pitch)
a8 = -sind(aa_RF)*cosd(Pitch)*cosd(Roll)
a9 = cosd(aa_RF)*sind(Pitch)
a10 = a2+a3
a11 = a4+a5
a12 = a6+a7
a13 = a8+a9
a14 = D_skf_LF + D_interference_LF
a15 = ((AR_FF/(AR_FF+2)*k2_LF)**2)/(pi*e*AR_FF)*q_LF*S_FF
a16 = q_LF*0.5/((F_h_LS*2)**2)*h_LS*b_FF*((AR_FF/(AR_FF+2)*k2_LF)**2)
a17 = q_LF*(S_FF)
a18 = D_skf_RF + D_interference_RF
a19 = ((AR_FF/(AR_FF+2)*k2_RF)**2)/(pi*e*AR_FF)*q_RF*S_FF + q_RF*0.5/((F_h_RS*2)**2)*h_RS*b_FF*((AR_FF/(AR_FF+2)*k2_RF)**2)
a20 = q_RF*(S_FF)
###########################
L = -Y_LS*LS_CoP - Z_LF*dist_FF_cm_y-Y_RS*RS_CoP + Z_RF*dist_FF_cm_y-Y_C*C_CoP - Y_Prop*(dist_ReF_cm_z-0.2) 
N = X_LF*dist_FF_cm_y+X_LS*dist_FF_cm_y-X_RF*dist_FF_cm_y-X_RS*dist_FF_cm_y-Y_C*dist_C_cm_x-Y_Prop*dist_C_cm_x
(u_phi + GyrTerm)*(Ixx*Izz - Ixz**2) = Izz*(-Y_LS*LS_CoP-Y_RS*RS_CoP - Y_C*C_CoP - Y_Prop*(dist_ReF_cm_z-0.2) ) - Izz*Z_LF*dist_FF_cm_y + Izz*Z_RF*dist_FF_cm_y - Ixz*(X_LS*dist_FF_cm_y-X_RS*dist_FF_cm_y+Y_LS*dist_FF_cm_x+Y_RS*dist_FF_cm_x-Y_C*dist_C_cm_x-Y_Prop*dist_C_cm_x) - Ixz*X_LF*dist_FF_cm_y + Ixz*X_RF*dist_FF_cm_y
a21 = a22 - Izz*dist_FF_cm_y*Z_LF + Izz*Z_RF*dist_FF_cm_y + a23 - Ixz*X_LF*dist_FF_cm_y + Ixz*X_RF*dist_FF_cm_y
a21 = a22 - Izz*dist_FF_cm_y*(L_b_LF[2] + D_b_LF[2])+ Izz*dist_FF_cm_y*(L_b_RF[2] + D_b_RF[2]) +a23 - Ixz*dist_FF_cm_y*(L_b_LF[0] + D_b_LF[0]) + Ixz*dist_FF_cm_y*(L_b_RF[0] + D_b_RF[0])
a21-a22-a23 = (a24+a25)*cl_LF + (a26+a27)*cl_RF + Izz*dist_FF_cm_y*sind(aa_LF)*D_LF - Izz*dist_FF_cm_y*sind(aa_RF)*D_RF + Ixz*dist_FF_cm_y*cosd(aa_LF)*D_LF - Ixz*dist_FF_cm_y*cosd(aa_RF)*D_RF
a28 = a29*cl_LF + a30*cl_RF + (a31+a32)*D_LF + (a33+a34)*D_RF
a28 = a29*cl_LF + a30*cl_RF + a35*D_LF + a36*D_RF
a28 = a29*cl_LF + a30*cl_RF + a35*(a14 + (a15+a16)*cl_LF**2 + a17*cdp_LF) + a36*(a18 + a19*cl_RF**2 + a20*cdp_RF)
a28 = a29*cl_LF + a30*cl_RF + a35*a14 + a35*(a15+a16)*cl_LF**2 + a36*a19*cl_RF**2 + a36*a18 + a35*a17*cdp_LF + a36*a20*cdp_RF
a28 - a35*a14 -a36*a18 = a29*cl_LF + a30*cl_RF + a35*(a15+a16)*cl_LF**2 + a36*a19*cl_RF**2 + a35*a17*cdp_LF + a36*a20*cdp_RF 2º equação

#constants
a21 = (u_phi + GyrTerm)*(Ixx*Izz - Ixz**2)
a22 = Izz*(-Y_LS*LS_CoP-Y_RS*RS_CoP - Y_C*C_CoP - Y_Prop*(dist_ReF_cm_z-0.2) )
a23 = - Ixz*(X_LS*dist_FF_cm_y-X_RS*dist_FF_cm_y+Y_LS*dist_FF_cm_x+Y_RS*dist_FF_cm_x-Y_C*dist_C_cm_x-Y_Prop*dist_C_cm_x)
a24 = Izz*dist_FF_cm_y*AR_FF/(AR_FF+2)*k2_LF*q_LF*S_FF*cosd(aa_LF)
a25 =  -Ixz*dist_FF_cm_y*AR_FF/(AR_FF+2)*k2_LF*q_LF*S_FF*sind(aa_LF)
a26 = -Izz*dist_FF_cm_y*AR_FF/(AR_FF+2)*k2_RF*q_RF*S_FF*cosd(aa_RF)
a27 = Ixz*dist_FF_cm_y*AR_FF/(AR_FF+2)*k2_F*q_RF*S_FF*sind(aa_RF)
a28 = a21 - a22 - a23
a29 = a24 + a25
a30 = a26 + a27
a31 = Izz*dist_FF_cm_y*sind(aa_LF)
a32 = Ixz*dist_FF_cm_y*cosd(aa_LF)
a33 = -Izz*dist_FF_cm_y*sind(aa_RF)
a34 = -Ixz*dist_FF_cm_y*cosd(aa_RF)
a35 = a31 + a32
a36 = a33 + a34
a37 = a28 - a35*a14 -a36*a18
a38 = a35*(a15+a16)
a39 = a36*a19
a40 = a35*a17
a41 = a36*a20
a42 = a1 - a12*a14 - a13*a18
a43 = a12*(a15+a16)
a44 = a12*a17
a45 = a13*a19
a46 = a13*a20

## SOLVING
a37 = a29*cl_LF + a30*cl_RF + a38*cl_LF**2 + a39*cl_RF**2 + a40*cdp_LF + a41*cdp_RF
a42 = a10*cl_LF + a11*cl_RF + a43*cl_LF**2 + a45*cl_RF**2 + a44*cdp_LF + a46*cdp_RF

## IN ORDER OF AA_LF and AA_RF
cl_LF = cl1*aa_LF + cl0
cl_RF = cl1*aa_RF + cl0
cdp_LF = cd1*aa_LF + cd0
cdp_RF = cl1*aa_RF + cd0

cl_LF² = cl1²*aa_LF² + 2*cl0*aa_LF + cl0²
cl_RF² = cl1²*aa_RF² + 2*cl0*aa_RF + cl0²

## SUBSTITUTING

a29*cl_LF = a29*cl1*aa_LF + a29*cl0
a30*cl_RF = a30*cl1*aa_RF + a30*cl0
a38*cl_LF**2 = a38*cl1²*aa_LF² + a38*2*cl0*aa_LF + a38*cl0²
a39*cl_RF**2 = a39*cl1²*aa_RF² + a39*2*cl0*aa_RF + a39*cl0²
a40*cdp_LF = a40*cd1*aa_LF + a40*cd0
a41*cdp_RF = a41*cd1*aa_RF + a41*cd0

contribuition aa_LF: a29*cl1*aa_LF + a29*cl0 + a38*cl1²*aa_LF² + a38*2*cl0*aa_LF + a38*cl0² + a40*cd1*aa_LF + a40*cd0
contribuition aa_RF: a30*cl1*aa_RF + a30*cl0 + a39*cl1²*aa_RF² + a39*2*cl0*aa_RF + a39*cl0² + a41*cd1*aa_RF + a41*cd0

contribution aa_LF: a38*cl1²*aa_LF² + (a29*cl1 + 2*a38*cl0 + a40*cd1)*aa_LF + a29*cl0 + a40*cd0 + a38*cl0²
contribution aa_RF: a39*cl1²*aa_RF² + (a30*cl1 + 2*a39*cl0 + a41*cd1)*aa_RF + a30*cl0 + a41*cd0 + a39*cl0²

a47 = a38*cl1²
a48 = a29*cl1 + 2*a38*cl0 + a40*cd1
a49 = a29*cl0 + a40*cd0 + a38*cl0²
a50 = a39*cl1²
a51 = a30*cl1 + 2*a39*cl0 + a41*cd1
a52 = a30*cl0 + a41*cd0 + a39*cl0²

## second equation 29 -> 10 , 30 -> 11 , 38 -> 43 , 39 -> 45 , 40 -> 44 , 41 -> 46
a37 = a29*cl_LF + a30*cl_RF + a38*cl_LF**2 + a39*cl_RF**2 + a40*cdp_LF + a41*cdp_RF
a42 = a10*cl_LF + a11*cl_RF + a43*cl_LF**2 + a45*cl_RF**2 + a44*cdp_LF + a46*cdp_RF

a53 = a43*cl1²
a54 = a10*cl1 + 2*a43*cl0 + a44*cd1
a55 = a10*cl0 + a44*cd0 + a43*cl0²
a56 = a45*cl1²
a57 = a11*cl1 + 2*a45*cl0 + a46*cd1
a58 = a11*cl0 + a46*cd0 + a45*cl0²

## final system
a37 = a47*aa_LF² + a48*aa_LF + a49 + a50*aa_RF² + a51*aa_RF + a52
a42 = a53*aa_LF² + a54*aa_LF + a55 + a56*aa_RF² + a57*aa_RF + a58

### Cubic
cl_LF = cl3*aa_LF**3 + cl2*aa_LF**2 + cl1*aa_LF + cl0
cl_RF = cl3*aa_RF**3 + cl2*aa_RF**2 + cl1*aa_RF + cl0
cdp_LF = cd3*aa_LF**3 + cd2*aa_LF**2 + cd1*aa_LF + cd0
cdp_RF = cd3*aa_LF**3 + cd2*aa_LF**2 + cd1*aa_LF + cd0

cl_LF² = (cl3*aa_LF**3 + cl2*aa_LF**2 + cl1*aa_LF + cl0)(cl3*aa_LF**3 + cl2*aa_LF**2 + cl1*aa_LF + cl0)
cl_RF² = (cl3*aa_RF**3 + cl2*aa_RF**2 + cl1*aa_RF + cl0)(cl3*aa_RF**3 + cl2*aa_RF**2 + cl1*aa_RF + cl0)

     = (cl3**2*aa_LF**6 + cl3*cl2*aa_LF**5 + cl3*cl1*aa_LF**4 + cl3*cl0*aa_LF**3
      + cl2*cl3*aa_LF**5 + cl2**2*aa_LF**4 + cl2*cl1*aa_LF**3 + cl2*cl0*aa_LF**2
      + cl1*cl3*aa_LF**4 + cl1*cl2*aa_LF**3 + cl1**2*aa_LF**2 + cl1*cl0*aa_LF
      + cl0*cl3*aa_LF**3 + cl0*cl2*aa_LF**2 + cl0*cl1*aa_LF + cl0**2)
     
cl_LF**2    = cl3**2*aa_LF**6 + 2*cl2*cl3*aa_LF**5 + (cl2**2 + 2*cl1*cl3)*aa_LF**4 + (cl3*cl0 + 2*cl2*cl1 + cl0*cl3)*aa_LF**3 + (2*cl2*cl0 + cl1**2)*aa_LF**2 + 2*cl1*cl0*aa_LF + cl0**2
cl_RF**2    = cl3**2*aa_RF**6 + 2*cl2*cl3*aa_RF**5 + (cl2**2 + 2*cl1*cl3)*aa_RF**4 + (cl3*cl0 + 2*cl2*cl1 + cl0*cl3)*aa_RF**3 + (2*cl2*cl0 + cl1**2)*aa_RF**2 + 2*cl1*cl0*aa_RF + cl0**2

cl4 = cl0**2
cl5 = 2*cl1*cl0
cl6 = 2*cl2*cl0 + cl1**2
cl7 = cl3*cl0 + 2*cl2*cl1 + cl0*cl3
cl8 = cl2**2 + 2*cl1*cl3
cl9 = 2*cl2*cl3
cl10 = cl3**2

cl_LF**2 = cl10*aa_LF**6 + cl9*aa_LF**5 + cl8*aa_LF**4 + cl7*aa_LF**3 + cl6*aa_LF**2 + cl5*aa_LF + cl4
cl_RF**2 = cl10*aa_RF**6 + cl9*aa_RF**5 + cl8*aa_RF**4 + cl7*aa_RF**3 + cl6*aa_RF**2 + cl5*aa_RF + cl4


a37 = a29*cl_LF + a30*cl_RF + a38*cl_LF**2 + a39*cl_RF**2 + a40*cdp_LF + a41*cdp_RF

 contribution LF:
 a29*cl_LF = a29*cl3*aa_LF**3 + a29*cl2*aa_LF**2 + a29*cl1*aa_LF + a29*cl0
 a38*cl_LF**2 = a38*cl10*aa_LF**6 + a38*cl9*aa_LF**5 + a38*cl8*aa_LF**4 + a38*cl7*aa_LF**3 + a38*cl6*aa_LF**2 + a38*cl5*aa_LF + a38*cl4
 a40*cdp_LF = a40*cd3*aa_LF**3 + a40*cd2*aa_LF**2 + a40*cd1*aa_LF + a40*cd0

 a47 = a38*cl10
 a48 = a38*cl9
 a49 = a38*cl8
 a50 = a29*cl3 + a38*cl7 + a40*cd3
 a51 = a29*cl2 + a38*cl6 + a40*cd2
 a52 = a29*cl1 + a38*cl5 + a40*cd1
 a53 = a29*cl0 + a38*cl4 + a40*cd0

contribution RF:
 a30*cl_RF = a30*cl3*aa_RF**3 + a30*cl2*aa_RF**2 + a30*cl1*aa_RF + a30*cl0
 a39*cl_RF**2 = a39*cl10*aa_RF**6 + a39*cl9*aa_RF**5 + a39*cl8*aa_RF**4 + a39*cl7*aa_RF**3 + a39*cl6*aa_RF**2 + a39*cl5*aa_RF + a39*cl4
 a41*cdp_RF = a41*cd3*aa_RF**3 + a41*cd2*aa_RF**2 + a41*cd1*aa_RF + a41*cd0

 a54 = a39*cl10
 a55 = a39*cl9
 a56 = a39*cl8
 a57 = a30*cl3 + a39*cl7 + a41*cd3
 a58 = a30*cl2 + a39*cl6 + a41*cd2
 a59 = a30*cl1 + a39*cl5 + a41*cd1
 a60 = a30*cl0 + a39*cl4 + a41*cd0

 a37 = a47*aa_LF**6 + a48*aa_LF**5 + a49*aa_LF**4 + a50*aa_LF**3 + a51*aa_LF**2 + a52*aa_LF + a53 + a54*aa_RF**6 + a55*aa_RF**5 + a56*aa_RF**4 + a57*aa_RF**3 + a58*aa_RF**2 + a59*aa_RF + a60

 a61 = a53 + a60

a37 = a47*aa_LF**6 + a48*aa_LF**5 + a49*aa_LF**4 + a50*aa_LF**3 + a51*aa_LF**2 + a52*aa_LF + a54*aa_RF**6 + a55*aa_RF**5 + a56*aa_RF**4 + a57*aa_RF**3 + a58*aa_RF**2 + a59*aa_RF + a61
 
## second equation 29 -> 10 , 30 -> 11 , 38 -> 43 , 39 -> 45 , 40 -> 44 , 41 -> 46
a62 = a43*cl10
a63 = a43*cl9
a64 = a43*cl8
a65 = a10*cl3 + a43*cl7 + a44*cd3
a66 = a10*cl2 + a43*cl6 + a44*cd2
a67 = a10*cl1 + a43*cl5 + a44*cd1
a68 = a10*cl0 + a43*cl4 + a44*cd0

a69 = a45*cl10
a70 = a45*cl9
a71 = a45*cl8
a72 = a11*cl3 + a45*cl7 + a46*cd3
a73 = a11*cl2 + a45*cl6 + a46*cd2
a74 = a11*cl1 + a45*cl5 + a46*cd1
a75 = a11*cl0 + a45*cl4 + a46*cd0

a76 = a68 + a75

a42 = a62*aa_LF**6 + a63*aa_LF**5 + a64*aa_LF**4 + a65*aa_LF**3 + a66*aa_LF**2 + a67*aa_LF + a69*aa_RF**6 + a70*aa_RF**5 + a71*aa_RF**4 + a72*aa_RF**3 + a73*aa_RF**2 + a74*aa_RF + a76

a77 = a61 - a37
a78 = a76 - a42

## final

f1 = a47*aa_LF**6 + a48*aa_LF**5 + a49*aa_LF**4 + a50*aa_LF**3 + a51*aa_LF**2 + a52*aa_LF + a54*aa_RF**6 + a55*aa_RF**5 + a56*aa_RF**4 + a57*aa_RF**3 + a58*aa_RF**2 + a59*aa_RF + a77 = 0
f2 = a62*aa_LF**6 + a63*aa_LF**5 + a64*aa_LF**4 + a65*aa_LF**3 + a66*aa_LF**2 + a67*aa_LF + a69*aa_RF**6 + a70*aa_RF**5 + a71*aa_RF**4 + a72*aa_RF**3 + a73*aa_RF**2 + a74*aa_RF + a78 = 0

grad(f1,aa_LF) = 6*a47*aa_LF**5 + 5*a48*aa_LF**4 + 4*a49*aa_LF**3 + 3*a50*aa_LF**2 + 2*a51*aa_LF + a52 
grad(f1,aa_RF) = 6*a54*aa_RF**5 + 5*a55*aa_RF**4 + 4*a56*aa_RF**3 + 3*a57*aa_RF**2 + 2*a58*aa_RF + a59
grad(f2,aa_LF) = 6*a62*aa_LF**5 + 5*a63*aa_LF**4 + 4*a64*aa_LF**3 + 3*a65*aa_LF**2 + 2*a66*aa_LF + a67 
grad(f2,aa_RF) = 6*a69*aa_RF**5 + 5*a70*aa_RF**4 + 4*a71*aa_RF**3 + 3*a72*aa_RF**2 + 2*a73*aa_RF + a74

hess(f1, aa_LF**2) = 30*a47*aa_LF**4 + 20*a48*aa_LF**3 + 12*a49*aa_LF**2 + 6*a50*aa_LF + 2*a51 
hess(f1, aa_LF*aa_RF) = 0
hess(f1, aa_RF**2) = 30*a47*aa_LF**4 + 20*a48*aa_LF**3 + 12*a49*aa_LF**2 + 6*a50*aa_LF + 2*a51 

##Square to have one function

f = (a*x**6 + b*x**5 + c*x**4 + d*x**3 + e*x**2 + f*x + g)*(a*x**6 + b*x**5 + c*x**4 + d*x**3 + e*x**2 + f*x + g)
= a²x**12 + abx**11 + acx**10 + adx**9 + aex**8 + afx**7 + agx**6
+ bax**11 + b²x**10 + bcx**9 + bdx**8 + bex**7 + bfx**6 + bgx**5
+ cax**10 + cbx**9 + c²x**8 + cdx**7 + cex**6 + cfx**5 + cgx**4
+ da**9 + dbx**8 + dcx**7 + d²x**6 + dex**5 + dfx**4 + dgx**3
+ ea**8 + ebx**7 + ecx**6 + edx**5 + eex**4 + efx**3 + egx**2
+ fa**7 + fbx**6 + fcx**5 + fdx**4 + fex**3 + f²x**2 + fgx
+ ga**6 + gbx**5 + gcx**4 + gdx**3 + gex**2 + gfx + g²

## SOLVING splitting equations
a37 = a29*cl_LF + a30*cl_RF + a38*cl_LF**2 + a39*cl_RF**2 + a40*cdp_LF + a41*cdp_RF
a42 = a10*cl_LF + a11*cl_RF + a43*cl_LF**2 + a45*cl_RF**2 + a44*cdp_LF + a46*cdp_RF

a37 = (a29+a30)*cl_com + (a38+a39)*cl_sqr_com + (a40+a41)*cdp_com

-> alfa_com
#Aprox (cl_LF + cl_RF)/2 = cl_com
#Aprox (CL_sqr_LF + cl_sqr_RF)/2 = cl_sqr_com

a42 = a10*cl_LF + a11*(2*cl_com-cl_LF) + a43*cl_sqr_LF + a45(2*cl_sqr_com - cl_sqr_LF) + a44*cdp_LF + a46*(2*cdp_com - cdp_LF)
a42 = (a10-a11)*cl_LF + (a43-a45)cl_sqr_LF + (a44-a46)cdp_LF + 2*a11*cl_com + 2*a45*cl_sqr_com + 2*a46*cdp_com

a47 = a10-a11
a48 = a43-a45
a49 = a44-a46
a50 = 2*a11*cl_com + 2*a45*cl_sqr_com + 2*a46*cdp_com

a42 = a47*cl_LF + a48*cl_sqr_LF + a49*cdp_LF + a50

cl_LF = cl3*aa_LF**3 + cl2*aa_LF**2 + cl1*aa_LF + cl0
cl_sqr_LF = clsq3*aa_LF**3 + clsq2*aa_LF**2 + clsq1*aa_LF + clsq0
cdp_LF = cd3*aa_LF**3 + cd2*aa_LF**2 + cd1*aa_LF + cd0

a42 = (cl3*a47 + clsq3*a48 + cd3*a49)*aa_LF**3 + (cl2*a47 + clsq2*a48 + cd2*a49)*aa_LF**2 + (cl1*a47 + clsq1*a48 + cd1*a49)*aa_LF + (cl0*a47 + clsq0*a48 + cd0*a49) + a50

a51 = (cl3*a47 + clsq3*a48 + cd3*a49)
a52 = (cl2*a47 + clsq2*a48 + cd2*a49)
a53 = (cl1*a47 + clsq1*a48 + cd1*a49)
a54 = (cl0*a47 + clsq0*a48 + cd0*a49) + a50

a42 = a51*aa_LF**3 + a52*aa_LF**2 + a53*aa_LF + a54

## com case

a37 = (a29+a30)*cl_com + (a38+a39)*cl_sqr_com + (a40+a41)*cdp_com

a55 = a29+a30
a56 = a38+a39
a57 = a40+a41

a37 = a55*cl_com + a56*cl_sqr_com + a57*cdp_com

a58 = (cl3*a55 + clsq3*a56 + cd3*a57)
a59 = (cl2*a55 + clsq2*a56 + cd2*a57)
a60 = (cl1*a55 + clsq1*a56 + cd1*a57)
a61 = (cl0*a55 + clsq0*a56 + cd0*a57) 

a37 = a58*aa_com**3 + a59*aa_com**2 + a60*aa_com + a61

## Final functions
a62 = a54 - a42
a63 = a61 - a37

cl_RF = 2*cl_com-cl_LF



a58*aa_com**3 + a59*aa_com**2 + a60*aa_com + a63 = 0
a51*aa_LF**3 + a52*aa_LF**2 + a53*aa_LF + a62 = 0

a64 = cl0 - 2*cl_com + cl_LF

cl3*aa_RF**3 + cl2*aa_RF**2 + cl1*aa_RF + a64 = 0

## Find minimum

aa_com = get_minimum(a58, a59, a60, a63)

cl_com = cl(aa_com)
cd_com = cd(aa_com)
clsqr_com = clsqr_com(aa_com)
a47 = a10-a11
a48 = a43-a45
a49 = a44-a46
a50 = 2*a11*cl_com + 2*a45*cl_sqr_com + 2*a46*cdp_com
a51 = (cl3*a47 + clsq3*a48 + cd3*a49)
a52 = (cl2*a47 + clsq2*a48 + cd2*a49)
a53 = (cl1*a47 + clsq1*a48 + cd1*a49)
a54 = (cl0*a47 + clsq0*a48 + cd0*a49) + a50
a62 = a54 - a42
aa_LF = get_minimum(a51, a52, a53, a62)

cl_LF = cl(aa_LF)
a64 = cl0 - 2*cl_com + cl_LF
aa_RF = get_minimum(cl3, cl2, cl1, a64)




def get_minimum(c3,c2,c1,c0):



     #p = abs(c3*x**3 + c2*x**2 + c1*x + c0)
     #if p > 0: p = c3*x**3 + c2*x**2 + c1*x + c0 else p = -(c3*x**3 + c2*x**2 + c1*x + c0)

     a = 3*c3
     b = 2*c2
     c = c1

     x1 = (-b + mt.sqrt(b**2 - 4*a*c))/(2*a)
     x2 = (-b - mt.sqrt(b**2 - 4*a*c))/(2*a)

     sol = x1

     p1 = abs(c3*x1**3 + c2*x1**2 + c1*x1 + c0)
     p2 = abs(c3*x2**3 + c2*x2**2 + c1*x2 + c0)

     if(p2 < p1):
          sol = x2

     return sol







