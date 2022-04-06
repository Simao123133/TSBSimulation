import math as mt
import numpy as np
from scipy.optimize import minimize, root, least_squares


def get_aa(Y, params, u_z, u_phi):

    a37, a29, a30, a38, a39, a40, a41, a42, a10, a11, a43, a45, a44, a46, Re_LF, Re_RF, d_aa_LF, d_aa_RF = constants(Y, params, u_z, u_phi)

    cl3, cl2, cl1, cl0 = cl_constants((Re_LF+Re_RF)/2)
    clsq3, clsq2, clsq1, clsq0 = clsq_constants((Re_LF+Re_RF)/2)
    cd3, cd2, cd1, cd0 = cd_constants((Re_LF+Re_RF)/2)

    ## Find minimum
    a55 = a10+a11
    a56 = a43+a45
    a57 = a44+a46
    a58 = (cl3*a55 + clsq3*a56 + cd3*a57)
    a59 = (cl2*a55 + clsq2*a56 + cd2*a57)
    a60 = (cl1*a55 + clsq1*a56 + cd1*a57)
    a61 = (cl0*a55 + clsq0*a56 + cd0*a57) 
    a63 = a61 - a42
    
    aa_com = get_minimum(a58, a59, a60, a63)

    cl_com = p3(aa_com, cl3, cl2, cl1, cl0)
    cd_com = p3(aa_com, cd3, cd2, cd1, cd0)
    clsqr_com = p3(aa_com, clsq3, clsq2, clsq1, clsq0)
    a47 = a29-a30 
    a48 = a38-a39 
    a49 = a40-a41
    a50 = 2*a30*cl_com + 2*a39*clsqr_com + 2*a41*cd_com
    a51 = (cl3*a47 + clsq3*a48 + cd3*a49)
    a52 = (cl2*a47 + clsq2*a48 + cd2*a49)
    a53 = (cl1*a47 + clsq1*a48 + cd1*a49)
    a54 = (cl0*a47 + clsq0*a48 + cd0*a49) + a50
    a62 = a54 - a37
    aa_LF = get_minimum(a51, a52, a53, a62)

    cl_LF = p3(aa_LF, cl3, cl2, cl1, cl0)
    a64 = cl0 - 2*cl_com + cl_LF
    aa_RF = get_minimum(cl3, cl2, cl1, a64)

    act_LF = aa_LF - d_aa_LF
    act_RF = aa_RF - d_aa_RF

    return act_LF, act_RF

def constants(states, parameters, u_z, u_phi):

    def cosd(ang):
        return np.cos(ang*np.pi/180)
    def sind(ang):
        return np.sin(ang*np.pi/180)
    def tand(ang):
        return np.tan(ang*np.pi/180)
    def asind(ang):
        return 180/np.pi*np.arcsin(ang)
    def atand(ang):
        return 180/np.pi*np.arctan(ang)    

    #FULL MODEL 

    #states
    u=states[0]
    w=states[1]
    qq=states[2]
    Pitch=states[3]
    z_cm=states[4]
    v=states[5]
    p=states[6]
    r=states[7]
    Roll=states[8]
    Yaw=states[9]

    #inputs
    Act_rear = parameters[0]
    Rudder = parameters[1]
    Rpms_motor = parameters[2]

    #Constants
    pi=3.14159265358979
    g = 9.81 #gravity (m/s**2)
    rho= 1025.00 # Water density(kg/m**3)   <-I
    mu= 1.14E-3 # Coef. Viscosidade Dinâmica para calular Reynolds(Pa*s)   <-I
    m = 180 #Weight of boat (Kg) <-I
    dist_cm_Us_x = 3.150   #distance between cm and US
    dist_S_Us_x = 1.850
    dist_S_cm_x = 1.3
    dist_C_Us_x = 5.520
    dist_FF_cm_y = 1.2/2 #front foils' y distance divided by 2
    dist_C_cm_x = 2.367 #distance from the column to cm in the x axis
    dist_FF_cm_x = 1.30156 #distance from the front foils to cm in the x axis
    dist_ReF_cm_z = 1.05 #distance from the rear foil to cm in the z axis
    dist_FF_cm_z = 0.85 #distance from the front foils to cm in the z axis
    Exposed_C = 0.776 #max column exposed size (the rest is covered)

    # Height at each point (z is from US)
    z = (z_cm - tand(Pitch)*dist_cm_Us_x)/cosd(Pitch)        
    z_S_no_roll = z+dist_S_Us_x*tand(Pitch) #Perpendicular distance to water of the struts and column
    z_C_no_roll = z+dist_C_Us_x*tand(Pitch) 
    z_LS = z_S_no_roll - tand(Roll)*dist_FF_cm_y
    z_RS = z_S_no_roll + tand(Roll)*dist_FF_cm_y

    # Depth of struts and column
    h_C = dist_ReF_cm_z + z_C_no_roll
    h_LS = dist_FF_cm_z + z_LS
    h_RS = dist_FF_cm_z + z_RS

    if h_C > Exposed_C:
        h_C = Exposed_C

    if h_LS > dist_FF_cm_z:
        h_LS = dist_FF_cm_z

    if h_RS > dist_FF_cm_z:
        h_RS = dist_FF_cm_z


    Iyy = 503        # boat's moment of inertia in kg.m**2 for pitch
    Ixx = 31         # moment of inertia for roll
    Izz = 300
    Ixz = 16.5

    #strut and column lift center depende on heave
    C_CoP = h_C/2 - z_C_no_roll #column center of pressure
    LS_CoP = h_LS/2 - z_LS
    RS_CoP = h_RS/2 - z_RS

    #Inertial angular speeds
    Yaw_dot = (qq*sind(Roll)/(cosd(Roll)*cosd(Roll)*cosd(Pitch))+r/(cosd(Roll)*cosd(Pitch)))/(1+sind(Roll)*sind(Roll)/(cosd(Roll)*cosd(Roll)))
    Pitch_dot = qq/cosd(Roll) - Yaw_dot*cosd(Pitch)*sind(Roll)/cosd(Roll)
    Roll_dot = p + Yaw_dot*sind(Pitch)
    W = pi/180*np.array([Roll_dot,Pitch_dot,Yaw_dot])

    #speed at each point - rotational kinematics Vp = Vo + WxRp
    #Left strut
    R_LS = np.array([dist_S_cm_x,-dist_FF_cm_y,LS_CoP])
    V_LS = np.array([u,v,w]) + np.cross(W,R_LS)
    #Right Strut
    R_RS = np.array([dist_S_cm_x,dist_FF_cm_y,RS_CoP])
    V_RS = np.array([u,v,w]) + np.cross(W,R_RS)
    #Column
    R_C = np.array([-dist_C_cm_x,0,C_CoP])
    V_C = np.array([u,v,w]) + np.cross(W,R_C)
    #Left foil
    R_LF = np.array([dist_S_cm_x,-dist_FF_cm_y,dist_FF_cm_z])
    V_LF = np.array([u,v,w]) + np.cross(W,R_LF)
    #Right foil
    R_RF = np.array([dist_S_cm_x,dist_FF_cm_y,dist_FF_cm_z])
    V_RF = np.array([u,v,w]) + np.cross(W,R_RF)
    #Rear foil
    R_RE = np.array([-dist_C_cm_x,0,dist_ReF_cm_z])
    V_ReF = np.array([u,v,w]) + np.cross(W,R_RE)

    #angle of attack at each point
    bb_C = asind(V_C[1]/V_C[0])
    bb_LS = atand(V_LS[1]/V_LS[0])
    bb_RS = atand(V_RS[1]/V_RS[0])
    aa_LF = atand(V_LF[2]/V_LF[0])
    aa_RF = atand(V_RF[2]/V_RF[0])
    aa_ReF = atand(V_ReF[2]/V_ReF[0])
    eff_aa_ReF = Act_rear + aa_ReF

    #Speed vector's magnitude at each point
    speed_cm = np.linalg.norm(np.array([u,v,w]))
    speed_LF = np.linalg.norm(V_LF)
    speed_RF = np.linalg.norm(V_RF)
    speed_ReF = np.linalg.norm(V_ReF)
    speed_C = np.linalg.norm(V_C)
    speed_LS = np.linalg.norm(V_LS)
    speed_RS = np.linalg.norm(V_RS)

    #Dynamic pressure at each point
    q_cm = 0.5*rho*speed_cm**2#Dynamic Pressure (Kg m**-1 s**-2)
    q_LF = 0.5*rho*speed_LF**2
    q_RF = 0.5*rho*speed_RF**2
    q_ReF = 0.5*rho*speed_ReF**2
    q_C = 0.5*rho*speed_C**2
    q_LS = 0.5*rho*speed_LS**2
    q_RS = 0.5*rho*speed_RS**2

    #Geometria foil frente(x2)
    c_FF=0.08 # Corda foil frente(m)
    t_FF=0.12 # Espessura maxima perfil foil frente (n63412: t=0.12c)
    b_FF=0.82 # Envergadura foil frente (m)
    S_FF=pi*b_FF*c_FF/4 #Area planform foil frente (m**2) --  0.04773 valor dos foils atuais
    SA_FF=S_FF*2 # Area superficie foil frente (assumir 2x planform area) (m**2)
    AR_FF=b_FF**2/S_FF #Aspect ratio foil frente

    def cordfun(x):
        c_FF=0.08 # Corda foil frente(m)
        b_FF=0.82 # Envergadura foil frente (m)
        return c_FF**2 *(x-(x**3/(b_FF/2)**2)/3)

    MAC_FF =  2/S_FF*(cordfun(b_FF/2)-cordfun(0))
    c_FF=MAC_FF #Chord according to MAC
    k2_LF = (16*(h_LS/MAC_FF)**2+2)/(16*(h_LS/MAC_FF)**2+1)
    k2_RF = (16*(h_RS/MAC_FF)**2+2)/(16*(h_RS/MAC_FF)**2+1)
    T_FF=c_FF*t_FF # Thickness front foil= chord*t

    #Geometria foil trás
    c_ReF=0.08 # Corda foil trás (m)
    t_ReF=0.12 # Espessura maxima foil tras (n63412: t=0.12c)
    b_ReF=0.82 # Envergadura foil tras (m)
    S_ReF=pi*b_ReF*c_ReF/4 #Area planform foil tras (m**2) --- 0.04584 valor dos foils atuais
    SA_ReF=S_ReF*2 # Area superficie foil tras (assumir 2x planform area) (m**2)
    AR_ReF=b_ReF**2/S_ReF #Aspect ratio foil tras

    def cordfun1(x):
        c_ReF=0.08 # Corda foil trás (m)
        b_ReF=0.82 # Envergadura foil tras (m)
        return c_ReF**2 *(x-(x**3/(b_ReF/2)**2)/3)

    MAC_ReF = 2/S_ReF*(cordfun1(b_ReF/2)-cordfun1(0))
    c_ReF=MAC_ReF #Chord according to MAC
    k2_ReF = (16*(h_C/MAC_ReF)**2+2)/(16*(h_C/MAC_FF)**2+1)
    T_ReF=c_ReF*t_ReF# Thickness rear foil= chord*t

    #Struts geometry
    c_S= 0.13 # Corda strut (m)
    t_S=0.12 # Espessura maxima perfil strut (n0012: t=0.12c)
    T_S=c_S*t_S # Thickness strut= chord*t
    #Left Strut
    S_LS=c_S*h_LS #area planform strut (m**2)
    SA_LS=2*S_LS # Area superficie strut (SOLIDWORKS) (m**2)
    sub_LS=h_LS/c_S #aspect ratio da seccao submergida strut
    #Right Strut
    S_RS=c_S*h_RS #area planform strut (m**2)
    SA_RS=2*S_RS # Area superficie strut (SOLIDWORKS) (m**2)
    sub_RS=h_RS/c_S #aspect ratio da seccao submergida strut

    #Column geometry
    c_C=0.17 # corda coluna (m)
    t_C=0.18 # Espessura maxima perfil coluna (n0012: t=0.12c)
    S_C=c_C*h_C #area planform coluna (m**2)
    SA_C=2*S_C # Area superficie coluna (m**2)
    sub_C=h_C/c_C #aspect ratio da seccao submergida coluna
    T_C=c_C*t_C# Thickness coluna= chord*t

    #Reynolds numbers
    Re_LF=(rho)*(speed_LF)*(c_FF)/(mu) #Num Reynolds foil frente
    Re_RF=(rho)*(speed_RF)*(c_FF)/(mu) #Num Reynolds foil frente
    Re_ReF=(rho)*(speed_ReF)*(c_ReF)/(mu) #Num Reynold foil tras
    Re_LS=(rho)*(speed_LS)*(c_S)/(mu) #Num Reynolds strut
    Re_RS=(rho)*(speed_RS)*(c_S)/(mu) #Num Reynolds strut
    Re_C=(rho)*(speed_C)*(c_C)/(mu) #Num Reynolds coluna

    # Determining cd of the foils, struts and prop column

    def cdp_S(x,y):
        #strut pressure drag. Column assumed to have same pressure drag as
        #strut
        p00 =     0.01162
        p10 =  -7.938e-07
        p01 =  -5.664e-05
        p20 =   0.0001378
        p11 =   3.735e-08
        p02 =   1.241e-07
        p30 =  -1.257e-07
        p21 =   -4.13e-07
        p12 =  -9.269e-11
        p03 =  -1.231e-10
        p40 =   1.604e-06
        p31 =    4.88e-10
        p22 =   4.056e-10
        p13 =   6.432e-14
        p04 =   5.515e-14
        p50 =  -1.201e-09
        p41 =  -1.096e-09
        p32 =   -1.95e-13
        p23 =   -1.04e-13
        p14 =   -1.36e-17
        p05 =   -9.13e-18

        return (p00 + p10*x + p01*y + p20*x**2 + p11*x*y + p02*y**2 + p30*x**3 + p21*x**2*y
        + p12*x*y**2 + p03*y**3 + p40*x**4 + p31*x**3*y + p22*x**2*y**2 
        + p13*x*y**3 + p04*y**4 + p50*x**5 + p41*x**4*y + p32*x**3*y**2 
        + p23*x**2*y**3 + p14*x*y**4 + p05*y**5)

    CDp_LS = cdp_S(bb_LS,Re_LS/1000)
    CDp_RS = cdp_S(bb_RS,Re_RS/1000)
    CDp_C = cdp_S(bb_C+Rudder,Re_C/1000)

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

    cl_ReF = cl(eff_aa_ReF,Re_ReF/1000)
    CL_3D_ReF =cl_ReF*(AR_ReF/(AR_ReF+2))*k2_ReF

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

    cdp_ReF = cdpf(eff_aa_ReF,Re_ReF/1000)


    #Lift (XFRL5)
    L_ReF=(CL_3D_ReF)*q_ReF*S_ReF

    k_LS=1/(13+(13/(sub_LS)**2)+(26/(sub_LS))) # coeficiente de lift que toma em conta a aspect ratio da seccao submergida da strut
    k_RS=1/(13+(13/(sub_RS)**2)+(26/(sub_RS))) # coeficiente de lift que toma em conta a aspect ratio da seccao submergida da strut
    k_C=1/(13+(13/(sub_C)**2)+(26/(sub_C))) # coeficiente de lift que toma em conta a aspect ratio da seccao submergida da coluna
    L_LS = -1/2*rho*(speed_LS**2)*S_LS*k_LS*bb_LS
    L_RS = -1/2*rho*(speed_RS**2)*S_RS*k_RS*bb_RS
    L_C = -1/2*rho*(speed_C**2)*S_C*k_C*(bb_C+Rudder)

    #Induced Drag (XFLR5)
    e=1 #oswald efficiency factor (assumir 1 para dist. lift puramente eliptica)
    C_Dind_ReF=(CL_3D_ReF)**2/(pi*e*(AR_ReF))
    D_induced_ReF=(C_Dind_ReF)*q_ReF*(S_ReF) #Drag induzido (N)

    # Skin Friction Drag

    Re_crit=500000  #Como Rule Of Thumb, assume-se que a transição se dá pontualemente a Re=500 000- Anderson's Fundamentals of Aerodynamics (pg 903)
    x_crit_LF=Re_crit*mu/(rho*speed_LF) #calcula o ponto onde se dá a transição na asa retângular
    x_crit_RF=Re_crit*mu/(rho*speed_RF)
    x_crit_ReF=Re_crit*mu/(rho*speed_ReF)
    x_crit_LS=Re_crit*mu/(rho*speed_LS)
    x_crit_RS=Re_crit*mu/(rho*speed_RS)
    x_crit_C=Re_crit*mu/(rho*speed_C)
    #x_c_crit_n63412= #  <-M
    #x_c_crit_n0012=0.3 #     <-M
    #VARIAVEIS DE INPUT VAZIAS PARA FORÇAR A
    #TRANSIÇAO NUM PONTO QUE NAO RE 500 000

    #Left foil

    C_skf_ff_lam=1.328/(Re_crit)**0.5 #Coef skin friction esc laminar- Hoerner’s Fluid-Dynamic Drag (pg 2-5)-Foil frente
    S_ff_lam=S_FF*x_crit_LF #Area de escoamento laminar- corresponde a zona 1-lam
    D_skf_ff_lam=q_LF*C_skf_ff_lam*S_ff_lam #Drag skf foil frente_flat plate_zona 1-lam
    C_skf_ff_tur_12=0.455/(np.log10(Re_LF)**2.58) #Coef skin friction esc turbulento- Hoerner’s Fluid-Dynamic Drag (pg 2-5)-fim zona 2-tur
    C_skf_ff_tur_1=0.455/(np.log10(Re_crit)**2.58) #coef skf ff esc turbulento- fim zona 1-trans
    D_skf_ff_tur_12=q_LF*S_FF*C_skf_ff_tur_12 #drag skf ff_flat plate_esc turbulento-fim zona 2-tur
    D_skf_ff_tur_1=q_LF*S_ff_lam*C_skf_ff_tur_1#drag skf ff_flat plate_esc turbulento- fim zona 1-trans
    D_skf_ff_tur=D_skf_ff_tur_12-D_skf_ff_tur_1#drag skf ff_flat plate_esc turbulento total
    D_skf_ff_tot_1=D_skf_ff_lam +D_skf_ff_tur #drag skf ff_flat plate_total
    CD_skf_ff=D_skf_ff_tot_1/(q_LF*S_FF) #coeficiente de drag de skin friction ff final
    D_skf_LF=1/2*(CD_skf_ff)*(rho)*(speed_LF)**2*(SA_FF) #Drag de skin friction final ff

    #Right foil
    C_skf_ff_lam=1.328/(Re_crit)**0.5 #Coef skin friction esc laminar- Hoerner’s Fluid-Dynamic Drag (pg 2-5)-Foil frente
    S_ff_lam=S_FF*x_crit_RF #Area de escoamento laminar- corresponde a zona 1-lam
    D_skf_ff_lam=q_RF*C_skf_ff_lam*S_ff_lam #Drag skf foil frente_flat plate_zona 1-lam
    C_skf_ff_tur_12=0.455/(np.log10(Re_RF)**2.58) #Coef skin friction esc turbulento- Hoerner’s Fluid-Dynamic Drag (pg 2-5)-fim zona 2-tur
    C_skf_ff_tur_1=0.455/(np.log10(Re_crit)**2.58) #coef skf ff esc turbulento- fim zona 1-trans
    D_skf_ff_tur_12=q_RF*S_FF*C_skf_ff_tur_12 #drag skf ff_flat plate_esc turbulento-fim zona 2-tur
    D_skf_ff_tur_1=q_RF*S_ff_lam*C_skf_ff_tur_1#drag skf ff_flat plate_esc turbulento- fim zona 1-trans
    D_skf_ff_tur=D_skf_ff_tur_12-D_skf_ff_tur_1#drag skf ff_flat plate_esc turbulento total
    D_skf_ff_tot_1=D_skf_ff_lam +D_skf_ff_tur #drag skf ff_flat plate_total
    CD_skf_ff=D_skf_ff_tot_1/(q_RF*S_FF) #coeficiente de drag de skin friction ff final
    D_skf_RF=1/2*(CD_skf_ff)*(rho)*(speed_RF)**2*(SA_FF) #Drag de skin friction final ff

    #Rear foil
    C_skf_ft_lam=1.328/(Re_crit)**0.5#Coef skin friction esc laminar- Hoerner’s Fluid-Dynamic Drag (pg 2-5)-Foil tras
    S_ft_lam=S_ReF*x_crit_ReF #Area de escoamento laminar- corresponde a zona 1-lam
    D_skf_ft_lam=q_ReF*C_skf_ft_lam*S_ft_lam#Drag skf foil tras_flat plate_zona 1-lam
    C_skf_ft_tur_12=0.455/(np.log10(Re_ReF)**2.58) #Coef skin friction esc turbulento- Hoerner’s Fluid-Dynamic Drag (pg 2-5)-fim zona 2-tur
    C_skf_ft_tur_1=0.455/(np.log10(Re_crit)**2.58) #coef skf ft esc turbulento- fim zona 1-trans
    D_skf_ft_tur_12=q_ReF*S_ReF*C_skf_ft_tur_12#drag skf ft_flat plate_esc turbulento-fim zona 2-tur
    D_skf_ft_tur_1=q_ReF*S_ft_lam*C_skf_ft_tur_1#drag skf ft_flat plate_esc turbulento- fim zona 1-trans
    D_skf_ft_tur=D_skf_ft_tur_12-D_skf_ft_tur_1#drag skf ft_flat plate_esc turbulento total
    D_skf_ft_tot_1=D_skf_ft_lam +D_skf_ft_tur #drag skf ft_flat plate_total
    CD_skf_ft=D_skf_ft_tot_1/(q_ReF*S_ReF)#coeficiente de drag de skin friction ft final
    D_ftskt=1/2*(CD_skf_ft)*(rho)*(speed_ReF)**2*(SA_ReF)#Drag de skin friction final ft



    #Strut Esquerda
    C_skf_s_lam=1.328/(Re_crit)**0.5#Coef skin friction esc laminar- Hoerner’s Fluid-Dynamic Drag (pg 2-5)-strut
    S_s_lam_left=S_LS*x_crit_LS #Area de escoamento laminar- corresponde a zona 1-lam
    D_skf_s_lam_left=q_LS*C_skf_s_lam*S_s_lam_left#Drag skf strut_flat plate_zona 1-lam
    C_skf_s_tur_12=0.455/(np.log10(Re_LS)**2.58)#Coef skin friction esc turbulento- Hoerner’s Fluid-Dynamic Drag (pg 2-5)-fim zona 2-tur
    C_skf_s_tur_1=0.455/(np.log10(Re_crit)**2.58)#coef skf s esc turbulento- fim zona 1-trans
    D_skf_s_tur_12_left=q_LS*S_LS*C_skf_s_tur_12#drag skf s_flat plate_esc turbulento-fim zona 2-tur
    D_skf_s_tur_1_left=q_LS*S_s_lam_left*C_skf_s_tur_1#drag skf s_flat plate_esc turbulento- fim zona 1-trans
    D_skf_s_tur_left=D_skf_s_tur_12_left-D_skf_s_tur_1_left #drag skf s_flat plate_esc turbulento total
    D_skf_s_tot_1_left=D_skf_s_lam_left +D_skf_s_tur_left #drag skf s_flat plate_total
    CD_skf_s_left=D_skf_s_tot_1_left/(q_LS*S_LS)#coeficiente de drag de skin friction s final
    D_strutskf_left=1/2*(CD_skf_s_left)*(rho)*(speed_LS)**2*(SA_LS)#Drag de skin friction final strut

    #Strut Direita
    C_skf_s_lam=1.328/(Re_crit)**0.5#Coef skin friction esc laminar- Hoerner’s Fluid-Dynamic Drag (pg 2-5)-strut
    S_s_lam_right=S_RS*x_crit_RS #Area de escoamento laminar- corresponde a zona 1-lam
    D_skf_s_lam_right=q_RS*C_skf_s_lam*S_s_lam_right#Drag skf strut_flat plate_zona 1-lam
    C_skf_s_tur_12=0.455/(np.log10(Re_RS)**2.58)#Coef skin friction esc turbulento- Hoerner’s Fluid-Dynamic Drag (pg 2-5)-fim zona 2-tur
    C_skf_s_tur_1=0.455/(np.log10(Re_crit)**2.58)#coef skf s esc turbulento- fim zona 1-trans
    D_skf_s_tur_12_right=q_RS*S_RS*C_skf_s_tur_12#drag skf s_flat plate_esc turbulento-fim zona 2-tur
    D_skf_s_tur_1_right=q_RS*S_s_lam_right*C_skf_s_tur_1#drag skf s_flat plate_esc turbulento- fim zona 1-trans
    D_skf_s_tur_right=D_skf_s_tur_12_right-D_skf_s_tur_1_right #drag skf s_flat plate_esc turbulento total
    D_skf_s_tot_1_right=D_skf_s_lam_right +D_skf_s_tur_right #drag skf s_flat plate_total
    CD_skf_s_right=D_skf_s_tot_1_right/(q_RS*S_RS)#coeficiente de drag de skin friction s final
    D_strutskf_right=1/2*(CD_skf_s_right)*(rho)*(speed_RS)**2*(SA_RS)#Drag de skin friction final strut


    #Column
    C_skf_c_lam=1.328/(Re_crit)**0.5#Coef skin friction esc laminar- Hoerner’s Fluid-Dynamic Drag (pg 2-5)-coluna
    S_c_lam=S_C*x_crit_C#Area de escoamento laminar- corresponde a zona 1-lam
    D_skf_c_lam=q_C*C_skf_c_lam*S_c_lam#Drag skf coluna_flat plate_zona 1-lam
    C_skf_c_tur_lt=0.455/(np.log10(Re_C)**2.58)#Coef skin friction esc turbulento- Hoerner’s Fluid-Dynamic Drag (pg 2-5)-fim zona 2-tur
    C_skf_c_tur_l=0.455/(np.log10(Re_crit)**2.58)#coef skf c esc turbulento- fim zona 1-trans
    D_skf_c_tur_lt=q_C*S_C*C_skf_c_tur_lt#drag skf c_flat plate_esc turbulento-fim zona 2-tur
    D_skf_c_tur_l=q_C*S_c_lam*C_skf_c_tur_l#drag skf c_flat plate_esc turbulento- fim zona 1-trans
    D_skf_c_tur=D_skf_c_tur_lt-D_skf_c_tur_l#drag skf c_flat plate_esc turbulento total
    D_skf_c_tot_1=D_skf_c_lam +D_skf_c_tur#drag skf c_flat plate_total
    CD_skf_c=D_skf_c_tot_1/(q_C*S_C)#coeficiente de drag de skin friction c final
    D_skfc=1/2*(CD_skf_c)*(rho)*(speed_C)**2*(SA_C) #Drag de skin friction final coluna

    #Form Drag
    D_form_ReF=(cdp_ReF)*q_ReF*(S_ReF)
    D_form_LS=(CDp_LS)*q_LS*(S_LS)
    D_form_RS=(CDp_RS)*q_RS*(S_RS)
    D_form_C=(CDp_C)*q_C*(2*S_C)

    #Interference Drag
    C_Dint_FF=17*(t_FF)**2-.05 #Coef Drag interferencia- Hoerner’s Fluid Dynamic Drag (pg 8-1)
    C_Dint_ReF=17*(t_ReF)**2-.05 #Coef Drag interferencia- Hoerner’s Fluid Dynamic Drag (pg 8-1)
    D_interference_LF=q_LF*(C_Dint_FF)*((T_FF+T_S)/2)**2#Drag Interferencia (N)
    D_interference_RF=q_RF*(C_Dint_FF)*((T_FF+T_S)/2)**2#Drag Interferencia (N)
    D_interference_ReF=q_ReF*(C_Dint_ReF)*((T_ReF+T_C)/2)**2#Drag Interferencia (N)

    #Spray Drag
    F_h_LS= speed_LS/(g*h_LS)**(0.5)#Numero de Froude a profundidade da strut
    F_h_RS= speed_RS/(g*h_RS)**(0.5)#Numero de Froude a profundidade da strut
    F_h_C= speed_C/(g*h_C)**(0.5)#Numero de Froude a profundidade da coluna
    #Froude numbers ~3<
    CDs=0.24 #coeficiente empirico de spray @Hoerner_Fluid-Dynamic Drag (10-13)
    Ds_LS=q_LS*CDs*(t_S*c_S)**2#Drag de spray para uma strut
    Ds_RS=q_RS*CDs*(t_S*c_S)**2#Drag de spray para uma strut
    Ds_C=q_ReF*CDs*(t_C*c_C)**2#Drag de spray para a coluna

    #Wave Drag 
    CDw_C=0.5*(CL_3D_ReF/(F_h_C*2))**2
    Dw_C=q_C*CDw_C*h_C*b_ReF

    #Total Drag
    D_ReF=D_induced_ReF+D_ftskt+D_form_ReF+D_interference_ReF+Dw_C
    D_LS=D_strutskf_left+D_form_LS+Ds_LS
    D_RS=D_strutskf_right+D_form_RS+Ds_RS
    D_C=D_skfc+D_form_C+Ds_C

    #Hull
    CD_H = 0.006# CFD Data
    zp = -z_cm*100
    if (zp < 32 + 5/22):
        Hull_area = -0.2167*zp+6.98365
    else:
        Hull_area = 0

    D_H = 0.5*CD_H*rho*Hull_area*speed_cm**2
    V_HULL = 0.001*zp**2-0.064*zp+1.03025
    if (zp >= 31):
        V_HULL = 0

    BUOYANCY = V_HULL*rho*g

    #calculate Thrust based on rpms
    Rpms_propeller = Rpms_motor/2.85
    diam = 0.25

    js = 60*speed_cm/(Rpms_propeller*diam)

    def ktfun(x):
        #Interpolation coefficients of thrust
        p1 =     -0.2264
        p2 =      0.8165
        p3 =     -0.9315
        p4 =      0.2713
        p5 =   -0.009195
        p6 =   -0.005493
        p7 =      0.1617
        return p1*x**6 + p2*x**5 + p3*x**4 + p4*x**3 + p5*x**2 + p6*x + p7

    kt =ktfun(js)

    def Efffun(x):
        #Interpolation coefficients of proppeler's efficiency
        k1 =      -170.6
        k2 =       30.69
        k3 =      -50.99
        k4 =       440.2
        k5 =      -1.744
        q1 =      -201.3
        q2 =       93.96
        q3 =       150.8
        q4 =      -151.9
        q5 =       398.2
        return ((k1*x**4 + k2*x**3 + k3*x**2 + k4*x + k5) /
            (x**5 + q1*x**4 + q2*x**3 + q3*x**2 + q4*x + q5))

    Eff = Efffun(js)

    Thrust = rho*(Rpms_propeller/60)**2*diam**4*kt
    P_propeller = Thrust*speed_cm/Eff

    #Calculate Power usage
    P_mec=P_propeller/0.96
    Torque=P_mec/(Rpms_motor*2*pi/60)

    #Inertial to body frame matrix
    R = np.array([[cosd(Pitch)*cosd(Yaw), cosd(Pitch)*sind(Yaw), -sind(Pitch)],
        [sind(Roll)*sind(Pitch)*cosd(Yaw)-cosd(Roll)*sind(Yaw), sind(Roll)*sind(Pitch)*sind(Yaw)+cosd(Roll)*cosd(Yaw), sind(Roll)*cosd(Pitch)],
        [cosd(Roll)*sind(Pitch)*cosd(Yaw)+sind(Roll)*sind(Yaw), cosd(Roll)*sind(Pitch)*sind(Yaw)-sind(Roll)*cosd(Yaw), cosd(Roll)*cosd(Pitch)]])

    #Inertial forces
    W_inertial = m*g*np.array([0,0,1])
    BUOY_inertial = BUOYANCY*np.array([0,0,-1])

    #Body frame forces
    L_b_ReF = L_ReF*np.array([sind(aa_ReF),0,-cosd(aa_ReF)])
    L_b_LS = L_LS*np.array([-sind(bb_LS),cosd(bb_LS),0])
    L_b_RS = L_RS*np.array([-sind(bb_RS),cosd(bb_RS),0])
    L_b_C = L_C*np.array([-sind(bb_C),cosd(bb_C),0])
    D_b_ReF = D_ReF*np.array([-cosd(aa_ReF),0,-sind(aa_ReF)])
    D_b_LS = D_LS*np.array([-cosd(bb_LS),-sind(bb_LS),0])
    D_b_RS = D_RS*np.array([-cosd(bb_RS),-sind(bb_RS),0])
    D_b_C = D_C*np.array([-cosd(bb_C),-sind(bb_C),0])
    D_b_H = D_H*np.array([-u/speed_cm,-v/speed_cm,-w/speed_cm])

    W_b = R@W_inertial
    BUOY_b = R@BUOY_inertial

    #Lift_Strut_left  = -1/2*rho*(speed)**2*Ss_left*0.1121*bb
    #Lift_Strut_right  = -1/2*rho*(speed)**2*Ss_right*0.1121*bb
    #Lift_Coluna  = -1/2*rho*(speed**2)*Sc*0.1121*(bb_rear+Rudder)

    #Sum of body frame forces in each location to calculate momenta
    X_ReF = L_b_ReF[0]+D_b_ReF[0]
    Y_ReF =  0
    Z_ReF = L_b_ReF[2]+D_b_ReF[2]
    X_LS = L_b_LS[0] + D_b_LS[0]
    Y_LS = L_b_LS[1] + D_b_LS[1]
    Z_LS = 0
    X_RS = L_b_RS[0] + D_b_RS[0]
    Y_RS = L_b_RS[1] + D_b_RS[1]
    Z_RS = 0
    X_C = L_b_C[0] + D_b_C[0]
    Y_C = L_b_C[1] + D_b_C[1]
    Z_C = 0
    Y_Prop = -Thrust*sind(Rudder)
    X_Prop = Thrust*cosd(Rudder)
        
    J = np.array([ [Ixx, 0, Ixz] , [0, Iyy , 0], [Ixz, 0, Izz] ])
    w = np.array([p,qq,r])*pi/180
    GyrTerm = np.linalg.inv(J)@np.cross(w,J@w)

    Y = Y_ReF + Y_LS + Y_RS + Y_C + Y_Prop + W_b[1] + BUOY_b[1] + D_b_H[1]

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
    a15 = ((AR_FF/(AR_FF+2)*k2_LF)**2)/(pi*e*(AR_FF))*q_LF*S_FF
    a16 = q_LF*0.5/((F_h_LS*2)**2)*h_LS*b_FF*((AR_FF/(AR_FF+2)*k2_LF)**2)
    a17 = q_LF*(S_FF)
    a18 = D_skf_RF + D_interference_RF
    a19 = ((AR_FF/(AR_FF+2)*k2_RF)**2)/(pi*e*(AR_FF))*q_RF*S_FF + q_RF*0.5/((F_h_RS*2)**2)*h_RS*b_FF*((AR_FF/(AR_FF+2)*k2_RF)**2)
    a20 = q_RF*(S_FF)
    a21 = (u_phi + GyrTerm[0])*(Ixx*Izz - Ixz**2)
    a22 = Izz*(-Y_LS*LS_CoP-Y_RS*RS_CoP - Y_C*C_CoP - Y_Prop*(dist_ReF_cm_z-0.2) )
    a23 = - Ixz*(X_LS*dist_FF_cm_y-X_RS*dist_FF_cm_y+Y_LS*dist_FF_cm_x+Y_RS*dist_FF_cm_x-Y_C*dist_C_cm_x-Y_Prop*dist_C_cm_x)
    a24 = Izz*dist_FF_cm_y*AR_FF/(AR_FF+2)*k2_LF*q_LF*S_FF*cosd(aa_LF)
    a25 =  -Ixz*dist_FF_cm_y*AR_FF/(AR_FF+2)*k2_LF*q_LF*S_FF*sind(aa_LF)
    a26 = -Izz*dist_FF_cm_y*AR_FF/(AR_FF+2)*k2_RF*q_RF*S_FF*cosd(aa_RF)
    a27 = Ixz*dist_FF_cm_y*AR_FF/(AR_FF+2)*k2_RF*q_RF*S_FF*sind(aa_RF)
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

    ###########################

    return a37, a29, a30, a38, a39, a40, a41, a42, a10, a11, a43, a45, a44, a46, Re_LF, Re_RF, aa_LF, aa_RF

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

    p = np.poly1d([c3, c2, c1, c0])

    x = -3
    for i in range(1,5):
        x = x - (c3*x**3+c2*x**2+c1*x+c0)/(a*x**2 + b*x + c)

    #def f(x,c3,c2,c1,c0):

        #return abs(c3*x**3 + c2*x**2 + c1*x + c0)

    #min = minimize(f, args = (c3,c2,c1,c0), x0=-1.55)

    return x

def cl_constants(Re):

    Re_LF = Re/1000

    p00 = 0.059320480930810174
    p10 = 0.13106498910132924
    p01 = 0.001499763421768745
    p20 = -0.0006103181410530192
    p11 = -0.0001128310042421913
    p02 = -2.7940023788168525e-06
    p30 = -0.00039240959649213165
    p21 = -5.402562318825203e-08
    p12 = 1.5450015808831076e-07
    p03 = 2.428005401349985e-09
    p31 = 6.021581789343582e-07
    p22 = -3.892020443343621e-09
    p13 = -7.221007926942857e-11
    p04 = -9.966049389256757e-13
    p41 = 1.0
    p32 = -2.6569952741295344e-10
    p23 = 2.2670065015109916e-12
    p14 = 9.7510028025199e-15
    p05 = 1.555016012253181e-16

    cl0 = p00 + p01*Re_LF + p02*Re_LF**2 + p03*Re_LF**3 + p04*Re_LF**4 + p05*Re_LF**5
    cl1 = p11*Re_LF + p12*Re_LF**2 + p13*Re_LF**3 + p14*Re_LF**4 + p10
    cl2 = p21*Re_LF + p22*Re_LF**2 + p23*Re_LF**3 + p20
    cl3 = p31*Re_LF + p32*Re_LF**2 + + p30

    return cl3, cl2, cl1, cl0

def cd_constants(Re):

    Re_LF = Re/1000

    p00 = 0.010300392869687475
    p10 = -0.00048176078975881555
    p01 = -5.712750436466757e-05
    p20 = 0.00026343810839221925
    p11 = 2.5604038784960047e-06
    p02 = 8.998879806473492e-08
    p30 = 8.175657100873038e-06
    p21 = 7.697207229383863e-07
    p12 = -2.11500349374705e-09
    p03 = -7.445739399129038e-11
    p31 = -9.815018505790247e-08
    p22 = -2.6039954027120664e-10
    p13 = 1.2130023576195772e-12
    p04 = 3.001760671936499e-14
    p41 = 1.0
    p32 = 1.913998778001434e-11
    p23 = 3.917985844337982e-14
    p14 = -2.5990037923286277e-16
    p05 = -4.681217559756671e-18

    cd0 = p00 + p01*Re_LF + p02*Re_LF**2 + p03*Re_LF**3 + p04*Re_LF**4 + p05*Re_LF**5
    cd1 = p11*Re_LF + p12*Re_LF**2 + p13*Re_LF**3 + p14*Re_LF**4 + p10
    cd2 = p21*Re_LF + p22*Re_LF**2 + p23*Re_LF**3 + p20
    cd3 = p31*Re_LF + p32*Re_LF**2 + + p30

    return cd3, cd2, cd1, cd0

def clsq_constants(Re):

    Re_LF = Re/1000

    p00 = 0.007781094029274539
    p10 = 0.04087970593500313
    p01 = 0.0003791002829874271
    p20 = 0.015173181882318807
    p11 = 0.00023217503697906588
    p02 = 6.664608108458437e-08
    p30 = -0.001110947027172005
    p21 = -2.2590811517247053e-05
    p12 = -4.1451454674740763e-07
    p03 = -1.1450348930619397e-09
    p31 = 1.1624579410125995e-06
    p22 = 2.6008638506434562e-08
    p13 = 2.748778509390584e-10
    p04 = 1.2018831107581225e-12
    p41 = 1.0
    p32 = -5.676551192730255e-10
    p23 = -9.091960835043724e-12
    p14 = -5.467193414476767e-14
    p05 = -3.751085100619171e-16

    clsq0 = p00 + p01*Re_LF + p02*Re_LF**2 + p03*Re_LF**3 + p04*Re_LF**4 + p05*Re_LF**5
    clsq1 = p11*Re_LF + p12*Re_LF**2 + p13*Re_LF**3 + p14*Re_LF**4 + p10
    clsq2 = p21*Re_LF + p22*Re_LF**2 + p23*Re_LF**3 + p20
    clsq3 = p31*Re_LF + p32*Re_LF**2 + + p30

    return clsq3, clsq2, clsq1, clsq0

def p3(x, c3, c2, c1, c0):

    val = c3*x**3 + c2*x**2 + c1*x + c0

    return val

