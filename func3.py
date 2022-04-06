import numpy as np

def func_func3(x,a47, a48, a49, a50, a51, a52, a54, a55, a56, a57, a58, a59, a77, a62, a63, a64, a65, a66, a67, a69, a70, a71, a72, a73, a74, a78, Re_LF, Re_RF, aa_LF, aa_RF):

    aa_LF = x[0] + aa_LF
    aa_RF = x[1] + aa_RF

    f1 = a47*aa_LF**6 + a48*aa_LF**5 + a49*aa_LF**4 + a50*aa_LF**3 + a51*aa_LF**2 + a52*aa_LF + a54*aa_RF**6 + a55*aa_RF**5 + a56*aa_RF**4 + a57*aa_RF**3 + a58*aa_RF**2 + a59*aa_RF + a77 
    f2 = a62*aa_LF**6 + a63*aa_LF**5 + a64*aa_LF**4 + a65*aa_LF**3 + a66*aa_LF**2 + a67*aa_LF + a69*aa_RF**6 + a70*aa_RF**5 + a71*aa_RF**4 + a72*aa_RF**3 + a73*aa_RF**2 + a74*aa_RF + a78

    return f1**2 + f2**2

def gradient(x,a47, a48, a49, a50, a51, a52, a54, a55, a56, a57, a58, a59, a77, a62, a63, a64, a65, a66, a67, a69, a70, a71, a72, a73, a74, a78, Re_LF, Re_RF, aa_LF, aa_RF):

    aa_LF = x[0] + aa_LF
    aa_RF = x[1] + aa_RF

    f1 = a47*aa_LF**6 + a48*aa_LF**5 + a49*aa_LF**4 + a50*aa_LF**3 + a51*aa_LF**2 + a52*aa_LF + a54*aa_RF**6 + a55*aa_RF**5 + a56*aa_RF**4 + a57*aa_RF**3 + a58*aa_RF**2 + a59*aa_RF + a77 
    f2 = a62*aa_LF**6 + a63*aa_LF**5 + a64*aa_LF**4 + a65*aa_LF**3 + a66*aa_LF**2 + a67*aa_LF + a69*aa_RF**6 + a70*aa_RF**5 + a71*aa_RF**4 + a72*aa_RF**3 + a73*aa_RF**2 + a74*aa_RF + a78

    f1_aa_LF = 6*a47*aa_LF**5 + 5*a48*aa_LF**4 + 4*a49*aa_LF**3 + 3*a50*aa_LF**2 + 2*a51*aa_LF + a52 
    f1_aa_RF = 6*a54*aa_RF**5 + 5*a55*aa_RF**4 + 4*a56*aa_RF**3 + 3*a57*aa_RF**2 + 2*a58*aa_RF + a59
    f2_aa_LF = 6*a62*aa_LF**5 + 5*a63*aa_LF**4 + 4*a64*aa_LF**3 + 3*a65*aa_LF**2 + 2*a66*aa_LF + a67 
    f2_aa_RF = 6*a69*aa_RF**5 + 5*a70*aa_RF**4 + 4*a71*aa_RF**3 + 3*a72*aa_RF**2 + 2*a73*aa_RF + a74

    gradf_aa_LF = 2*f1_aa_LF*f1 + 2*f2_aa_LF*f2
    gradf_aa_RF = 2*f1_aa_RF*f1 + 2*f2_aa_RF*f2

    return gradf_aa_LF, gradf_aa_RF

def hessian(x,a47, a48, a49, a50, a51, a52, a54, a55, a56, a57, a58, a59, a77, a62, a63, a64, a65, a66, a67, a69, a70, a71, a72, a73, a74, a78, Re_LF, Re_RF, aa_LF, aa_RF):

    aa_LF = x[0] + aa_LF
    aa_RF = x[1] + aa_RF

    f1 = a47*aa_LF**6 + a48*aa_LF**5 + a49*aa_LF**4 + a50*aa_LF**3 + a51*aa_LF**2 + a52*aa_LF + a54*aa_RF**6 + a55*aa_RF**5 + a56*aa_RF**4 + a57*aa_RF**3 + a58*aa_RF**2 + a59*aa_RF + a77 
    f2 = a62*aa_LF**6 + a63*aa_LF**5 + a64*aa_LF**4 + a65*aa_LF**3 + a66*aa_LF**2 + a67*aa_LF + a69*aa_RF**6 + a70*aa_RF**5 + a71*aa_RF**4 + a72*aa_RF**3 + a73*aa_RF**2 + a74*aa_RF + a78

    f1_aa_LF = 6*a47*aa_LF**5 + 5*a48*aa_LF**4 + 4*a49*aa_LF**3 + 3*a50*aa_LF**2 + 2*a51*aa_LF + a52 
    f1_aa_RF = 6*a54*aa_RF**5 + 5*a55*aa_RF**4 + 4*a56*aa_RF**3 + 3*a57*aa_RF**2 + 2*a58*aa_RF + a59
    f2_aa_LF = 6*a62*aa_LF**5 + 5*a63*aa_LF**4 + 4*a64*aa_LF**3 + 3*a65*aa_LF**2 + 2*a66*aa_LF + a67 
    f2_aa_RF = 6*a69*aa_RF**5 + 5*a70*aa_RF**4 + 4*a71*aa_RF**3 + 3*a72*aa_RF**2 + 2*a73*aa_RF + a74

    gradf_aa_LF = 2*f1_aa_LF*f1 + 2*f2_aa_LF*f2
    gradf_aa_RF = 2*f1_aa_RF*f1 + 2*f2_aa_RF*f2

    f1_aa_LF_2 = 30*a47*aa_LF**4 + 20*a48*aa_LF**3 + 12*a49*aa_LF**2 + 6*a50*aa_LF + 2*a51 
    f1_aa_RF_2 = 30*a54*aa_RF**4 + 20*a55*aa_RF**3 + 12*a56*aa_RF**2 + 6*a57*aa_RF + 2*a58
    f2_aa_LF_2 = 30*a62*aa_LF**4 + 20*a63*aa_LF**3 + 12*a64*aa_LF**2 + 6*a65*aa_LF + 2*a66 
    f2_aa_RF_2 = 30*a69*aa_RF**4 + 20*a70*aa_RF**3 + 12*a71*aa_RF**2 + 6*a72*aa_RF + 2*a73


    hess_aa_LF_2 = 2*f1_aa_LF_2*f1 + 2*f1_aa_LF**2 + 2*f2_aa_LF_2*f2 + 2*f2_aa_LF**2 
    hess_aa_LF_RF = 2*f1_aa_LF*f1_aa_RF + 2*f2_aa_LF*f2_aa_RF
    hess_aa_RF_2 = 2*f1_aa_RF_2*f1 + 2*f1_aa_RF**2 + 2*f2_aa_RF_2*f2 + 2*f2_aa_RF**2
      

    return [[hess_aa_LF_2, hess_aa_LF_RF],[hess_aa_LF_RF, hess_aa_RF_2]]