from calculate_constants import constants
from derivatives import gradient, hessian
from func import func_func
from scipy.optimize import minimize
import numpy as np
import time
from SMC_MODEL import SMC_MODEL_func
from Nelder_mead import nelder_mead

a37, a29, a30, a38, a39, a40, a41, a42, a10, a11, a43, a45, a44, a46, Re_LF, Re_RF, aa_LF, aa_RF = constants([7, 0, 0, 0, -0.5, 0, 0, 0, 0, 0], [2.02, 0, 4780], 0, 0)

start = time.time()
res = minimize(func_func, np.array([1, 1]),args = (a37, a29, a30, a38, a39, a40, a41, a42, a10, a11, a43, a45, a44, a46, Re_LF, Re_RF, aa_LF, aa_RF), method='nelder-mead', bounds = {(-3,14), (-3,14)})
print(time.time()-start)
print(res)


