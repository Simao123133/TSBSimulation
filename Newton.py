from calculate_constants import constants
from derivatives import gradient, hessian
from func_newton import func
from scipy.optimize import minimize
import numpy as np
import time


def Newton(f,df, args):

    
