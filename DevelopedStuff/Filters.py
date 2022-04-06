import numpy as np
import scipy.io
from scipy.io import loadmat 
import matplotlib.pyplot as plt
import pprint   
from tabulate import tabulate

class median_filter(object):
    
    def __init__(self, initial_state, window):
        self.window = window
        self.state = []
        for i in range(0, window):
            self.state.append(initial_state)

    def update(self, new):
        
        aux = self.state.copy()
        self.state[0] = new
        for i in range(0, self.window-1):
            self.state[i+1] = aux[i]

        self.median = np.median(self.state)


        return self.median

class states_filter(object):

    def __init__(self, initial_states, window):

        self.filter = []
        for i in range(0,len(initial_states)):
            self.filter.append(median_filter(initial_states[i], window[i]))

    def update(self, new):

        output = []

        for i in range(0,len(new)):
            output.append(self.filter[i].update(new[i]))

        return output

class Kalman_filter(object):
    
    def __init__(self, initial_states, A, B, H, Q, R):

        self.A = np.array(A)
        self.B = np.array(B)
        self.H = np.array(H)
        self.x0 = initial_states
        self.Q = np.array(Q)
        self.R = np.array(R)

    def update(self, uk, zk):

        self.xk_pred = self.A@self.xk_corr + self.B@uk
        self.Pk_pred = self.A@self.Pk_corr@np.transpose(self.A) + self.Q 

        self.K = self.Pk_pred@np.transpose(self.H)@np.linalg.inv(self.H@self.Pk_pred@np.transpose(self.H)+self.R)
        self.xk_corr = self.xk_pred + self.K@(zk-self.H@self.xk_pred)
        self.Pk_corr = (np.identity(len(self.xk_pred))-self.K@self.H)@self.Pk_pred

        return self.xk_corr

# states : z, zdot
#zdot = Vz
#zdotdot = Accz
# A = [0]
# 









