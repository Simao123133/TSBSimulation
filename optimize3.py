from calculate_constants3 import constants3
from func3 import func_func3
from func3 import gradient
from func3 import hessian
from scipy.optimize import minimize, root, least_squares
import numpy as np
import time
from func_newton import func_split, gradient_split
from Simplified_sys import get_aa
import matplotlib.pyplot as plt
import random 
from scipy.signal import medfilt

Y =[7, 0.1, 0, 0, -0.5, 0, 0, 0, 0, 0]
params = [2.02, -5, 4780]
u_z = 0
u_phi = 0

a47, a48, a49, a50, a51, a52, a54, a55, a56, a57, a58, a59, a77, a62, a63, a64, a65, a66, a67, a69, a70, a71, a72, a73, a74, a78, Re_LF, Re_RF, aa_LF, aa_RF = constants3(Y, params, u_z, u_phi)

start = time.time()
sol1_low = np.zeros(100)
sol1_med = np.zeros(100)
sol1_mean = np.zeros(100)
sol1_cool = np.zeros(100)
Y_sensor = np.array([7, 0, 0, 0, -0.5, 0, 0, 0, 0, 0])
Y_filtered = np.array([7, 0, 0, 0, -0.5, 0, 0, 0, 0, 0])
Y_data = np.zeros([120, 10])
Y_mean = np.array([7, 0, 0, 0, -0.5, 0, 0, 0, 0, 0])
Y_med = np.array([7, 0, 0, 0, -0.5, 0, 0, 0, 0, 0])

def median(lst): return np.median(np.array(lst))
def mean(lst): return sum(lst)/len(lst)

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
            self.state[i+1] = aux[i])

        self.median = median(self.state)


        return self.median

class states_filter(object):
    def __init__(self, initial_states, window):

        self.filter = []
        for i in range(0,len(initial_states)):
            self.filter.append(median_filter(initial_states[i], window[i]))

    def update(self, new):

        output = []

        for i in range(0,len(new)):
            print(i)
            output.append(self.filter[i].update(new[i]))

        return output
    

a47, a48, a49, a50, a51, a52, a54, a55, a56, a57, a58, a59, a77, a62, a63, a64, a65, a66, a67, a69, a70, a71, a72, a73, a74, a78, Re_LF, Re_RF, aa_LF, aa_RF = constants3(Y, params, u_z, u_phi)
res = root(func_split, [1.55, 1.55], args=(a47, a48, a49, a50, a51, a52, a54, a55, a56, a57, a58, a59, a77, a62, a63, a64, a65, a66, a67, a69, a70, a71, a72, a73, a74, a78, Re_LF, Re_RF, aa_LF, aa_RF), method='hybr')
sol1 = res.x[0]

for i in range(0,30):
    Y_data[i,0] = Y[0] + random.gauss(0,0.05)
    Y_data[i,1] = Y[1] + random.gauss(0,0.05)
    Y_data[i,2] = Y[2] + random.gauss(0,0.035)
    Y_data[i,3] = Y[3] + random.gauss(0,0.2)
    Y_data[i,4] = Y[4] + random.gauss(0,0.001)
    Y_data[i,5] = Y[5] + random.gauss(0,0.05)
    Y_data[i,6] = Y[6] + random.gauss(0,0.035)
    Y_data[i,7] = Y[7] + random.gauss(0,0.035)
    Y_data[i,8] = Y[8] + random.gauss(0,0.2)
    Y_data[i,9] = Y[9] + random.gauss(0,0.5)

filter_c = 0.95

u_win = 30
w_win = 30
q_win = 30
pitch_win = 30
z_win = 30
v_win = 30
p_win = 30
r_win = 30
roll_win = 30
yaw_win = 30

windows = [u_win, w_win, q_win, pitch_win, z_win, v_win, p_win, r_win, roll_win, yaw_win]

filters = states_filter(Y, windows)

for i in range(0,100):
    Y_sensor[0] = Y[0] + random.gauss(0,0.05)
    Y_sensor[1] = Y[1] + random.gauss(0,0.05)
    Y_sensor[2] = Y[2] + random.gauss(0,0.035)
    Y_sensor[3] = Y[3] + random.gauss(0,0.2)
    Y_sensor[4] = Y[4] + random.gauss(0,0.001)
    Y_sensor[5] = Y[5] + random.gauss(0,0.05)
    Y_sensor[6] = Y[6] + random.gauss(0,0.035)
    Y_sensor[7] = Y[7] + random.gauss(0,0.035)
    Y_sensor[8] = Y[8] + random.gauss(0,0.2)
    Y_sensor[9] = Y[9] + random.gauss(0,0.5)

    Y_made = filters.update(Y_sensor)

    Y_low = np.multiply(1-filter_c,Y_filtered) + np.multiply(filter_c,np.array(Y_sensor[0:10]))

    for j in range(0,10):
        Y_mean[j] = mean(Y_data[i:i+30,j])
        Y_med[j] = median(Y_data[i:i+30,j])

    a47, a48, a49, a50, a51, a52, a54, a55, a56, a57, a58, a59, a77, a62, a63, a64, a65, a66, a67, a69, a70, a71, a72, a73, a74, a78, Re_LF, Re_RF, aa_LF, aa_RF = constants3(Y_low, params, u_z, u_phi)
    res = root(func_split, [1.55, 1.55], args=(a47, a48, a49, a50, a51, a52, a54, a55, a56, a57, a58, a59, a77, a62, a63, a64, a65, a66, a67, a69, a70, a71, a72, a73, a74, a78, Re_LF, Re_RF, aa_LF, aa_RF), method='hybr')
    sol1_low[i] = res.x[0]
 

    a47, a48, a49, a50, a51, a52, a54, a55, a56, a57, a58, a59, a77, a62, a63, a64, a65, a66, a67, a69, a70, a71, a72, a73, a74, a78, Re_LF, Re_RF, aa_LF, aa_RF = constants3(Y_med, params, u_z, u_phi)
    res = root(func_split, [1.55, 1.55], args=(a47, a48, a49, a50, a51, a52, a54, a55, a56, a57, a58, a59, a77, a62, a63, a64, a65, a66, a67, a69, a70, a71, a72, a73, a74, a78, Re_LF, Re_RF, aa_LF, aa_RF), method='hybr')
    sol1_med[i] = res.x[0]


    a47, a48, a49, a50, a51, a52, a54, a55, a56, a57, a58, a59, a77, a62, a63, a64, a65, a66, a67, a69, a70, a71, a72, a73, a74, a78, Re_LF, Re_RF, aa_LF, aa_RF = constants3(Y_mean, params, u_z, u_phi)
    res = root(func_split, [1.55, 1.55], args=(a47, a48, a49, a50, a51, a52, a54, a55, a56, a57, a58, a59, a77, a62, a63, a64, a65, a66, a67, a69, a70, a71, a72, a73, a74, a78, Re_LF, Re_RF, aa_LF, aa_RF), method='hybr')
    sol1_mean[i] = res.x[0]

    a47, a48, a49, a50, a51, a52, a54, a55, a56, a57, a58, a59, a77, a62, a63, a64, a65, a66, a67, a69, a70, a71, a72, a73, a74, a78, Re_LF, Re_RF, aa_LF, aa_RF = constants3(Y_made, params, u_z, u_phi)
    res = root(func_split, [1.55, 1.55], args=(a47, a48, a49, a50, a51, a52, a54, a55, a56, a57, a58, a59, a77, a62, a63, a64, a65, a66, a67, a69, a70, a71, a72, a73, a74, a78, Re_LF, Re_RF, aa_LF, aa_RF), method='hybr')
    sol1_cool[i] = res.x[0]
    

a_list = list(range(0, 100))
plt.plot(a_list,sol1_cool, a_list, sol1_med)

MSE_median = sum((sol1 - sol1_med)**2)/len(sol1_med)
MSE_mean = sum((sol1 - sol1_mean)**2)/len(sol1_med)
MSE_low = sum((sol1 - sol1_low)**2)/len(sol1_med)
print(MSE_median)
print(MSE_mean)
print(MSE_low)

plt.show()



print(res)  
print(time.time()-start)


start = time.time()
res = get_aa(Y, params, u_z, u_phi)
print(res)
print(time.time()-start)