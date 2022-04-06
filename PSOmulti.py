import numpy as np
from SMC_MODEL import SMC_MODEL_func
import time
from calculate_constants import constants
from func import func_func
from func_newton import func

def PSO_func(ref, states, params, ub, lb):
    # Define the details of the table design problem
    nVar = 2
    # Define the PSO's paramters
    noP = 3
    maxIter = 40
    wMax = 0.9
    wMin = 0.4
    c1 = 1
    c2 = 1
    vMax = (ub - lb)* 0.2
    vMin  = -vMax
    # The PSO algorithm

    # Initialize the particles

    Swarm_Particles_X = (ub-lb) * np.random.rand(noP,nVar) + lb
    Swarm_Particles_V = np.zeros((noP, nVar))
    Swarm_Particles_PBEST_X = np.zeros((noP,nVar))
    Swarm_Particles_PBEST_O = np.ones(noP)*np.inf
    Swarm_Particles_O = np.ones(noP)

    Swarm_GBEST_X = np.zeros(nVar)
    Swarm_GBEST_O = np.inf
    t = 1

    a37, a28, a29, a38, a39, a40, a41, a42, a10, a11, a43, a45, a44, a46, Re_LF, Re_RF = constants([7, 0, 0, 0, -0.5, 0, 0, 0, 0, 0], [2.02, 0, 4780], 0, 1)

    # Main loop
    while t<=maxIter:
        # Calculate the objective value
        for k in range(0, noP):
            #u_z, u_phi = SMC_MODEL_func(Swarm_Particles_X[k,:], states, params)
            #Swarm_Particles_O[k] = (ref[0] - u_z)**2+(ref[1] - u_phi)**2
            Swarm_Particles_O[k]  = func_func(Swarm_Particles_X[k,:], a37, a28, a29, a38, a39, a40, a41, a42, a10, a11, a43, a45, a44, a46, Re_LF, Re_RF)
        
        index = np.nonzero(Swarm_Particles_O - Swarm_Particles_PBEST_O < 0)
 
        Swarm_Particles_PBEST_X[index,:] = Swarm_Particles_X[index,:]
        
        
        Swarm_Particles_PBEST_O[index] = Swarm_Particles_O[index]
        
        Omin = Swarm_Particles_O.min(0)

        Pmin = np.argmin(Swarm_Particles_O)
  
        aux = Swarm_GBEST_O
        if Omin < Swarm_GBEST_O:
            Swarm_GBEST_X = Swarm_Particles_X[Pmin,:]
            Swarm_GBEST_O = Omin
        
        # Update the X and V vectors
        w = wMax - t * ((wMax - wMin) / maxIter)
        
        Swarm_Particles_V = (w*Swarm_Particles_V + c1*np.random.rand(noP,nVar)*(Swarm_Particles_PBEST_X - Swarm_Particles_X) 
            + c2 * np.random.rand(noP,nVar) * (Swarm_GBEST_X - Swarm_Particles_X))
        
        # Check velocities
        
        np.putmask(Swarm_Particles_V, Swarm_Particles_V > vMax, vMax)
        np.putmask(Swarm_Particles_V, Swarm_Particles_V < vMin, vMin)

        Swarm_Particles_X = Swarm_Particles_X + Swarm_Particles_V
        

        # Check positions

        np.putmask(Swarm_Particles_X, Swarm_Particles_X > ub, ub)
        np.putmask(Swarm_Particles_X, Swarm_Particles_X < lb, lb)
        
        t = t + 1 
    
    return Swarm_GBEST_X[0], Swarm_GBEST_X[1], Swarm_GBEST_O

start = time.time()
print(PSO_func([0,0], [7,0,0,0,-0.5,0,0,0,0,0], [2.02, 0, 4780], np.array([14,14]), np.array([-3,-3])))
print(time.time()-start)

