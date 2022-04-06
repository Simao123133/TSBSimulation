import numpy as np

class SMC_controller(object):
    def __init__(self, c, d, rho, dt, init):
        self.c = c
        self.d = d
        self.rho = rho
        self.x = 0
        self.xdot = 0
        self.u = 0
        self.dt = dt
        self.integral = init

    def compute(self, target, x, xdot):
        self.error = target - x
        self.error_dot = - xdot
        self.u_eq = self.c*self.error_dot
        self.sigma = self.c*self.error + self.error_dot
        self.integral += self.rho*np.sign(self.sigma)*self.dt
        self.u_sw = self.d*np.sqrt(np.absolute(self.sigma))*np.sign(self.sigma) + self.integral
        self.u = self.u_eq + self.u_sw

        return self.u

class PID_controller(object):
    def __init__(self, Kp, Kd, Ki, init_state, dt):
        self.Kp = Kp
        self.Kd = Kd
        self.Ki = Ki
        self.setpoint = 0
        self.error = 0
        self.u = 0
        self.dt = dt
        self.error_int = init_state/Ki

    def compute(self, target, x, xdot):
        self.error = target - x
        self.error_dot = -xdot
        self.error_int += self.error*self.dt
        self.Pu = self.error*self.Kp
        self.Du = self.error_dot*self.Kd
        self.Iu = self.error_int*self.Ki
        self.u = self.Pu + self.Du + self.Iu

        return self.u