from Trim_boat import trim_boat, linearize_boat, print_ss_matrix
import numpy as np
from ReadData import get_data, get_signal
from scipy.io import loadmat 
import matplotlib.pyplot as plt

def cosd(ang):
    return np.cos(ang/180*np.pi)

def sind(ang):
    return np.sin(ang/180*np.pi)

def tand(ang):
    return np.tan(ang/180*np.pi)

sys, xeq, ueq = trim_boat(7, -0.4)

A, B = linearize_boat(sys, xeq, ueq)

B = np.delete(B, [2, 4], axis = 1)

states = ["u", "w", "q", "Pitch", "z", "v", "p", "r", "Roll", "yaw"]
inputs = ["LF", "RF", "Rpms"]

print_ss_matrix(B, "B", states, inputs)

data = loadmat("mont.mat")

Velx = get_signal("_velX_", data)
Vely = get_signal("_velY_", data)
Velz = get_signal("_velZ_", data)
Rollv = get_signal("_roll_01", data)
Pitchv = get_signal("_pitch_00", data)
Yawv = get_signal("_yaw_02", data)
z = get_signal("_dist_01", data)
aa_LFv = get_signal("_left_angle_02", data)
aa_RFv = get_signal("_right_angle_03", data)
Rpmsv = get_signal("_M1_RPM_01", data)

start_time =  11327.413
end_time = 11359.066

z_idx = np.argmin(np.abs(start_time - z.time))
z_idx_f = np.argmin(np.abs(end_time - z.time))

time = z.time[z_idx]
idx = 0

Vel_x = np.zeros(z_idx_f - z_idx)
Vel_y = np.zeros(z_idx_f - z_idx)
Vel_z = np.zeros(z_idx_f - z_idx)
us = np.zeros(z_idx_f - z_idx)
vs = np.zeros(z_idx_f - z_idx)
ws = np.zeros(z_idx_f - z_idx)
dRolls = np.zeros(z_idx_f - z_idx)
dPitchs = np.zeros(z_idx_f - z_idx)
dYaws = np.zeros(z_idx_f - z_idx)
Rolls = np.zeros(z_idx_f - z_idx)
Pitchs = np.zeros(z_idx_f - z_idx)
Yaws = np.zeros(z_idx_f - z_idx)
zs = np.zeros(z_idx_f - z_idx)
aa_LFs = np.zeros(z_idx_f - z_idx)
aa_RFs = np.zeros(z_idx_f - z_idx)
Rpmss = np.zeros(z_idx_f - z_idx)

while time <= end_time - 0.1:

    time = z.time[z_idx]
    V_idx = np.argmin(np.abs(time - Velx.time))
    eul_idx = np.argmin(np.abs(time - Rollv.time))

    Roll = Rollv.value[eul_idx]
    Pitch = Pitchv.value[eul_idx]
    Yaw = Yawv.value[eul_idx]

    vx = Velx.value[V_idx]
    vy = Vely.value[V_idx]
    vz = Velz.value[V_idx]

    Vel_x[idx] = vx
    Vel_y[idx] = vy
    Vel_z[idx] = vz

    # Inertial to body frame matrix
    R = np.array([[cosd(Pitch)*cosd(Yaw), cosd(Pitch)*sind(Yaw), -sind(Pitch)],
        [sind(Roll)*sind(Pitch)*cosd(Yaw)-cosd(Roll)*sind(Yaw), sind(Roll)*sind(Pitch)*sind(Yaw)+cosd(Roll)*cosd(Yaw), sind(Roll)*cosd(Pitch)],
        [cosd(Roll)*sind(Pitch)*cosd(Yaw)+sind(Roll)*sind(Yaw), cosd(Roll)*sind(Pitch)*sind(Yaw)-sind(Roll)*cosd(Yaw), cosd(Roll)*cosd(Pitch)]])

    Vb = R@np.array([vx, vy, vz])

    us[idx] = Vb[0]
    vs[idx] = Vb[1]
    ws[idx] = Vb[2]

    Rolls[idx] = Roll
    Pitchs[idx] = Pitch
    Yaws[idx] = Yaw

    zs[idx] = z.value[z_idx]

    Rpms_idx = np.argmin(np.abs(time - Rpmsv.time))
    Rpmss[idx] = Rpmsv.value[Rpms_idx]

    aa_LF_idx = np.argmin(np.abs(time - aa_LFv.time))
    aa_LFs[idx] = aa_LFv.value[aa_LF_idx]

    aa_RF_idx = np.argmin(np.abs(time - aa_RFv.time))
    aa_RFs[idx] = aa_LFv.value[aa_RF_idx]

    idx = idx + 1
    z_idx = z_idx + 1

dRolls[0:-1] = np.diff(Rolls)/0.1
dPitchs[0:-1] = np.diff(Pitchs)/0.1
dYaws[0:-1] = np.diff(Yaws)/0.1

dRolls[-1] = (Rollv.value[eul_idx+1] - Rollv.value[eul_idx])/0.1
dPitchs[-1] = (Pitchv.value[eul_idx+1] - Pitchv.value[eul_idx])/0.1
dYaws[-1] = (Yawv.value[eul_idx+1] - Yawv.value[eul_idx])/0.1

plt.plot(aa_LFs)
plt.show()


