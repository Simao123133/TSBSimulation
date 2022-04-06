import control as ct
from Boat_Model import boat_update_trim
from tabulate import tabulate

def print_ss_matrix(m, name, states, inputs):

    if(name == "A"):
        print("\nA matrix:\n")
        print(tabulate(m, states, tablefmt="fancy_grid", showindex=states, numalign="center", floatfmt=".2f"))
    elif(name == "B"):
        print("\nB matrix:\n")
        print(tabulate(m, inputs, tablefmt="fancy_grid", showindex=states, numalign="center", floatfmt=".2f"))
    else:
        print("Select name A or B")

def trim_boat(speed, height):

    states = ["u", "w", "q", "Pitch", "z", "v", "p", "r", "Roll", "Yaw"]
    inputs = ["LF", "RF", "ReF", "Rpms", "Rudder"]
    outputs = states

    io_sys = ct.NonlinearIOSystem(updfcn = boat_update_trim, outfcn = None, inputs=inputs, outputs=outputs, states=states, name = "Boat")

    x0 = [speed, 0, 0, 0, height, 0, 0, 0, 0, 0]
    u0 = [1.91, 1.91, 1.6, 4, 0]
    y0 = [speed, 0, 0, 0, height, 0, 0, 0, 0, 0]

    iu = [4]
    iy = [0, 3, 4, 8]

    xeq, ueq = ct.find_eqpt(io_sys, x0=x0, y0=y0, u0=u0, iu = iu, iy=iy)

    if type(xeq) == type(None):
        print("Could not find equilibrium point")
    else:
        #Print the states
        print("\nEquilibrium states:\n")
        for idx, val in enumerate(xeq):
            if idx > 0 and idx < 9:
                print(states[idx] + ": " + "{:.2e}".format(val))
            else:
                print(states[idx] + ": " + str(val))

        #Print the inputs
        print("\nEquilibrium inputs:\n")
        for idx, val in enumerate(ueq):
            print(inputs[idx] + ": " + str(round(val, 2)))

    return io_sys, xeq, ueq

def linearize_boat(io_sys, xeq, ueq):

    states = ["u [m/s]", "w [m/s]", "q [º/s]", "Pitch [º]", "z [m]", "v [m/s]", "p [º/s]", "r [º/s]", "Roll [º]", "Yaw [º]"]
    inputs = ["Left foil [º]", "Right foil [º]", "Rear foil [º]", "Motor Rpms [10^3 Rpms]", "Rudder [º]"]
    outputs = states

    ss_sys = io_sys.linearize(xeq, ueq)

    print_ss_matrix(ss_sys.A, "A", states, inputs)

    print_ss_matrix(ss_sys.B, "B", states, inputs)

    return ss_sys.A, ss_sys.B


