import sympy as sym

Iyy = 503        # boat's moment of inertia in kg.m**2 for pitch
Ixx = 31         # moment of inertia for roll
Izz = 300
Ixz = 16.5

Ixx, Iyy, Izz, Ixz = sym.symbols('Ixx Iyy Izz, Ixz')
J = sym.Matrix([[Ixx, 0, Ixz], [0, Iyy, 0], [Ixz, 0, Izz]])
InvJ = J.inv()

params = [1,2,3,4,5]

print(params[0:])