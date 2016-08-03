import sys
sys.path.append('/home/tim/Dropbox/EigenDocumenten/Doctoraat/SitControl/Python/CasADi3.0/') #add spline toolbox to path 
from spline import *
import spline_extra as spl_extra
import numpy as np
import casadi as cas
import matplotlib.pyplot as plt
import time
plt.interactive(True)

from AMPLexport import AMPLexport

# Find the shortest path from start to end, avoid circular obstacle in the middle
# There is no velocity or time info which is important here, we just look for the shortest path
# We will test what the influence of the initial guess is and look what different solvers give us.

veh_r =  0.1  #enclosing circle of vehicle

pos_veh0 = [-2. , -2.]    
goal     = [2. , 2.]

Nm_c = 1 #amount of circular obstacles
pos_obj0_c   = {0:np.array([[0],[0.1]])} #midpoint
vel_obj0_c   = {0:np.array([[0],[0]])}     #velocity
r_obj        = {0:0.5}

#Borders
xU =  3
xL = -3
yU =  3
yL = -3

# Parameterization
deg   = 3
n     = 2
knots = np.r_[np.zeros(deg), np.linspace(0, 1, n), np.ones(deg)]
B     = BSplineBasis(knots, deg)       #basis for motion trajectory
Blin  = BSplineBasis([0, 0, 1, 1], 1); #basis for prediction

# The optimization variables
X = cas.MX.sym('x', 2 * len(B))
cx = X[:len(B)]
cy = X[len(B):2*len(B)]

# The parameters
P = cas.MX.sym('p', 4 + 4*Nm_c)
P0_veh = P[0:2]
V0_veh = P[2:4]  # True velocity, remember to divide by T!
P0_obj_c = []
V0_obj_c = []
for i in range(Nm_c):
    P0_obj_c.append(P[4+i*4:6+i*4])
    V0_obj_c.append(P[6+i*4:8+i*4])

# initial parameter values
vel_veh0 = [0,0] #initially at rest
p0 = np.array([pos_veh0[0], pos_veh0[1], vel_veh0[0], vel_veh0[1]])
for i in range(Nm_c):
    p0 = np.append(p0, pos_obj0_c[i])
    p0 = np.append(p0, vel_obj0_c[i])

# The relevant splines
x =   BSpline(B, cx)
y =   BSpline(B, cy)
dx =  x.derivative(1)
dy =  y.derivative(1)
ddx = x.derivative(2)
ddy = y.derivative(2)

#Constraints
con = []
#-(ddx**2 + ddy**2) + (amax * T**2)**2,
#-(dx**2 + dy**2)   + (vmax * T)**2] #acceleration and velocity constraints

#Add circle constraints
for i in range(Nm_c):    
    # Velocity is supposed to be 0
    obj_pred = [BSpline(Blin, cas.vertcat(P0_obj_c[i][0], P0_obj_c[i][0])),
                BSpline(Blin, cas.vertcat(P0_obj_c[i][1], P0_obj_c[i][1]))] #make prediction for obstacle i    
    con.extend([(x - 0)*(x - 0) + (y - 0)*(y - 0) - (r_obj[i] + veh_r)**2])
    
#Add border constraints  
# con.extend([ xU - veh_r - x])
# con.extend([-xL - veh_r + x])
# con.extend([ yU - veh_r - y])
# con.extend([-yL - veh_r + y])
    
#Translate to constraints on the coeffs
con = cas.vertcat(*[c.coeffs for c in con])

#Equality constraints
con = cas.vertcat(con, x(0) - P0_veh[0], y(0) - P0_veh[1],
                  x(1), y(1)) #position constraints

# Shortest path objective
obj = spl_extra.definite_integral(dx*dx, 0, 1) + spl_extra.definite_integral(dy*dy, 0, 1)
#obj=T

nlp = {'x':X, 'p':P, 'f':obj, 'g':con}
options = {}
options['ipopt.linear_solver'] = 'ma57'
options['ipopt.tol'] = 1e-4
options['ipopt.warm_start_init_point'] = 'yes'
options['ipopt.max_iter'] = 3000
solver = cas.nlpsol('solver','ipopt', nlp, options)

lbg = np.r_[np.zeros(con.size1()-4),  0, 0, goal[0], goal[1]]
ubg = np.r_[np.inf * np.ones(con.size1()-4),  0, 0, goal[0], goal[1]]
lbx = -100 * np.ones(X.size1())
ubx = 100 * np.ones(X.size1())

solver_input = {}
solver_input['lbx'] = lbx
solver_input['ubx'] = ubx
solver_input['lbg'] = lbg
solver_input['ubg'] = ubg

plt.figure()
plt.axis(1.1*np.array([xL, xU, yL, yU]))#('auto')
plt.hold(True)
plt.show()
plt.plot(goal[0], goal[1], 'rx')
plt.plot([xU , xU , xL , xL , xU],[yL , yU , yU , yL , yL],'r-') #plot border

s = np.linspace(0, 1, 101) 

Ts = 0.05 #0.025 #sample time
x0 = np.array(np.zeros(X.size1())) #initial guess
# Adapt initial guess here, study several cases
#x0 = np.array([-2, 0.5, 1., 2, -2, 1.07805, 0.627224, 2])

solver_input['x0'] = x0
solver_input['p'] = p0

t_solve = []
# Do the loop
running = True

AMPLexport(nlp,solver_input)

bla

while running:
    solver_input['x0'] = x0
    solver_input['p'] = p0
    solver_output = solver(**solver_input)

    t_solve.append(solver.stats()['t_proc_mainloop'])

    X = solver_output['x']
    cx = X[:len(B)]
    cy = X[len(B):2*len(B)]
    print cx
    print cy
    x = BSpline(B, cx)
    y = BSpline(B, cy)
    dx = x.derivative(1)
    dy = y.derivative(1)
    Tprev = p0[-1] #get Tprev out of initial parameters
    obstacles_c = []
    for i in range(Nm_c):
        obstacles_c.append(plt.plot(r_obj[i] * np.sin(2 * np.pi * s) + pos_obj0_c[i][0], r_obj[i] * np.cos(2 * np.pi * s) + pos_obj0_c[i][1], 'b'))
    vehicle_circle = plt.plot(veh_r * np.sin(2 * np.pi * s) + p0[0], veh_r * np.cos(2 * np.pi * s) + p0[1], 'g')
    plt.plot(p0[0], p0[1], 'kx')
    vehicle_path = plt.plot(x(s), y(s), 'g')
    
    plt.draw()
    time.sleep(10)
    running = False
