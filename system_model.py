"""
Created for P3-Project -- Aalborg University

Authors: MATTEK 4.211
"""

import numpy as np
from scipy import signal
from matplotlib import pyplot as plt


# =============================================================================
# Constants
# =============================================================================

#Physical values
M           = 6.28 #Mass of cart
m           = 0.200 #Mass of pendulum
l           = 0.3235 #Lenght of pendulum arm
g           = 9.82 #Acceleration due to gravity
mu          = 0.052 #Friction Coefficient
F_c         = -(g*M)*mu #Coloumb Force
displ       = 0 #x-position to stabilise system

#Inital values
x_0         = 0 #Start postion of cart
theta_0     = 0.23 #Start angle of pendulum
x_dot_0     = 0 #Start velocity of cart
theta_dot_0 = 0 #Start angular velocity of pendulum

#RK4 parameters
t_start     = 0.0 #Start time
t_stop      = 10 #End time
N           = 50000 #Number of steps

#Stabilisation parameters
s = 0
v = 1
P           = [-5*v-s,-12*v-s,-13*v-s,-14*v-s] #Eigenvalues used for stabilisation
#P           = [-0.97592931, -0.74017246, -0.88919959, -0.69022008]

# =============================================================================
# RK4 Implementation
# =============================================================================

t_step = (t_stop - t_start)/float(N) #Step size
t_arr = np.linspace(t_start, t_stop, N+1) #Time array for plotting

def rk4(f, t, x, h):
    k1 = f(t, x)
    k2 = f(t + 0.5*h, x + 0.5*h*k1)
    k3 = f(t + 0.5*h, x + 0.5*h*k2)
    k4 = f(t + h, x + h*k3)
    xp = x + h*(k1 + 2.0*(k2 + k3) + k4)/6.0
    return xp, t + h

#Empty arrays to be filled
X1 = np.zeros(N + 1) #Cart positions
X2 = np.zeros(N + 1) #Pendulum angles
X3 = np.zeros(N + 1) #Cart velocities
X4 = np.zeros(N + 1) #Pendulum velocities

#Inserting initial values into arrays
X1[0] = x_0
X2[0] = theta_0
X3[0] = x_dot_0
X4[0] = theta_dot_0


# =============================================================================
# System Definition
# =============================================================================

#Defining the matricies of the system
A = np.array([[0,0,1,0],[0,0,0,1],[0,m*g/M,F_c/M,0],[0,g*(m+M)/(l*M),F_c/(l*M),0]])
B = np.array([[0,0,-1/M,-1/(l*M)]]).T

#Calculating the gain matrix
K = signal.place_poles(A,B,np.array(P)).gain_matrix
print("Lambda =",P)
print("K = (",K[0][0], ",", K[0][1], ",", K[0][2], ",", K[0][3], ")")

#Function definition describing the system
def fun(t, z):
    z_dot = (A - B @ K) @ (z-np.array([displ,0,0,0]))
    return z_dot
    
t = 0.0
for k in range(N):
    Xp, t = rk4(fun, t, np.array([X1[k], X2[k], X3[k], X4[k]]), t_step)
    X1[k+1] = Xp[0]
    X2[k+1] = Xp[1]
    X3[k+1] = Xp[2]
    X4[k+1] = Xp[3]

# =============================================================================
# Plotting setup
# =============================================================================

cart_pos = X1 #Cart positions
pend_ang = X2 #Pendulum angles
cart_vel = X3 #Cart velocities
pend_vel = X4 #Pendulum velocities

#The variables to put on each axis
x_vars = [t_arr, t_arr]
y_vars = [cart_pos, pend_ang]
# =============================================================================
# 
# plt.figure(figsize=(16, 4*len(y_vars)))
# for i in range(len(x_vars)):
#     
#     nr_plots = len(x_vars)
#     x_var = x_vars[i]
#     y_var = y_vars[i]
#     x_label = ""
#     y_label = ""
#     
#     if np.array_equal(x_var, t_arr):
#         x_label = "Time [s]"
#     elif np.array_equal(x_var, cart_pos):
#         x_label = "Position [m]"
#     elif np.array_equal(x_var, pend_ang):
#         x_label = "Angle [rads]"
#     elif np.array_equal(x_var, cart_vel):
#         x_label = "Velocity [m/s]"
#     elif np.array_equal(x_var, pend_vel):
#         x_label = "Angular Velocity [rads/s]"
#     
#     if np.array_equal(y_var, t_arr):
#         y_label = "Time [s]"
#     elif np.array_equal(y_var, cart_pos):
#         y_label = "Position [m]"
#     elif np.array_equal(y_var, pend_ang):
#         y_label = "Angle [rads]"
#     elif np.array_equal(y_var, cart_vel):
#         y_label = "Velocity [m/s]"
#     elif np.array_equal(y_var, pend_vel):
#         y_label = "Angular Velocity [rads/s]"
# 
#     plt.subplot(nr_plots, 1, i+1, xlabel=x_label, ylabel=y_label)
#     plt.plot(x_var, y_var)
#     plt.axhline(y=0, color="r")
#     plt.axvline(x=0, color ="r")
#     
#     if i < 1:
#         plt.title("$\lambda$ = {}".format(P), fontsize="xx-large")
#     
#     plt.grid()
# plt.show()
# =============================================================================
#plt.savefig("$\lambda$ = {}.png".format(P))









