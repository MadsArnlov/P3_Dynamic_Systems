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
M           = 5.273 #Mass of cart
m           = 0.250 #Mass of pendulum
l             = 0.334 #Lenght of pendulum arm
g            = -9.82 #Acceleration due to gravity
mu          = 0.052 #Friction Coefficient
F_c         = -g*mu #Coloumb Force
displ      = 0 #x-position to stabilise system

#Inital values
x_0         = 0 #Start postion of cart
theta_0     = np.pi/6 #Start angle of pendulum
x_dot_0     = 0 #Start velocity of cart
theta_dot_0 = 0 #Start angular velocity of pendulum

#RK4 parameters
t_start     = 0.0 #Start time
t_stop      = 10 #End time
N           = 5000 #Number of steps


# =============================================================================
# RK4 Implementation
# =============================================================================

t_step = (t_stop - t_start)/float(N) #Step size
t_arr = np.linspace(t_start, t_stop, N+1) #Time array for plotting

def rk4(f, t, x, h, K):
    k1 = f(t, x, K)
    k2 = f(t + 0.5*h, x + 0.5*h*k1, K)
    k3 = f(t + 0.5*h, x + 0.5*h*k2, K)
    k4 = f(t + h, x + h*k3, K)
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


#Function definition describing the system
def fun(t, z, K):

    z_dot = (A - B @ K) @ (z-np.array([displ,0,0,0]))
    return z_dot



def runrk4(K):
    t = 0.0
    for k in range(N):
        Xp, t = rk4(fun, t, np.array([X1[k], X2[k], X3[k], X4[k]]), t_step, K)
        X1[k+1] = Xp[0]
        X2[k+1] = Xp[1]
        X3[k+1] = Xp[2]
        X4[k+1] = Xp[3]

# =============================================================================
# Plotting setup
# =============================================================================

A = np.array([[0,0,1,0],[0,0,0,1],[0,m*g/M,F_c/M,0],[0,g*(m+M)/(l*M),F_c/(l*M),0]])
B = np.array([[0,0,-1/M,-1/(l*M)]]).T
name = []
x_label = ""
y_label = ""
plt.figure(figsize=(16,8))

def plotfig(P, x_var, y_var):
    for i in P:
        name.append(i)
        K = signal.place_poles(A,B,np.array(i)).gain_matrix
        print(K)
        runrk4(K)
        plt.plot(x_var, y_var)
        global x_label
        global y_label
    
    if np.array_equal(x_var, t_arr):
        x_label = "Time [s]"
    elif np.array_equal(x_var, cart_pos):
        x_label = "Position [m]"
    elif np.array_equal(x_var, pend_ang):
        x_label = "Angle [rads]"
    elif np.array_equal(x_var, cart_vel):
        x_label = "Velocity [m/s]"
    elif np.array_equal(x_var, pend_vel):
        x_label = "Angular Velocity [rads/s]"
    
    if np.array_equal(y_var, t_arr):
        y_label = "Time [s]"
    elif np.array_equal(y_var, cart_pos):
        y_label = "Position [m]"
    elif np.array_equal(y_var, pend_ang):
        y_label = "Angle [rads]"
    elif np.array_equal(y_var, cart_vel):
        y_label = "Velocity [m/s]"
    elif np.array_equal(y_var, pend_vel):
        y_label = "Angular Velocity [rads/s]"
    

cart_pos = X1 #Cart positions
pend_ang = X2 #Pendulum angles
cart_vel = X3 #Cart velocities
pend_vel = X4 #Pendulum velocities

# Eigenvalues
P = [[-4,-5,-6,-7]]
# =============================================================================
# P = [[-1,-2,-3,-1300],
#      [-2,-4,-3,-1300]]
# =============================================================================
# =============================================================================
# P = [
#      [-1, -2, -3, -4],
#      [-4, -5, -6, -7],
#      [-8, -9, -10, -11], #Umildbart den bedste
#      [-12, -13, -15, -16],
#      [-17, -18, -19, -20],
#      [-21, -22, -23, -24]
#      ]
# 
# =============================================================================
plotfig(P, t_arr, cart_pos)
#plt.axis([0,11,-10,0])
plt.axhline(y=0, color="r")
plt.axvline(x=0, color="r")
plt.legend(name, fontsize="xx-large")
plt.xlabel(x_label, fontsize = "xx-large")
plt.ylabel(y_label, fontsize = "xx-large")
plt.tick_params(labelsize = "xx-large")
plt.show()





