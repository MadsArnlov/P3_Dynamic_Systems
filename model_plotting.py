# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 11:01:24 2019

@author: Andreas
"""

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
displ       = 0.445 #x-position to stabilise system

#Inital values
x_0         = 0.445 #Start postion of cart
theta_0     = 0.14765 #Start angle of pendulum
x_dot_0     = 0 #Start velocity of cart
theta_dot_0 = 0 #Start angular velocity of pendulum

#RK4 parameters
t_start     = 0.0 #Start time
t_stop      = 5 #End time
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
plt.figure(figsize=(16,12))

def plotfig(P):
    for i in P:
        name.append(i)
        K = signal.place_poles(A,B,np.array(i)).gain_matrix
        runrk4(K)
        plt.subplot(2,1,1)
        plt.xlabel("Time [s]", fontsize=14)
        plt.ylabel("Position [m]", fontsize=14)
        plt.plot(t_arr, cart_pos)
        plt.ylim(0,0.9) #Position
        plt.tick_params(labelsize = "xx-large")
        plt.legend(name, fontsize="xx-large", loc="upper right")
        plt.xlim(0,5)
        plt.subplot(2,1,2)
        plt.xlabel("Time [s]", fontsize=14)
        plt.ylabel("Angle [rad]", fontsize=14)
        plt.plot(t_arr, pend_ang)
        plt.ylim(-0.15,0.15) #Angle
        plt.tick_params(labelsize = "xx-large")
        plt.legend(name, fontsize="xx-large", loc="upper right")
    

cart_pos = X1 #Cart positions
pend_ang = X2 #Pendulum angles
cart_vel = X3 #Cart velocities
pend_vel = X4 #Pendulum velocities

#Generate Eigenval
P1p = []
P2p = []
P3p = []

P1s = []
P2s = []
P3s = []
s = 1
for i in range(0,4):
    p = i
    P1p.append([-3*s-p,-4*s-p,-5*s-p,-6*s-p])
    P2p.append([-3*s-p,-5*s-p,-7*s-p,-9*s-p])
    P3p.append([-3*s-p,-6*s-p,-9*s-p,-12*s-p])
    
for j in range(0,4):
    p = 0
    P1s.append([-3*s-p,-4*s-p,-5*s-p,-6*s-p])
    P2s.append([-3*s-p,-5*s-p,-7*s-p,-9*s-p])
    P3s.append([-3*s-p,-6*s-p,-9*s-p,-12*s-p])
    s += 0.5

# Eigenvalues
plotfig(P3s)
#plt.savefig("Figures_of_models/Model_Ref3_Scale.png")





