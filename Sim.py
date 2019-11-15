# -*- coding: utf-8 -*-
"""
Created on Fri Nov 15 13:28:47 2019

@author: Andreas
"""

import numpy as np
from scipy import signal
from matplotlib import pyplot as plt

M = 5
m = 0.251
l = 0.334
k = 0.003
g = 9.82

A = np.array([[0,0,1,0],[0,0,0,1],[0,m*g/M,k/M,0],[0,g*(m+M)/(l*M),k/(l*M),0]])
B = np.array([[0,0,-1/M,-1/(l*M)]]).T

P = np.array([-1,-2,-3,-4])

K = signal.place_poles(A,B,P).gain_matrix



#z = np.array([x, theta, x_dot, theta_dot])

#z_dot = (A - B @ K) @ z

# =============================================================================
# RK4 Implementation
# =============================================================================

def rk4(f, t, x, h):
    k1 = f(t, x)
    k2 = f(t + 0.5*h, x + 0.5*h*k1)
    k3 = f(t + 0.5*h, x + 0.5*h*k2)
    k4 = f(t + h, x + h*k3)
    xp = x + h*(k1 + 2.0*(k2 + k3) + k4)/6.0
    return xp, t + h

t_start = 0.0
t_stop   = 20
N = 5000
t_step = (t_stop - t_start)/float(N)

def fun(t, z):
    z_dot = (A - B @ K) @ z
    return z_dot
    
x_0 = 0
theta_0 = np.pi/3
x_dot_0 = 0
theta_dot_0 = 0

X1 = np.zeros(N + 1)
X2 = np.zeros(N + 1)
X3 = np.zeros(N + 1)
X4 = np.zeros(N + 1)

X1[0] = x_0
X2[0] = theta_0
X3[0] = x_dot_0
X4[0] = theta_dot_0

# Time variable
t = 0.0

for k in range(N):
    Xp, t = rk4(fun, t, np.array([X1[k], X2[k], X3[k], X4[k]]), t_step)
    X1[k+1] = Xp[0]
    X2[k+1] = Xp[1]
    X3[k+1] = Xp[2]
    X4[k+1] = Xp[3]

plt.figure(figsize=(16, 9))
# =============================================================================
# plt.subplot(3,1,1)
# #plt.plot(X1, X3, 'k-', [x_0], [x_dot_0], 'kd')
# #plt.plot(X2, X4, 'k-', [theta_0], [theta_dot_0], 'kd')
# plt.plot(X2, X1, 'k-', [theta_0], [x_0], 'kd')
# plt.ylabel("Position", fontsize=14)
# plt.xlabel("Vinkel", fontsize=14)
# plt.yticks([-np.pi/2,0,np.pi/2])
# 
# plt.subplot(3,1,2)
# plt.plot(X3, X4, 'k-', [x_dot_0], [theta_dot_0], 'kd')
# plt.xlabel("Hastighed", fontsize=14)
# plt.ylabel("Vinkelhasighed", fontsize=14)
# 
# =============================================================================
plt.subplot(2,1,1)
plt.plot(X1)
plt.subplot(2,1,2)
plt.plot(X2)















