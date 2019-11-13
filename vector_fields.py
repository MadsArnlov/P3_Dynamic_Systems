# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 16:19:35 2019

@author: Lars
"""

import matplotlib.pyplot as plt
import numpy as np

# Interval - Axes
x = np.linspace(-4, 4, 31)
y = np.linspace(-3, 3, 16)

X, Y = np.meshgrid(x, y)

g = 9.82
M = 1.1
m = 0.175
l_p = 0.282

dy = (g*M + g*m)/(M*l_p)*X
dx = Y
L = (dx**2 + dy**2)**0.5
dx = dx/L
dy = dy/L

# one step RK4 adapted to system
def rk4(f, t, x, h):
    k1 = f(t, x)
    k2 = f(t + 0.5*h, x + 0.5*h*k1)
    k3 = f(t + 0.5*h, x + 0.5*h*k2)
    k4 = f(t + h, x + h*k3)
    xp = x + h*(k1 + 2.0*(k2 + k3) + k4)/6.0
    return xp, t + h

t_start = 0.0
t_stop   = 1.0
N = 500
t_step = (t_stop - t_start)/float(N)

def fun(t, x):
    return np.array([x[1],
                      (g*(M+m))/(M*l_p)*x[0]])

x10 = 0#-0.52605243180
x20 = 0#1.6087558040

X1 = np.zeros(N + 1)
X2 = np.zeros(N + 1)

X1[0] = x10
X2[0] = x20

# Time variable
t = 0.0

for k in range(N):
    Xp, t = rk4(fun, t, np.array([X1[k], X2[k]]), t_step)
    X1[k+1] = Xp[0]
    X2[k+1] = Xp[1]

plt.figure(figsize=(16, 9))
plt.quiver(X, Y, dx, dy, color='red', headlength=5)
#plt.quiver(Z3, Z4, dz3, dz4, color='red', headlength=5)
plt.plot(X1, X2, 'k-', [x10], [x20], 'kd')
plt.savefig("direction_field.png")
plt.show()