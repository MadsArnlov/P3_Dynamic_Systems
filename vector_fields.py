# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 16:19:35 2019

@author: Lars
"""

import matplotlib.pyplot as plt
import numpy as np

# Interval - Axes
x = np.linspace(-2.5*np.pi, 2.5*np.pi, 31)
y = np.linspace(-2, 2, 16)

X, Y = np.meshgrid(x, y)

a = 0.0

dy = -np.sin(X) - a*Y
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
t_stop   = 10.0
N = 500
t_step = (t_stop - t_start)/float(N)

def fun(t, x):
    return np.array([x[1],
                      -1*np.sin(x[0]) - a*x[1]])

x10 = -np.pi*1.0
x20 = 0.0

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
plt.plot(X1, X2, 'k-', [x10], [x20], 'kd')
plt.savefig("direction_field.png")
plt.show()