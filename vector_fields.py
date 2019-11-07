# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 16:19:35 2019

@author: Lars
"""

import matplotlib.pyplot as plt
import numpy as np

# Interval - Axes
x = np.linspace(-4, 4, 31)
y = np.linspace(-10, 10, 16)

#z1 = np.linspace(0, 5, 31)
#z2 = np.linspace(-2.5*np.pi, 2.5*np.pi, 31)
#z3 = np.linspace(-2, 2, 16)
#z4 = np.linspace(-2*np.pi, 2*np.pi, 16)

X, Y = np.meshgrid(x, y)
#Z1, Z2, Z3, Z4 = np.meshgrid(z1, z2, z3, z4)
a = 0.0

g, M, m, l_p = 9.82, 10, 7, 1.5

dy = (g*M)/(m*l_p)*X
dx = Y
L = (dx**2 + dy**2)**0.5
dx = dx/L
dy = dy/L

#m, M, l, g = 0.5, 2, 1, 9.82
#
#dz1 = Z3
#dz2 = Z4
#dz3 = (m*g*np.cos(Z2)*np.sin(Z2)-m*l*np.sin(Z2)*Z4**2)/(m+M - m*np.cos(Z2)**2)
#dz4 = (m*g*np.cos(Z2)**2*np.sin(Z2)-m*l*np.cos(Z2)*np.sin(Z2)*Z4**2)/(m*l+M*l-m*l*np.cos(Z2)**2)+g/l*np.sin(Z2)
#L = (dz1**2 + dz2**2 + dz3**2 + dz4**2)**0.5
#dz1 = dz1/L
#dz2 = dz2/L
#dz3 = dz3/L
#dz4 = dz4/L


# one step RK4 adapted to system
def rk4(f, t, x, h):
    k1 = f(t, x)
    k2 = f(t + 0.5*h, x + 0.5*h*k1)
    k3 = f(t + 0.5*h, x + 0.5*h*k2)
    k4 = f(t + h, x + h*k3)
    xp = x + h*(k1 + 2.0*(k2 + k3) + k4)/6.0
    return xp, t + h

t_start = 0.0
t_stop   = 2.0
N = 500
t_step = (t_stop - t_start)/float(N)

def fun(t, x):
    return np.array([x[1],
                      (g*M)/(m*l_p)*x[0]])

x10 = -0.52605243180*4
x20 = 1.6087558040*4

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