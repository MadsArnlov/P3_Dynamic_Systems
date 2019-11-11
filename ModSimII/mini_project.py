# -*- coding: utf-8 -*-
"""
Created on Wed Nov  6 09:04:13 2019

@author: arnlo
"""

import numpy as np
import matplotlib.pyplot as plt
#from plotting import plot
from scipy.signal import savgol_filter
"""
data = np.array([[1,0],
                 [1,2],
                 [2,3]])

A = np.array([[1,1],
             [1,1],
             [1,2]])

b = np.array([0,2,3])


Q, R = np.linalg.qr(A)

b_hat = Q.T@b

x = np.linalg.solve(R,b_hat)

a = x[1]
b = x[0]

t = np.linspace(0, 2, 100)

y = lambda t: x[0] + x[1]*t

plot("{} + {}*x".format(b, a), 0, 2, 'x', 'y')

plt.plot(t, y(t))
plt.plot(data[0][0], data[0][1], 'kx', data[1][0], data[1][1], 'kx',
         data[2][0], data[2][1], 'kx')
plt.show()
"""


def midlingsfilter(x):
    newValue = sum(M) + x    
    for i in range(len(M)-1):
        M[i] = M[i+1]
    M[-1] = x
    return newValue/(len(M)+1)


file = np.loadtxt("data1_2019.txt", delimiter=",")
file = file[:16000]

m, g, l_p = 0.175, 9.82, 0.282
c = m*g*l_p
M = np.zeros(60)

samplingtime = 0.005
w = np.zeros(len(file))
dtw = np.zeros_like(w)

for k in range(len(file)-1):
    w[k] = (file[k+1][1] - file[k][1])/samplingtime

for i in range(len(w)):
    w[i] = midlingsfilter(w[i])

for k in range(len(file)-1):
    dtw[k] = (w[k+1] - w[k])/samplingtime

for i in range(len(dtw)):
    dtw[i] = midlingsfilter(dtw[i])

w = savgol_filter(w,151,3)
dtw = savgol_filter(dtw,151,3)

A = np.array([file[:, 0], -np.sign(w), -w, -dtw])
A = A.T
b = c*np.sin(file[:, 1])
Q, R = np.linalg.qr(A)
b_hat = Q.T @ b
x = np.linalg.solve(R, b_hat)
Km, f, alpha, J = x
residual = b - (A @ x)
orthogonal = [A[:, i] @ residual for i in range(4)]

plt.figure(figsize=(16, 9))
plt.subplot(3, 1, 1)
plt.plot(file[:, 1])
plt.subplot(3, 1, 2)
plt.plot(np.arcsin((A @ x)/c))
plt.subplot(3, 1, 3)
plt.plot(file[:, 1] - np.arcsin((A @ x)/c))
plt.show()

print(max(file[:, 1] - np.arcsin((A @ x)/k)))