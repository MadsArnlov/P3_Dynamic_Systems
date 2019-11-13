# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt


def midlingsfilter(x):
    newValue = sum(M) + x
    for i in range(len(M)-1):
        M[i] = M[i+1]
    M[-1] = x
    return newValue/(len(M)+1)


file = np.loadtxt("data1_2019.txt", delimiter=",")
file = file[:13490]

m, g, l_p = 0.175, 9.82, 0.282
c = m*g*l_p
M = np.zeros(60)

t_s = 0.00667
w = np.zeros(len(file))
dtw = np.zeros_like(w)

for k in range(len(file)-1):
    w[k] = (file[k+1][1] - file[k][1])/t_s

for i in range(len(w)):
    w[i] = midlingsfilter(w[i])

for k in range(len(file)-1):
    dtw[k] = (w[k+1] - w[k])/t_s

for i in range(len(dtw)):
    dtw[i] = midlingsfilter(dtw[i])

A = np.array([file[:, 0], -np.sign(w), -w, -dtw])
A = A.T
b = c*np.sin(file[:, 1])
Q, R = np.linalg.qr(A)
b_hat = Q.T @ b
x = np.linalg.solve(R, b_hat)
Km, f, alpha, J = x
residual = b - (A @ x)

t = np.linspace(0, len(file)*t_s, len(file))

plt.figure(figsize=(16, 4.5))
plt.plot(t, file[:, 1], 'b-', label="Measured")
plt.plot(t, np.arcsin((A @ x)/c), 'orange', label="Simulated")
plt.ylabel("Angle [radians]", fontsize=14)
plt.xlabel("Time [s]", fontsize=14)
plt.legend()
plt.show()
