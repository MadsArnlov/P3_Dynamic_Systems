# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 10:45:04 2019

@author: Andreas
"""

import numpy as np
import matplotlib.pyplot as plt


# Exercise 1

A = np.array([[4,0],[0,2],[1,1]])
b = np.array([2,0,11])

A_pse = np.linalg.pinv(A)

QR = A
Q, R = np.linalg.qr(QR) # Doing it manually we would use Gram-Schmitt

Qb = Q.T @ b

sol = np.linalg.solve(R, Qb)

x_points = []
y_points = []
for i in range(len(A[:,0])):
    x_points.append(A[i][0])
    y_points.append(A[i][1])

plt.plot(x_points, y_points, "kx",)

def f(x):
    y = sol[0]*x + sol[1]*x**2
    
    return y

x = np.linspace(0,5, 100)
plt.plot(x,f(x))


r = b - (A @ sol)

for i in range(len(A[:,0])-1):
    print("The dot product is {}".format(np.dot(r, A[:,i])))
plt.show()
    

# Exercise 2

# =============================================================================
# def f2(omega, N): # Omega = Vector, N = Scalar
#     
#     A1 = []
#     A2 = []
#     
#     for n in range(0, N):
#         row1 = []
#         row2 = []
#         
#         for l in range(0,len(omega)):
#             row1.append(np.cos(n*omega[l]))
#             row2.append(np.cos(n*omega[l]))
#             
#         A1.append(row1)
#         A2.append(row2)
#     
#     A1 = np.array(A1)      
#     A2 = np.array(A2)
#     
#     return A1, A2
# 
# 
# omega = np.array([0.1, 0.5, 0.6])
# 
# A1, A2 = f2(omega, 100)
# A = np.concatenate((A1, A2), axis = 1)
# 
# # Generating
# w = np.random.normal(0,1,100)
# theta = [1, 0, 0.5, 0, 0.15, 0.1]
# y = A @ theta + w
# 
# y = y - w
# 
# 
# 
# # Least Squares
# Q, R = np.linalg.qr(A)
# 
# Qy = Q.T @ y
# 
# sol = np.linalg.solve(R, Qy)
# 
# r = y - (A @ sol)
#  
# for i in range(len(A[0,:])-1):
#     print("The dot product is {}".format(np.dot(r, A[:,i])))
# 
# 
# x = np.arange(0,100)
# 
# def poly(x):
#     
#     newpoly = sol[0]*x + sol[1]*x**2 + sol[2]*x**3 + sol[3]*x**4 +sol[4]*x**5 + sol[5]*x**6
# 
#     return newpoly
# 
# plt.plot(x,y+w)
# #plt.plot(x,poly(x))
# =============================================================================
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        