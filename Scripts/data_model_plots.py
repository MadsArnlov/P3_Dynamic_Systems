# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 09:33:51 2019

Comparison of model and data

@author: Andreas
"""

import numpy as np
from scipy import signal
from matplotlib import pyplot as plt

# =============================================================================
# DATA
# =============================================================================

def loadData(eigen, change, state = "pos"):

    #Load data
    eig1 = eigen[0]
    eig2 = eigen[1]
    eig3 = eigen[2]
    eig4 = eigen[3]
    folder = "m_0.200_eig_{}_{}_{}_{}".format(eig1, eig2, eig3, eig4)
    if change == "displacement":
        Data_file1=np.loadtxt("Test_Recordings/"+folder+"/template.txt", delimiter=",")
        Data_file2=np.loadtxt("Test_Recordings/"+folder+"/displacement_-1.txt", delimiter=",")
        Data_file3=np.loadtxt("Test_Recordings/"+folder+"/displacement_-2.txt", delimiter=",")
        Data_file4=np.loadtxt("Test_Recordings/"+folder+"/displacement_-3.txt", delimiter=",")
    elif change == "scale":
        Data_file1=np.loadtxt("Test_Recordings/"+folder+"/template.txt", delimiter=",")
        Data_file2=np.loadtxt("Test_Recordings/"+folder+"/scale_1.5.txt", delimiter=",")
        Data_file3=np.loadtxt("Test_Recordings/"+folder+"/scale_2.txt", delimiter=",")
        Data_file4=np.loadtxt("Test_Recordings/"+folder+"/scale_2.5.txt", delimiter=",")
    
    #Data for user pushing cart
    Data_file_curve=np.loadtxt("Test_Recordings/curveball.txt",delimiter=",")
    
    #Filling the correct data into arrays
    if state == "pos":
        cart_pos1=Data_file1[:,0]
        cart_pos2=Data_file2[:,0]
        cart_pos3=Data_file3[:,0]
        cart_pos4=Data_file4[:,0]
        pos = [cart_pos1, cart_pos2, cart_pos3, cart_pos4]
    
        pend_ang1=Data_file1[:,1]
        pend_ang2=Data_file2[:,1]
        pend_ang3=Data_file3[:,1]
        pend_ang4=Data_file4[:,1]
        ang = [pend_ang1, pend_ang2, pend_ang3, pend_ang4]
        return pos, ang
    
    elif state == "vel":
        cart_vel1=Data_file1[:,2]
        cart_vel2=Data_file2[:,2]
        cart_vel3=Data_file3[:,2]
        cart_vel4=Data_file4[:,2]
        c_vel = [cart_vel1, cart_vel2, cart_vel3, cart_vel4]
    
        pend_ang_vel1=Data_file1[:,3]
        pend_ang_vel2=Data_file2[:,3]
        pend_ang_vel3=Data_file3[:,3]
        pend_ang_vel4=Data_file4[:,3]
        p_vel = [pend_ang_vel1, pend_ang_vel2, pend_ang_vel3, pend_ang_vel4]
        return c_vel, p_vel
    
    elif state == "curve":
        cart_pos_curve=Data_file_curve[:,0]
        pend_ang_curve=Data_file_curve[:,1]
        cart_vel_curve=Data_file_curve[:,2]
        pend_ang_vel_curve=Data_file_curve[:,3]
        return cart_pos_curve, pend_ang_curve, cart_vel_curve, pend_ang_vel_curve


# =============================================================================
# MODEL
# =============================================================================

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
stab        = 0.445 #x-position to stabilise system

#Inital values
x_0         = 0.445 #Start postion of cart
theta_0     = 0.14765 #Start angle of pendulum
x_dot_0     = 0 #Start velocity of cart
theta_dot_0 = 0 #Start angular velocity of pendulum

#RK4 parameters
t_start     = 0.0 #Start time
t_stop      = 20 #End time
N           = 5000 #Number of steps


# =============================================================================
# RK4 Implementation
# =============================================================================

t_step = (t_stop - t_start)/float(N) #Step size

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

def runRK4():
    t = 0.0
    for k in range(N):
        Xp, t = rk4(fun, t, np.array([X1[k], X2[k], X3[k], X4[k]]), t_step)
        X1[k+1] = Xp[0]
        X2[k+1] = Xp[1]
        X3[k+1] = Xp[2]
        X4[k+1] = Xp[3]

# =============================================================================
# System Definition
# =============================================================================

#Defining the matricies of the system
A = np.array([[0,0,1,0],[0,0,0,1],[0,m*g/M,F_c/M,0],[0,g*(m+M)/(l*M),F_c/(l*M),0]])
B = np.array([[0,0,-1/M,-1/(l*M)]]).T

#Function definition describing the system
def fun(t, z):
    z_dot = (A - B @ K) @ (z-np.array([stab,0,0,0]))
    return z_dot


# =============================================================================
# Plotting setup
# =============================================================================


##Change values here##
# =============================================================================
eigen_template = [-3,-6,-9,-12] #Graps both position and angle
s = 1 #Scale
p = -2 #Displacement
change = "displacement" #Displacement or Scale
# =============================================================================

#Calculating the gain matrix and running RK4
eigen_model = [eigen_template[0]*s+p,eigen_template[1]*s+p,eigen_template[2]*s+p,eigen_template[3]*s+p]
K = signal.place_poles(A,B,np.array(eigen_model)).gain_matrix
runRK4()
model_cart_pos = X1 #Cart positions
model_pend_ang = X2 #Pendulum angles
model_cart_vel = X3 #Cart velocities
model_pend_vel = X4 #Pendulum velocities

#Loading and cropping data
data1, data2 = loadData(eigen_template, change)
n1 = len(data1[0])
n2 = len(data1[1])
n3 = len(data1[2])
n4 = len(data1[3])

#Sampling
sampling_time=0.00667
sampling_frequency=1/sampling_time
N_data=max(n1,n2,n3,n4)
N_20=int(sampling_frequency*20)

#Creating time arrays for plotting
t_arr_model = np.linspace(t_start, t_stop, N+1)
t_arr_data = np.linspace(0,sampling_time*N_data,N_data)


#Plotting parameters
plt.figure(figsize=(16,12))
plt.subplot(2,1,1)
plt.xlim(0,20)
if change == "displacement":
    plt.title("Template: {}, Displacement: {}".format(eigen_template,p), fontsize=20)
elif change == "scale":
    plt.title("Template: {}, Scaling: {}".format(eigen_template,s), fontsize=20)
plt.plot(t_arr_model, model_cart_pos) #Model
if p == 0 and s == 1:
    plt.plot(t_arr_data[:min(n1,N_20)], data1[0][:min(n1,N_20)]) #Template
elif p == -1 or s == 1.5:
    plt.plot(t_arr_data[:min(n2,N_20)], data1[1][:min(n2,N_20)]) #Displacement +1 or Scale 1.5
elif p == -2 or s == 2:
    plt.plot(t_arr_data[:min(n3,N_20)], data1[2][:min(n3,N_20)]) #Displacement +2 or Scale 2
elif p == -3 or s == 2.5:
    plt.plot(t_arr_data[:min(n4,N_20)], data1[3][:min(n4,N_20)]) #Displacement +3 or Scale 2.5

plt.xlabel("Time [s]",fontsize=14)
plt.ylabel("Position [m]",fontsize=14)
plt.ylim(0,0.9)
plt.legend(["Model", "Data"])

plt.subplot(2,1,2)

plt.plot(t_arr_model, model_pend_ang) #Model
if p == 0 and s == 1:
    plt.plot(t_arr_data[:min(n1,N_20)], data2[0][:min(n1,N_20)]) #Template
elif p == -1 or s == 1.5:
    plt.plot(t_arr_data[:min(n2,N_20)], data2[1][:min(n2,N_20)]) #Displacement +1 or Scale 1.5
elif p == -2 or s == 2:
    plt.plot(t_arr_data[:min(n3,N_20)], data2[2][:min(n3,N_20)]) #Displacement +2 or Scale 2
elif p == -3 or s == 2.5:
    plt.plot(t_arr_data[:min(n4,N_20)], data2[3][:min(n4,N_20)]) #Displacement +3 or Scale 2.5

plt.xlabel("Time [s]",fontsize=14)
plt.ylabel("Angle [rad]",fontsize=14)
plt.ylim(-0.15,0.15)
plt.legend(["Model", "Data"])


