# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 08:22:56 2019

@author: bergl
"""

import numpy as np
import matplotlib.pyplot as plt

folder = "m_0.200_eig_-3_-4_-5_-6"
Data_file1=np.loadtxt("Test_Recordings/"+folder+"/template.txt", delimiter=",")
Data_file2=np.loadtxt("Test_Recordings/"+folder+"/displacement_-1.txt", delimiter=",")
Data_file3=np.loadtxt("Test_Recordings/"+folder+"/displacement_-2.txt", delimiter=",")
Data_file4=np.loadtxt("Test_Recordings/"+folder+"/displacement_-3.txt", delimiter=",")


Data_file_curve=np.loadtxt("Test_Recordings/curveball.txt",delimiter=",")

cart_pos_curve=Data_file_curve[:,0]
pend_ang_curve=Data_file_curve[:,1]
cart_vel_curve=Data_file_curve[:,2]
pend_ang_vel_curve=Data_file_curve[:,3]
input_current_curve=Data_file_curve[:,4]


cart_pos1=Data_file1[:,0]
pend_ang1=Data_file1[:,1]
cart_vel1=Data_file1[:,2]
pend_ang_vel1=Data_file1[:,3]
input_current1=Data_file1[:,4]

cart_pos2=Data_file2[:,0]
pend_ang2=Data_file2[:,1]
cart_vel2=Data_file2[:,2]
pend_ang_vel2=Data_file2[:,3]
input_current2=Data_file2[:,4]

cart_pos3=Data_file3[:,0]
pend_ang3=Data_file3[:,1]
cart_vel3=Data_file3[:,2]
pend_ang_vel3=Data_file3[:,3]
input_current3=Data_file3[:,4]

cart_pos4=Data_file4[:,0]
pend_ang4=Data_file4[:,1]
cart_vel4=Data_file4[:,2]
pend_ang_vel4=Data_file4[:,3]
input_current4=Data_file4[:,4]

sampling_time=0.00667
sampling_frequency=1/sampling_time


n1=len(Data_file1)
n2=len(Data_file2)
n3=len(Data_file3)
n4=len(Data_file4)
N=max(n1,n2,n3,n4)
N_20=int(sampling_frequency*20)

t=np.linspace(0,sampling_time*N,N)


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
    p=0
    P1s.append([-3*s-p,-4*s-p,-5*s-p,-6*s-p])
    P2s.append([-3*s-p,-5*s-p,-7*s-p,-9*s-p])
    P3s.append([-3*s-p,-6*s-p,-9*s-p,-12*s-p])
    s += 0.5



plt.figure(figsize=(16,12))
plt.subplot(2,1,1)
plt.plot(t[:min(n1,N_20)],cart_pos1[:min(n1,N_20)])
plt.plot(t[:min(n2,N_20)],cart_pos2[:min(n2,N_20)])
plt.plot(t[:min(n3,N_20)],cart_pos3[:min(n3,N_20)])
plt.plot(t[:min(n4,N_20)],cart_pos4[:min(n4,N_20)])
plt.legend(P1p,loc="upper right")
plt.ylim(0,0.9)
plt.xlim(left=0,right=20)
plt.xlabel("Time [s]",fontsize=14)
plt.ylabel("Position [m]",fontsize=14)
plt.subplot(2,1,2)
plt.plot(t[:min(n1,N_20)],pend_ang1[:min(n1,N_20)])
plt.plot(t[:min(n2,N_20)],pend_ang2[:min(n2,N_20)])
plt.plot(t[:min(n3,N_20)],pend_ang3[:min(n3,N_20)])
plt.plot(t[:min(n4,N_20)],pend_ang4[:min(n4,N_20)])
plt.legend(P1p)
plt.ylim(-0.15,0.15)
plt.xlim(left=0,right=20)
plt.xlabel("Time [s]",fontsize=14)
plt.ylabel("Angle [rad]",fontsize=14)
plt.savefig("Test_Recordings/"+folder+"/Displacement_plot.png")
#plt.show()
#plt.plot(cart_pos4[:min(n3,N_20)],pend_ang4[:min(n3,N_20)],) 

cart_ac=[]
for i in range(len(cart_vel4)-1):
    ac=(cart_vel4[i+1]-cart_vel4[i])/sampling_time
    cart_ac.append(ac)