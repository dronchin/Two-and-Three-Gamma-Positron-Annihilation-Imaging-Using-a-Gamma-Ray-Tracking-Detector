# -*- coding: utf-8 -*-
"""
Created on Tue May  9 12:20:21 2017

@author: Nicolas D

The purpose of this program is to create 2d data for a future 2 gamma 
reconstruction program. It will simulate lines of response on write to a text
file the 2 points of decay.
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def point_on_detector(V, POd=np.array([0,0,0]), r=190):
    '''This function takes in vector and the point of decay. It outputs the point that would be detected.'''
    V = V/100
    Q = np.array([0,0,0])
    #uses the quadratic formula to find the intersect on the circle
    a = np.dot(V,V)
    b = 2* np.dot(V, POd)
    c = np.dot(POd,POd) + np.dot(Q,Q) - 2* np.dot(POd,Q) - r**2
    disc = b**2 - 4*a*c
    if disc > 0: #doulbe checks to stop potential errors
        #solve the quadratic equation for the points of interection
        t1 = (-b + np.sqrt(disc))/(2*a)
        t2 = (-b - np.sqrt(disc))/(2*a)
        
        detect1 = POd + t1*V
        detect2 = POd + t2*V
        
        return [detect1, detect2]
    else:
        print("Sorry, point wouldn't have been detected")

def makepoint(centx, centy, width, height, weightmult=1000, detector_radius=190, dist=100):
    global CURfile, X_all, Y_all, Z_all
    r_max = (width/2)**2
    volume = height*np.pi*(r_max)
    with open(CURfile,'a') as wroteto:
        for i in range(int(volume*weightmult)):
            R = np.random.uniform(0,r_max)
            
            th = np.random.uniform(0, 2*np.pi)
            z = np.random.uniform(-height/2, height/2)
            
            x = np.sqrt(R)*np.cos(th) + centx
            y = np.sqrt(R)*np.sin(th) + centy
            
            X_all.append(x)
            Y_all.append(y)
            Z_all.append(z)
            
            POdecay = np.array([x,y,z])
        
            rand_thet = np.random.uniform(0, 2*np.pi) #(np.pi/3, 2*np.pi/3)
            rand_phi = np.random.uniform(0, 2*np.pi)
            
            initialdirection = np.array([dist*np.sin(rand_thet)*np.cos(rand_phi),
                                     dist*np.sin(rand_thet)*np.sin(rand_phi),
                                     dist*np.cos(rand_thet)])
            detect1, detect2 = point_on_detector(initialdirection, POdecay, detector_radius)
        
            writie = '{} {} {} {} {} {}  \n'.format(detect1[0], detect1[1], detect1[2], 
                                                detect2[0], detect2[1], detect2[2])
            wroteto.write(writie)      
            
        
    
# =============================================================================
# MAIN PROGRAM
# =============================================================================
#current_file = '2gamma_fakedata.txt'
#with open(current_file, 'w'): pass


phantom_master = [[6.4952, 3.75, 2.5, 12.7], [11.4952, 3.75, 2.5, 12.7], [16.4952, 3.75, 2.5, 12.7], #2.5 first line
           [8.9952, 8.0801, 2.5, 12.7], [13.9952, 8.0801, 2.5, 12.7], #2.5 second line
           [11.4952, 12.4103, 2.5, 12.7], #2.5 third line
           
           [6.4952, -3.75, 2.0, 12.7], [10.4952, -3.75, 2.0, 12.7], [14.4952, -3.75, 2.0, 12.7], [18.4952, -3.75, 2.0, 12.7], #2.0 first line
           [8.4952, -7.2141, 2.0, 12.7], [12.4952, -7.2141, 2.0, 12.7], [16.4952, -7.2141, 2.0, 12.7], #2.0 second line
           [10.4952, -10.6782, 2.0, 12.7], [14.4952, -10.6782, 2.0, 12.7], #2.0 third line
           [12.4952, -14.1423, 2.0, 12.7], #2.0 fourth line
           
           [0, -7.5, 1.5, 12.7], #1.5 first line
           [-1.5, -10.0981, 1.5, 12.7], [-1.5 + 3, -10.0981, 1.5, 12.7], #1.5 second line
           [-3, -12.6962, 1.5, 12.7], [0, -12.6962, 1.5, 12.7], [3, -12.6962, 1.5, 12.7], #1.5 third line
           [-4.5, -15.2942, 1.5, 12.7], [-4.5 + 3, -15.2942, 1.5, 12.7], [-4.5 + 3*2, -15.2942, 1.5, 12.7], [-4.5 + 3*3, -15.2942, 1.5, 12.7], #1.5 fourth line
           [-6.0, -17.8923, 1.5, 12.7], [-3.0, -17.8923, 1.5, 12.7], [0., -17.8923, 1.5, 12.7], [3.0, -17.8923, 1.5, 12.7], [6.0, -17.8923, 1.5, 12.7], #1.5 fifth line
           
           [-6.4952, -3.75, 1.25, 12.7], [-6.4952 - 2.5, -3.75, 1.25, 12.7], [-6.4952 - 2.5*2, -3.75, 1.25, 12.7], [-6.4952 - 2.5*3, -3.75, 1.25, 12.7], [-6.4952 - 2.5*4, -3.75, 1.25, 12.7], [-6.4952 - 2.5*5, -3.75, 1.25, 12.7], #1.25 first line
           [-7.7452, -5.9151, 1.25, 12.7], [-7.7452 - 2.5, -5.9151, 1.25, 12.7], [-7.7452 - 2.5*2, -5.9151, 1.25, 12.7], [-7.7452 - 2.5*3, -5.9151, 1.25, 12.7], [-7.7452 - 2.5*4, -5.9151, 1.25, 12.7], #1.25 second line
           [-8.9952 - 2.5*0, -8.0801, 1.25, 12.7], [-8.9952 - 2.5*1, -8.0801, 1.25, 12.7], [-8.9952 - 2.5*2, -8.0801, 1.25, 12.7], [-8.9952 - 2.5*3, -8.0801, 1.25, 12.7], #1.25 third line
           [-10.2452 - 2.5*0, -10.2452, 1.25, 12.7], [-10.2452 - 2.5*1, -10.2452, 1.25, 12.7], [-10.2452 - 2.5*2, -10.2452, 1.25, 12.7], #1.25 fourth line
           [-11.4592 - 2.5*0, -12.4103, 1.25, 12.7], [-11.4592 - 2.5*1, -12.4103, 1.25, 12.7], #1.25 fifth line
           [-12.7452, -14.5753, 1.25, 12.7], #1.25 sixth line
           
           [-12.4952, 14.1423, 1.0, 12.7], #1.0 first line
           [-11.4952 - 2*0, 12.4103, 1.0, 12.7], [-11.4952 - 2*1, 12.4103, 1.0, 12.7], #1.0 second line
           [-10.4952 - 2*0, 10.6782, 1.0, 12.7], [-10.4952 - 2*1, 10.6782, 1.0, 12.7], [-10.4952 - 2*2, 10.6782, 1.0, 12.7], #1.0 thrid line
           [-9.4952 - 2*0, 8.9462, 1.0, 12.7], [-9.4952 - 2*1, 8.9462, 1.0, 12.7], [-9.4952 - 2*2, 8.9462, 1.0, 12.7], [-9.4952 - 2*3, 8.9462, 1.0, 12.7], #1.0 fourth line
           [-8.4952 - 2*0, 7.2141, 1.0, 12.7], [-8.4952 - 2*1, 7.2141, 1.0, 12.7], [-8.4952 - 2*2, 7.2141, 1.0, 12.7], [-8.4952 - 2*3, 7.2141, 1.0, 12.7], [-8.4952 - 2*4, 7.2141, 1.0, 12.7], #1.0 fifth line
           [-7.4952 - 2*0, 5.4821, 1.0, 12.7], [-7.4952 - 2*1, 5.4821, 1.0, 12.7], [-7.4952 - 2*2, 5.4821, 1.0, 12.7], [-7.4952 - 2*3, 5.4821, 1.0, 12.7], [-7.4952 - 2*4, 5.4821, 1.0, 12.7], [-7.4952 - 2*5, 5.4821, 1.0, 12.7], #1.0 sixth line
           [-6.4952 - 2*0, 3.75, 1.0, 12.7], [-6.4952 - 2*1, 3.75, 1.0, 12.7], [-6.4952 - 2*2, 3.75, 1.0, 12.7], [-6.4952 - 2*3, 3.75, 1.0, 12.7], [-6.4952 - 2*4, 3.75, 1.0, 12.7], [-6.4952 - 2*5, 3.75, 1.0, 12.7], [-6.4952 - 2*6, 3.75, 1.0, 12.7], #1.0 seventh line
           
           [-5.6 + 1.6*0, 17.1995, 1.0, 12.7], [-5.6 + 1.6*1, 17.1995, 1.0, 12.7], [-5.6 + 1.6*2, 17.1995, 1.0, 12.7], [-5.6 + 1.6*3, 17.1995, 1.0, 12.7], [-5.6 + 1.6*4, 17.1995, 1.0, 12.7], [-5.6 + 1.6*5, 17.1995, 1.0, 12.7], [-5.6 + 1.6*6, 17.1995, 1.0, 12.7], [-5.6 + 1.6*7, 17.1995, 1.0, 12.7], #0.8 first line
           [-4.8 + 1.6*0, 15.8138, 1.0, 12.7], [-4.8 + 1.6*1, 15.8138, 1.0, 12.7], [-4.8 + 1.6*2, 15.8138, 1.0, 12.7], [-4.8 + 1.6*3, 15.8138, 1.0, 12.7], [-4.8 + 1.6*4, 15.8138, 1.0, 12.7], [-4.8 + 1.6*5, 15.8138, 1.0, 12.7], [-4.8 + 1.6*6, 15.8138, 1.0, 12.7], #0.8 second line
           [-4.0 + 1.6*0, 14.4282, 1.0, 12.7], [-4.0 + 1.6*1, 14.4282, 1.0, 12.7], [-4.0 + 1.6*2, 14.4282, 1.0, 12.7], [-4.0 + 1.6*3, 14.4282, 1.0, 12.7], [-4.0 + 1.6*4, 14.4282, 1.0, 12.7], [-4.0 + 1.6*5, 14.4282, 1.0, 12.7], #0.8 third line
           [-3.2 + 1.6*0, 13.0426, 1.0, 12.7], [-3.2 + 1.6*1, 13.0426, 1.0, 12.7], [-3.2 + 1.6*2, 13.0426, 1.0, 12.7], [-3.2 + 1.6*3, 13.0426, 1.0, 12.7], [-3.2 + 1.6*4, 13.0426, 1.0, 12.7], #0.8 fourth line
           [-2.4 + 1.6*0, 11.6569, 1.0, 12.7], [-2.4 + 1.6*1, 11.6569, 1.0, 12.7], [-2.4 + 1.6*2, 11.6569, 1.0, 12.7], [-2.4 + 1.6*3, 11.6569, 1.0, 12.7], #0.8 fifth line
           [-1.6 + 1.6*0, 10.2713, 1.0, 12.7], [-1.6 + 1.6*1, 10.2713, 1.0, 12.7], [-1.6 + 1.6*2, 10.2713, 1.0, 12.7], #0.8 sixth line
           [-0.8 + 1.6*0, 8.8856, 1.0, 12.7], [-0.8 + 1.6*1, 8.8856, 1.0, 12.7], #0.8 seventh line
           [0, 7.5, 1.0, 12.7], #0.8 eight line
           ]

CURfile = '2gamma_3Dfakedata.txt'
with open(CURfile, 'w'): pass

#opens up a textfile and writes all the final data    

X_all = []
Y_all = []
Z_all = []


for point in phantom_master:
    x, y, width, height = point
    makepoint(x, y, width, height)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.scatter(X_all,Y_all,Z_all)
plt.show()


        
        

#    