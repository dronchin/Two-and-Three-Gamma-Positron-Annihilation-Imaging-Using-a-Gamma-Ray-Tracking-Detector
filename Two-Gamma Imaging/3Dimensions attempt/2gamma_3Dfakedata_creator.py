# -*- coding: utf-8 -*-
"""
Created on Tue May  9 12:20:21 2017

@author: Nicolas D

The purpose of this program is to create 2d data for a future 2 gamma 
reconstruction program. It will simulate lines of response on write to a text
file the 2 points of decay.
"""

import numpy as np

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

def makepoint(x,y,z, std, amount, dist=100, detector_radius=190):
    for i in range(amount):
        xer = np.random.normal(x, std)
        yer = np.random.normal(y, std)
        zer = np.random.normal(z, std) 
        POdecay = np.array([xer,yer,zer])
        #POdecay = np.array([10,10,0])
        
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

#opens up a textfile and writes all the final data    
with open('2gamma_3Dfakedata.txt','w') as wroteto:
    makepoint(0,0,0,2,50000)
    makepoint(10,0,0,2,50000)
    makepoint(50,0,0,2,50000)
    makepoint(15,15,15,3,112500)
    makepoint(30,30,30,4,200000)
    makepoint(45,45,45,5,312500)
    
#    dist = 100
#    detector_radius = 190
#    
#    for i in range(50000):
#        x = np.random.normal(50, 2)
#        y = np.random.normal(50, 2)
#        z = np.random.normal(50, 2) 
#        POdecay = np.array([x,y,z])
#        #POdecay = np.array([10,10,0])
#        
#        
#        rand_thet = np.random.uniform(np.pi/3, 2*np.pi/3)
#        rand_phi = np.random.uniform(0, 2*np.pi)
#        
#        initialdirection = np.array([dist*np.sin(rand_thet)*np.cos(rand_phi),
#                                     dist*np.sin(rand_thet)*np.sin(rand_phi),
#                                     dist*np.cos(rand_thet)])
#        detect1, detect2 = point_on_detector(initialdirection, POdecay, detector_radius)
#        
#        writie = '{} {} {} {} {} {}  \n'.format(detect1[0], detect1[1], detect1[2], 
#                                                detect2[0], detect2[1], detect2[2],)
#        wroteto.write(writie)
#        
#        zlim = 40    
#        if detect1[2] < zlim and detect2[2] < zlim:
#            if detect1[2] > -zlim and detect2[2] > -zlim:
#                writie = '{} {} {} {} {} {}  \n'.format(detect1[0], detect1[1], detect1[2], 
#                                                        detect2[0], detect2[1], detect2[2],)
#                wroteto.write(writie)
        
        
        

#    