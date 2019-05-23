# -*- coding: utf-8 -*-
"""
Created on Tue May 30 14:28:39 2017

@author: dronchin
"""

import numpy as np
import winsound
import timeit


def unpackdata(initial_points_data):
    ''' This function needs to take a line of data and change it into usable 
    arrays
    '''
    detect1 = np.array(initial_points_data[0:3])
    detect2 = np.array(initial_points_data[3:])
    return detect1, detect2

def slope_angles(vec):
    '''Function that takes in a vector and defines it based on the polar 
    coordinate angles theta and phi while disregarding the length of the vector
    Input:
        vec = any vector in 3D
    Output:
        theta = the angle from the z-axis
        phi = the angle from the x axis to the vector's projection onto the x,y
            plane
    '''
    theta = np.arccos(np.dot(np.array([0,0,1]), vec)/np.linalg.norm(vec))
    proj = vec - np.dot(vec, np.array([0,0,1]))*np.array([0,0,1])
    phi = np.arccos(np.dot(np.array([1,0,0]), proj)/np.linalg.norm(proj))
    return theta, phi



starttime = timeit.default_timer()

#Erases any txt files already in place
ThetaMain_bins = []
ThetaMain_indexed = []
ThetaMain_bin_number = 50
for i in np.linspace(np.pi/3, 2*np.pi/3 - np.pi/150 , num=ThetaMain_bin_number): #0, np.pi - np.pi/ThetaMain_bin_number , num=ThetaMain_bin_number): #
    with open("C:/data/thetabins/bin_for_{0:.4f}.txt".format(i), "w") as writeto:
        ThetaMain_bins.append("{0:.4f}".format(i))
        ThetaMain_indexed.append(i)
            
        

data_values = []
with open("2gamma_3Dfakedata.txt", "r") as data:
    for line in data:
        #unpacks the data and stores them in a list
        numbers_float = map(float, line.split())
        initial_points_data = list(numbers_float)
        detect1, detect2 = unpackdata(initial_points_data) #saves data as detect1 and dectect2
        
        if detect1[1] > detect2[1]:
            point1 = detect1
            point2 = detect2
        else:
            point1 = detect2
            point2 = detect1
    
        vector = point1 - point2
        theta, phi = slope_angles(vector)
        vectornorm = vector/np.linalg.norm(vector)

        data_values.append([theta, phi, point1, vectornorm])
                
for i in range(len(ThetaMain_indexed)):
    with open("C:/data/thetabins/bin_for_{}.txt".format(ThetaMain_bins[i]), "w" ) as wroteto:
        correct = []
        for data in data_values:
            #check if it is correct
            if data[0] > ThetaMain_indexed[i] and data[0] < ThetaMain_indexed[i] + np.pi/150:
                correct.append(data)
        for data in correct:
            theta, phi, detect1, vectornorm = data
            writie = '{} {} {} {} {} {} {} {}\n'.format(theta, phi, detect1[0], detect1[1], detect1[2], vectornorm[0], vectornorm[1], vectornorm[2])
            wroteto.write(writie)    


stoptime = timeit.default_timer()
print("\n run info")
print("time taken: ", stoptime - starttime)  


winsound.Beep(800,1000)  



#