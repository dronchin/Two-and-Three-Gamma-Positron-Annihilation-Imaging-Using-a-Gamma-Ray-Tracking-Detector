# -*- coding: utf-8 -*-
"""
Created on Wed May 31 13:36:45 2017

@author: dronchin
"""

import numpy as np
import winsound
from os import listdir
import matplotlib.pyplot as plt
import timeit
import sys

#==============================================================================
# Functions
#==============================================================================
def unpackdata2(initial_points_data):
    ''' This function needs to take a line of data and change it into usable 
    arrays
    '''
    theta = initial_points_data[0]
    phi = initial_points_data[1]
    equ0 = np.array(initial_points_data[2:5])
    equ1 = np.array(initial_points_data[5:8])
    return theta, phi, equ0, equ1

def RotateZY(pointToRotate, phi, theta):
    '''This function will rotate a point/vector around the Z-axis phi amount, 
    and then rotate it around the y-axis theta amount.'''
    pointToRotate = np.asarray(pointToRotate)
    phi = -phi
    rotate_matrix_Z = np.array([[np.cos(phi), -np.sin(phi), 0],
                                [np.sin(phi), np.cos(phi) , 0],
                                [0          , 0           , 1]])
    theta = np.pi/2 - theta
    rotate_matrix_Y = np.array([[np.cos(theta) , 0, np.sin(theta)],
                                [0             , 1, 0            ],
                                [-np.sin(theta), 0, np.cos(theta)]])
    rotated1 = np.matmul(rotate_matrix_Z, pointToRotate)
    rotated2 = np.matmul(rotate_matrix_Y, rotated1)
    return rotated2

def intersect2D(equ0, equ1, planetheta, planephi):
    '''This function will find the intersection of a line on a plane that is 
    angled. The intersection on the angled plane isn't useful though so the 
    intersection point is rotated so that it lies in a vertical plane like 
    there is a new basis.
    Inputs:
        equ0 = if <x0,y0,z0> + <a,b,c>*t is the equation for a line, then equ0 = <x0,y0,z0>
        equ1 = <a,b,c>
        planetheta = the current theta that is being sampled from the data
        planephi = the current phi that is being sampled from the data
    Outputs:
        ys = the distance in the x,y plane
        zs = the distance in the z direction of the rotated plane
    '''
    n = np.array([np.sin(planetheta)*np.cos(planephi),
                  np.sin(planetheta)*np.sin(planephi),
                  np.cos(planetheta)]) #calculates the normal vector
    #solves for the time part of the intersection between a line and plane
    t = -np.dot(n, equ0)/ np.dot(n, equ1)
    #uses time to solve for the intersection
    x = equ0[0] + equ1[0]*t 
    y = equ0[1] + equ1[1]*t
    z = equ0[2] + equ1[2]*t
    #rotates around Z and then Y
    xs, ys, zs = RotateZY([x,y,z], planephi, planetheta)
    return ys, zs

def update_progress(progress):
    '''Displays or updates a console progress bar. Accepts a float between 0 
    and 1. Any int will be converted to a float.
    A value under 0 represents a 'halt'.
    A value at 1 or bigger represents 100%.
    '''
    barLength = 10 # Modify this to change the length of the progress bar
    status = ""
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
        status = "error: progress var must be float\r\n"
    if progress < 0:
        progress = 0
        status = "Halt...\r\n"
    if progress >= 1:
        progress = 1
        status = "Done...\r\n"
    block = int(round(barLength*progress))
    text = "\rPercent: [{0}] {1}% {2}".format( "#"*block + "-"*(barLength-block), round(progress*100, 1), status)
    sys.stdout.write(text)
    sys.stdout.flush()

#==============================================================================
# Main program
#==============================================================================

starttime = timeit.default_timer() #for diagnostic purposes

#finds all theta bin file names
fileslist = list(listdir("C:/data/thetabins"))

PhiMain_bin_number = 90 #Set by user.If changed double check nothing else changes0
Philist = [[] for _ in range(PhiMain_bin_number)] #prepares all Philists

thetalistofphilists = [] #master list that holds the other lists

##FINDS THE INTERCEPT FOR ALL DATA AND SAVES IN SIMULAR PHI BINS


for progress, file in enumerate(fileslist): #picks a theta file
    thetastring = file
#    print(thetastring)
    planetheta = float(thetastring[8:14]) #gets the plane's theta out of the file name
    
    Philist = [[] for _ in range(PhiMain_bin_number)] #prepares all Philists
    
    #opens theta file
    with open("C:/data/thetabins/{}".format(file), "r") as data:
        for line in data:
            #unpacks the data and stores them in correct variables
            numbers_float = map(float, line.split())
            initial_points_data = list(numbers_float)     
            theta, phi, equ0, equ1 = unpackdata2(initial_points_data)
            
            #Use theta, phi, equ0, equ1 to find 2D intersect
            xdis, ydis = intersect2D(equ0, equ1, theta, phi)
            
            #picks the correct part of list for folder
            for i in range(PhiMain_bin_number):
                if phi >= i*(np.pi/PhiMain_bin_number) and phi <= (i+1)*(np.pi/PhiMain_bin_number):
                    Philist[i].append([xdis, ydis, phi]) #saves data in the list
                
    thetalistofphilists.append(Philist) #saves the Philist in the master list
    update_progress((progress+1)/len(fileslist)) #progess update


##2D HISTOGRAM EACH PHI BIN IN EACH THETA BIN
for num, Philist in enumerate(thetalistofphilists): #opens up each Philist 
    sinogram3d = [] #saves the total 3d sinogram
#    if num == 3: 
#        break
    for philistnum, current in enumerate(Philist): #looks at the group of points for a single plane
        current = np.array(current)
        if len(current) == 0: #checks for empty array
            hist = np.zeros((360,360)) #189,189
            sinogram3d.append(hist) #appends empty 2d 0array
            continue
        
        xcol = current[:,0] #cuts all x values
        ycol = current[:,1] #cuts all y values
        #2d histogram groups points close to eachother in a grid
        hist, _, _ = np.histogram2d(xcol, ycol, bins=[np.arange(-180,181,1), 
                                    np.arange(-180,181,1)], normed=False)
        
#        plt.imshow(hist)
#        plt.show()
        hist = hist.T #transposes the 2d plane
        hist = hist[::-1,:] #flips vertically as well
#        print(hist.shape, type(hist))
        sinogram3d.append(hist)
        
#        #optional visualization (732.1 sec with pics/ 87.7 sec without pics) 
#        plt.imshow(hist, cmap="Greys")  #hist[::-1,:], cmap="Greys") #need to invert to get correct looking picture 
#        plt.title(np.mean(current[:,2]))
#        plt.show()
#        print("")   
    #saves the 3d sinograms in a folder with an incremented numbmer
    sinogram3d = np.asarray(sinogram3d)
    np.save("C:/data/sinograms3d/{name}{num:03d}".format(
            name="sinogram3d", num=(num + 1)), sinogram3d) 
    update_progress((num+1)/len(thetalistofphilists))    

#runtime info
stoptime = timeit.default_timer()
print("\n run info")
print("time taken: ", stoptime - starttime) 

print("done!")
winsound.Beep(800,1000)  

#