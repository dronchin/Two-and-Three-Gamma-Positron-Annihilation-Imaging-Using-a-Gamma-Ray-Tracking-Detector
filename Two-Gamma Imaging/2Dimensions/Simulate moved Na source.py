# -*- coding: utf-8 -*-
"""
Created on Mon Jan 14 15:14:07 2019

@author: dronchin
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d
import timeit

starttime = timeit.default_timer()

def unpackdata(initial_points_data):
    ''' This function needs to take a line of data and change it into usable 
    varibales
    '''
    detect1 = initial_points_data[0:3]
    detect2 = initial_points_data[3:]
    return detect1, detect2

def point_on_detector(V, POd=np.array([0,0,0]), r=180):
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
        
def getsource(filename):
    ''' function opens up a csv spreadsheet to retreve a nxn grid 
    representing the source
    ''' 
    source = []
    total = 0
    with open(filename,'r') as data:
        for line in data:
            #unpacks the data and stores them in a list
            numbers_float = map(float, line.split(","))
            number, counts = list(numbers_float)
            source.append(counts)
            total += counts
    
    source = np.asarray(source)
    source /= total
    
    return source

def getlocation(i, n, centx, centy, centz, pixsize):
    ''' calculates the decay position of the positronium based on what pixel was selected from the 
    autoradiograph data '''
    row = i//n - n/2
    col = i%n - n/2
    xpos = centx + row*pixsize + np.random.normal(0,scale=pixsize/3)
    ypos = centy + col*pixsize + np.random.normal(0,scale=pixsize/3)
    zpos = centz + np.random.normal(0,scale=pixsize/3)
    return np.array([xpos,ypos,zpos])
        
def getangle(D):
    phi = np.arctan2(D[1], D[0])
    theta = np.arctan2(np.sqrt(D[0]**2 + D[1]**2), D[2])
    return theta, phi

source = getsource("old Na 16x16.csv")
plt.imshow(source.reshape((20,20)))
plt.show()


zdir = 5
#not shifted
Allpoints = [[-33, -0.5, zdir], #12 oclock #35
             [38, 1, zdir], #6 oclock
             [5, -34 , zdir], #3 oclock #33
             [2, 35, zdir],#9 oclock
             [3, 7, zdir]] #center #-8
#shifted down by 20 mm             
Allpoints = [[-33+15, -0.5, zdir], #12 oclock #35
             [38+15, 1, zdir], #6 oclock
             [5+15, -34 , zdir], #3 oclock #33
             [2+15, 35, zdir],#9 oclock
             [3+15, 7, zdir]] #center #-8
             
             
             
             
hist = np.load("hist_90Deg_Gent4.npy")
thcenter = np.load("thcenter_90Deg.npy")
phicenter = np.load("phicenter_90Deg.npy")

plt.imshow(hist.T, origin='lower')
plt.show()

pointsaver = []

#Create test function to sample the 2dhistogram of the PDF(gretina angular distribution in this case)
PDFtest = interp2d(phicenter, thcenter, hist.T)
#ideal detector angular distribution
#test = 1*np.exp(-(np.pi/2 - randtheta)**2/(2*((np.pi/6)/2.355)**2)) 

current_file = "C:/data/simulate_moved_Na_downshifted_15.txt"
datasize = 78000    #78000#2000000


with open(current_file,'w') as wroteto:
    cnt = 0
    for pos in Allpoints:
        accepted = 0
        cnt += 1
        print(cnt)
        while accepted < datasize:
            if accepted%1000 == 0:
                print(cnt + round(accepted/datasize, 3))
                #just to get it to not show up again (gives 0.1% less counts than expected)
                accepted += 1
                
            #randomize initial decay direction
            randtheta = np.random.uniform(0, np.pi)
            randphi = np.random.uniform(-np.pi, np.pi)
            initialdirection = np.array([10*np.sin(randtheta)*np.cos(randphi), 
                                                 10*np.sin(randtheta)*np.sin(randphi), 
                                                 10*np.cos(randtheta)])
            
            #POdecay is based on a probabilistic choice of source geometry and center of the real source
            index = int(np.random.choice(len(source), 1, p=source))
            POdecay = getlocation(index, 20, pos[0], pos[1], pos[2], 16*0.025)
            
            #Detect where both photons hit in the detector
            detect1, detect2 = point_on_detector(initialdirection, POdecay, 190)
            
            #determine agnles from origin instead of POdecay
            theta_orig, phi_orig = getangle(detect1)
            
            #check first point to PDF
            if np.random.random() < PDFtest(phi_orig, theta_orig):
                #other detection is oposite direction
                theta_orig2, phi_orig2 = getangle(detect2)
                
                #check second point to PDF
                if np.random.random() < PDFtest(phi_orig2, theta_orig2):
                    accepted += 1
                    pointsaver.append([theta_orig, phi_orig])
                    pointsaver.append([theta_orig2, phi_orig2])
                    #if both events pass, accept the simulated point/write to file                 
                    writie = '{} {} {} {} {} {}  \n'.format(detect1[0], detect1[1], -detect1[2], 
                                                            detect2[0], detect2[1], -detect2[2],)
                    wroteto.write(writie)


#plot final angular distribution
pointsaver = np.asarray(pointsaver)

plt.hist2d(pointsaver[:,1], pointsaver[:,0], bins=[np.arange(-np.pi,np.pi,0.01), np.arange(0,np.pi,0.01)], normed=True)
plt.colorbar()
plt.title("Angular Distribution in Downshifted Simulation")
plt.xlabel("Phi (rad)")
plt.ylabel("Theta (rad)")
plt.savefig("Angular Distribution in Downshifted Simulation.png", dpi=1200, format='png')
plt.show()

# timer info
stoptime = timeit.default_timer()
print("time taken: ", stoptime - starttime)  

#