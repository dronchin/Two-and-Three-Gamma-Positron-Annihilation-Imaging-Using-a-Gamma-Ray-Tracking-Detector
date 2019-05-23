# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 14:45:28 2019

@author: dronchin
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d
import timeit

starttime = timeit.default_timer()


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
        
        return detect1
    else:
        print(V, disc)
        print("Sorry, point wouldn't have been detected")

def R(theta, u): #theta = amount of angle rotated, u = given axis (unit vect)
    return [[np.cos(theta) + u[0]**2 * (1-np.cos(theta)), 
             u[0] * u[1] * (1-np.cos(theta)) - u[2] * np.sin(theta), 
             u[0] * u[2] * (1 - np.cos(theta)) + u[1] * np.sin(theta)],
            [u[0] * u[1] * (1-np.cos(theta)) + u[2] * np.sin(theta),
             np.cos(theta) + u[1]**2 * (1-np.cos(theta)),
             u[1] * u[2] * (1 - np.cos(theta)) - u[0] * np.sin(theta)],
            [u[0] * u[2] * (1-np.cos(theta)) - u[1] * np.sin(theta),
             u[1] * u[2] * (1-np.cos(theta)) + u[0] * np.sin(theta),
             np.cos(theta) + u[2]**2 * (1-np.cos(theta))]]
#poitToRotate = point being rotated, around = axis of rotation, 
#theata=amount of angle rotated
def Rotate(pointToRotate, around, theta): 
    '''takes in an initial point, a vector to rotate it arond, and then the angle of how much it rotates.
    This uses a 3d rotational function defined as R()'''
    rotate_matrix = np.array(R(theta, around))
    rotated = np.matmul(rotate_matrix,pointToRotate)
    return rotated

def angles(e1,e2,e3):
    '''finds the angles of the unique trianlge created by the lengths of the energies
    Inputs = (e1,e2,e3) where e1 = energy1, e2 = energy2, and e3 = energy3
    '''
    alpha = np.arccos((e2**2-e1**2-e3**2)/(-2*e1*e3))
    beta = np.arccos((e3**2-e1**2-e2**2)/(-2*e1*e2))
    gamma = np.arccos((e1**2-e2**2-e3**2)/(-2*e2*e3))

    return [alpha,beta,gamma]

def getangle(D):
    phi = np.arctan2(D[1], D[0])
    theta = np.arctan2(np.sqrt(D[0]**2 + D[1]**2), D[2])
    return theta, phi

def create_detection_locations(PO_decay, energy):
    # used to find random starting angle
    rand_phi = np.random.uniform(-np.pi, np.pi) 
    #pick theta smarter
    rand_theta = np.random.uniform(0, np.pi)
    vec1 = np.array([np.sin(rand_theta)*np.cos(rand_phi),
                     np.sin(rand_theta)*np.sin(rand_phi),
                     np.cos(rand_theta)]) #intial vector is in random position

    e1, e2, e3 = energy
    alpha, beta, gamma = angles(e1, e2, e3)
    
    #find the angle to rotate around, has to be normal to random theta picked
    ntheta = rand_theta - np.pi/2
    if ntheta < 0:
        ntheta *= -1
        nphi = rand_phi + np.pi
    else:
        nphi = rand_phi
    #create normal vector
    norm = np.array([np.sin(ntheta)*np.cos(nphi),
                     np.sin(ntheta)*np.sin(nphi),
                     np.cos(ntheta)])
    
    #rotate the initial vector around the norm to find the other 2 vectors
    vec2 = Rotate(vec1, norm, (alpha+gamma))
    vec3 = Rotate(vec1, norm, (-(beta+gamma)))
    
    #finds the points of detection based on the vectors that point to where there should be a detection
    detect1 = point_on_detector(vec1, PO_decay)
    detect2 = point_on_detector(vec2, PO_decay)
    detect3 = point_on_detector(vec3, PO_decay)

    return [detect1,detect2,detect3]

def find_perfect_point_position(PO_decay, PDFtest, energy):
    '''This function takes in a radius and angle, then calculates where the detections would be on the detector.
    This does not assume the first decay is vertical and will therefore rotate where the detection points are.
    This assumes that the decay splits the energies perfectly.
    It inputs r and theta which correspond to polar coordinates for where the point of decay is.'''

#    PO_decay = np.array([0,0,0])
    
    # Monte Carlo to determine 3gamma that makes sense 
    # give limited num of tries so to not get stuck
    for i in range(100000):
        
        # find actual angles based on location
        det1, det2, det3 = create_detection_locations(PO_decay, energy)
        
        # calculate the angle made from the center. cent_angle = [theta, phi]
        cent_angle1 = getangle(det1)
        cent_angle2 = getangle(det2)
        cent_angle3 = getangle(det3)
        
        # need to calculate the angle between the detections because of later cuts
        # project down to xy plane
        xydet1 = np.array([det1[0],det1[1]])
        xydet2 = np.array([det2[0],det2[1]])
        xydet3 = np.array([det3[0],det3[1]])
        # use formula for cosine to calcluate angle between
        anglebetween1 = np.arccos(np.dot(xydet1,xydet2) / (np.linalg.norm(xydet1) * np.linalg.norm(xydet2)))
        anglebetween2 = np.arccos(np.dot(xydet2,xydet3) / (np.linalg.norm(xydet2) * np.linalg.norm(xydet3)))
        anglebetween3 = np.arccos(np.dot(xydet1,xydet3) / (np.linalg.norm(xydet1) * np.linalg.norm(xydet3)))

        # Dirk made cuts to angles where angles are >20 and <160 degrees
        if anglebetween1 < np.radians(20) or anglebetween2 > np.radians(160):
            continue
        if anglebetween2 < np.radians(20) or anglebetween2 > np.radians(160):
            continue
        if anglebetween3 < np.radians(20) or anglebetween3 > np.radians(160):
            continue
        
        # check is decay is accepted based on MC methods
        if np.random.rand() < PDFtest(cent_angle1[1], cent_angle1[0]):
            if np.random.rand() < PDFtest(cent_angle2[1], cent_angle2[0]):
                if np.random.rand() < PDFtest(cent_angle3[1], cent_angle3[0]):
                    # Decay is accepted
                    return 1, [det1, det2, det3]

    return 0, [0,0,0]

            

def cross_section(e1, e2):
    ''' A function that uses the cross section for positronium formation to 
    return the probability that an event with e1 = energy1 and e2 = energy2.
    The energies are scaled between the range [0,1) where 0.5 = 511keV. 
    
    The equation for the cross section comes from Kamińska, et.al. (2016) 
    https://doi.org/10.1140/epjc/s10052-016-4294-3
    '''
    # reject if the energies are 0 because it blows up (rare case in np.random.unifrom())
    if e1==0 or e2==0:
        return 0
    # reject because e1 + e2 > 0.5 has to be true as e3 !> 0.5
    elif e2 < (0.5 - e1) :
        return 0
    # return the cross section of the energies/4 where 4 is the maximum value
    # the cross section outputs with e1,e2 <= 0.5
    else:
        return ((e1 + e2 - 0.5)**2/(e1**2*e2**2))/4

def getlocation(i, n, centx, centy, centz, pixsize):
    row = i//n - n/2
    col = i%n - n/2
    xpos = centx + row*pixsize + np.random.normal(0,scale=pixsize/3)
    ypos = centy + col*pixsize + np.random.normal(0,scale=pixsize/3)
    zpos = centz + np.random.normal(0,scale=pixsize/3)
    return np.array([xpos,ypos,zpos])

##########################################################################
# MAIN PROGRAM
##########################################################################

source = getsource("magic powder 8x8.csv")
plt.imshow(source.reshape((20,20)))
plt.grid(False)
plt.show()

current_file = 'sim_source2_magicpowder.txt'
#with open('C:/data/3gamma data sets/{}'.format(current_file), 'w'): pass


hist = np.load("hist_config2.npy")
thcenter = np.load("thcenter_config2.npy")
phicenter = np.load("phicenter_config2.npy")

PDFtest = interp2d(phicenter, thcenter, hist.T)

plt.imshow(hist.T, origin="lower", cmap="Greys")
plt.show()

datasize = 500

all_x = []
all_y =[]
all_z = []
accepted_points = []
count = 0

with open('C:/data/3gamma data sets/{}'.format(current_file),'a') as wroteto:
    Energies = []
    for i in range(datasize): # Note: Does not return the same number of events
        e1, e2 = np.random.uniform(0,0.5,2) # Randomize e1, e2  between [0,0.5)
        prob = cross_section(e1,e2) # Returns the probability based on the cross section
        if np.random.rand() < prob: # Monte Carlo to determine if we accept
            # If accepted, add them to the list of all of the energies
            Energies.append([e1, e2,(1-e1-e2)])
    for energy in Energies:
        e1,e2,e3 = energy
        maxE = max(energy)
        #approximate energy cut to data
        if maxE > 0.5:
            continue
        
        # pick loc within source
        # uniform cylinder source
#        r_max = (10.0)**2
#        r_pick = np.random.uniform(0,r_max)
#        th = np.random.uniform(0, 2*np.pi)
#        x = np.sqrt(r_pick)*np.cos(th)
#        y = np.sqrt(r_pick)*np.sin(th)
#        z = np.random.uniform(-5,5)
#        pos = np.array([x,y,z])
        # realistic autoradiograph source
        pos = [4,0,0]
        index = int(np.random.choice(len(source), 1, p=source))
        PO_decay = getlocation(index, 20, pos[0], pos[1], pos[2], 8*0.05)
        # use a real point source
#        PO_decay = np.array([0,0,0])
                
        # determine angles based on energies
        alpha,beta,gamma = angles(e1,e2,e3)
        th1 = alpha + beta
        th2 = beta + gamma
        th3 = alpha + gamma
        
        # Use loc and angles to pic angle detected at
        ind, ans = find_perfect_point_position(PO_decay, PDFtest, energy)
        if ind == 1:
            d1,d2,d3 = ans
        else:
            continue
        
        accepted_points.append(PO_decay)

        # Point is accepted
        e1 *= 1022
        e2 *= 1022
        e3 *= 1022
        
        count += 1

        print(count)
        all_x.append(d1[0])
        all_x.append(d2[0])
        all_x.append(d3[0])
        all_y.append(d1[1])
        all_y.append(d2[1])
        all_y.append(d3[1])
        all_z.append(d1[2])
        all_z.append(d2[2])
        all_z.append(d3[2])
        
        accepted_points.append([d1,d2,d3])
        
        
        writie = '{} {} {} {} {} {} {} {} {} {} {} {} \n'.format(d1[0],d1[1],d1[2],e1,
                                                d2[0],d2[1],d2[2],e2,
                                                d3[0],d3[1],d3[2],e3)
        wroteto.write(writie)


accepted_points = np.asarray(accepted_points)

def Heatmapit(data1, data2, name, step=7.5, lim=240):
    xmin = -lim #xr.min()
    xmax = lim #xr.max()
    ymin = -lim #yr.min()
    ymax = lim #yr.max()
    xbin = np.arange(xmin,xmax+step,step)
    ybin = np.arange(ymin,ymax+step,step)
    print(xbin)
 
    plt.hist2d(data1, data2, bins=[xbin,ybin],cmap='Greys')
    plt.colorbar() #ticks=[0,1,2,3])
    plt.grid(color='k', alpha=0.15, linestyle='--')
    plt.axes().set_aspect('equal')
    plt.xlabel('x (mm)')
    plt.ylabel('z (mm)')
    plt.title(name)
    
#    plt.savefig("{}.png".format(name), dpi=1200, format='png' )
    plt.show()

Heatmapit(all_x,all_y, "simulated three-gamma detection positions (XY)")
Heatmapit(all_x,all_z, "simulated three-gamma detection positions (XZ)")
Heatmapit(accepted_points[:,0], accepted_points[:,1], "accepted points map", 0.2, 15)