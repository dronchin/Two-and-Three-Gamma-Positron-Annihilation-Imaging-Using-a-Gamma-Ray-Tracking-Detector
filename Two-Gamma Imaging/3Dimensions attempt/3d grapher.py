# -*- coding: utf-8 -*-
"""
Created on Wed Jul 12 09:46:57 2017

@author: dronchin
"""

import numpy as np
import timeit
from math import ceil 
import winsound, sys
from os import listdir


#==============================================================================
# Functions
#==============================================================================

def get_line_3d(start, end):
    ''' Gets the points for the voxels that form a line in 3 dimensional space.
    Parameters    
    ----------
    start: array_like, cordinates for the starting point of the line in 3D
    end: array_like, cordinates for the ending point of the line in 3D
    Returns
    -------
    points: list, a list of points giving information about which voxels the
        line will go through. This is useful in conjuction with put_on_grid_3d
    '''
    p = start
    d = end-start
    N= int(round(max(abs(d)), 0)) #Number of ponts needed
    points = [] #empty list for saving points
    #rounds each point and then makes it an integer
    points.append(np.round(p, decimals= 0).astype(int))
    if N != 0: #only one point if N=0
        s = d/N #vector of change per interval
        for i in range(N):
            p = p+s #new point in float
            #rounds each point and then makes it an integer
            points.append(np.round(p, decimals= 0).astype(int))            
    return points #returns a list of points

def put_on_grid_3d(grid, points, value):
    ''' Puts the weight of a specific pixel from a 2d picture across all of the
    voxels in 3d on a specific grid.
    Parameters    
    ----------
    grid: 3D array_like, the current 3d array that represents a 3d picture.
    points: list of arrays, a list of arrays that contain information on the 
    points in space that make up a line.
    value: float, the value from the original picture that will be set for each
    point on the grid.
    Returns
    -------
    grid: the updated grid with the new values changed.
    '''
    for point in points:       
        x,y,z = point
        #!!! maybe average here in the future!!!
        #updates the grid's values at indexes. takes into account points have origin
        grid[grid_origin[0]+x][grid_origin[1]+y][grid_origin[2]+z] += value           
    return grid

def RotateYZ(pointToRotate, phi, theta):
    '''This function will rotate a point/vector around the y-axis theta amount, 
    and then rotate it around the Z-axis phi amount.'''
    pointToRotate = np.asarray(pointToRotate)
    theta = theta - np.pi/2 #want to return theta to pi/2
    #rotation matrix round Y axis
    rotate_matrix_Y = np.array([[np.cos(theta) , 0, np.sin(theta)],
                                [0             , 1, 0            ],
                                [-np.sin(theta), 0, np.cos(theta)]])
    #rotation matrix round Z axis
    rotate_matrix_Z = np.array([[np.cos(phi), -np.sin(phi), 0],
                                [np.sin(phi), np.cos(phi) , 0],
                                [0          , 0           , 1]])
    #multiplies point by both matrices
    rotated1 = np.matmul(rotate_matrix_Y, pointToRotate)
    rotated2 = np.matmul(rotate_matrix_Z, rotated1)
    return rotated2

def linepoints(phi, theta, radius, a, b, origin):
    '''Finds the starting point and end point of a line based on the 
    orientation of the picture and what element is being examined.
    Parameters    
    ----------
    phi, theta: float, phi or theta angle that the picture was taken at. 
    radius: float, distance out from the center to the center of the picture
    a, b: int, idexes for the position on the picture
    origin: array, set up relation between indexes and position on grid
    Returns
    -------
    start, end: 3x1 arrays, these are the starting and ending positions for a 
        line
    '''
    norm = np.array([np.sin(theta)*np.cos(phi),np.sin(theta)*np.sin(phi), np.cos(theta)])
    R = norm * radius
    C = np.array([0, -(origin[1] - b - 0.5), (origin[0] - a - 0.5)])
    Q = RotateYZ(C, phi, theta)
    start = R + Q
    #math for end points
    t = -np.dot(start,norm)/np.dot(norm,norm) #how many normal vectors the start is from the plane
    end = start+ norm*2*t #the mirrored point is twice the distance to the plane
    return start, end

def set_grid(sinogram3d):
    '''Takes in the initial 3d signogram and estimates what the grid size needs
    to be through padding
    Returns
    -------
    grid: array, a 3x3 array full of zeros
    radius: float, estimate on the radius based on the initial signogram
    '''
    gridshape = sinogram3d.shape
    print(gridshape)
    radius = max(gridshape[1],gridshape[2])//2 #raidus is max dimension of the picture
    grid = np.zeros([radius,radius,radius]) #grid of zeros
    a = ceil(radius*2) #padding on each side
    padwidth = [[a,a],[a,a],[a,a]]
    grid = np.pad(grid, padwidth, mode='constant', constant_values=0) #pads the grid
    return grid, radius

def unpad(grid, radius):
    '''takes in the padded grid and cutts it back to the unpadded size
    '''
    beg = ceil(radius)
    end = grid.shape
    #print(beg, end[0]-beg, end[1]-beg, end[2]-beg)
    #grid = grid[beg:end[0]-beg][beg:end[1]-beg][beg:end[2]-beg]
    grid = grid[beg:(end[0]-beg), beg:(end[1]-beg), beg:(end[2]-beg)]
    return grid

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
    text = "\rPercent: [{0}] {1}% {2}".format( "#"*block + "-"*(barLength-block), round(progress*100, 4), status)
    sys.stdout.write(text)
    sys.stdout.flush()

#==============================================================================
# Main Program
#==============================================================================
starttime = timeit.default_timer() #starts timer

#looks to see what files are in the folder
fileslist = list(listdir("C:/data/sinograms3d"))
 
#takes first file and sets grid based off size
sinogram3d_1st = np.load("C:/data/sinograms3d/{}".format(fileslist[1]))
grid, orig_radius = set_grid(sinogram3d_1st) #sets the grid
grid_origin = (grid.shape[0] // 2, grid.shape[1] // 2, grid.shape[2] // 2)

#sets the bin amount numbers
ThetaMain_bin_number = 50
PhiMain_bin_number = 90
#reclaculates the thetabin intervals
thetabins = np.linspace(np.pi/3, 2*np.pi/3 , num=ThetaMain_bin_number+1)
#thetabins = np.linspace(0, np.pi, num=ThetaMain_bin_number+1)
    

for filenum, filename in enumerate(fileslist): #goes through all files in folder
    #opens up the 3d sinogram
    sinogram3d = np.load("C:/data/sinograms3d/{}".format(filename))
    theta = (thetabins[filenum] + thetabins[filenum+1])/2 #picks the middle of the theta bin
    
    for phinum,picture in enumerate(sinogram3d): #loops through each 2d hist in each sinogram
        #get phi from the middle of the bin and then gets the pictures origin
        
        
        
        
        #======================================================================
        # added np.pi/4 for rotation purposes        
        #======================================================================
        phi = (phinum + 0.5)*(np.pi/PhiMain_bin_number) #+ np.pi/4
        
        
        
        
        pic_origin = (picture.shape[0] / 2, picture.shape[1] / 2)
        
        #Loops through each cell of the picture while keeping track of which one its at
        for a, line in enumerate(picture):
            for b, weight in enumerate(line):
                #skips the space if it has no weight
                if weight == 0:
                    continue
                if weight != 0 :
                    #finds the start and end points for a line that square follows
                    start, end = linepoints(phi, theta, orig_radius, a, b, pic_origin)
                    #gets a set of points that the line follows. The points represent a voxel in the grid
                    points = get_line_3d(start, end)
                    #uses the points to change the weight of the voxels
                    grid = put_on_grid_3d(grid, points, weight)
                    
                    
    #updates the progress to show the program is running
    update_progress((filenum+1)/(len(fileslist)+1) )
    
#unpads the data to remove extraneous zeros
grid = unpad(grid, orig_radius)
#saves the grid to be used in next program
np.save("C:/data/grid", grid)

#end of run info
update_progress(1)
stoptime = timeit.default_timer()
print("\n run info")
print("time taken: ", stoptime - starttime)  
winsound.Beep(800,1000)  


##Extra
#dataline = np.reshape(grid, grid.size) #puts data in a 1D array




#