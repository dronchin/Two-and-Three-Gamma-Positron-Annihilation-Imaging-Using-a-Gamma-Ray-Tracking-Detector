# -*- coding: utf-8 -*-
"""
Created on Mon Oct 22 16:17:38 2018

@author: dronchin
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage.filters import gaussian_filter
from scipy import optimize

file = '3Phot_001_output.txt'

x = []
y = []
z = []

# import all of the 3 gamma reconstructions
with open(file,'r') as data:
    for line in data:
        #unpacks the data and stores them in a list
        numbers_float = map(float, line.split())
        point = list(numbers_float)
          
        x.append(point[0])
        y.append(point[1])
        z.append(point[2])
xr = np.array(x)
yr = np.array(y)
zr = np.array(z)

# import the autoradio data
xpos = []
intensity = []
with open("magicpowder source.csv",'r') as data:
    next(data)
    for line in data:
        #unpacks the data and stores them in a list
        numbers_float = map(float, line.split(","))
        initial_points_data = list(numbers_float)
        xpos.append(initial_points_data[0]*0.05)
        
        intensity.append(np.average(initial_points_data[1:3]))
        
xpos = xpos[::20]
xpos = xpos[0:len(xpos)]
intensity = intensity[::20]
intensity = intensity[0:len(intensity)]
plt.plot(xpos, intensity)
plt.title("autoradio line")
#plt.plot(xpos[::20],intensity2[::20])
plt.show()


length = len(xpos)


xmin = -12
ymin = -17

xbin = np.arange(xmin, xmin+length+1,1)
ybin = np.arange(ymin, ymin+length+1,1)    

grid, _, _ = np.histogram2d(x, y, bins=[xbin,ybin])
grid = grid.T
plt.plot(19,19, "go")
plt.imshow(grid, cmap="inferno") #grid[120:220,100:200])
plt.grid(True)
plt.title("heat map")
plt.show()
grid = np.asarray(grid)



#

filtered = gaussian_filter(intensity, sigma= 3)
plt.plot(xpos, filtered)
plt.title("smeared autoradio line sig=3")
plt.show()


origin = [grid.shape[0]//2, grid.shape[1]//2]

x = np.linspace(0,length,length+1)
xdirection = []
xline = []
for j in range(grid.shape[0]):
    averager = []
    for i in range(5):
        averager.append(grid[origin[0]-5+i,j])
    xdirection.append(np.average(averager))
    xline.append(int(x[j]))

plt.plot(xline,xdirection, "b.")
plt.show()

y = np.linspace(0,length,length+1)
ydirection = []
yline = []
for j in range(grid.shape[1]):
    averager = []
    for i in range(5):
        averager.append(grid[j,origin[0]-xmin-16+i])
    ydirection.append(np.average(averager))
    yline.append(int(y[j]))

plt.plot(yline,ydirection, "b.")
plt.show()

##xline = np.linspace(0,length,length)
##xdirection = np.zeros_like(grid[0])
##for line in grid:
##    
##    xdirection += line
##xdirection = xdirection[origin[0]-len(intensity)//2: origin[0]+len(intensity)//2]
##xline = range(len(xdirection))
#
#ydirection = ydirection[origin[1]-len(intensity)//2: origin[1]+len(intensity)//2]
#yline = range(len(ydirection))
#plt.plot(ydirection)
#plt.show()
#
def gaussian(x, amp, sig, back):
    global intensity
    
    if sig < 0:
        return line*100000000
    f = gaussian_filter(intensity, sigma= np.abs(sig))
    
#    print("amp, sig, back", amp, sig, back)
#    print(x)
    return amp*f[x] + back  #runtime warning when it gets close to the value of zero

    
init_vals = [1, 4, 0.1]  # for [amp, sig, back]
best_vals, covar = optimize.curve_fit(gaussian, xline, xdirection, p0=init_vals)

stdx = np.sqrt(best_vals[1])
FWHMx = 2.355*stdx
reg = 6
cent = best_vals[1]

plt.plot(xline, xdirection, "bo")
plt.plot(xline, gaussian(xline, best_vals[0],best_vals[1], best_vals[2]), "r")
plt.title("X FWHM fit")
plt.show()

init_vals = [1, 4, 0.1]  # for [amp, sig, back]
best_vals, covar = optimize.curve_fit(gaussian, yline, ydirection, p0=init_vals)

stdy = np.sqrt(best_vals[1])
FWHMy = 2.355*stdy
reg = 6
cent = best_vals[1]

plt.plot(yline, ydirection, "bo")
plt.plot(yline, gaussian(yline, best_vals[0],best_vals[1], best_vals[2]), "r")
plt.title("Y FWHM fit")
plt.show()

step = 1
print("standard Deviation x: ", stdx*step)
print("FWHMx: ", FWHMx*step)
print("standard Deviation y: ", stdy*step)
print("FWHMy: ", FWHMy*step)

##
##grid, _, _ = np.histogram2d(x, z, bins=[xbin,zbin])
##plt.imshow(grid, cmap="inferno") #grid[120:220,100:200])
##plt.grid(False)
##plt.show()
##
##xline = np.linspace(0,length,length)
##xdirection = np.zeros_like(grid[0])
##for line in grid:
##    xdirection += line
##plt.plot(xdirection)
##plt.show()
##
##
##init_vals = [1, 15, 1]  # for [amp, cen, wid, back]
##best_vals, covar = optimize.curve_fit(gaussian, xline, xdirection, p0=init_vals)
##
##stdx = np.sqrt(best_vals[2]/2)
##FWHMx = 2.355*stdx
##reg = 6
##cent = best_vals[1]
##
##plt.plot(xline, xdirection, "bo")
##plt.plot(xline, gaussian(xline, best_vals[0],best_vals[1], best_vals[2]), "r")
##plt.show()
## 
##print("standard Deviation: ", stdx*step)
##print("FWHM: ", FWHMx*step)
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
##