# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize




# increase x to go down
# increase y to go right

pos = [145,179] #12oclock
#pos = [215,179] #6oclock
#pos = [183,213] #3oclock
#pos = [181,144] #9oclock
#pos = [181,172] #center


#pos = [175,175] #z axis center
#pos = [196, 175] #z axis downstream
#pos = [175,172] #z axis 3oclock
#pos = [175, 138] #z axis 6oclock
#pos = [175,177] #z axis 9oclock
#pos = [172,210] #z axis 12oclock


reconstruction = np.load("reconstruction2.npy")
plt.imshow(reconstruction)
plt.plot(pos[1], pos[0], "r.")
plt.plot(179,179, "b.")
plt.show()


pad = 12
length = pad*2
step = 1

grid = reconstruction[pos[0]-pad:pos[0]+pad, pos[1]-pad:pos[1]+pad]
np.save("picofpoint.npy", grid)
plt.imshow(grid)
plt.plot([pad],[pad],"ro")
plt.show()

origin = [grid.shape[0]//2, grid.shape[1]//2]


x = np.linspace(0,length,length)
xdirection = []
xline = []
for i in range(3):
    print(i)
    for j in range(grid.shape[0]):
        xdirection.append(grid[origin[0]-1+i,j])
        xline.append(x[j])
y = np.linspace(0,length,length)
ydirection = []
yline = []
for i in range(3):
    print(i)
    for j in range(grid.shape[0]):
        ydirection.append(grid[j,origin[1]-1+i])
        yline.append(y[j])
        
def gaussian(x, amp, cen, wid, back): #,back):
    return amp * np.exp(-(x-cen)**2 / wid) + back #runtime warning when it gets close to the value of zero

## X Direction
init_vals = [1, 15, 1, 0.1]  # for [amp, cen, wid, back]
best_vals, covar = optimize.curve_fit(gaussian, xline, xdirection, p0=init_vals)

stdx = np.sqrt(best_vals[2]/2) 
FWHMx = 2.355*stdx
reg = 6
cent = best_vals[1]

plt.plot(xline, xdirection, "bo")
plt.plot(x, gaussian(y, best_vals[0],best_vals[1], best_vals[2], best_vals[3]), "r")
plt.show()

## Y Direction
init_vals = [1, 15, 1, 0.3]  # for [amp, cen, wid, back]
best_vals, covar = optimize.curve_fit(gaussian, yline, ydirection, p0=init_vals)

stdy = np.sqrt(best_vals[2]/2)
FWHMy = 2.355*stdy
reg = 6
cent = best_vals[1]
print(best_vals[3])

plt.plot(yline, ydirection, "bo")
plt.plot(y, gaussian(y, best_vals[0],best_vals[1], best_vals[2], best_vals[3]), "r")
plt.show()
 
print("standard Deviation: ", stdx*step)
print("FWHMx: ", FWHMx*step)
print("standard Deviation: ", stdy*step)
print("FWHMy: ", FWHMy*step)


#