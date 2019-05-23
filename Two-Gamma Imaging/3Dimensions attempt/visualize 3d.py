# -*- coding: utf-8 -*-
"""
Created on Mon Jul 17 14:20:47 2017

@author: dronchin
This program takes in a 3d array created from 3d back projectin and 
"""

from pyqtgraph.Qt import QtCore, QtGui
import pyqtgraph.opengl as gl
import pyqtgraph as pg
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from OpenGL.GL import *


#initializing the application window
app = QtGui.QApplication([])
w = gl.GLViewWidget()
w.opts['distance'] = 200
w.show()
w.setWindowTitle('3D Back Projection')
w.setBackgroundColor('k')

#adding a grid to determine distance to the window
b = gl.GLBoxItem()
w.addItem(b)
g = gl.GLGridItem()
#g = Grid(color=(0,0,0), size=None) 
g.setSize(100,100,100)
g.scale(5, 5, 5)
w.addItem(g)

slider = 0.45
def unpad(grid, radius):
    '''takes in the padded grid and cutts it back to the unpadded size
    '''
    beg = int(np.ceil(radius))
    end = grid.shape
    grid = grid[beg:(end[0]-beg), beg:(end[1]-beg), beg:(end[2]-beg)]
    return grid


# DATA ########################################################################
#make sure data is saved on ssd or "c:/" for fastest speeds
grid = np.load("C:/data/grid.npy")
print("working")
gridshape = grid.shape
origin = [gridshape[0]//2, gridshape[1]//2, gridshape[2]//2]

grid = grid.transpose((1,0,2))
#grid = grid[:,:,::-1]

#shift = [102, 102, -16]
#
#boxwidth = 100
#
##(origin[0] - shift[0]) - boxwidth//2  
#
#
#grid = grid[(origin[0] + shift[0]) - boxwidth//2:(origin[0] + shift[0]) + boxwidth//2,
#     (origin[1] + shift[1]) - boxwidth//2:(origin[1] + shift[1]) + boxwidth//2,
#     (origin[2] + shift[2]) - boxwidth//2:(origin[2] + shift[2]) + boxwidth//2]

#grid = grid[(origin[0] + shift[0]) - 30:(origin[0] + shift[0]) + 30,
#     (origin[1] + shift[1]) - 50:(origin[1] + shift[1]) + 50,
#     (origin[2] + shift[2]) - 30:(origin[2] + shift[2]) + 30]
##
data = np.copy(grid) #makes sure we don't change the data
###
###
#xpos = len(grid)//2
#ypos = len(grid[0])//2
#zpos = len(grid[0][0])//2
#xshift = 102
#yshift = 102
#zshift = -16
#for i in range(grid.shape[0]):
#    data[xpos+xshift,i,zpos+zshift] = grid.max()
#    data[i,ypos+yshift,zpos+zshift] = grid.max()
#    data[xpos+xshift,ypos+yshift,i] = grid.max()
#    
#line = []
#xpos = len(grid)//2
#ypos = len(grid[0])//2
#zpos = len(grid[0][0])//2
#xshift = -5
#yshift = 10
#zshift = -2
#for i in range(grid.shape[0]):
#    grid[xpos+xshift,i,zpos+zshift] = grid.max()



#data = np.copy(grid) #makes sure we don't change the data
#grid = [] #saves memory


#data = unpad(data, 70)

#data = data[111:171,111:171,175:205]

#data = data[111:171,111:171,111:171]

dataline = np.reshape(data, data.size) #puts data in a 1D array
dataline = dataline / dataline.max() #rescales data to [0,1]
dataline = (1/(1+np.exp(-10*(dataline - slider)))) #optional data scaling

#adds color to each point in the 1d array.
cmap = cm.plasma #other useful color maps: (magma, plasma,viridis, hot, afmhot)
colors = cmap(dataline) #each point is replaced with [r,g,b,a]
colors = colors*255 #scales [0,1] to [0,255]
colors = np.reshape(colors, (data.shape[0],data.shape[1],data.shape[2],4 )) #changes to (x,y,z,4)

alpha = 255*(1/(1+np.exp(-100*(dataline - slider)))) #exagerates high values
colors[...,3] = alpha.reshape(data.shape) #turns into a 1d array to copy values over

##adds axis to edge of box
#colors[:, 0, 0] = [255,0,0,100] #xaxis = red
#colors[49:51,0,0] = [255,255,0,255]
#colors[99:101,0,0] = [255,255,0,255]
#colors[149:151,0,0] = [255,255,0,255]
#colors[199:201,0,0] = [255,255,0,255]
#
#colors[0, :, 0] = [0,255,0,100] #yaxis = green
#colors[0, 0, :] = [0,0,255,100] #zaxis = blue

###############################################################################
#creates 3d volume item
v = gl.GLVolumeItem(data = colors, glOptions='translucent')  #translucent or additive

#brings 3d graph to the center
translator = data.shape
v.translate(-translator[0]//2, -translator[1]//2, -translator[2]//2 - 8)# 95 )

##empties lists to save memory
#alpha = []
#data = []
#dataline = []

w.addItem(v)
ax = gl.GLAxisItem()
w.addItem(ax)


QtGui.QApplication.instance().exec_()

### Start Qt event loop unless running in interactive mode.
#if __name__ == '__main__':
#    import sys
#    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
#        QtGui.QApplication.instance().exec_()
