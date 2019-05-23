# -*- coding: utf-8 -*-
"""
Created on Tue Jul 18 11:38:10 2017

@author: dronchin
"""
from pyqtgraph.Qt import QtCore, QtGui
import pyqtgraph.opengl as gl
import numpy as np
from matplotlib import cm

#initializes the application window
app = QtGui.QApplication([])
w = gl.GLViewWidget()
w.opts['distance'] = 80 #starting view dist
w.show()
w.setWindowTitle('3D Back Projection with slider')

#adds a grid to the application
b = gl.GLBoxItem()
w.addItem(b)
g = gl.GLGridItem()
g.setSize(100,100,100)
g.scale(1, 1, 1)
w.addItem(g)

#loads data set
grid = np.load("C:/data/grid.npy")
data = np.copy(grid) #makes sure we don't change the data
grid = [] #saves memory\

#trims further to reduce calculations and speed up updating 
#print(data.shape)
beg = 200
end = 500
data = data[beg:(end-beg), beg:(end-beg), beg:(end-beg)]

#sets slider variables. slider is the start and sliderchange is the difference
#between each update
slider = 0.2
sliderchange = 0.01

#all data manipulation
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

#creates 3d volume item
v = gl.GLVolumeItem(data=colors, glOptions='translucent')  #translucent or additive

#brings 3d graph to the center
translator = data.shape
v.translate(-translator[0]//2, -translator[1]//2, -translator[2]//2) # - 50 )

w.addItem(v)
ax = gl.GLAxisItem()
w.addItem(ax)


def update():
    global slider, sliderchange
    if slider < 0.7: #creates errors with the exp values if >1    
        slider += sliderchange
    else:
        slider=0.2
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
    
    v.setData(data=colors)

#sets the updating function to work    
t = QtCore.QTimer()
t.timeout.connect(update)
t.start(100) #the time interval in ms 

#starts the application
QtGui.QApplication.instance().exec_()

