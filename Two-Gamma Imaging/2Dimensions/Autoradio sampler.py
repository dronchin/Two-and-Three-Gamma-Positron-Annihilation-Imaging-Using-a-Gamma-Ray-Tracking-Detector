# -*- coding: utf-8 -*-
"""
Created on Thu Jan 17 14:14:56 2019

@author: dronchin
"""

import numpy as np
import matplotlib.pyplot as plt


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

def getlocation(i, n, center, pixsize):
    row = i//n - n/2
    col = i%n - n/2
    xpos = center[0] + row*pixsize +np.random.normal(0,scale=pixsize/2)
    ypos = center[1] + col*pixsize+np.random.normal(0,scale=pixsize/2)
    
    return np.array([xpos,ypos,0])
    

source = getsource("old Na 16x16.csv")
plt.imshow(source.reshape((20,20)))
plt.show()


dist = np.zeros(20*20)
pointpos = [100,-50]
alldecays = []
for i in range(10000):
    index = np.random.choice(len(source), 1, p=source)
    dist[index] += 1
    POdecay = getlocation(index, 20, pointpos, 16*0.025)
    alldecays.append(POdecay)

plt.imshow(dist.reshape((20,20)))
plt.show()

alldecays = np.asarray(alldecays)
plt.plot(alldecays[:,0],alldecays[:,1],"b,")
plt.xlim(95,105)
plt.ylim(-45,-55)
plt.show()