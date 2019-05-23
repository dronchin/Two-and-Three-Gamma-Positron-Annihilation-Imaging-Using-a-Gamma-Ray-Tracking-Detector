# -*- coding: utf-8 -*-
"""
Created on Thu May 23 11:23:08 2019

@author: dronchin
This file shows an example of adding plotting functionality to the main ThreeGammaImaging class 
through class inheritance. Here a 3D plot method is added.
"""

from ThreeGammaSolver import ThreeGammaEvent, ThreeGammaImaging
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#create a new clas that inherits the methods of ThreeGammaImaging but adds a new method.
class ThreeGamma3DPlot(ThreeGammaImaging):
    def plot3D(self):
        '''
        '''
        x = []
        y = []
        z = []
        for ans in self.succeeded:
            x.append(ans.PO_decay[0])
            y.append(ans.PO_decay[1])
            z.append(ans.PO_decay[2])
        
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(x,y,z, c='red', marker='.')
        
        plt.title('3 gamma reconstruction', fontsize=22)
        ax.set_xlabel('X pos (mm)')
        ax.xaxis.label.set_size(15)
        ax.set_ylabel('Y pos (mm)')
        ax.yaxis.label.set_size(15)
        ax.set_zlabel('Z pos (mm)')
        ax.zaxis.label.set_size(15)
        
        #ax.view_init(90, 0)
        #plt.savefig("wrong1022Na_002_output_3d.png", dpi=1200, format='png' )
        plt.show() 


file4 = '3phot_001_output.txt'

Reconstructed = ThreeGamma3DPlot()
Reconstructed.import_previous_run('C:/data/3gamma data sets/{}'.format(file4))

Reconstructed.plot3D()







#