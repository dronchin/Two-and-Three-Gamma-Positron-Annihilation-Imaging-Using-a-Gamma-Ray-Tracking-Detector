 # -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 22:12:33 2017

@author: Nicolas Dronchi
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
from IPython.display import set_matplotlib_formats
from scipy import optimize
from matplotlib.colors import LogNorm

set_matplotlib_formats('pdf', 'svg') # this makes higher quality plots than the default (PNG)

plt.rcParams["figure.figsize"] = [10,5]

plot_type = 'heat1'


file1 = 'output_good.txt'
file2 = 'output_all_events.txt'

file3 = '3Phot_001_output.txt'
file4 = 'd:/Documents/School/Research/python 3 files/3Phot_002_output.txt'
file5 = 'd:/Documents/School/Research/python 3 files/3Phot_003_output.txt'

file6 = 'ThreePhotCo_001_output.txt'
file7 = 'ThreePhotCo_002_output.txt'
file8 = 'wrong1022Na_001_output.txt'
file9 = 'wrong1022Na_002_output.txt'

file10 = '3Phot_001_output_innercircle.txt'

file11 = 'output_good.txt'
file12 = 'output_all_events.txt'
file13 = 'wrong1022Na_001_output_with_orig_points.txt'
file14 = 'wrong1022Na_002_output_with_orig_points.txt'
file15 = '3Phot_001_output_with_orig_points.txt'
file16 = 'sim_uniformsource_with_orig_points.txt'
file17 = 'sim_source2_magicpowder_with_orig_points.txt'
file18 = 'testing_file_out.txt'

x = []
y = []
z = []
datapoints = []
good_energies =[]
with open('C:/data/3gamma data sets/{}'.format(file3),'r') as data:
    for line in data:
        #unpacks the data and stores them in a list
        numbers_float = map(float, line.split())
        point = list(numbers_float)
        
        # [x, y, z, det1x, det1y, det1z, det1E, det2x, det2y, det2z, det2E, det3x, det3y, det3z, det3E]
        
#        if point[2] > 10 or point[2] < -10:
#            continue
        x.append(point[0])# + np.random.normal(0,scale=5.4))
        y.append(point[1])# + np.random.normal(0,scale=4.8))
        z.append(point[2])# + np.random.normal(0,scale=3))
        
#        datapoints.append(point)
#        
#        r = np.sqrt(point[0]**2 + point[1]**2 + point[2]**2)
#        x.append(point[3])
#        y.append(point[4])
#        z.append(point[5])
#        
#        x.append(point[7])
#        y.append(point[8])
#        z.append(point[9])
#        
#        x.append(point[11])
#        y.append(point[12])
#        z.append(point[13])
#        good_energies.append(point[6] + point[10] + point[14]-0.001)
#        
        
#plt.hist(good_energies1, alpha = 0.3, bins = 20)
#plt.hist(good_energies2, alpha = 0.3, bins = 20)
##plt.hist(good_energies3, alpha = 0.3, bins = 20)
##plt.show()
#xa = []
#ya = []
#za = []
#xb = []
#yb = []
#zb = [] 
#bad_energies = []
##datapoints2 = []
#
#with open('C:/data/3gamma data sets/{}'.format(file13),'r') as data:
#    for line in data:
#        #unpacks the data and stores them in a list
#        numbers_float = map(float, line.split())
#        point = list(numbers_float)
#         
#        r = np.sqrt(point[0]**2 + point[1]**2 + point[2]**2)
#        if r > 20:
#            xa.append(point[3])
#            ya.append(point[4])
#            za.append(point[5])
#            
#            xa.append(point[7])
#            ya.append(point[8])
#            za.append(point[9])
#            
#            xa.append(point[11])
#            ya.append(point[12])
#            za.append(point[13])
#        if r < 20:
#            xb.append(point[3])
#            yb.append(point[4])
#            zb.append(point[5])
#            
#            xb.append(point[7])
#            yb.append(point[8])
#            zb.append(point[9])
#            
#            xb.append(point[11])
#            yb.append(point[12])
#            zb.append(point[13])
#        bad_energies.append(point[6] + point[10] + point[14]-0.001)
#
##        detects.append(detect1)
#        detects.append(detect2)
#        detects.append(detect3)
#        energies.append(energy1)
#        energies.append(energy2)
#        energies.append(energy3)
#        datapoints2.append(point)
        

if plot_type == 'gretinaheat':
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
    Heatmapit(x, z, "all from good dataset")
    Heatmapit(xa, za, "outter wrong energy")
    Heatmapit(xb, zb, "inner wrong energy")
    
    plt.hist(good_energies, bins=np.arange(1017,1030))
    plt.show()
    plt.hist(bad_energies, bins=np.arange(1027,1035))
    plt.show()
    

if plot_type == '3d': 
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
    
if plot_type == 'heat1':
    
    lim = 20
    step = 1
    
    xr = np.array(x)
    yr = np.array(y)
    zr = np.array(z)
    xmin = -lim #xr.min()
    xmax = lim #xr.max()
    ymin = -lim #yr.min()
    ymax = lim #yr.max()
    xbin = np.arange(xmin,xmax+step,step)
    ybin = np.arange(ymin,ymax+step,step)
 
    plt.hist2d(x,y,bins=[xbin,ybin],cmap='Greys')
    plt.colorbar() #ticks=[0,1,2,3])
    plt.grid(color='k', alpha=0.15, linestyle='--')
    plt.axes().set_aspect('equal')
    plt.xlabel('x (mm)')
    plt.ylabel('z (mm)')
    plt.title('Error Added before Reconstruction')

#    plt.savefig("Error_added_before_Reconstruction.png", dpi=1200, format='png' )

    plt.show()
    
if plot_type == 'heat2':
    
    lim = 220
    step = 7.5
    
    xr = np.array(xa)
    yr = np.array(ya)
    zr = np.array(za)
    xmin = -lim #xr.min()
    xmax = lim #xr.max()
    ymin = -lim #yr.min()
    ymax = lim #yr.max()
    xbin = np.arange(xmin,xmax+step,step)
    ybin = np.arange(ymin,ymax+step,step)
    print(xbin)
 
    plt.hist2d(xa,ya,bins=[xbin,ybin],cmap='Greys')
    plt.colorbar() #ticks=[0,1,2,3])
    plt.grid(color='k', alpha=0.15, linestyle='--')
    plt.axes().set_aspect('equal')
    plt.xlabel('x (mm)')
    plt.ylabel('y (mm)')
    plt.title('Wrong energy Na Reconstruction')
    
#    plt.savefig("WrongNaE_002.png", dpi=1200, format='png' )

    plt.show()
    
    suspicious = []
    bad = []
    count = 0
    for i in range(len(xr)):
        r = np.sqrt(xr[i]**2 + yr[i]**2 + zr[i]**2)
        if r < 20:
            suspicious.append(datapoints2[i])
            count += 1
        else:
            bad.append(datapoints2[i])
    print("center: ", count)
    
    detects = []
    sus_energies1 = []
    sus_energies2 = []
    for point in suspicious:
        detect1 = point[3:6]
        energy1 = point[6]
        detect2 = point[7:10]
        energy2 = point[10]
        detect3 = point[11:14]
        energy3 = point[14]

        sus_energies1.append(energy1)
        sus_energies2.append(energy2)
    
    bad_energies1 = []
    bad_energies2 = []
    for point in bad:
        detect1 = point[3:6]
        energy1 = point[6]
        detect2 = point[7:10]
        energy2 = point[10]
        detect3 = point[11:14]
        energy3 = point[14]

        bad_energies1.append(energy1)
        bad_energies2.append(energy2)
    
    
    plt.hist2d(bad_energies1,bad_energies2, bins=[np.linspace(0,511,30), np.linspace(0,511,30)], cmap='Greys')
    plt.plot(sus_energies1, sus_energies2, "ro")
    plt.show()
    
    plt.hist2d(good_energies1,good_energies2, bins=[np.linspace(0,511,30), np.linspace(0,511,30)], cmap='Greys')
    plt.plot(sus_energies1, sus_energies2, "ro")
    plt.show()
    
if plot_type == 'heat3':
    
    lim = 240
    step = 7.5
    
    xr = np.array(x)
    yr = np.array(y)
    zr = np.array(z)
    xmin = -lim #xr.min()
    xmax = lim #xr.max()
    ymin = -lim #yr.min()
    ymax = lim #yr.max()
    xbin = np.arange(xmin,xmax+step,step)
    ybin = np.arange(ymin,ymax+step,step)
 
    plt.hist2d(x,y,bins=[xbin,ybin],cmap='Greys')
    plt.colorbar() #ticks=[0,1,2,3])
    plt.grid(color='k', alpha=0.15, linestyle='--')
    plt.axes().set_aspect('equal')
    plt.xlabel('x (mm)')
    plt.ylabel('y (mm)')
    plt.title('XY Projection of Source 2 Three-gamma Events')
    
#    plt.savefig("Fig10-XYSource2.png", dpi=1200, format='png' )

    plt.show()    

if plot_type == 'kde':
    xr = np.array(x)
    yr = np.array(y)
    zr = np.array(z)
    xmin = xr.min()/10
    xmax = xr.max()/10
    ymin = yr.min()/10
    ymax = yr.max()/10
    zmin = zr.min()/10
    zmax = zr.max()/10
    
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(6, 9), sharex=False, sharey=False) #sets up the sublots 3,1 is needed

    cmap = sns.cubehelix_palette(start=0, light=1, as_cmap=True) #sets up the color palet

    sns.kdeplot(xr, yr, cmap=cmap, shade=True, cut=5, ax=ax1) #uses seaborn to plot the kernel density in 2d
    ax1.set(xlim=(xmin, xmax), ylim=(ymin, ymax))
    ax1.set(xlabel='X (mm)', ylabel='Y (mm)')
    sns.kdeplot(yr, zr, cmap=cmap, shade=True, cut=5, ax=ax2)
    ax2.set(xlim=(ymin, ymax), ylim=(zmin, zmax))
    ax2.set(xlabel='Y (mm)', ylabel='Z (mm)')
    sns.kdeplot(xr, zr, cmap=cmap, shade=True, cut=5, ax=ax3)
    ax3.set(xlim=(xmin, xmax), ylim=(zmin, zmax))
    ax3.set(xlabel='X (mm)', ylabel='Z (mm)')
    
    plt.savefig('wrong1022Na_002_output_kde.png', dpi=1200, format='png')


if plot_type == "quantify-z":
    lim = 30
    step = 2
    length = int(lim/step*2)
    
    xr = np.array(x)
    yr = np.array(y)
    zr = np.array(z)
    xmin = -lim #xr.min()
    xmax = lim #xr.max()
    zmin = -lim #yr.min()
    zmax = lim #yr.max()
    xbin = np.arange(xmin+3,xmax+step+3,step)
    zbin = np.arange(zmin,zmax+step,step)    
    
    grid, _, _ = np.histogram2d(x, z, bins=[xbin,zbin])
    grid = grid.T
    plt.imshow(grid, cmap="inferno") #grid[120:220,100:200])
    plt.grid(False)
    
    xpos = 15
    zpos = 15
    
    x = range(length)
    z = [zpos for i in range(length)]
    plt.plot(x,z)
    z = range(length)
    x = [xpos for i in range(length)]
    plt.plot(x,z)
    plt.show()
    
    
    x = np.linspace(0,length,length)
    xdata = np.zeros(length)
    linex = []
    rangesval = []
    for i in range(3):
        for j in range(grid.shape[0]):
            linex.append(grid[ypos-1+i,j])
            rangesval.append(x[j])

    z = np.linspace(0,length,length)
    zdata = np.zeros(length)
    linez = []
    zrangesval = []
    for i in range(3):
        for j in range(grid.shape[1]):
            linez.append(grid[j,xpos-1+i])
            zrangesval.append(z[j])
            
    def gaussian(x, amp, cen, wid): #,back):
        return amp * np.exp(-(x-cen)**2 / wid) #+back #runtime warning when it gets close to the value of zero
       
    init_vals = [1, 15, 1]  # for [amp, cen, wid, back]
    best_vals, covar = optimize.curve_fit(gaussian, rangesval, linex, p0=init_vals)
   
    stdx = np.sqrt(best_vals[2]/2)
    FWHMx = 2.355*stdx
    reg = 6
    cent = best_vals[1]
   
    plt.plot(rangesval, linex, "bo")
    plt.plot(x, gaussian(x, best_vals[0],best_vals[1], best_vals[2]), "r")
    plt.show()

    print("standard Deviation: ", stdx*step)
    print("FWHM: ", FWHMx*step)

    init_vals = [1, 15, 1]  # for [amp, cen, wid, back]
    best_vals, covar = optimize.curve_fit(gaussian, zrangesval, linez, p0=init_vals)
    
    stdz = np.sqrt(best_vals[2]/2)
    FWHMz = 2.355*stdz
    reg = 6
    cent = best_vals[1]
    
    plt.plot(zrangesval, linez, "bo")
    plt.plot(z, gaussian(z, best_vals[0],best_vals[1], best_vals[2]), "r")
    plt.show()
 
    print("standard Deviation: ", stdz*step)
    print("FWHM: ", FWHMz*step)

if plot_type == "quantify":
    lim = 30
    step = 2
    length = int(lim/step*2)
    
    xr = np.array(x)
    yr = np.array(y)
    zr = np.array(z)
    xmin = -lim #xr.min()
    xmax = lim #xr.max()
    ymin = -lim #yr.min()
    ymax = lim #yr.max()
    xbin = np.arange(xmin+2,xmax+step+2,step)
    ybin = np.arange(ymin-2,ymax+step-2,step)    
    
    grid, _, _ = np.histogram2d(x, y, bins=[xbin,ybin])
    grid = grid.T
    plt.imshow(grid, cmap="inferno") #grid[120:220,100:200])
    plt.grid(False)
    
    xpos = 15
    ypos = 15
    
    x = range(length)
    y = [ypos for i in range(length)]
    plt.plot(x,y)
    y = range(length)
    x = [xpos for i in range(length)]
    plt.plot(x,y)
    plt.show()
    
    
    x = np.linspace(0,length,length)
    xdata = np.zeros(length)
    linex = []
    rangesval = []
    for i in range(3):
        for j in range(grid.shape[0]):
            linex.append(grid[ypos-1+i,j])
            rangesval.append(x[j])
    xdata = xdata/3
#        for j in range(grid.shape[0]):
#            linex.append(grid[xpos-1+i,j])
#            rangesval.append(x[j])
            
#    plt.plot(x, xdata, "bo")
#    plt.show()
    y = np.linspace(0,length,length)
    ydata = np.zeros(length)
    liney = []
    yrangesval = []
    for i in range(3):
        for j in range(grid.shape[1]):
            liney.append(grid[j,xpos-1+i])
            yrangesval.append(y[j])
            
    def gaussian(x, amp, cen, wid): #,back):
        return amp * np.exp(-(x-cen)**2 / wid) #+back #runtime warning when it gets close to the value of zero
       
    init_vals = [1, 15, 1]  # for [amp, cen, wid, back]
    best_vals, covar = optimize.curve_fit(gaussian, rangesval, linex, p0=init_vals)
   
    stdx = np.sqrt(best_vals[2]/2)
    FWHMx = 2.355*stdx
    reg = 6
    cent = best_vals[1]
   
    plt.plot(rangesval, linex, "bo")
    plt.plot(x, gaussian(x, best_vals[0],best_vals[1], best_vals[2]), "r")
    plt.show()

    print("standard Deviation: ", stdx*step)
    print("FWHM: ", FWHMx*step)

    init_vals = [1, 15, 1]  # for [amp, cen, wid, back]
    best_vals, covar = optimize.curve_fit(gaussian, yrangesval, liney, p0=init_vals)
    
    stdy = np.sqrt(best_vals[2]/2)
    FWHMy = 2.355*stdy
    reg = 6
    cent = best_vals[1]
    
    plt.plot(yrangesval, liney, "bo")
    plt.plot(y, gaussian(y, best_vals[0],best_vals[1], best_vals[2]), "r")
    plt.show()
 
    print("standard Deviation: ", stdy*step)
    print("FWHM: ", FWHMy*step)








#
