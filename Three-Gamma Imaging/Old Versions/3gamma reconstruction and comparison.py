# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 19:02:29 2017

@author: Nicolas Dronchi

Program for 3-gamma reconstruction.
The objective is to take in detection and energy data points and solve for the 
inital position that the positronim was in before it decays.
The data was taken at the NSCL using the Gretina detector. 
"""

from math import sin,cos,degrees,pi
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
import matplotlib.ticker as plticker


#takes in the data and stores them as variables
#inital data in (x,y,z) with energies *3
def initial_points_old(x1,y1,z1,e1,x2,y2,z2,e2,x3,y3,z3,e3): 
    global detect1,detect2,detect3,energy1,energy2,energy3

    #energies divided by 1000
    detect1 = np.array([x1,y1,z1])
    energy1 = float(e1)/100
    detect2 = np.array([x2,y2,z2])
    energy2 = float(e2)/100
    detect3 = np.array([x3,y3,z3])
    energy3 = float(e3)/100
    

def initial_points(data): 
    ''' takes in the data and stores them as variables
    data[i] -> 0-x1;1-y1;2-z1;3-e1;4-x2;5-y2;6-z2;7-e2;8-x3;9-y3;10-z3;11-e3
    '''

    #energies divided by 1000
    detect1 = np.array(data[0:3])
    energy1 = data[3]
    detect2 = np.array(data[4:7])
    energy2 = data[7]
    detect3 = np.array(data[8:11])
    energy3 = data[11]
    
    return [detect1,detect2,detect3,energy1,energy2,energy3]
    

def angles(e1,e2,e3):
    '''finds the angles of the unique trianlge created by the lengths of the energies
    Inputs = (e1,e2,e3) where e1 = energy1, e2 = energy2, and e3 = energy3
    '''
    #global alpha,beta,gamma
    alpha = np.arccos((e2**2-e1**2-e3**2)/(-2*e1*e3))
    beta = np.arccos((e3**2-e1**2-e2**2)/(-2*e1*e2))
    gamma = np.arccos((e1**2-e2**2-e3**2)/(-2*e2*e3))
    #prints out the 3 angles (Rad, Deg)
    return [alpha,beta,gamma]

    
#finds a normal vector to the plane that the detections lie in
def findnorm():
    #global vect12,vect13,vect23,vectnorm
    vect12 = detect2-detect1
    vect13 = detect3-detect1
    vect23 = detect3-detect2

    vectnorm = .005*np.cross(vect12,vect13)
    return [vect12,vect13,vect23,vectnorm]
   
    
#these functions are used to rotate any point around a given axis as long as
#it's a unit vector
def R(theta, u): #theta = amount of angle rotated, u = given axis (unit vect)
    return [[cos(theta) + u[0]**2 * (1-cos(theta)), 
             u[0] * u[1] * (1-cos(theta)) - u[2] * sin(theta), 
             u[0] * u[2] * (1 - cos(theta)) + u[1] * sin(theta)],
            [u[0] * u[1] * (1-cos(theta)) + u[2] * sin(theta),
             cos(theta) + u[1]**2 * (1-cos(theta)),
             u[1] * u[2] * (1 - cos(theta)) - u[0] * sin(theta)],
            [u[0] * u[2] * (1-cos(theta)) - u[1] * sin(theta),
             u[1] * u[2] * (1-cos(theta)) + u[0] * sin(theta),
             cos(theta) + u[2]**2 * (1-cos(theta))]]
#poitToRotate = point being rotated, around = axis of rotation, 
#theata=amount of angle rotated
def Rotate(pointToRotate, around, theta): 
    rotate_matrix = np.array(R(theta, around))
    rotated = np.matmul(rotate_matrix,pointToRotate)
    return rotated


#used to find the thrid points of the circle based off of the angles
#takes in 2 positions =A,B and uses angle= angle of triangle inscibed in circle
def thirdpoint(A,B,angle): 
    #length of line PR, useful for trig
    LlinePR = np.sqrt((B[0]-A[0])**2 + (B[1]-A[1])**2 + (B[2]-A[2])**2)

    #length of line PQ, the hypotenus for the right triangle formed by spliting a isosceles triangle in half
    LlinePQ = (LlinePR/2)/sin(angle/2)

    #norm of vector PQ
    ax1 = (1/LlinePR)*(B-A)
    
    #rotate the norm of vecPQ, rotate it 
    vectnorm_unit = np.array(vectnorm)/np.sqrt(vectnorm[0]**2+vectnorm[1]**2+vectnorm[2]**2)
    ax2 = np.array(Rotate(ax1,vectnorm_unit,(pi-(angle))/2))
    
    #translates point back to it's position
    p3 = A+LlinePQ*ax2
    p3_list.append(p3)
    

#changes the coordinate system that the points are in and saves a way to go
#back and forth by multiplying a matrix. New basis = B, old basis = E
def changebasisE_B():
    #global detect1_E,detect2_E,detect3_E,detect1_B,detect2_B,detect3_B,basis_B_inversed,basis_B
    #translation from original detect points to the origin
    detect1_E = detect1-detect1
    detect2_E = detect2-detect1
    detect3_E = detect3-detect1
    detect4_E = vectnorm

    #calculations using QR factorization in order to create an orthonormal basis
    basis_B = np.linalg.qr([[detect2_E[0],detect3_E[0],detect4_E[0]],
                               [detect2_E[1],detect3_E[1],detect4_E[1]],
                               [detect2_E[2],detect3_E[2],detect4_E[2]]])
    
    #the inverse of the O-N basis is the transition matrix for points
    basis_B_inversed = np.linalg.inv(basis_B[0])
    #print("ortho normal basis: ",basis_B[0])
    #print("converts between basis E and B: ",basis_B_inversed)

    #converts detect points in basis E to the basis B
    detect1_B = np.dot(basis_B_inversed, detect1_E) 
    detect2_B = np.dot(basis_B_inversed, detect2_E) 
    detect3_B = np.dot(basis_B_inversed, detect3_E) 
    #print("3 points in new basis B:", detect1_B,detect2_B,detect3_B)
    return detect1_E,detect2_E,detect3_E,detect1_B,detect2_B,detect3_B,basis_B_inversed,basis_B
    
    

#converts all of the 3rd points into the new basis labeled B   
def p3_B():
    p3_list[0] -= detect1
    p3_list[1] -= detect1
    p3_list[2] -= detect1
    #doting the points on the basis_B_inversed gives you the points in B
    p3_list[0] = np.dot(basis_B_inversed, p3_list[0])
    p3_list[1] = np.dot(basis_B_inversed, p3_list[1])
    p3_list[2] = np.dot(basis_B_inversed, p3_list[2])


#Takes in 3 points A,p3,B and calculates the equation of a circle. Saves them
#in the list called DEF. x^2 + y^2 + Dx + Ey + F = 0
def cirlce(A,p3,B):
    #takes a matrix and dots it to find the solution
    DEFcoef = np.array([[A[0],A[1],1],
                        [p3[0],p3[1],1],
                        [B[0],B[1],1]])
    DEFans = np.array([-(A[0])**2-(A[1])**2, 
                       -(p3[0])**2-(p3[1])**2,
                        -(B[0])**2-(B[1])**2])
    DEF = np.linalg.solve(DEFcoef, DEFans)

    #stores the values for later use
    cirpara.append(DEF)    
    
    
#finds the intersection point of the 2 circles and tests which one is the point
#where all 3 interesct. function takes in the circle parameters for 2 circles    
def cirintersect(D1,E1,F1,D2,E2,F2):
    #solves for a,b,c in the quadratic formula
    QuadA = ((-E1+E2)/(D1-D2))**2 + 1
    QuadB = 2*((-E1+E2)/(D1-D2) * (-F1+F2)/(D1-D2))+ D1*(-E1+E2)/(D1-D2) + E1
    QuadC = ((-F1+F2)/(D1-D2))**2 + D1*(-F1+F2)/(D1-D2) + F1
    
    #the disc is the part of the quadratic formula under the sqrt
    disc = QuadB**2 - 4*QuadA*QuadC
    #tests for imaginary numbers
    if disc > 0:
        #solves for x,y of the quadratic
        QuadansY1 = (-QuadB + np.sqrt(disc))/(2*QuadA)
        QuadansY2 = (-QuadB - np.sqrt(disc))/(2*QuadA)
        QuadansX1= (-F1+F2)/(D1-D2) + (-E1+E2)/(D1-D2)*QuadansY1
        QuadansX2= (-F1+F2)/(D1-D2) + (-E1+E2)/(D1-D2)*QuadansY2
        
        radius_srd = -cirpara[0][2]+(cirpara[0][0]/2)**2+(cirpara[0][1]/2)**2
        radius_srd -= .01* radius_srd
         
        radius_test1 = (QuadansX1 + ((cirpara[0][0]) / 2))**2 + (QuadansY1 + (cirpara[0][1] / 2))**2
        radius_test2 = (QuadansX2 + ((cirpara[0][0]) / 2))**2 + (QuadansY2 + (cirpara[0][1] / 2))**2       
        
        if radius_test1 < radius_srd:
            ans = (QuadansX1,QuadansY1,0)
            ans_list.append(ans)
            
            #print("Basis B ans:      ", ans)
        if radius_test2 < radius_srd:
            ans = (QuadansX2,QuadansY2,0)
            ans_list.append(ans)
            
            #print("Basis B ans:      ", ans)
        
    else:
        print("Sorry, error occurred with circle equations", end=' ')




#takes the answer given in the basis B and converts it back to the basis E    
def changebasisB_E():
    global success_num, failed_num0, errors
    if len(ans_list)>2:
        success_num += 1
        
        if (np.allclose(ans_list[0], ans_list[1]) and np.allclose(ans_list[1], ans_list[2])):
            pass
        else:
            errors += 1
        
        #averages the 3 answers
        ans_B = (np.asarray(ans_list[0])
                +np.asarray(ans_list[1])
                +np.asarray(ans_list[2]))/3
        ans_E = np.dot(basis_B[0], ans_B) + detect1
        return ans_E
    
    else:
        failed_num0 += 1
        return [50,50,50]

def MSE(A,B):
    total = (np.square(A - B)).mean(axis=None)
    return total
    
def gaussian(x, amp, cen, wid): #,back):
    return amp * np.exp(-(x-cen)**2 / wid) #+back #runtime warning when it gets close to the value of zero

def FWHMcalc(grid, step=1):
    
    length = grid.shape[1]
    x = np.linspace(0,length,length)
    xdata = np.zeros(length)
    for i in range(3):
        xdata += grid[length//2 - 1 + i]
    xdata = xdata/3

       
    init_vals = [0.005, length/2, 0.1]  # for [amp, cen, wid, back]
    best_vals, covar = optimize.curve_fit(gaussian, x, xdata, p0=init_vals)
   
    stdx = np.sqrt(best_vals[2]/2)
    FWHMx = 2.355*stdx
   
    plt.plot(x, xdata, "bo")
    plt.plot(x, gaussian(x, best_vals[0],best_vals[1], best_vals[2]), "r")
    plt.show()

    return FWHMx*step

def radialdensity(rlist, bins):
    r_min = 0
    dens = []
    total_count = 0
    #loop through dr
    for r_max in bins:
        count = 0
        #find all points within dr 
        for r_cur in rlist:
            if (r_cur > r_min and r_cur < r_max):
                count += 1
                total_count += 1
        dens.append(count/(r_max**2 - r_min**2))
        r_min = r_max

    dens = np.asarray(dens)/float(total_count)
    print(total_count)    
    return dens

#==============================================================================
# MAIN PROGRAM
#==============================================================================

set_matplotlib_formats('pdf', 'svg') # this makes higher quality plots than the default (PNG)

plt.rcParams["figure.figsize"] = [10,5]



x_original = []
y_original = []
z_original = []
r_original = []
with open('C:/data/3gamma data sets/3Phot_001_output.txt','r') as data:
    for line in data:
        #unpacks the data and stores them in a list
        numbers_float = map(float, line.split())
        point = list(numbers_float)
        
        x_original.append(point[0])
        y_original.append(point[1])
        z_original.append(point[2])
        r_original.append(np.sqrt(point[0]**2 + (point[1]-3)**2))

lim = 20
step = 1

xmin = -lim #xr.min()
xmax = lim #xr.max()
ymin = -lim #yr.min()
ymax = lim #yr.max()
xbin = np.arange(xmin,xmax+step,step)
ybin = np.arange(ymin,ymax+step,step)

real_image, xedges, yedges = np.histogram2d(x_original, y_original, bins=[xbin,ybin])
plt.imshow(real_image.T)
plt.colorbar() #ticks=[0,1,2,3])
plt.grid(color='k', alpha=0.15, linestyle='--')
plt.axes().set_aspect('equal')
plt.xlabel('x (mm)')
plt.ylabel('y (mm)')
plt.title('Real image')
plt.show()

fwhmx_list = []
fwhmy_list = []

error_list = []
sigmaiter = np.arange(3,3.2,1)
#opens up ThreePhotFile_00_theta with all the data and reads line by line
for sigma in sigmaiter:
    final_ans = []
    initial = []
    success_num = 0
    failed_num0 = 0
    errors = 0
    #sim_source2_magicpowder
    #sim_source2_magicpowder_10mm
    with open('C:/data/3gamma data sets/sim_source2_magicpowder.txt','r') as data:
        for line in data:
            #unpacks the data and stores them in a list
            numbers_float = map(float, line.split())
            initial_points_data = list(numbers_float)
            
            
#            if np.random.random() < 0.3:
#                r_sqr = np.sqrt(np.random.uniform(0,15**2))
#                angle = np.random.uniform(0,np.pi*2)
#                z = np.random.uniform(-5,5)
#                
#                answer = [3+r_sqr*np.cos(angle), r_sqr*np.sin(angle), z]
#                final_ans.append(answer)
#                continue
#            
            #inital 3 lists for storing later values (all 3 reset for each data set)
            #stores values form the function: circle
            cirpara = []
            #stores values from the function: thirdpoint
            p3_list = []
            #stores values from the function: circleintersect
            ans_list = []
            
            
            detect1,detect2,detect3,energy1,energy2,energy3 = initial_points(initial_points_data)
            
            # add on error for sigma in the loop as well as 
#            if np.random.random() < 0.5:
##                detect1 += np.random.uniform(-18,18,3)
##                detect2 += np.random.uniform(-18,18,3)
##                detect3 += np.random.uniform(-18,18,3)
##                
#                detect1 += np.random.normal(0, scale=10, size=3)
#                detect2 += np.random.normal(0, scale=10, size=3)
#                detect3 += np.random.normal(0, scale=10, size=3)
#            else:
            detect1 += np.random.normal(0, scale=sigma, size=3)# + 1/4*np.random.normal(0, scale=sigma*4, size=3)
            detect2 += np.random.normal(0, scale=sigma, size=3)# + 1/4*np.random.normal(0, scale=sigma*4, size=3)
            detect3 += np.random.normal(0, scale=sigma, size=3)# + 1/4*np.random.normal(0, scale=sigma*4, size=3)
            # Add energy resolution error on the detections
            energy1 += np.random.normal(0, scale=0.85)
            energy2 += np.random.normal(0, scale=0.85)
            energy3 += np.random.normal(0, scale=0.85)
            
            alpha,beta,gamma = angles(energy1,energy2,energy3)
            vect12,vect13,vect23,vectnorm = findnorm()
            
            thirdpoint(detect1,detect2,(alpha+gamma))
            thirdpoint(detect2,detect3,(alpha+beta))
            thirdpoint(detect3,detect1,(beta+gamma))
    
            detect1_E,detect2_E,detect3_E,detect1_B,detect2_B,detect3_B,basis_B_inversed,basis_B = changebasisE_B()
            p3_B()
            #print("p3_list: ",p3_list)
            
            cirlce(detect1_B,detect2_B,detect3_B)
            cirlce(detect1_B,p3_list[0],detect2_B)
            cirlce(detect2_B,p3_list[1],detect3_B)
            cirlce(detect1_B,p3_list[2],detect3_B)
    
            cirintersect(cirpara[1][0],cirpara[1][1],cirpara[1][2],
                         cirpara[2][0],cirpara[2][1],cirpara[2][2])
            cirintersect(cirpara[1][0],cirpara[1][1],cirpara[1][2],
                         cirpara[3][0],cirpara[3][1],cirpara[3][2])
            cirintersect(cirpara[2][0],cirpara[2][1],cirpara[2][2],
                         cirpara[3][0],cirpara[3][1],cirpara[3][2]) 
            
            answer = changebasisB_E()
            final_ans.append(answer)
            initial.append(initial_points_data)

    
    # Work with the current list of answers here
    
    final_ans = np.asarray(final_ans)
    x = final_ans[:,0] 
    y = final_ans[:,1]
    

 
    sim_image, xedges, yedges = np.histogram2d(x, y, bins=[xbin,ybin])
    
    diff_image = sim_image - real_image
    
    
    plt.hist2d(x,y,bins=[xbin,ybin],cmap='Greys')
    plt.colorbar() #ticks=[0,1,2,3])
    plt.grid(color='k', alpha=0.15, linestyle='--')
    plt.axes().set_aspect('equal')
    plt.xlabel('x (mm)')
    plt.ylabel('y (mm)')
    plt.title(r'Simulation of Three-gamma Events $\sigma_{pos} = 3$')

    plt.savefig("Simulation_sig_3mm.png", dpi=1200, format='png' )

    plt.show()
    
#    fig,ax=plt.subplots()
#
#    im1 = ax.imshow(sim_image.T, cmap='Greys')
#    f.colorbar(im1)
##    plt.axes().set_aspect('equal')
##    plt.xlabel('x (mm)')
##    plt.ylabel('y (mm)')
##    plt.title("Simulation at sig={}".format(np.round(sigma,1)))
#
#    loc = plticker.MultipleLocator(base=5)
#    ax.xaxis.set_major_locator(loc)
#    ax.yaxis.set_major_locator(loc)
#    ax.grid(which='major', color='k', alpha=0.15, linestyle='--')
#
#
#    plt.show()

#    plt.imshow(real_image.T, cmap='Greys')
#    plt.colorbar()
#    plt.grid(color='k', alpha=0.15, linestyle='--')
#    plt.axes().set_aspect('equal')
#    plt.xlabel('x (mm)')
#    plt.ylabel('y (mm)')
#    plt.title('Real image')
#    plt.show()
#    
#    f, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
#    f.set_figwidth(8)
#    
#    ax1.imshow(sim_image.T, cmap='Greys')
#    ax1.set_aspect('equal')
#    ax1.set_title("Simulation at sig={}".format(np.round(sigma,1)))
#    ax1.set_xlabel('x (mm)')
#    ax1.set_ylabel('y (mm)')
#    ax1.colorbar()
#    
#    ax2.imshow(real_image.T, cmap='Greys')
#    ax2.set_aspect('equal')
#    ax2.set_xlabel('x (mm)')
#    ax2.set_ylabel('y (mm)')
#    ax2.set_title("Source 2 reconstructions")
#    
#    ax3.imshow(diff_image.T)
#    ax3.set_aspect('equal')
#    ax3.set_xlabel('x (mm)')
#    ax3.set_ylabel('y (mm)')
#    ax3.set_title("Difference at sig={}".format(np.round(sigma,1)))
#    

    
    # looking at mse between the plots
    diff = MSE(sim_image, real_image)
    error_list.append(diff)

    #look at radial density plots
    r = []
    for ans in final_ans:
        r.append(ans[0]**2 + ans[1]**2)
    
    radial_bin = np.arange(2,40,2)

    density_original = radialdensity(r_original, radial_bin)
    density = radialdensity(r, radial_bin)
    
    plt.plot(radial_bin, density_original, label='original')
    plt.plot(radial_bin, density, label='sim')
#    plt.ylim(top = 0.02)
    plt.legend(loc='best')
    plt.show()    
    
    fwhmx = FWHMcalc(sim_image.T)
    print(fwhmx)
    fwhmx_list.append(fwhmx)
    
    fwhmy = FWHMcalc(sim_image)
    print(fwhmy)
    fwhmy_list.append(fwhmy)


fwhmx_list = np.asarray(fwhmx_list)
fwhmy_list = np.asarray(fwhmy_list)
z = np.polyfit(sigmaiter*2.355, fwhmx_list, 1)
p = np.poly1d(z)

print(z)

plt.plot(sigmaiter*2.355, fwhmx_list)
plt.plot(sigmaiter*2.355, fwhmy_list)
plt.plot(sigmaiter*2.355, p(sigmaiter*2.355))
plt.show()

#print("Average fwhmx: ", np.average(fwhmx_list), " sigma = ", np.average(fwhmx_list)/2.355)
#print("Average fwhmy: ", np.average(fwhmy_list), " sigma = ", np.average(fwhmy_list)/2.355)

#    
#
error_list = error_list/max(error_list)
plt.plot(sigmaiter, error_list)    
plt.show()
    
##    
