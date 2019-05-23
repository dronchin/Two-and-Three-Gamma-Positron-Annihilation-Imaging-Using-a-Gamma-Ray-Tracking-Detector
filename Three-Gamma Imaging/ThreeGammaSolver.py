# -*- coding: utf-8 -*-
"""
Created on Wed May 22 10:56:03 2019

@author: dronchin
Restructuring of the file "3gamma reconstruction 4-8" in a way that makes the program more versitile. Here each three
gamma point will be an object. This will still complete the primary goal of solving decay positions from thee-gamma
positronium annihilation events.
"""

import numpy as np
import matplotlib.pyplot as plt
from IPython.display import set_matplotlib_formats

class ThreeGammaEvent:
    ''' data object used for holding all information on a single ThreeGammaEvent. Also contains the methods invloved
    with solving for the point of decay. '''
    def __init__(self, data_line):
        #unpacks the data and stores them in a list
        numbers_float = map(float, data_line.split())
        parsed_data = list(numbers_float)
        
        #in the case of raw data
        if len(parsed_data) == 12:
            #sets initial values for the three photon decays
            self.det_pos1 = np.array(parsed_data[0:3])
            self.energy1 = parsed_data[3]
            self.det_pos2 = np.array(parsed_data[4:7])
            self.energy2 = parsed_data[7]
            self.det_pos3 = np.array(parsed_data[8:11])
            self.energy3 = parsed_data[11]
            
            self.three_answers = []
            self.PO_decay = np.array([180.0, 180.0, 180.0])
            
            self.good_reconstruction = True
        
        #in the case of opening from a file
        if len(parsed_data) == 3:
            self.PO_decay = np.array(parsed_data)
            self.good_reconstruction = True
            
        if len(parsed_data) == 15:
            self.PO_decay = np.array(parsed_data[0:3])
            self.good_reconstruction = True
            
            #sets initial values for the three photon decays
            self.det_pos1 = np.array(parsed_data[3:6])
            self.energy1 = parsed_data[6]
            self.det_pos2 = np.array(parsed_data[7:10])
            self.energy2 = parsed_data[10]
            self.det_pos3 = np.array(parsed_data[11:14])
            self.energy3 = parsed_data[14]
                
    def angles(self):
        ''' finds the angles of the unique trianlge created by the lengths of the energies Inputs = (e1,e2,e3) 
        where e1 = energy1, e2 = energy2, and e3 = energy3 '''
        self.alpha = np.arccos((self.energy2**2 - self.energy1**2 - self.energy3**2)/(-2*self.energy1*self.energy3))
        self.beta = np.arccos((self.energy3**2 - self.energy1**2 - self.energy2**2)/(-2*self.energy1*self.energy2))
        self.gamma = np.arccos((self.energy1**2 - self.energy2**2 - self.energy3**2)/(-2*self.energy2*self.energy3))
    
    def find_norm(self):
        ''' finds a normal vector to the plane that the detections lie in '''
        vect12 = self.det_pos2 - self.det_pos1
        vect13 = self.det_pos3 - self.det_pos1
        
        #cross product to find normal vector, then transformed into a unit vector
        self.norm = .005*np.cross(vect12,vect13)
        self.norm = self.norm / np.sqrt(self.norm[0]**2 + self.norm[1]**2 + self.norm[2]**2)
        

    def Rotate(self, pointToRotate, u, theta):
        ''' produces a matrix used to rotate an angle theta (amount of angle rotated) around a vector u (unit vector).
        This transformation vector is then applied to the pointToRotate vector producing a resulting point'''
        #build matrix
        rotate_matrix = np.array([[np.cos(theta) + u[0]**2 * (1 - np.cos(theta)), 
                 u[0] * u[1] * (1 - np.cos(theta)) - u[2] * np.sin(theta), 
                 u[0] * u[2] * (1 - np.cos(theta)) + u[1] * np.sin(theta)],
                [u[0] * u[1] * (1 - np.cos(theta)) + u[2] * np.sin(theta),
                 np.cos(theta) + u[1]**2 * (1 - np.cos(theta)),
                 u[1] * u[2] * (1 - np.cos(theta)) - u[0] * np.sin(theta)],
                [u[0] * u[2] * (1 - np.cos(theta)) - u[1] * np.sin(theta),
                 u[1] * u[2] * (1 - np.cos(theta)) + u[0] * np.sin(theta),
                 np.cos(theta) + u[2]**2 * (1 - np.cos(theta))]])
        #multiply to find the rotated vector
        rotated = np.matmul(rotate_matrix, pointToRotate)
        return rotated
    
    def find_third_point(self, A, B, angle):
        ''' Every set of two points and an angle can be used to make a circle where the angle is the inscribed angle. 
        To eventually find the equation of this circle, finding a third point that lays on the circle is useful. Here
        we pick a vector between the two detections, then rotate this vector around the vector normal to the detection
        plane of the three detections. How much it is rotated is dependent on the inscribed angle of the circle. '''
        
        #length of line PR, useful for trig
        LlinePR = np.sqrt((B[0]-A[0])**2 + (B[1]-A[1])**2 + (B[2]-A[2])**2)
    
        #length of line PQ, the hypotenus for the right triangle formed by spliting a isosceles triangle in half
        LlinePQ = (LlinePR/2)/np.sin(angle/2)
    
        #norm of vector PQ
        ax1 = (1/LlinePR)*(B-A)
        
        #rotate the norm of vecPQ, rotate it 
        ax2 = np.array(self.Rotate(ax1, self.norm, (np.pi-(angle))/2))
        
        #translates point back to it's position
        p3 = A + LlinePQ*ax2
        
        return p3
        
    def changebasisE_B(self):
        ''' changes the coordinate system that the points are in and saves a way to go back and forth by multiplying 
        a matrix. New basis = B, old basis = E
        '''
        #translation from original detect points to the origin
        detect1_E = self.det_pos1 - self.det_pos1
        detect2_E = self.det_pos2 - self.det_pos1
        detect3_E = self.det_pos3 - self.det_pos1
        detect4_E = self.norm
    
        #QR factorization is used to create an orthonormal basis
        self.basis_B, _ = np.linalg.qr([[detect2_E[0],detect3_E[0],detect4_E[0]],
                                   [detect2_E[1],detect3_E[1],detect4_E[1]],
                                   [detect2_E[2],detect3_E[2],detect4_E[2]]])
        
        #the inverse of the O-N basis is the transition matrix for points
        basis_B_inversed = np.linalg.inv(self.basis_B)
    
        #converts detect points in basis E to the basis B
        self.det_pos1_B = np.dot(basis_B_inversed, detect1_E) 
        self.det_pos2_B = np.dot(basis_B_inversed, detect2_E) 
        self.det_pos3_B = np.dot(basis_B_inversed, detect3_E) 
        
    def cirlce(self, A,p3,B):
        ''' Takes in 3 points A, p3, B and returns parameters for the equation of a circle. Saves these parameters 
        for the future use of solving for the intersections. 
        x^2 + y^2 + Dx + Ey + F = 0 --> Dx + Ey + F = x^2 + y^2 '''
        #builds matricies 
        #DEFcoef * [D, E, F] = DEFans
        DEFcoef = np.array([[A[0],A[1],1],
                            [p3[0],p3[1],1],
                            [B[0],B[1],1]])
        DEFans = np.array([-(A[0])**2-(A[1])**2, 
                           -(p3[0])**2-(p3[1])**2,
                            -(B[0])**2-(B[1])**2])
        DEF = np.linalg.solve(DEFcoef, DEFans)
    
        return DEF   
        
    def cirintersect(self, detector, DEF1, DEF2):
        ''' Finds the intersection point of the 2 circles and tests which one is the point where all 3 interesct.
        function takes in the circle parameters for 2 circles. '''
        #solves for a,b,c in the quadratic formula
        D1 = DEF1[0]
        E1 = DEF1[1]
        F1 = DEF1[2]
        D2 = DEF2[0]
        E2 = DEF2[1]
        F2 = DEF2[2]
        
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
            
            #calculates the radius of the three detected points
            radius_srd = -detector[2] + (detector[0]/2)**2 + (detector[1]/2)**2
            radius_srd -= .01* radius_srd
             
            radius_test1 = (QuadansX1 + ((detector[0]) / 2))**2 + (QuadansY1 + (detector[1] / 2))**2
            radius_test2 = (QuadansX2 + ((detector[0]) / 2))**2 + (QuadansY2 + (detector[1] / 2))**2       

            #takes the point that doesn't lay on the edge of the detection ring
            if radius_test1 < radius_srd:
                ans = (QuadansX1,QuadansY1,0)
                return ans
            elif radius_test2 < radius_srd:
                ans = (QuadansX2,QuadansY2,0)
                return ans
            else:
                self.good_reconstruction = False
                return np.array([50,50,50])
                            
        else:
            self.good_reconstruction = False
            return np.array([50,50,50])

    def changebasisB_E(self, ans1, ans2, ans3):
        ''' takes the three answers given in the basis B and converts it back to the standard basis. 
        Also averages the three values which should all be within floating point error anyways. '''
        #make sure answers are all the same
        if (np.allclose(ans1, ans2) and np.allclose(ans2, ans3)):        
            #averages the 3 answers
            ans_B = (np.asarray(ans1) + np.asarray(ans2) + np.asarray(ans3))/3
            ans_E = np.dot(self.basis_B, ans_B) + self.det_pos1
            self.PO_decay = ans_E
        
        #if the three answers aren't the same, just make it a failed point
        else:
            self.good_reconstruction = False
            self.PO_decay = np.NAN
        
    def solvepoint(self):
        ''' do the work to solve the point '''
        self.angles()
        self.find_norm()  
        
        thrid_point1 = self.find_third_point(self.det_pos1, self.det_pos2, (self.alpha + self.gamma))
        thrid_point2 = self.find_third_point(self.det_pos2, self.det_pos3, (self.alpha + self.beta))
        thrid_point3 = self.find_third_point(self.det_pos3, self.det_pos1, (self.beta + self.gamma))
        
        self.changebasisE_B()
        thrid_point1_B = thrid_point1 - self.det_pos1
        thrid_point1_B = np.dot(np.linalg.inv(self.basis_B), thrid_point1_B)
        thrid_point2_B = thrid_point2 - self.det_pos1
        thrid_point2_B = np.dot(np.linalg.inv(self.basis_B), thrid_point2_B)
        thrid_point3_B = thrid_point3 - self.det_pos1
        thrid_point3_B = np.dot(np.linalg.inv(self.basis_B), thrid_point3_B)
        
        DEF0 = self.cirlce(self.det_pos1_B, self.det_pos2_B, self.det_pos3_B)
        DEF1 = self.cirlce(self.det_pos1_B, thrid_point1_B, self.det_pos2_B)
        DEF2 = self.cirlce(self.det_pos2_B, thrid_point2_B, self.det_pos3_B)
        DEF3 = self.cirlce(self.det_pos1_B, thrid_point3_B, self.det_pos3_B)

        ans1 = self.cirintersect(DEF0, DEF1, DEF2)
        ans2 = self.cirintersect(DEF0, DEF1, DEF3)
        ans3 = self.cirintersect(DEF0, DEF2, DEF3)
        
        if self.good_reconstruction:            
            self.changebasisB_E(ans1, ans2, ans3)
    
        
class ThreeGammaImaging:
    def __init__(self):
        self.all_three_gamma_points = []
        self.succeeded = []
        self.failed = []
    
    def import_raw_data(self, filename):
        ''' When called, opens up the file specified by filename and creates the ThreeGammaEvent objests for each line
        in the file '''
        #opens specified file
        with open(filename, 'r') as data_file:
            #runs through each line in the file
            for data_line in data_file:
                #puts in line of data for all ThreeGammaEvents
                self.all_three_gamma_points.append(ThreeGammaEvent(data_line))
    
    def import_previous_run(self, filename):
            ''' When called, opens up the file specified by filename which has only the solutions to the three gamma
            decays '''
            #opens specified file
            with open(filename, 'r') as data_file:
                #runs through each line in the file
                for data_line in data_file:
                    #puts in line of data for all ThreeGammaEvents
                    self.succeeded.append(ThreeGammaEvent(data_line))     
        
    def solve(self):
        ''' steps through each ThreeGammaEvent object and runs the solvepoint() method so each object is correctly
        solved. '''
        #step through each object in list
        for ans in self.all_three_gamma_points:
            #solves for the point of decay
            ans.solvepoint()
            #decides whether to keep it or not
            if ans.good_reconstruction:
                self.succeeded.append(ans)
            else:
                self.failed.append(ans)
        self.succeeded = np.asarray(self.succeeded)
        self.failed = np.asarray(self.failed)
            
    def write_to_file(self, filename):
        ''' opens up a file and and writes all the final points of decay data as x y z'''
        with open(filename,'w') as wroteto:
            for ans in self.succeeded:
                writie = '{} {} {}\n'.format(ans.PO_decay[0], ans.PO_decay[1], ans.PO_decay[2])
                wroteto.write(writie)
    
    def heatmapplot(self, lim, step, name, save=False, xstring="x (mm)", ystring="y (mm)"):    
        set_matplotlib_formats('pdf', 'svg') # this makes higher quality plots than the default (PNG)
        plt.rcParams["figure.figsize"] = [10,5]
        #seperate x and y from all points of decay
        x = []
        y = []
        for ans in self.succeeded:
            x.append(ans.PO_decay[0])
            y.append(ans.PO_decay[1])
        #use lim and step to set bins
        xbin = np.arange(-lim, lim+step, step)
        ybin = np.arange(-lim, lim+step, step)
        
        #plotting commands
        plt.hist2d(x, y, bins=[xbin,ybin],cmap='Greys')
        plt.colorbar() #ticks=[0,1,2,3])
        plt.grid(color='k', alpha=0.15, linestyle='--')
        plt.axes().set_aspect('equal')
        plt.xlabel(xstring)
        plt.ylabel(ystring)
        plt.title(name)
        
        if save:
            plt.savefig("{}.png".format(name), dpi=1200, format='png' )
        plt.show()



if __name__ == '__main__':
    Reconstructed = ThreeGammaImaging()
    Reconstructed.import_data('C:/data/3gamma data sets/3Phot_001.txt')
    Reconstructed.solve()
    Reconstructed.heatmapplot(20,1, "testfile", False)

