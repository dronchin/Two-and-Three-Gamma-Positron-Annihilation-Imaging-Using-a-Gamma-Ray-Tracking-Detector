# -*- coding: utf-8 -*-
"""
Created on Tue May  9 13:48:58 2017

@author: Nicolas D

The purpose of this program is to take in the detection data and preform a 
radon transformation on it to construct a sinogram. After the radon 
transformation, an inverse radon transformation is applied to reconstruct what 
the picture is using back projection. The sinogram and reconstructed image are 
then displayed.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft, ifft, fftfreq
from warnings import warn
from skimage.transform import warp
from scipy.interpolate import interp1d
from IPython import display
import winsound
import timeit

starttime = timeit.default_timer()


def unpackdata(initial_points_data):
    ''' This function needs to take a line of data and change it into usable 
    varibales
    '''
    detect1 = initial_points_data[0:3]
    detect2 = initial_points_data[3:]
    
#    #take x-z projection now
#    z1 = detect1[2]
#    z2 = detect2[2]
#    detect1[1] = z1
#    detect2[1] = z2
#    detect1[2] = 0
#    detect2[2] = 0
    return detect1, detect2

def adderror(det1, det2, sigma):
    det1 += np.random.normal(0,scale=sigma, size=3)
    det2 += np.random.normal(0,scale=sigma, size=3)
    if np.random.random() < 0.05:
        det1 += np.random.normal(0,scale=50, size=3)
        det2 += np.random.normal(0,scale=50, size=3)
    return det1, det2

def unpackdata_RD(initial_points_data):
    ''' This function needs to take a line of data and change it into usable 
    varibales
    '''
    detect1 = initial_points_data[0:3]
    detect2 = initial_points_data[4:7]
    
#    #rotate detections by pi/4
#    x = detect1[0]
#    y = detect1[1]
#    detect1[0] = x*np.cos(np.pi/4) - y*np.sin(np.pi/4)
#    detect1[1] = x*np.sin(np.pi/4) + y*np.cos(np.pi/4)
##    
#    x = detect2[0]
#    y = detect2[1]
#    detect2[0] = x*np.cos(np.pi/4) - y*np.sin(np.pi/4)
#    detect2[1] = x*np.sin(np.pi/4) + y*np.cos(np.pi/4)
    return detect1, detect2

def unpackdata_noY(initial_points_data):
    ''' This function needs to take a line of data and change it into usable 
    varibales
    '''
    detect1 = initial_points_data[0:3]
    detect2 = initial_points_data[4:7]
    x1 = detect1[0]
    x2 = detect2[0]
    z1 = detect1[2]
    z2 = detect2[2]
    detect1[1] = x1
    detect2[1] = x2
    detect1[0] = z1 
    detect2[0] = z2
    
    #print(detect1, detect2)
    return detect1, detect2, z1, z2

def intersect(phi, theta, b1):
    '''This function's goal is to solve the intersection of the axis line and 
    the LOR for a sinlge angle projection.
    Inputs:
        phi = angle for axis
        theta = slope anlge for LOR
        b1 = y intercept of the slope
    Output:
        s = distance on the axis. + is right, - is left
    '''
    m2 = np.tan(phi)
    m1 = np.tan(theta)
    
    x = b1/(m2-m1)
    y = m2*x
    
    #print("x", x, "y",y)
    
    if y >= 0:    
        s = np.sqrt(x**2 + y**2)
    else:
        s = -np.sqrt(x**2 + y**2)
    
    return s
    
def _sinogram_circle_to_square(sinogram):
    '''This function is used internally inside of the iradon function.
    Its job is to take data for a sinogram given in a circle and pad the edges 
    so that it will make a square final picture later on.
    '''
    diagonal = int(np.ceil(np.sqrt(2) * sinogram.shape[0]))
    pad = diagonal - sinogram.shape[0]
    old_center = sinogram.shape[0] // 2
    new_center = diagonal // 2
    pad_before = new_center - old_center
    pad_width = ((pad_before, pad - pad_before), (0, 0))
    return np.pad(sinogram, pad_width, mode='constant', constant_values=0.0)


def iradon(radon_image, theta=None, output_size=None,
           filter="ramp", interpolation="linear", circle=None):
    """
    Inverse radon transform.
    Reconstruct an image from the radon transform, using the filtered
    back projection algorithm.
    Parameters
    
    adapted from scikit-image>transform>radon_transform sourcecode found at:
    https://github.com/scikit-image/scikit-image/blob/master/skimage/transform/radon_transform.py#L127
    ----------
    radon_image : array_like, dtype=float
        Image containing radon transform (sinogram). Each column of
        the image corresponds to a projection along a different angle. The
        tomography rotation axis should lie at the pixel index
        ``radon_image.shape[0] // 2`` along the 0th dimension of
        ``radon_image``.
    theta : array_like, dtype=float, optional
        Reconstruction angles (in degrees). Default: m angles evenly spaced
        between 0 and 180 (if the shape of `radon_image` is (N, M)).
    output_size : int
        Number of rows and columns in the reconstruction.
    filter : str, optional (default ramp)
        Filter used in frequency domain filtering. Ramp filter used by default.
        Filters available: ramp, shepp-logan, cosine, hamming, hann.
        Assign None to use no filter.
    interpolation : str, optional (default 'linear')
        Interpolation method used in reconstruction. Methods available:
        'linear', 'nearest', and 'cubic' ('cubic' is slow).
    circle : boolean, optional
        Assume the reconstructed image is zero outside the inscribed circle.
        Also changes the default output_size to match the behaviour of
        ``radon`` called with ``circle=True``.
        The default behavior (None) is equivalent to False.
    Returns
    -------
    reconstructed : ndarray
        Reconstructed image. The rotation axis will be located in the pixel
        with indices
        ``(reconstructed.shape[0] // 2, reconstructed.shape[1] // 2)``.
    References
    ----------
    .. [1] AC Kak, M Slaney, "Principles of Computerized Tomographic
           Imaging", IEEE Press 1988.
    .. [2] B.R. Ramesh, N. Srinivasa, K. Rajgopal, "An Algorithm for Computing
           the Discrete Radon Transform With Some Applications", Proceedings of
           the Fourth IEEE Region 10 International Conference, TENCON '89, 1989
    Notes
    -----
    It applies the Fourier slice theorem to reconstruct an image by
    multiplying the frequency domain of the filter with the FFT of the
    projection data. This algorithm is called filtered back projection.
    """
    if radon_image.ndim != 2:
        raise ValueError('The input image must be 2-D')
    if theta is None:
        m, n = radon_image.shape
        theta = np.linspace(0, 180, n, endpoint=False)
    else:
        theta = np.asarray(theta)
    if len(theta) != radon_image.shape[1]:
        raise ValueError("The given ``theta`` does not match the number of "
                         "projections in ``radon_image``.")
    interpolation_types = ('linear', 'nearest', 'cubic')
    if interpolation not in interpolation_types:
        raise ValueError("Unknown interpolation: %s" % interpolation)
    if not output_size:
        # If output size not specified, estimate from input radon image
        if circle:
            output_size = radon_image.shape[0]
        else:
            output_size = int(np.floor(np.sqrt((radon_image.shape[0]) ** 2
                                               / 2.0)))
    if circle is None:
        warn('The default of `circle` in `skimage.transform.iradon` '
             'will change to `True` in version 0.15.')
        circle = False
    if circle:
        radon_image = _sinogram_circle_to_square(radon_image)
        #print(radon_image.shape)
        #plt.imshow(radon_image, cmap="Greys")
        #plt.show()

    th = (np.pi / 180.0) * theta
    #print(th)
    # resize image to next power of two (but no less than 64) for
    # Fourier analysis; speeds up Fourier and lessens artifacts
    projection_size_padded = \
        max(64, int(2 ** np.ceil(np.log2(2 * radon_image.shape[0]))))

    pad_width = ((0, projection_size_padded - radon_image.shape[0]), (0, 0))
    img = np.pad(radon_image, pad_width, mode='constant', constant_values=0)

    # Construct the Fourier filter
    f = fftfreq(projection_size_padded).reshape(-1, 1)   # digital frequency
    omega = 2 * np.pi * f                                # angular frequency
    fourier_filter = 2 * np.abs(f)                       # ramp filter
    if filter == "ramp":
        pass
    elif filter == "shepp-logan":
        # Start from first element to avoid divide by zero
        fourier_filter[1:] = fourier_filter[1:] * np.sin(omega[1:]) / omega[1:]
    elif filter == "cosine":
        fourier_filter *= np.cos(omega)
    elif filter == "hamming":
        fourier_filter *= (0.54 + 0.46 * np.cos(omega / 2))
    elif filter == "hann":
        fourier_filter *= (1 + np.cos(omega / 2)) / 2
    elif filter is None:
        fourier_filter[:] = 1
    else:
        raise ValueError("Unknown filter: %s" % filter)
    # Apply filter in Fourier domain
    projection = fft(img, axis=0) * fourier_filter
    
    radon_filtered = np.real(ifft(projection, axis=0))

    # Resize filtered image back to original size
    radon_filtered = radon_filtered[:radon_image.shape[0], :]
    reconstructed = np.zeros((output_size, output_size))
    # Determine the center of the projections (= center of sinogram)
    mid_index = radon_image.shape[0] // 2

    [X, Y] = np.mgrid[0:output_size, 0:output_size]
    xpr = X - int(output_size) // 2
    ypr = Y - int(output_size) // 2

    # Reconstruct image by interpolation
    for i in range(len(theta)):
        t = ypr * np.cos(th[i]) - xpr * np.sin(th[i])
        x = np.arange(radon_filtered.shape[0]) - mid_index
        if interpolation == 'linear':
            backprojected = np.interp(t, x, radon_filtered[:, i],
                                      left=0, right=0)
            
        else:
            interpolant = interp1d(x, radon_filtered[:, i], kind=interpolation,
                                   bounds_error=False, fill_value=0)
            backprojected = interpolant(t)
        reconstructed += backprojected
        
    if circle:
        radius = output_size // 2
        reconstruction_circle = (xpr ** 2 + ypr ** 2) <= radius ** 2
        reconstructed[~reconstruction_circle] = 0.
    print("it's working")
    return reconstructed * np.pi / (2 * len(th))



slopelist = [] #stores the slope values
filelist = ['detector datasets/PET_03_r53.dat', 'detector datasets/PET_06_r52.dat',
            'detector datasets/PET_09_r54.dat', 'detector datasets/PET_12_r50.dat',
            'detector datasets/PET_center_r48.dat'] #'PET_06_r51.dat', 
filelist2 = ['detector datasets/PET_center_r48.dat']
filelist3 = ['PET_center_r48.dat']
            
filelist4 = ['Run0007.dat']            
filelist5 = ['testfile_ringdata.txt']    

filelist6 = ['PET_run05.txt']
filelist7 = ['PET_run05.txt','PET_run07.txt','PET_run09.txt','PET_run12.txt','PET_run13.txt']#,'PET_run15.txt','PET_run23.txt','PET_run27.txt']
filelist8 = ['PET_run13.txt']#,'PET_run15.txt','PET_run23.txt','PET_run27.txt']
filelist9 = ['PET_zeq15mm_run15.txt']

filelist10 = ['PET_run23.txt', 'PET_run25.txt', 'PET_run27.txt']

filelist11 = ['PET_zeq15mm_run05.txt','PET_zeq15mm_run07.txt','PET_zeq15mm_run09.txt','PET_zeq15mm_run12.txt','PET_zeq15mm_run13.txt']
filelist12 = ['PET_zeq15mm_run15.txt']
filelist13 = ['PET_zeq15mm_run23.txt','PET_zeq15mm_run25.txt','PET_zeq15mm_run27.txt']
filelist14 = ['PET_run25.txt']
filelist15 = ['C:/data/detector datasets/PET_powder.asc']

filelist16 = ["PET_03_r53.dat", "PET_06_r52.dat",
             "PET_09_r54.dat", "PET_12_r50.dat", "PET_center_r48.dat"]
filelist17 = ["PET_03_r53.dat"]
filelist18 = ["PET_center_r48.dat", "PET_dwnstr_e49.dat"]
filelist19 = ["PET_03_r53.dat", "PET_06_r52.dat", "PET_09_r54.dat", "PET_12_r50.dat", "PET_center_r48.dat", "PET_dwnstr_e49.dat"]

filelist20 = ["simulate_moved_Na.txt"]
filelist21 = ["simulate posdepend/position_dep_0.txt","simulate posdepend/position_dep_1.txt",
              "simulate posdepend/position_dep_2.txt","simulate posdepend/position_dep_3.txt",
              "simulate posdepend/position_dep_4.txt","simulate posdepend/position_dep_5.txt",
              "simulate posdepend/position_dep_6.txt","simulate posdepend/position_dep_7.txt",
              "simulate posdepend/position_dep_8.txt","simulate posdepend/position_dep_9.txt",
              "simulate posdepend/position_dep_10.txt","simulate posdepend/position_dep_11.txt",
              "simulate posdepend/position_dep_12.txt"]

filelist22 = ["simulate posdepend/position_dep_0.txt","simulate posdepend/position_dep_5.txt"]
filelist23 = ["simulate_moved_Na_downshifted_15.txt"]
filelist24 = ["simulate_point_source.txt"]
filelist25 = ['simulate_phantom_temp.txt']

zrange = []    
for filenum, filename in enumerate(filelist25):
    print(filenum)
    with open("C:/data/{}".format(filename, "r")) as data:
        for count,line in enumerate(data):
#            if count > 78733:
#                print("cut", count)
#                break
            #unpacks the data and stores them in a list
            numbers_float = map(float, line.split())
            initial_points_data = list(numbers_float)
            #print(initial_points_data)
            detect1, detect2 = unpackdata(initial_points_data) #saves data as detect1 and dectect2
            detect1, detect2 = adderror(detect1, detect2, 3)

            
            if detect2[0] - detect1[0] == 0:
                continue
                        
            m = (detect2[1] - detect1[1])/(detect2[0] - detect1[0]) #finds the y-intercept of the line
            b = detect1[1] - m*detect1[0] #calculates the slope in dy/dx form
            
            slope = np.arctan(m) #calculates the slope in angles form
            
            slopelist.append([slope,b])
        print(count)

      

phistep = np.pi/100
slist = [] #stores all calculated s vlaues in a 1d list
angles = [] #stores all matched angles to the slist number


slopelist = np.array(slopelist)
thetas = slopelist[:,0]
plt.hist(thetas, bins = np.arange(-np.pi/2,np.pi/2,0.1), normed = True)
plt.show()


#Need to serach and sort through the angles and group them by their angles
for phi in np.arange(0, np.pi + phistep, phistep):
    normphi = phi - np.pi/2 #need to be parallel to the normal angle
#    print(phi)
    for theta,b1 in slopelist:
        if (theta > normphi and theta < normphi+phistep):
            #uses interesect function to find the distance on the new rotated axis
            s = intersect(phi, theta, b1) 
            #print("s", s)
            angles.append(phi)
            slist.append(s)
            
#
#stoptime = timeit.default_timer()
#print("time taken: ", stoptime - starttime)  
#
#
#A histogram will group together all the s values and weight them accordingly to form a sinogram 
heatmap, xedges, yedges = np.histogram2d(slist, angles,bins=[np.arange(-60,61,1),np.linspace(0,np.pi,100)], normed=True)
##uses the iradon function to turn a sinogram into a reconstructed image using back projection
reconstruction = iradon(heatmap, circle=True, filter="ramp")

reconstruction -= reconstruction.min()
reconstruction = reconstruction/reconstruction.max()
#reconstruction = reconstruction.T
reconstruction = reconstruction[::-1, ...]


#np.save("reconstruction_upstream.npy", reconstruction)

#
##plots the sinogram and reconstructed image
#plt.figure(figsize=(8.5,11))
##plt.subplot(2,1,1)
#plt.imshow(heatmap, cmap="Greys")
##plt.subplot(2,1,2)
plt.style.use('default')

plt.imshow(reconstruction, cmap="Greys", origin='lower')
plt.colorbar()
plt.xlabel('X (mm)')
plt.ylabel('Y (mm)')
plt.grid(color='w', linestyle='--')
plt.title('Simulated Phantom')
#
#plt.imshow(heatmap) 
#plt.savefig("Sim_phantom_1.png", dpi=1200, format='png') 
plt.show()

winsound.Beep(800,1000)  
    
#        