# -*- coding: utf-8 -*-
"""
Created on Wed Jan 14 00:39:13 2015

Plot 2D RDF

@author: Ben
"""

import numpy as np
import skimage.io
import skimage.feature
import skimage.morphology
import os
import matplotlib.pyplot as plt
import time
from sys import stdout

from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection
from matplotlib import rcParams

from paircorrelation import PairCorrelationFunction_2D



##################### BACKGROUND REMOVAL FUNCTIONS #############################


def make_HanningWindow(image, pos_x, pos_y, r):
    """
    Creates a mask by making a 2D Hanning window of radius r, about position (pos_x, pos_y)

    Note that this is quite different from the Hanning window function used to remove low frequency noise!
    This one is much simpler (and maybe doesn't really need to be a function...)
    """
    r = int(r)
    mask = np.zeros(np.shape(image))
    x_max, y_max = np.shape(image)

    hanning_1D = np.hanning(2*r+1)
    hanningWindow = np.outer(hanning_1D,hanning_1D)

    for i in range(0,2*r+1):
        for j in range(0,2*r+1):
            x = pos_y - r + i    # I did a weird thing and flipped these somehow
            y = pos_x - r + j    # but it works fine like this, so....yeah.
            if (x >= 0 and x < x_max and y >= 0 and y < y_max):
                mask[x][y]= hanningWindow[i][j]

    return mask


def get_spacing_blobs(image, Hann=True, plot=False):
    """
    Finds the lowest periodicity and SL orientation of an image (i.e. particle spacing in a lattice)
    by taking an FFT and finding lowest spots with a blob detecting algorithm.

    Accepts:
    image: ndarray containing image
    Hann: boolean value - if True, applies a Hanning window before taking an fft
    plot: booliean value - if true, plots output, displaying the identified lowest peak

    Returns:
    spacing: float, spacing between adjecent particles, in real space, in pixels
    nearest_peak: ndarray with (x,y) coords of peak
    sigma: std. dev. of blob which found nearest peak
    """
    if Hann:
        fft=np.fft.fftshift(np.fft.fft2(image * np.outer(np.hanning(np.shape(image)[0]),np.hanning(np.shape(image)[0]))))
    else:
        fft=np.fft.fftshift(np.fft.fft2(image))     # Takes the FFT, shifts central spot to center
    im_fft=np.log(np.abs(fft))                      # For displaying the FFT, need log(abs(fft))

    thresh=2
    sigRat = 1.6
    print "Finding blobs with Difference of Gaussian method.  Parameters are:"
    print "Threshold = {}".format(repr(thresh))
    print "Sigma ratio = {}".format(repr(sigRat))
    stdout.flush()
    blobs = skimage.feature.blob_dog(im_fft,min_sigma=1.0,max_sigma=20.0,sigma_ratio=sigRat,threshold=thresh,overlap=0.4)
    # blobs = (n,3) ndarray
    # each index is (y,x,sigma) for a blob
    # play with threshold, sigmas/ratios a bit...

    x_max,y_max=np.shape(fft)
    nearest_peak_dist=float(max((x_max,y_max)))
    center = [x_max/2,y_max/2]

    for blob in blobs:
        dist = np.sqrt( (blob[1]-center[0])**2 + (blob[0]-center[1])**2  )
        if ( dist > 2 and dist < nearest_peak_dist ):
            nearest_peak_dist = dist
            nearest_peak=np.array([blob[1],blob[0]])
            sigma = blob[2]
    spacing = min(np.shape(image))/nearest_peak_dist

    if plot:
        fig, ax = plt.subplots()
        skimage.io.imshow(im_fft)     # Note that the filtered center will appear bright b/c we take a log(0)
        plt.plot(blobs[:, 1], blobs[:, 0], 'b.')
        plt.plot(nearest_peak[0],nearest_peak[1],'ro')
        plt.axis([int(center[0]-3*nearest_peak_dist),int(center[0]+3*nearest_peak_dist),int(center[1]-3*nearest_peak_dist),int(center[1]+3*nearest_peak_dist)])
        plt.show()


    return spacing, nearest_peak, sigma


def get_spacing_peaks(image, numPeaks=10, minDist=3, Hann=True, plot=False):
    """
    Finds the lowest periodicity and SL orientation of an image (i.e. particle spacing in a lattice)
    by taking an FFT and finding lowest spots by finding local maxima.

    Accepts:
    image: ndarray containing image
    numPeaks: int, the number of peaks to find
    minDist: int, minimum distance (in pixels) between adjacent peaks
    Hann: boolean value - if True, applies a Hanning window before taking an fft
    plot: booliean value - if true, plots output, displaying the identified lowest peak

    Returns:
    spacing: float, spacing between adjecent particles, in real space, in pixels
    nearest_peak: ndarray with (x,y) coords of peak
    """
    if Hann:
        fft=np.fft.fftshift(np.fft.fft2(image * np.outer(np.hanning(np.shape(image)[0]),np.hanning(np.shape(image)[0]))))
    else:
        fft=np.fft.fftshift(np.fft.fft2(image))     # Takes the FFT, shifts central spot to center
    im_fft=np.log(np.abs(fft))                      # For displaying the FFT, need log(abs(fft))

    # Find FFT spots - niave approach, via maxima
    fft_maxima=skimage.feature.peak_local_max(im_fft, min_distance=minDist, num_peaks=numPeaks)

    x_max,y_max=np.shape(fft)
    nearest_peak_dist=float(max((x_max,y_max)))
    center = [x_max/2,y_max/2]

    for maximum in fft_maxima:
        dist = np.sqrt( (maximum[0]-center[0])**2 + (maximum[1]-center[1])**2  )
        if ( dist > 2 and dist < nearest_peak_dist ):
            nearest_peak_dist = dist
            nearest_peak=np.array([maximum[1],maximum[0]])
    spacing = min(np.shape(image))/nearest_peak_dist

    nearest_peak[0] = nearest_peak[0] - center[0]
    nearest_peak[1] = nearest_peak[1] - center[1]

    if plot:
        fig, ax = plt.subplots()
        skimage.io.imshow(im_fft)     # Note that the filtered center will appear bright b/c we take a log(0)
        plt.plot(fft_maxima[:, 1], fft_maxima[:, 0], 'b.')
        plt.plot(nearest_peak[0],nearest_peak[1],'ro')
        plt.axis([int(center[0]-3*nearest_peak_dist),int(center[0]+3*nearest_peak_dist),int(center[1]-3*nearest_peak_dist),int(center[1]+3*nearest_peak_dist)])
        plt.show()

    return spacing, nearest_peak


def find_blobs(image, spacing, minSigma=1, maxSigma=10.0, sigRat=1.6, thresh=2, method="dog", plot=False):
    """
    Finds blobs using the selected algorithm.  Basically just a wrapper around the skimage functions, with added plotting utility.
    method variable must be "dog", "doh", or "log" (Diff. of Gaussians, Deter. of Hessian, and Log of Gaussian, respectively)

    See http://scikit-image.org/docs/dev/auto_examples/plot_blob.html, and http://scikit-image.org/docs/dev/api/skimage.feature.html#skimage.feature.blob_dog
    """
    sigma_0 = spacing/np.sqrt(8)

    if method=="dog":
        # Difference of Gaussian (DoG) blob detector
        print "Finding blobs with Difference of Gaussian method.  Parameters are:"
        print "Threshold = {}".format(repr(thresh))
        print "Sigma ratio = {}".format(repr(sigRat))
        stdout.flush()
        blobs = skimage.feature.blob_dog(image,min_sigma=minSigma,max_sigma=maxSigma,sigma_ratio=sigRat,threshold=thresh,overlap=0.4)
        ###### Blob detection algorithm may require dtype=float images to already be scaled to between -1 and 1?
        ###### Look into this........
    elif method=="doh":
        # Determinant of Hessian (DoH) blob detector
        thresh=0.002
        numSigma = 3
        print "Finding blobs with Determinant of Hessian method.  Parameters are:"
        print "Threshold = {}".format(repr(thresh))
        print "Number of sigma = {}".format(repr(numSigma))
        stdout.flush()
        blobs = skimage.feature.blob_doh(image,min_sigma=sigma_0*0.25,max_sigma=sigma_0*1.75,num_sigma=numSigma,threshold=thresh,overlap=0.4, log_scale=False)
    elif method=="log":
        # Laplacian of Gaussian (LoG) blob detector
        blobs = skimage.feature.blob_log(image,min_sigma=sigma_0*0.25,max_sigma=sigma_0*1.75,num_sigma=3,threshold=0.1,overlap=0.4)
    else:
        print "method = {}.  Must be set to 'dog', 'doh', or 'log'.".format(repr(method))

    if plot:
        fig, ax = plt.subplots()
        skimage.io.imshow(image)
        if len(blobs) != 0:
            plt.plot(blobs[:, 1], blobs[:, 0], 'r.')
        plt.show()
        fig.canvas.draw()

    return blobs

def get_NC_orientation(image, spacing, pos_x, pos_y, NCnum=0, plot=True):
    """
    Finds the orientation of a single NC by identifying the FFT spots from its atomic lattice.
    Returns the vector to the nearest FFT spot

    -crops around the NC of interest
    -applies a Hanning windows around the NC of interest
    -takes fft
    -finds local maxima
    -chooses local max which is the appropriate distance to be the atomic lattice
    """
    imageYsize,imageXsize = np.shape(image)

    # Find cropping window
    r = spacing/1.6
    minSize = 2*r
    powsOfTwo = np.array([2**i for i in range(1,11)])
    diff = powsOfTwo - minSize
    windowSize = powsOfTwo[np.array([1024 if i<=0 else i for i in diff]).argmin()]
    window = [pos_x - windowSize/2, pos_x + windowSize/2, pos_y - windowSize/2, pos_y + windowSize/2]
    if window[0] < 0:
        window[0],window[1] = 0, windowSize
    if window[1] >= imageXsize:
        window[0],window[1] = imageXsize-windowSize, imageXsize
    if window[2] < 0:
        window[2],window[3] = 0, windowSize
    if window[3] >= imageYsize:
        window[2],window[3] = imageYsize-windowSize, imageYsize

    croppedImage = image[window[2]:window[3], window[0]:window[1]]
    shiftedX, shiftedY = pos_x - window[0], pos_y - window[2]

    fft=np.fft.fftshift(np.fft.fft2(croppedImage*make_HanningWindow(croppedImage,shiftedX,shiftedY,r)))
    im_fft=np.log(np.abs(fft))                      # For displaying the FFT, need log(abs(fft))

    # Find FFT spots - niave approach, via maxima
    fft_maxima=skimage.feature.peak_local_max(im_fft, min_distance=5, num_peaks=7)

    x_max,y_max=np.shape(fft)
    nearest_peak_dist=float(max((x_max,y_max)))
    nearest_peak=np.array([0,0])
    center = [x_max/2,y_max/2]

    for maximum in fft_maxima:
        dist = np.sqrt( (maximum[0]-center[0])**2 + (maximum[1]-center[1])**2  )
#        if ( dist > 2 and dist < nearest_peak_dist ):
        if (dist > 38 and dist < 42):                       # Checked by hand - could also calculate this but it'd be annoying...
            nearest_peak_dist = dist
            nearest_peak=np.array([maximum[1],maximum[0]])

    nearest_peak[0] = nearest_peak[0] - center[0]
    nearest_peak[1] = nearest_peak[1] - center[1]

    if plot:
        fig1 = plt.figure(2)
        fig1.clf()
        fig1.suptitle("Orientation of NC {}".format(repr(NCnum)), fontsize=16)
        ax1 = fig1.add_subplot(211)
        ax1.imshow(image*make_HanningWindow(image,pos_x,pos_y,r), cmap='gray') 
        ax1.axis([int(pos_x - spacing/1.5), int(pos_x + spacing/1.5),int(pos_y - spacing/1.5), int(pos_y + spacing/1.5)])

        ax2 = fig1.add_subplot(212)
        ax2.imshow(im_fft, cmap='gray')     # Note that the filtered center will appear bright b/c we take a log(0)
        ax2.plot(fft_maxima[:, 1], fft_maxima[:, 0], 'b.')
        ax2.plot(nearest_peak[0]+center[0],nearest_peak[1]+center[1],'ro')
        ax2.axis([int(center[0]-3*nearest_peak_dist),int(center[0]+3*nearest_peak_dist),int(center[1]-3*nearest_peak_dist),int(center[1]+3*nearest_peak_dist)])
        plt.show()

        fig1.canvas.draw()

    return nearest_peak


def get_angle(AL_peak):
    """
    Finds the angle between a vector mod 90deg and the x-axis, formatted between -45 and 45 deg
    """

    A_vect = np.array(AL_peak)
    S_vect = np.array([1,0])

    if (A_vect[0] <= 0 and A_vect[1] > 0):
        A_vect[0], A_vect[1] = A_vect[1], -A_vect[0]
    elif (A_vect[0] < 0 and A_vect[1] <= 0):
        A_vect[0], A_vect[1] = -A_vect[0], -A_vect[1]
    elif (A_vect[0] >= 0 and A_vect[1] < 0):
        A_vect[0], A_vect[1] = -A_vect[1], A_vect[0]

    angle = np.degrees(np.arccos(np.dot(A_vect/np.linalg.norm(A_vect),S_vect/np.linalg.norm(S_vect))))
    if angle > 45:
        angle = angle - 90

    return angle

def plot_angles(image, spacing, blobs_DoG, showFFTs=True):
    NC_angles = []
    r = spacing/2.0
    patches = []
    angles = []
    FFTspotDistances = []
    for i in range(len(blobs_DoG)):
#    for i in range(50,60):      # For testing!
        x = blobs_DoG[i,1]
        y = blobs_DoG[i,0]
        print "Finding orientation of NC #{}".format(repr(i))
        stdout.flush()
        AL_peak = get_NC_orientation(image, spacing, x, y, NCnum=i, plot=showFFTs)
        angle = get_angle(AL_peak)
        print "Orientation of particle {} is {} degrees to x-axis".format(repr(i),repr(angle))
        stdout.flush()

        NC_angles.append([x,y,angle])
        circle = Circle((x,y),r)
        patches.append(circle)
        angles.append(angle)
        FFTspotDistances.append(np.linalg.norm(np.array(AL_peak)))

    p=PatchCollection(patches, alpha=0.5)
    p.set_array(np.array(angles))

    fig = plt.figure(3)
    ax = fig.add_subplot(111)
    ax.imshow(image, cmap="gray")
    ax.add_collection(p)
    plt.colorbar(p)

    plt.savefig("AL_v_SL_orientationMap.png")
    p.set_clim([-20, 5])
    plt.savefig("AL_v_SL_orientationMap2.png")
    plt.show()
    fig.canvas.draw()

    output = np.array(zip(blobs_DoG[:,1],blobs_DoG[:,0],angles, FFTspotDistances))
    np.save("AL_v_SL_centersAndAngles.npy",output) # Save output array as an npy file, if we wanna mess with fig params later...

    return output

def plotRDF(g,r):
    rcParams['figure.facecolor'] = 'w'  # Gets rid of gray background

    # No vertical lines

    fig1 = plt.figure()
    ax1 = fig1.add_subplot(1,1,1)

    # Plot g(r)
    plt.axhspan(0.9, 1.1, xmin=0, xmax=1, facecolor='y', alpha=0.5)
    #plt.axvspan(0, 11, ymin=0, ymax=1, facecolor='g',alpha=0.3)         # Color inset
    plt.axhline(y=1, color='k',linestyle='-')
    ax1.plot(r,g, color='r', linewidth=2)
    plt.axis([0,270,0,8])
    ax1.set_title('Radial Distribution Function', size=20)
    ax1.set_ylabel(r'$g(r)$', size=16)
    ax1.set_xlabel(r'$r$ (nm)', size=16)

    # Horizontal line at y=1


    # Add inset
    ax_inset = plt.axes([0.48,0.4,0.4,0.4], axisbg='w')
    #plt.axvspan(0, 11, ymin=0, ymax=1, facecolor='g',alpha=0.3)        # Color inset
    ax_inset.plot(r,g,'r-', linewidth=3)
    ax_inset.axis([5,11,0,4])
    #ax_inset.set_xlabel('r (nm)', size=12)
    ax_inset.set_xticks([6,8,10])
    ax_inset.set_yticks([0,4,8])

    return

############## TESTING ###########
"""
AL_peak = get_NC_orientation(image_all, 133.9, 2858, 3826, plot=True)

center_x,center_y = np.shape(image_all)[0]/2, np.shape(image_all)[1]/2,
angle = get_angle(AL_peak,SL_peak, center_x, center_y)

print "Relative orientation = {}".format(repr(angle))
#get_NC_orientation(image_all, spacing, blob[1], blob[0],spacing/1.6, plot=True)

"""

#patches = []

#centers = [[blobs_DoG[i][0],blobs_DoG[i][1]] for i in range(len(blobs_DoG))]
#vor = Voronoi(points=centers)


############## RUNNING ###########

if __name__ == "__main__":

    # read an image
    # keep a couple of sample images in this directory
    # one with and one without low freq. background?

    filename = "06c_115kx_100mm_16us_2048.png"
    image = skimage.io.imread(filename)

    t = time.localtime()
    time_stamp = "date={}_{}_{}_time={}h_{}min_{}s".format(t.tm_year,t.tm_mon,t.tm_mday,t.tm_hour,t.tm_min,t.tm_sec)

    if not os.path.exists(os.path.splitext(filename)[0]+"_outputs"):
        os.mkdir(os.path.splitext(filename)[0]+"_outputs")
    output_name=os.path.splitext(filename)[0]+"_outputs/"

    # Find SL location/spacing
    print "Getting SL spacing and orientation..."
    spacing, SL_peak = get_spacing_peaks(image, numPeaks=10, minDist=30, Hann=True, plot=True)
    print "Done.\nSL spacing = {} pixels".format(repr(spacing))
    stdout.flush()

    # Find positions of NCs.
    blob_method = "dog"
    print "Finding blobs via {} method...".format(repr(blob_method))
    stdout.flush()
    time_init = time.time()

    # Set blob detector parameters
    sigRat=1.6
    minSigma = spacing/(2*np.sqrt(2)*(sigRat**2))
    maxSigma = minSigma*(sigRat**2)
    thresh = .000001

    # Run blob detection.  Display results for checking...
    blobs_DoG = find_blobs(image,spacing, minSigma=minSigma, maxSigma=maxSigma, sigRat=sigRat, thresh=thresh, method=blob_method,plot=True)
    time_tot = time.time() - time_init
    print "Done.  Found {} blobs.  Process took {} seconds.".format(repr(len(blobs_DoG)), repr(int(time_tot)))
    if len(blobs_DoG) == 0:
        raise RuntimeError("Failed to find any blobs! Try altering threshold value.")
    np.savez(output_name+"blobs_DoG_"+time_stamp, blobs_DoG=blobs_DoG, spacing=spacing)

    maxDist = 20    # Distance to which we'll calculate the RDF, in particle diameters
    pixels_per_nm = 4096.0/600.0   # We'll need this if we want the RDF in nm...
    g,r,x,y = PairCorrelationFunction_2D(blobs_DoG[:,1],blobs_DoG[:,0],float(np.shape(image)[0]),maxDist*spacing,1)

    r_nm = r/pixels_per_nm
    r_particles = r/spacing

    plotRDF(g,r_nm)

    execfile("applause.py")

