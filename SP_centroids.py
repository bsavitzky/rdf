# -*- coding: utf-8 -*-
"""
Loads an image and a .npz file, and displays an overlay.
Used to determine the success of centroids.py or centroids_SP.py

inputs: image file, as a tiff (other formats compatible with skimage.io.load should work, i.e. png, etc)
        centroids file, as .npz
output: centroids file, as .npz

@author: BHS
"""

# Import Libraries
import numpy as np
import matplotlib.pyplot as plt
import os
import time
import argparse
import skimage.io
import scipy.optimize as opt

from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection

# Import local libraries
import mdscrape as md

# Parse command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("image_file")
parser.add_argument("centroid_file")
args = parser.parse_args()

# Handle IO
if not os.path.exists("outputs"):
    os.mkdir("outputs")
output_name="outputs/"+"SP_centroids"

# Get image, centroids, pixels
print "Loading image and centers and generating overlay."
image = skimage.io.imread(args.image_file)
centroids = np.load(args.centroid_file)
x = centroids['x']
y = centroids['y']
spacing = float(centroids['spacing'])

# Get metadata
metadata = md.extract_metadata(args.image_file)
fov_pixels = metadata['pixels']

# Check if a particle is on the edge of the image (within 1/2 SL spacing)
def on_edge(x0,y0,rad,image):
    xpixels = np.shape(image)[0]
    ypixels = np.shape(image)[1]
    xmin, xmax = x0-rad, x0+rad
    ymin, ymax = y0-rad, y0+rad
    if (xmin <= 0) or (xmax >= xpixels) or (ymin <= 0) or (ymax >= ypixels):
        return True
    else:
        return False

# Create small array about a given center (x0,y0) of width "size"
def filter_dot(x0,y0,size,image):
    """
    Inputs: x0,y0 - center of particle to filter about, in pixels
            size - size of box filter (i.e. SL spacing)
            image - original image, as ndarray
    Output: smallImage - small ndarray about area of interest
            xshift, yshift - shifts relative to original image
    """
    x0, y0 = int(x0), int(y0)
    rad = int(np.ceil(size/2.0))
    xpixels = np.shape(image)[0]
    ypixels = np.shape(image)[1]
    # Determine window max and min and make window array
    xmin, xmax = x0-rad, x0+rad
    ymin, ymax = y0-rad, y0+rad
    # Throw away data points on edge of image
    if (xmin < 0) or (xmax > xpixels) or (ymin < 0) or (ymax > ypixels):
        return None
    smallImage = image[ymin:ymax, xmin:xmax]
    return smallImage, xmin, ymin


# Convenience function to display a small window and particle center
def display(smallImage, xcenter, ycenter):
    fig, ax = plt.subplots()
    plt.imshow(smallImage, cmap="gray")
    plt.plot(xcenter,ycenter,'r.')
    plt.show()
    fig.canvas.draw()


# 2D Gaussian
# Returns result as a 1D array that can be passed to scipy.optimize.curve_fit
def gauss2d((x,y), amplitude, x0, y0, sigma_x, sigma_y, theta, offset):
    x0, y0 = float(x0), float(y0)
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + np.abs(amplitude)*np.exp( - (a*((x-x0)**2) + 2*b*(x-x0)*(y-y0) + c*((y-y0)**2) ))
    return g.ravel()


# Get subpixel center for a single spot with a 2D gaussian fit
# Edge particles should be filtered out before passing to xySP
def xySP(x0, y0, spacing, image, plot=False):
    # Get smaller image for faster processing
    smallImage, xshift, yshift = filter_dot(x0,y0,spacing,image)
    # Get new centers
    x0_smallIm, y0_smallIm = x0-xshift, y0-yshift
    # Define mesh for input values and initial guess
    xs,ys = np.meshgrid(range(np.shape(smallImage)[1]),range(np.shape(smallImage)[0]))
    initial_guess = (smallImage[y0_smallIm,x0_smallIm], y0_smallIm, x0_smallIm, spacing/4.0, spacing/4.0,0,0)
    # Set all values outside circle of radius spacing/2 to min value in a small image to account for adjacent particles
    baseline = smallImage.min()
    for i in range(np.shape(smallImage)[0]):
        for j in range(np.shape(smallImage)[1]):
            if (x0_smallIm - i)**2 + (y0_smallIm - j)**2 > (spacing/2.0)**2:
                smallImage[i,j] = baseline
    # Perform fit and pull out centers
    try:
        popt, pcov = opt.curve_fit(gauss2d, (xs,ys), smallImage.ravel(), p0=initial_guess)
    except RuntimeError:
        print "Particle could not be fit to a 2D gaussian.  Returning original centroid."
        return x0, y0, 1
    x_SP, y_SP = popt[2]+xshift, popt[1]+yshift
    # Plotting for troubleshooting
    if plot:
        data_fitted = gauss2d((xs,ys), *popt)
        fig,ax=plt.subplots(1,1)
        ax.imshow(smallImage)
        ax.contour(xs,ys,data_fitted.reshape(np.shape(smallImage)[0],np.shape(smallImage)[1]),8,colors='w')
        plt.show()
    return x_SP, y_SP, 0


# Iterate over all particles
print "Performing 2D gaussian fit to all peaks..."
time_init = time.time()
x_SP=[]
y_SP=[]
shift = spacing/2.0
unfit_particles = 0
edge_particles = 0
for i in range(len(x)):
    print "Fitting particle {} of {}".format(i, len(x))
    xcurr,ycurr = x[i],y[i]
    if not on_edge(xcurr,ycurr,shift,image):
        xSPcurr,ySPcurr, fit = xySP(xcurr,ycurr,spacing,image,plot=False)
        # Fitting may have shifted particle into edge region, so check again...
        if not on_edge(xSPcurr,ySPcurr,shift,image):
            x_SP.append(xSPcurr)
            y_SP.append(ySPcurr)
            unfit_particles += fit
        else:
            edge_particles += 1
    else:
        edge_particles += 1
time_tot = time.time() - time_init
# Shift centers
fov_pixels = fov_pixels - 2*shift
x_SP, y_SP = np.array(x_SP) - shift, np.array(y_SP) - shift
print "Done.\n{} total particles.".format(len(x))
print "{} edge particles were discarded.".format(edge_particles)
print "{} particles could not be fit and original pixel-accuracy centroid was used.".format(unfit_particles)
print "Process took {} seconds.".format(repr(int(time_tot)))
print "Shifting centers and FOV, and saving in {}.npz".format(output_name)
np.savez(output_name, x=x_SP, y=y_SP, spacing=spacing, fov_pixels=fov_pixels)









