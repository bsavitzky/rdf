# -*- coding: utf-8 -*-
"""
Loads an image and a .npz file, and displays an overlay.
Used to determine the success of centroids.py or centroids_SP.py

inputs: image file, as a tiff (other formats compatible with skimage.io.load should work, i.e. png, etc)
        centers file, as .npz
output: none

@author: BHS
"""

# Import Libraries
import numpy as np
import matplotlib.pyplot as plt
import os
import argparse
import skimage.io

from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection


# Parse command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("image_file")
parser.add_argument("centroid_file")
args = parser.parse_args()


# Handle IO
if not os.path.exists("outputs"):
    os.mkdir("outputs")
output_name="outputs/"+"centroids_display.pdf"


# Get image and centroids
print "Loading image and centers and generating overlay."
image = skimage.io.imread(args.image_file)
centroids = np.load(args.centroid_file)
x = centroids['x']
y = centroids['y']


# Create and display overlay
fig, ax = plt.subplots()
plt.imshow(image,cmap='gray')
if len(x) != 0:
    plt.plot(x, y, 'r.')
plt.show()
fig.canvas.draw()






