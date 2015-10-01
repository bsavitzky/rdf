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
from matplotlib import ion


# Handle I/O

# Parse command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("image_file")
parser.add_argument("centroid_file")
args = parser.parse_args()

# Output path
if not os.path.exists("outputs"):
    os.mkdir("outputs")
output_name="outputs/"+"centroids_edited"

# Get image and centroids
print "Loading image and centers and generating overlay."
image = skimage.io.imread(args.image_file)
centroids = np.load(args.centroid_file)
x = centroids['x']
y = centroids['y']
spacing = float(centroids['spacing'])

# Turn on mpl interactive mode
#ion()

# Define callback event
def onclick(event):
    print ("Button={}\nx={}\ny={}\nxdata={}\nydata={}".format(event.button,
           event.x, event.y, event.xdata, event.ydata))
    return None

# Create and display overlay
fig, ax = plt.subplots()
plt.imshow(image,cmap='gray')
if len(x) != 0:
    plt.plot(x, y, 'r.')
fig.canvas.draw()


### Add code here ###
cid = fig.canvas.mpl_connect('button_press_event',onclick)
fig.canvas.mpl_disconnect(cid)

plt.show()
#raw_input("Click on a spot!")

x_edited=x
y_edited=y

# Save as .npz
print "Saving as {}.npz".format(output_name)
np.savez(output_name, x=x_edited, y=y_edited, spacing=np.array(spacing))




