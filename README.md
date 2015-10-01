This file describes the data processing to be performed on the PbSe
QDS data acquired over large fields of view (600-1000 nm) with sufficient
pixel density to resolve QD centeroids.

The goal is to extract the RDF directly from real space, extract a 
self-consistent correlation length which described the length scale over
which long range order decays, and compare this number to the sample
thickness.

Ideally, this entire process will be non-interactive.  As of now, user
input is required.  See below.


########


Script 1: centroids.py
Description: Finds all of the QDs in an image and identifies their centroidswith subpixel resolution.
Usage:  python centroids.py filename.tif threshold
        May need to manually set minDist at line 475 if spacing is not
        extracted correctly
Input:  Image file, as a tiff
        Threshold value (~0.1-0.01 seems to work, takes ~5min)
Output: .npz file with all particle centers and the average spacing.
        Extract by:
            import numpy as np
            data = np.load("centroids.npz")
            x = data['x']   # np.ndarray
            y = data['y']   # np.ndarray
            spacing = float(data['spacing']     # float 


Script 2: view_centroids.py
Description: Loads image and .npz centers file, and displays overlay.  Use
to determine success of a run of centroids.py, or centroids_SP.py.
Usage:      python view_centroids.py image_file.tif centroid_file.npz
input:  image file, as a tiff
        centers file, as .npz
output: none


###### TODO: write this script ######
Script 3: edit_centroids.py
Description: Displays image with centers overlaid.  Allows erroneous
centers to be removed, and missed centers to be added manually. Be sure to
refine centers after with centroids_SP.
input:  image file, as a tiff
        centers file, as .npz
output: centers file, as .npz


Script 4: SP_centroids.py
Description: Refines a set of centroids with subpixel resolution.
Usage:      python view_centroids.py image_file.tif centroid_file.npz
Input:  image file, as tiff
        centers file, as .npz
output: centers file, as .npz


Script 5: rdf.py
Description: Generates the RDF for the data given.  Finds the best fit
ideal paracrystalline and projection paracrystalline RDF.
Note that rdf.py does not save plots of g(r) - these must be calculated 
separately from plot_rdf.py
Usage:      python rdf.py image_file.tif centroid_file.npz 
Input:  -.tif image file
        -.npz file with particle centers
        -rmax, maximum distance to calculate g(r) to
        -dr, step size in g(r) calculation
Output: .npz files of...
        -Experimental RDF
        -ideal paracrystalline RDF fit
        -projection RDF fit


Script 6: rdf_plot.py
Description: Makes plots of various rdfs (experimental, ideal, projection)
Input:  -.npz file with g(r) (exper, ideal, proj)


Script 7: corr_len.py
Description: Finds the correlation length for an experimental RDF
Input: .npz file with RDF
Output: -.npz file with correlation length and envelope curve used in fit
        -.pdf with RDF and envelope fit to find correlation length


