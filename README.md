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
Optionally highlights a given particle with the -pn flag.
Usage:  python view_centroids.py [-pn particleNumber] image_file.tif 
centroid_file.npz
input:  image file, as a tiff
        centers file, as .npz
output: none
options:    -pn #   Highlights particle number #


Script 3: edit_centroids.ipynb
Description: Displays image with centers overlaid.  Allows erroneous
centers to be removed, and missed centers to be added manually. Be sure to
refine centers after with centroids_SP.
input:  image file, as a tiff
        centers file, as .npz
output: centers file, as .npz


Script 4: SP_centroids.py
Description: Refines a set of centroids with subpixel resolution. Returns
pixelated centroid for any particles which cannot be fit to a gaussian.  
Throws out any particles on images egde, then shifts the positions and FOV 
accordingly in the output .npz file.
Usage:      python SP_centroids.py image_file.tif centroid_file.npz
Input:  image file, as tiff
        centers file, as .npz
output: centers file, as .npz


Script 5: rdf_exp.py
Description: Generates the RDF for the data given, using centroids.
Note that rdf.py does not save plots of g(r) - these must be calculated 
separately from plot_rdf.py.  Also note that 'fov_pixels' should be read
directly from the centroids file, as it will have been shifted while 
eliminating edger particles - do NOT use 'fov_pixels' from mdscape.py!
Usage:      python rdf.py image_file.tif centroid_file.npz dr 
Input:  -.tif image file
        -.npz file with particle centers
        -dr, step size in g(r) calculation, in particle diameters
Output: .npz files with:
        -Radii, in nm
        -Experimental RDF
        Grab data via:
            data = np.load('rdf_exp.npz')
            r_nm, g_exp = data['r_nm'],data['g_exp']


Script 6: rdf_models.py
Desciption: Calculates best fit paracrystalline and projection paracrystalline
rdfs based on g_exp
Usage:      python rdf_models.py exp_rdf_file.npz
Input:  -.npz file with r and g_exp
Output: -Radii, in nm
        -Experimental RDF
        -Ideal paracrystalline RDF
        -Projection RDF
        -best fit sigma, mu, and c
        Grab data via:
            data = np.load('rdf_exp.npz')
            r_nm, g_exp, g_ideal, g_proj, sigma_nm, mu_nm, c = data['r_nm'],data['g_exp'],...etc


Script 7: plot_rdf.py
Description: Makes plots of various rdfs (experimental, ideal, projection)
Input:  -.npz file with g(r) (exper, ideal, proj)


Script 8: corr_len.py
Description: Finds the correlation length for an experimental RDF
Input: .npz file with RDF
Output: -.npz file with correlation length and envelope curve used in fit
        -.pdf with RDF and envelope fit to find correlation length


