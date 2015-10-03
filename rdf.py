"""
Contains all functions which calculate RDFs. All are in 2D.
(1) sq_para - calculates the ideal paracrystalline rdf for a square lattice
(2) sq_para_projection - ideal paracrystalline rdf in projection for some thickness c
(3) experimental - calculates the experimental rdf from a set of particle centers
"""

# Import global libraries
import numpy as np
import warnings

# Paracrystal fit for square lattice...
def sq_para(r_data, sigma, mu):
    density = 1.0/(mu**2)
    g=np.zeros(len(r_data))
    sigma_uv = 0
    r_uv = 0
    uv_max =  int(max(r_data)/mu) + 2  # Plus 2 accounts for edges

    for u in range(-uv_max,uv_max):
        for v in range(-uv_max,uv_max):
            r_uv = np.sqrt(u**2 + v**2)*mu
            sigma_uv = np.sqrt(abs(u)+abs(v))*sigma
            if r_uv > 0:
                np.add( g[1::], (1.0/(r_data[1::]*sigma_uv)) * np.exp(-(r_data[1::]-r_uv)**2/(2*sigma_uv**2)), out=g[1::] )
    g = g/((2*np.pi)**(1.5)*density)

    return g


# Paracrtystal fit for square lattice with STEM projection through c layers
def sq_para_projection(r_data, sigma, mu, c):
    density = 1.0/(mu**2)
    g=np.zeros(len(r_data))
    sigma_uv = 0
    r_uv = 0
    uv_max =  int(max(r_data)/mu) + 2 # Plus 2 accounts for edges

    # Add in the projection term
    summ = 0
    for i in range(1,c+1):
        summ += (c - i)**2
    summ = summ/float(c**2)

    for u in range(-uv_max,uv_max):
        for v in range(-uv_max,uv_max):
            r_uv = np.sqrt(u**2 + v**2)*mu
            sigma_uv = np.sqrt(abs(u)+abs(v)+summ)*sigma
            if r_uv > 0:
                np.add( g[1::], (1.0/(r_data[1::]*sigma_uv)) * np.exp(-(r_data[1::]-r_uv)**2/(2*sigma_uv**2)), out=g[1::] )
    g = g/((2*np.pi)**(1.5)*density)

    return g


# Calculate the excess area that must be subtracted off an annulus with inner radius r and
# outer radius R with a center x from the box edge, where x<r.
def calc_excess(r, R, x):
    # Should be one line of code, but np.where appears to be written poorly
    # Mask and alternate r/R arrays prevent acrcos and sqrt from throwing errors
    mask = x<r
    r_mask = np.where(mask, r, x+1)   # (x+1) is a dummy, never appears in final output
    R_mask = np.where(mask, R, x+1)   #   but ensures arccos and sqrt won't break
    return np.where(mask, R_mask**2*np.arccos(x/R_mask) - r_mask**2*np.arccos(x/r_mask)
                        - x*(np.sqrt(R_mask**2-x**2) - np.sqrt(r_mask**2-x**2)), 0)


# Calculate overlap of x and y excess
def calc_xy_excess(r, R, x, y):
    mask = (x<r)*(y<r)  # i.e. and statement
    r_mask = np.where(mask, r, x+y+1)   # as in calc_excess, (x+y+1) is a dummy variable,
    R_mask = np.where(mask, R, x+y+1)   #   chosen such that arccos and sqrts are always valid
    return np.where(mask, 0.5*( R_mask**2*(np.arccos(x/R_mask)+np.arccos(y/R_mask)-np.pi/2)
                    - r_mask**2*(np.arccos(x/r_mask)+np.arccos(y/r_mask)-np.pi/2)
                    - x*(np.sqrt(R_mask**2-x**2) - np.sqrt(r_mask**2-x**2))
                    - y*(np.sqrt(R_mask**2-y**2) - np.sqrt(r_mask**2-y**2)) ), 0)


# Find members of a sorted list closest to a given value.
# Returns values bounding input value
def nearest_vals_array(arr, vals_ar):
    """
    Accepts:
        arr - a sorted array
        vals_arr - an array of numbers
    Returns:
        lower - an array with the shape of vals_ar, populated by the closest element of
                arr to each value in vals_ar that is less than that value
        upper - same idea
    """
    if vals_ar.min()<arr.min() or vals_ar.max()>arr.max():
        print "Error: value in vals_ar outside range of arr."
        return 0,0
    lower = np.empty_like(vals_ar)
    upper = np.empty_like(vals_ar)
    for i in range(len(vals_ar)):
        ind = (np.abs(arr - vals_ar[i])).argmin()
        if (arr[ind]-vals_ar[i]) <= 0 :
            lower[i],upper[i] = arr[ind], arr[ind+1]
        else:
            lower[i],upper[i] = arr[ind-1], arr[ind]
    return lower, upper


# Calculate experimental RDF from particle centers
def experimental(x,y,fov,dr):
    """
    Calculates the 2D radial distribution function from experimental data.
    Handles boundaries by calculating appropriate normalization area of each annulus near edges

    Inputs: x - ndarray of x coordinates of particle centers, in pixels
            y - ndarray of y coordinates, in pixels
            fov - float of bounding box side length, in pixels
            dr - binsize for rdf histogram, in pixels
    Output: g_average - ndarray of the radial distribution function
            radii - ndarray of the radius at the center of each bin used to calculate g(r)
    """
    # Number of particles in ring/area of ring/number of reference particles/number density
    # area of ring = pi*(r_outer**2 - r_inner**2)
    import warnings
    warnings.filterwarnings('error')

    rbins = np.arange(0., fov, dr)
    num_increments = len(rbins)
    g = np.zeros([len(x), num_increments-1])
    radii = np.zeros(num_increments)
    numberDensity = float(len(x))/fov**2

    # Compute pairwise correlation for each interior particle
    for i in range(len(x)):
        print "Analyzing particle {} of {}".format(i, len(x))
        d = np.sqrt((x[i]-x)**2 + (y[i]-y)**2)
        d_lower, d_upper = nearest_vals_array(np.arange(0, np.sqrt(2)*fov+dr, dr),d)
        d[i] = fov

        # Determine normalization, if (x[i],y[i]) is near an edge
        x_e, y_e = min(x[i], fov-x[i]), min(y[i], fov-y[i]) # Find distance to edge

        x_excess = calc_excess(d_lower, d_upper, x_e)
        y_excess = calc_excess(d_lower, d_upper, y_e)
        xy_excess = calc_xy_excess(d_lower, d_upper, x_e, y_e)

        annulus_areas = np.pi*( d_upper**2 - d_lower**2 )
        area_normalization = annulus_areas - x_excess - y_excess + xy_excess


        (result,bins) = np.histogram(d, bins=rbins, normed=False, weights=1.0/area_normalization)
        g[i,:] = result/numberDensity


    # Average g(r) for all interior particles and compute radii
    g_average = np.zeros(num_increments)
    for i in range(num_increments-1):
        radii[i] = (bins[i] + bins[i+1])/2.
        g_average[i] = np.average(g[:,i])

    return (g_average, radii)




if __name__=="__main__":
    # Import global libraries
    import argparse
    import os
    import skimage.io
    import time

    # Import local libraries
    import mdscrape as md

    # Get input data, deal with I/O
    parser=argparse.ArgumentParser()
    parser.add_argument("image_file")
    parser.add_argument("centroid_file")
    parser.add_argument("dr")
    args=parser.parse_args()

    if not os.path.exists("outputs"):
        os.mkdir("outputs")
    output_name="outputs/"+"rdf"

    print "Loading image and centers..."
    image = skimage.io.imread(args.image_file)
    centroids = np.load(args.centroid_file)
    x = centroids['x']
    y = centroids['y']
    spacing_pixels = centroids['spacing']

    print "Done. Loaded {} centroids.\nExtracting metadata...".format(len(x))
    metadata = md.extract_metadata(args.image_file)
    fov_nm = metadata['fov']
    fov_units = metadata['fov_units']
    fov_pixels = metadata['pixels']
    pixels_per_nm = float(fov_pixels)/float(fov_nm)

    # Perform rdf calculations
    print "Done. Calculating experimental RDF."
    time_init = time.time()
    g_exp,r_pixels = experimental(x,y,fov_pixels, float(args.dr)*spacing_pixels)
    r_nm = r_pixels / pixels_per_nm  # convert to nm

###### TODO: extract sigma (also possibly mu, c) properly here....!!! ######
###### Then save sigma (etc) in .npz file ##################################
    sigma_nm = 0.4
    mu_nm = 6.5
    c = 3

#    g_ideal = sq_para(r_pixels, sigma=sigma_nm*pixels_per_nm, mu=mu_nm*pixels_per_nm)
#    g_proj = sq_para_projection(r_pixels, sigma=sigma_nm*pixels_per_nm, mu=mu_nm*pixels_per_nm, c=3)
    g_ideal=0
    g_proj=0

    print "Done. Calculation took {} seconds.".format(time.time()-time_init)
    print "Saving as {}.npz".format(output_name)
    np.savez(output_name, r_nm=r_nm, g_exp=g_exp, g_ideal=g_ideal, g_proj=g_proj, sigma_nm=sigma_nm,mu_nm=mu_nm,c=c)

