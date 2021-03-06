"""
Calculates the RDF from centroid positions.  Handles boundary conditions such that all
particles can be used.
"""

# Import global libraries
import numpy as np

# Calculate the excess area that must be subtracted off an annulus with inner radius r and
# outer radius R with a center x from the box edge, where x<r.
def calc_excess(r, R, x):
    # Should be one line of code, but np.where appears to be written poorly
    # Mask and alternate r/R arrays prevent acrcos and sqrt from throwing errors
    mask = (x<r)
    r_mask = np.where(mask, r, x+1)   # (x+1) is a dummy, never appears in final output
    R_mask = np.where(mask, R, x+1)   #   but ensures arccos and sqrt won't break
#    print r_mask.min(),R_mask.min()
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

    rbins = np.arange(0., fov/2.0, dr)
    num_increments = len(rbins)
    g = np.zeros([len(x), num_increments-1])
    radii = np.zeros(num_increments-1)
    numberDensity = float(len(x))/fov**2

    # Compute pairwise correlation for each interior particle
    for i in range(len(x)):
#    for i in range(825,826):
        print "Analyzing particle {} of {}".format(i+1, len(x))
        d = np.sqrt((x[i]-x)**2 + (y[i]-y)**2)
        d_lower, d_upper = nearest_vals_array(np.arange(0, np.sqrt(2)*fov+dr, dr),d)
        d[i] = 0

        # Determine normalization, if (x[i],y[i]) is near an edge
        x_e, y_e = min(x[i], fov-x[i]), min(y[i], fov-y[i]) # Find distance to edge

#        print x_e, y_e

        x_excess = calc_excess(d_lower, d_upper, x_e)
        y_excess = calc_excess(d_lower, d_upper, y_e)
        xy_excess = calc_xy_excess(d_lower, d_upper, x_e, y_e)

        annulus_areas = np.pi*( d_upper**2 - d_lower**2 )
        area_normalization = annulus_areas - x_excess - y_excess + xy_excess

        (result,bins) = np.histogram(d, bins=rbins, normed=False, weights=1.0/area_normalization)
        g[i,:] = result/numberDensity

    # Average g(r) for all interior particles and compute radii
    g_average = np.zeros(num_increments-1)
    for i in range(num_increments-1):
        radii[i] = bins[i] + dr/2.0
        g_average[i] = np.average(g[:,i])

    return (g_average, radii)

def stretch_lattice(x, y, fov, stretch, angle):
    if stretch < 1:
        raise RunTimeError("Error. Stretch values must be >1.  {} entered.".format(stretch))
    print "Stretching by a factor of {} at an angle of {} degrees to the x-axis about the image center".format(stretch, angle)
    x_new=[]
    y_new=[]
    x_shifted, y_shifted = x-fov/2.0, y-fov/2.0
    x = ( x_shifted + np.sqrt(x_shifted**2+y_shifted**2) * 
          np.cos(np.angle(x_shifted+y_shifted*1j)-np.radians(angle)) *
          np.cos(np.radians(angle))*(stretch-1) + fov/2.0 )
    y = ( y_shifted + np.sqrt(x_shifted**2+y_shifted**2) * 
          np.cos(np.angle(x_shifted+y_shifted*1j)-np.radians(angle)) *
          np.sin(np.radians(angle))*(stretch-1)  + fov/2.0 )
    for i in range(len(x)):
        if x[i] < fov and x[i] >= 0 and y[i] < fov and y[i] >= 0:
            x_new.append(x[i])
            y_new.append(y[i])
    print "{} particles discarded".format(len(x) - len(x_new))
    return np.array(x_new), np.array(y_new)

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
    parser.add_argument("--dr", type=float, default=0.1, help="Bin size for rdf")
    parser.add_argument("-o","--output",type=str, default="outputs/", help="Directory for output files")
    parser.add_argument("-s","--stretch",type=float, nargs=2,
                        help="Amount to stretch lattice (>=1), and angle of stretch from x-axis, in degrees.")
    args=parser.parse_args()

    # Make output path with appropriate metadata
    if not os.path.exists(args.output):
        os.makedirs(args.output)
    if args.output[-1]!="/":
        output_name = args.output+"/"
    else:
        output_name = args.output
    output_name+="rdf_exp"
    if args.stretch:
       output_name += "_stretch={}_stretchangle={}".format(args.stretch[0], args.stretch[1])

    # Read data from file
    print "Loading image and centers..."
    image = skimage.io.imread(args.image_file)
    centroids = np.load(args.centroid_file)
    x = centroids['x']
    y = centroids['y']
    spacing_pixels = centroids['spacing']
    fov_pixels = centroids['fov_pixels']

    # Perform any stretches
    if args.stretch:
        x, y = stretch_lattice(x, y, fov_pixels, args.stretch[0], args.stretch[1])

    print "Done. Loaded {} centroids.\nExtracting metadata...".format(len(x))
    metadata = md.extract_metadata(args.image_file)
    fov_nm = metadata['fov']
    fov_units = metadata['fov_units']
    pixels_per_nm = float(fov_pixels)/float(fov_nm)

    # Perform rdf calculations
    print "Done. Calculating experimental RDF."
    time_init = time.time()
    g_exp,r_pixels = experimental(x,y,fov_pixels, float(args.dr)*spacing_pixels)
    r_nm = r_pixels / pixels_per_nm  # convert to nm

    print "Done. Calculation took {} seconds.".format(time.time()-time_init)
    print "Saving as {}.npz".format(output_name)
    np.savez(output_name, r_nm=r_nm, g_exp=g_exp)
