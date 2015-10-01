"""
Contains all functions which calculate RDFs. All are in 2D.
(1) sq_para - calculates the ideal paracrystalline rdf for a square lattice
(2) sq_para_projection - ideal paracrystalline rdf in projection for some thickness c
(3) experimental - calculates the experimental rdf from a set of particle centers
"""

# Import global libraries
import numpy as np


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


# Calculate experimental RDF from particle centers
def experimental(x,y,S,rMax,dr):
    """Compute the two-dimensional pair correlation function, also known
    as the radial distribution function, for a set of circular particles
    contained in a square region of a plane.  This simple function finds
    reference particles such that a circle of radius rMax drawn around the
    particle will fit entirely within the square, eliminating the need to
    compensate for edge effects.  If no such particles exist, an error is
    returned. Try a smaller rMax...or write some code to handle edge effects! ;)

    Arguments:
        x               an array of x positions of centers of particles
        y               an array of y positions of centers of particles
        S               length of each side of the square region of the plane
        rMax            outer diameter of largest annulus
        dr              increment for increasing radius of annulus

    Returns a tuple: (g, radii, interior_x, interior_y)
        g(r)            a numpy array containing the correlation function g(r)
        radii           a numpy array containing the radii of the
                        annuli used to compute g(r)
        interior_x      x coordinates of reference particles
        interior_y      y coordinates of reference particles

    Citation - this function pulled from:
    https://github.com/cfinch/Shocksolution_Examples/tree/master/PairCorrelation
    """
    from numpy import zeros, sqrt, where, pi, average, arange, histogram
    # Number of particles in ring/area of ring/number of reference particles/number density
    # area of ring = pi*(r_outer**2 - r_inner**2)

    # Find particles which are close enough to the box center that a circle of radius
    # rMax will not cross any edge of the box
    bools1 = x>1.1*rMax
    bools2 = x<(S-1.1*rMax)
    bools3 = y>rMax*1.1
    bools4 = y<(S-rMax*1.1)
    interior_indices, = where(bools1*bools2*bools3*bools4)
    num_interior_particles = len(interior_indices)

    print "Analyzing "+str(num_interior_particles)+" particles out of "+str(len(x))

    if num_interior_particles < 1:
        raise  RuntimeError ("No particles found for which a circle of radius rMax\
                will lie entirely within a square of side length S.  Decrease rMax\
                or increase the size of the square.")

    edges = arange(0., rMax+1.1*dr, dr)
    num_increments = len(edges)-1
    g = zeros([num_interior_particles, num_increments])
    radii = zeros(num_increments)
    numberDensity = float(len(x))/S**2

    # Compute pairwise correlation for each interior particle
    for p in range(num_interior_particles):
        index = interior_indices[p]
        d = sqrt((x[index]-x)**2 + (y[index]-y)**2)
        d[index] = 2*rMax

        (result,bins) = histogram(d, bins=edges, normed=False)
        g[p,:] = result/numberDensity

    # Average g(r) for all interior particles and compute radii
    g_average = zeros(num_increments)
    for i in range(num_increments):
        radii[i] = (edges[i] + edges[i+1])/2.
        rOuter = edges[i+1]
        rInner = edges[i]
        g_average[i] = average(g[:,i])/(pi*(rOuter**2 - rInner**2))

    return (g_average, radii, x[interior_indices], y[interior_indices])


if __name__=="__main__":
    # Import global libraries
    import argparse
    import os
    import skimage.io

    # Import local libraries
    import mdscrape as md

    # Get input data, deal with I/O
    parser=argparse.ArgumentParser()
    parser.add_argument("image_file")
    parser.add_argument("centroid_file")
    parser.add_argument("rmax")
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
    g_exp,r_pixels,x0,y0 = experimental(x,y,fov_pixels,float(args.rmax)*spacing_pixels, float(args.dr)*spacing_pixels)
    r_nm = r_pixels / pixels_per_nm  # convert to nm

###### TODO: extract sigma (also possibly mu, c) properly here....!!! ######
###### Then save sigma (etc) in .npz file ##################################
    sigma_nm = 0.4
    mu_nm = 6.5
    c = 3

    g_ideal = sq_para(r_pixels, sigma=sigma_nm*pixels_per_nm, mu=mu_nm*pixels_per_nm)
    g_proj = sq_para_projection(r_pixels, sigma=sigma_nm*pixels_per_nm, mu=mu_nm*pixels_per_nm, c=3)

    print "Saving as {}.npz".format(output_name)
    np.savez(output_name, r_nm=r_nm, g_exp=g_exp, g_ideal=g_ideal, g_proj=g_proj, sigma_nm=sigma_nm,mu_nm=mu_nm,c=c)

