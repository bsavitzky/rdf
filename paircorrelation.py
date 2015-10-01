"""
https://github.com/cfinch/Shocksolution_Examples/tree/master/PairCorrelation
This module contains routines for analyzing distributions of particles.  Currently,
this consists of:
    PairCorrelationFunction_2D
    PairCorrelationFunction_3D
    computePhi
        Returns True if the given point at (x,y) lies within any of the particles with
        centers given by (x_array, y_array) with radius "radius."  Otherwise, returns False.
    computePhi_PBC
    compute_AVF
"""

def PairCorrelationFunction_2D(x,y,S,rMax,dr):
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
    numberDensity = len(x)/S**2

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
####

def PairCorrelationFunction_3D(x,y,z,S,rMax,dr):
    """Compute the three-dimensional pair correlation function for a set of
    spherical particles contained in a cube with side length S.  This simple 
    function finds reference particles such that a sphere of radius rMax drawn
    around the particle will fit entirely within the cube, eliminating the need
    to compensate for edge effects.  If no such particles exist, an error is 
    returned.  Try a smaller rMax...or write some code to handle edge effects! ;) 
    
    Arguments:
        x               an array of x positions of centers of particles
        y               an array of y positions of centers of particles
        z               an array of z positions of centers of particles
        S               length of each side of the cube in space
        rMax            outer diameter of largest spherical shell
        dr              increment for increasing radius of spherical shell

    Returns a tuple: (g, radii, interior_x, interior_y, interior_z)
        g(r)            a numpy array containing the correlation function g(r)
        radii           a numpy array containing the radii of the
                        spherical shells used to compute g(r)
        interior_x      x coordinates of reference particles
        interior_y      y coordinates of reference particles
        interior_z      z coordinates of reference particles
    """
    from numpy import zeros, sqrt, where, pi, average, arange, histogram

    # Find particles which are close enough to the cube center that a sphere of radius
    # rMax will not cross any face of the cube
    bools1 = x>rMax
    bools2 = x<(S-rMax)
    bools3 = y>rMax
    bools4 = y<(S-rMax)
    bools5 = z>rMax
    bools6 = z<(S-rMax)

    interior_indices, = where(bools1*bools2*bools3*bools4*bools5*bools6)
    num_interior_particles = len(interior_indices)

    if num_interior_particles < 1:
        raise  RuntimeError ("No particles found for which a sphere of radius rMax\
                will lie entirely within a cube of side length S.  Decrease rMax\
                or increase the size of the cube.")

    edges = arange(0., rMax+1.1*dr, dr)
    num_increments = len(edges)-1
    g = zeros([num_interior_particles, num_increments])
    radii = zeros(num_increments)
    numberDensity = len(x)/S**3

    # Compute pairwise correlation for each interior particle
    for p in range(num_interior_particles):
        index = interior_indices[p]
        d = sqrt((x[index]-x)**2 + (y[index]-y)**2 + (z[index]-z)**2)
        d[index] = 2*rMax

        (result,bins) = histogram(d, bins=edges, normed=False)
        g[p,:] = result/numberDensity
        
    # Average g(r) for all interior particles and compute radii
    g_average = zeros(num_increments)
    for i in range(num_increments):
        radii[i] = (edges[i] + edges[i+1])/2.        
        rOuter = edges[i+1]
        rInner = edges[i]
        g_average[i] = average(g[:,i])/(4./3.*pi*(rOuter**3 - rInner**3))

    return (g_average, radii, x[interior_indices], y[interior_indices], z[interior_indices])
    # Number of particles in shell/total number of particles/volume of shell/number density
    # shell volume = 4/3*pi(r_outer**3-r_inner**3)
####

def computePhi(x_grid,y_grid,x_array,y_array,r):
    """Returns True if the given point at (x,y) lies within any of the particles with
    centers given by (x_array, y_array) with radius "radius."  Otherwise, returns False.
    """
    from scipy import weave

    numParticles = len(x_array)
    gridsize = len(x_grid)

    code = r"""
           double d2;
           double radius_squared = 4.0*pow(r,2);
           int avail = gridsize*gridsize;
           bool collision_not_found = true;
           int p = 0;

            for (int i=0; i<gridsize; i++) {
                for (int j=0; j<gridsize; j++) {
                    p = 0;
                    collision_not_found = true;
                    while ((p<numParticles) && collision_not_found) {
                        if ((pow(x_grid[i]-x_array[p], 2) + pow(y_grid[j]-y_array[p],2)) < radius_squared) {
                            avail = avail - 1;
                            collision_not_found = false;
                            }
                        p+=1;
                    }
                }
            }
           return_val = (float)avail/(float)(gridsize*gridsize);
    """
    return weave.inline(code, ['x_grid', 'y_grid', 'x_array', 'y_array', 'r', 'numParticles', 'gridsize'], \
        compiler='gcc')

def computePhi_PBC(x_grid,y_grid,x_array,y_array,r,width):
    """Returns True if the given point at (x,y) lies within any of the particles with
    centers given by (x_array, y_array) with radius "radius."  Otherwise, returns False.
    """
    from scipy import weave

    numParticles = len(x_array)
    gridsize = len(x_grid)

    code = r"""
           double d2;
           double radius_squared = 4.0*pow(r,2);
           int avail = gridsize*gridsize;
           bool collision_not_found = true;
           int p = 0;
           double x, y, x_p, y_p;

            for (int i=0; i<gridsize; i++) {
                x = x_grid[i];
                for (int j=0; j<gridsize; j++) {
                    y = y_grid[j];
                    p = 0;
                    collision_not_found = true;
                    while ((p<numParticles) && collision_not_found) {
                        // Handle periodic boundary conditions
                        // x sides
                        if ((x < 2*r) && (x_array[p] > width-2*r)) x_p = x_array[p]-width;
                        else if ((x > width-2*r) && (x_array[p] < 2*r)) x_p = width+x_array[p];
                        else x_p = x_array[p];
                        // y
                        if ((y < 2*r) && (y_array[p] > width-2*r)) y_p = y_array[p]-width;
                        else if ((y > width-2*r) && (y_array[p] < 2*r)) y_p = width+y_array[p];
                        else y_p = y_array[p];

                        if ((pow(x-x_p, 2) + pow(y-y_p,2)) < radius_squared) {
                            avail = avail - 1;
                            collision_not_found = false;
                        }
                        p+=1;
                    }
                }
            }
           return_val = (float)avail/(float)(gridsize*gridsize);
    """
    return weave.inline(code, ['x_grid', 'y_grid', 'x_array', 'y_array', 'r', 'width', 'numParticles', 'gridsize'], \
        compiler='gcc')

def compute_AVF(test_x, test_y, test_h, adsorbed_x, adsorbed_y, r, num_replicates, random_radius):
    """ Computes the available volume function as a function of distance from the surface.  This is done
    by checking whether a sphere at each of the grid points (test_x, test_y, test_h) overlaps with 
    any of the adsorbed spheres centered at (adsorbed_x, adsorbed_y, r).  The test particles are offset by a
    random vector for each replicate.  If image particles are required for periodic boundary conditions,
    it is assumed that these are already included in the adsorbed_x and adsorbed_y arrays.
    
    Assumes all spheres have radius "r".
    
    Arguments:
    test_x          array of x coordinates of test particles
    test_y          array of y coordinates of test particles
    test_h          array of z coordinates of test particles
    adsorbed_x      array of x coordinates of adsorbed particles
    adsorbed_y      array of y coordinates of adsorbed particles
    r               radius of all particles
    num_replicates  number of replicates to perform with different random offsets at each z level
    random_radius   each test particle will be randomly placed within this distance of its nominal position
    
    Returns the available volume fraction for each h value in test_h.
    """
    from compiled import compiled_functions

    return compiled_functions.compute_avf(len(test_x), test_x, test_y, len(test_h), test_h, len(adsorbed_x), \
            adsorbed_x, adsorbed_y, r, num_replicates, random_radius)

# Vector correlation function
def cov(u,v):
    from numpy import mean, sum
    return 1.0 / (len(u) - 1.0) * sum((u - mean(u)) * (v - mean(v)))

def vector_corr(w1, w2):
    """
    Returns the vector correlation coefficient.
    Arguments:
    w1  vector of x coordinates
    w2  vector of y coordinates

    References:
    A Proposed Definition for Vector Correlation in Geophysics: 
        Theory and Application
    Crosby, D. S.; Breaker, L. C.; Gemmill, W. H.
    Journal of Atmospheric and Oceanic Technology, vol. 10, issue 3, p. 355
    1993
    DOI: 10.1175/1520-0426(1993)010<0355:APDFVC>2.0.CO;2

	The Application of a Technique for Vector Correlation to Problems in 
        Meteorology and Oceanography
    Breaker, L. C.; Gemmill, W. H.; Crosby, D. S.
	Journal of Applied Meteorology, vol. 33, Issue 11, pp.1354-1365
	11/1994
    DOI: 10.1175/1520-0450(1994)033<1354:TAOATF>2.0.CO;2
    """
    u1 = w1[:,0]; v1 = w1[:,1]
    u2 = w2[:,0]; v2 = w2[:,1]

    f = cov(u1,u1)*(cov(u2,u2)*cov(v1,v2)**2 + cov(v2,v2)*cov(v1,u2)**2) \
            + cov(v1,v1)*(cov(u2,u2)*cov(u1,v2)**2 \
            + cov(v2,v2)*cov(u1,u2)**2) \
            + 2*(cov(u1,v1)*cov(u1,v2)*cov(v1,u2)*cov(u2,v2)) \
            + 2*(cov(u1,v1)*cov(u1,u2)*cov(v1,v2)*cov(u2,v2)) \
            - 2*(cov(u1,u1)*cov(v1,u2)*cov(v1,v2)*cov(u2,v2)) \
            - 2*(cov(v1,v1)*cov(u1,u2)*cov(u1,v2)*cov(u2,v2)) \
            - 2*(cov(u2,u2)*cov(u1,v1)*cov(u1,v2)*cov(v1,v2)) \
            - 2*(cov(v2,v2)*cov(u1,v1)*cov(u1,u2)*cov(v1,u2))

    g = (cov(u1,u1)*cov(v1,v1) \
            - cov(u1,v1)**2)*(cov(u2,u2)*cov(v2,v2)-cov(u2,v2)**2)

    return f/g

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Circle

def generate_test_grid(shape, period, circle_radius, radius_std, lattice_std, intensity_std ):

    x = np.linspace(circle_radius, shape[1]-circle_radius, (shape[1]-2*circle_radius)/period)
    y = np.linspace(circle_radius, shape[0]-circle_radius, (shape[0]-2*circle_radius)/period)

    grid_x, grid_y = np.meshgrid(x,y)

    grid_x = grid_x.reshape((-1,1))
    grid_y = grid_y.reshape((-1,1))

    if lattice_std > 0:
        np.add(grid_x,period*np.random.normal(loc=0, scale=lattice_std, size=grid_x.shape),out=grid_x)
        np.add(grid_y,period*np.random.normal(loc=0, scale=lattice_std, size=grid_y.shape),out=grid_y)

    cell_patches = []

    for point_index in range(len(grid_x)):
        cell_patches.append(Circle((grid_x[point_index],grid_y[point_index]),radius=circle_radius*radius_offset(radius_std),facecolor=(0,0,0),edgecolor='none',alpha=np.random.normal(loc=0.5, scale=intensity_std)))

    pc = PatchCollection(cell_patches,match_original=True)

    plt.figure(1)
    plt.gca().add_collection(pc)

    # set the limits for the plot
    # set the x axis range
    plt.gca().set_xlim(0, shape[1])

    # set the y-axis range and flip the y-axis
    plt.gca().set_ylim(0, shape[0])

    # save this plot to a file
    plt.gca().set_axis_off()
    plt.savefig('../test_data/input/test_image_10p_contrast.png',bbox_inches='tight')

def radius_offset(radius_std):
    if radius_std > 0:
        return np.random.normal(loc=1, scale=radius_std)
    else:
        return 1
        
def generate_disordered_rdf(a, b, std, max_distance, resolution):
    
    
    g = np.zeros(max_distance/resolution)
    r = np.linspace(0,max_distance,max_distance/resolution)
    
    # assumes a, b are primitive lattice vectors
    density = 1.0/(a[0]*b[1] - b[0]*a[1])
    
    a_mag = np.sqrt(a[0]**2 + a[1]**2)
    b_mag = np.sqrt(b[0]**2 + b[1]**2)
    max_peak_num = int(np.ceil(max_distance/a_mag)+1)

    dist = []
    width = std

    for u in range(-max_peak_num, max_peak_num):
        for v in range(-max_peak_num, max_peak_num):
            peak_r = np.sqrt((u*a[0] + v*b[0])**2 + (u*a[1] + v*b[1])**2)
            if peak_r > 0:
                np.add(g[1::],(2.0*np.pi*r[1::]*width)**(-1)*np.exp(-(r[1::]-peak_r)**2/(2*(width)**2)),out=g[1::])

    # normalize by density
    # not sure why the factor 4/10 comes in...
    g = g*4.0/(10.0*density)
    
    return r,g


def generate_paracrystal_rdf(a, b, std, max_distance, resolution):
    
    g = np.zeros(max_distance/resolution)
    r = np.linspace(0,max_distance,max_distance/resolution)
    
    # assumes a, b are primitive lattice vectors
    density = 1.0/(a[0]*b[1] - b[0]*a[1])
    
    a_mag = np.sqrt(a[0]**2 + a[1]**2)
    b_mag = np.sqrt(b[0]**2 + b[1]**2)
    max_peak_num = int(np.ceil(max_distance/a_mag)+1)

    width = 0

    for u in range(-max_peak_num, max_peak_num):
        for v in range(-max_peak_num, max_peak_num):
            peak_r = np.sqrt((u*a[0] + v*b[0])**2 + (u*a[1] + v*b[1])**2)
            if peak_r > 0:
                width = std * np.sqrt(np.sqrt((peak_r)**2 / ( a_mag * b_mag ) ))
                np.add(g[1::],(2.0*np.pi*r[1::]*width)**(-1)*np.exp(-(r[1::]-peak_r)**2/(2*(width)**2)),out=g[1::])

    # Use this code to plot the individual pair distribution functions for each lattice spacing in dist
    # plt.figure(1)
    # unique_dist,unique_count = np.unique(dist,return_index=False,return_inverse=False,return_counts=True)
    #
    # for peak,count in zip(unique_dist[1::],unique_count[1::]):
    #     peak_r = peak*period
    #     width = np.sqrt(peak)*std
    #     plt.plot(r,count*(2.0*np.pi*width)**(-1)*np.exp(-(r-peak_r)**2/(2*(width)**2))*(period**2)/(peak_r*2.5))
    # plt.xlabel('r')
    # plt.ylabel('$g(r_k)$')
    # plt.show()

    # normalize by density
    # not sure why the factor 4/10 comes in...
    g = g*4.0/(10.0*density)
    
    return r,g

# assume a,b lattice vectors are equal length
def fit_hex_paracrystal_rdf(r_data, std, r0):
    r,g = generate_paracrystal_rdf((r0,0.0), (-0.5*r0, r0 * np.sqrt(3.0)/2.0), std, r_data[-1]+(r_data[1]-r_data[0]), r_data[1]-r_data[0])
    return g
        
def fit_square_paracrystal_rdf(r_data, std, r0):
    r,g = generate_paracrystal_rdf((r0,0.0), (0.0, r0), std, r_data[-1]+(r_data[1]-r_data[0]), r_data[1]-r_data[0])
    return g
    
def fit_hex_disordered_rdf(r_data, std, r0):
    r,g = generate_disordered_rdf((r0,0.0), (-0.5*r0, r0 * np.sqrt(3.0)/2.0), std, r_data[-1]+(r_data[1]-r_data[0]), r_data[1]-r_data[0])
    return g
    
def fit_square_disordered_rdf(r_data, std, r0):
    r,g = generate_disordered_rdf((r0,0.0), (0.0, r0), std, r_data[-1]+(r_data[1]-r_data[0]), r_data[1]-r_data[0])
    return g

#generate_test_grid((1032,1376), period=20.0, circle_radius=10.0, radius_std=0.0, lattice_std=0.0, intensity_std=0.1)

# from scipy.optimize import curve_fit
#rh,gh = generate_paracrystal_rdf(a=(1.0,0.0), b=(-0.5*1.0,1.0*np.sqrt(3.0)/2.0), std=0.05, max_distance=50.0, resolution=0.01)
# r1,g1 = generate_paracrystal_rdf(a=(1.0,0.0), b=(0.0,1.0), std=0.05, max_distance=20.0, resolution=0.01)
# r2,g2 = generate_paracrystal_rdf(a=(1.0,0.0), b=(0.0,1.0), std=0.05, max_distance=20.0, resolution=0.01)
#r3,g3 = generate_paracrystal_rdf(a=(1.0,0.0), b=(0.0,1.0), std=0.05, max_distance=20.0, resolution=0.01)
# plt.plot(r1, g1, 'k-')
# plt.plot(r2, g2, 'b-')
#plt.plot(r3, g3, 'r-')
# plt.gca().legend()
#plt.show()

# g = g + 0.2*np.random.normal(size=len(r))
# popt, pcov = curve_fit(fit_paracrystal, r, g, p0=0.5)
# print(popt)
# print(pcov)
