"""
Calculates radial distribution function for paracrystal and projection paracrystal models.
Fits to experimental RDF curve.
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
    parser.add_argument("exp_rdf_file")
    args=parser.parse_args()

    if not os.path.exists("outputs"):
        os.mkdir("outputs")
    output_name="outputs/"+"rdf_all"

    print "Loading experimental radial distribution function file..."
    exp_rdf_files = np.load(args.exp_rdf_file)
    r_nm = centroids['r_nm']
    g_exp = centroids['g_exp']

    # Perform rdf calculations
    print "Done. Fitting paracrystalline and projection paracrystalline RDFs..."
    time_init = time.time()
    # TODO: figure out how to do this fit properly.....!
    sigma_nm = 0.4
    mu_nm = 6.5
    c = 3
    g_ideal = sq_para(r_pixels, sigma=sigma_nm*pixels_per_nm, mu=mu_nm*pixels_per_nm)
    g_proj = sq_para_projection(r_pixels, sigma=sigma_nm*pixels_per_nm, mu=mu_nm*pixels_per_nm, c=3)

    print "Done. Calculation took {} seconds.".format(time.time()-time_init)
    print "Saving as {}.npz".format(output_name)
    np.savez(output_name, r_nm=r_nm, g_exp=g_exp, g_ideal=g_ideal, g_proj=g_proj, sigma_nm=sigma_nm,mu_nm=mu_nm,c=c)

