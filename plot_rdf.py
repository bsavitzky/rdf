# Import libraries
import numpy as np
import matplotlib.pyplot as plt
import os
import argparse

# Import local libraries
import mdscrape

# Enable TeX rendering (slower but prettier)
#from matplotlib import rc
#rc('text', usetex=True)

# Deal with I/O
parser=argparse.ArgumentParser()
parser.add_argument("rdf_file")
args=parser.parse_args()

if not os.path.exists("outputs"):
    os.mkdir("outputs")
output_name="outputs/"+os.path.basename(args.rdf_file)

rdf = np.load(args.rdf_file)
r_nm = rdf['r_nm']
g_exp = rdf['g_exp']
g_exp[0]=0
#g_ideal = rdf['g_ideal']
#g_proj = rdf['g_proj']
#sigma_nm = float(rdf['sigma_nm'])
#mu_nm = float(rdf['mu_nm'])
#c = int(rdf['c'])

# Plot all_g's together
fig,ax = plt.subplots()
ax.plot(r_nm, g_exp,'k-',label=r"Experimental")
#ax.plot(r_nm, g_ideal, 'r-',label=r"Ideal Para, $\sigma = {}$ nm".format(sigma_nm))
#ax.plot(r_nm, g_proj, 'y-',label=r"Projection Para, $c = {}$ layers".format(c))
ax.set_xlabel(r"$r$ (nm)")
ax.set_ylabel(r"$g(r)$")
ax.set_title("Radial Distribution Functions")
ax.set_ylim(0,4.5)
ax.legend()
plt.savefig(output_name+".pdf")




