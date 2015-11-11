# Import libraries
import numpy as np
import matplotlib.pyplot as plt
import os
import argparse

# Enable TeX rendering (slower but prettier)
#from matplotlib import rc
#rc('text', usetex=True)

# Deal with I/O
parser=argparse.ArgumentParser()
parser.add_argument("rdf_file")
parser.add_argument("-o","--output",default="outputs/",help="Directory for output files.")
args=parser.parse_args()

if not os.path.exists(args.output):
    os.mkdir(args.output)
if args.output[-1]!="/":
    output_name=args.output+"/"
else:
    output_name=args.output
output_name+=os.path.splitext(os.path.basename(args.rdf_file))[0]

rdf = np.load(args.rdf_file)
r_nm = rdf['r_nm']
g_exp = rdf['g_exp']
g_exp[0]=0

# Plot
fig,ax = plt.subplots()
ax.plot(r_nm, g_exp,'k-',label=r"Experimental")
ax.set_xlabel(r"$r$ (nm)")
ax.set_ylabel(r"$g(r)$")
ax.set_title("Radial Distribution Functions")
ax.set_ylim(0,4.5)
ax.legend()
plt.savefig(output_name+".pdf")




