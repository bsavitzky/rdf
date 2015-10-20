#! /bin/bash

# Once centroids_edited has been obtained via user input, use run.py
#
# Operation:
# Performs subpixel refinement (SP_centroids.py)
# Calculates experimental radial distribution function (rdf_exp.py)
# Plots g_exp(r) (plt_rdf.py)

source import_vars.sh
python $DIR/SP_centroids.py $I $Ce
python $DIR/rdf_exp.py $I $Csp 0.1
python $DIR/plot_rdf.py $Gexp
