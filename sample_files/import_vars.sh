#! /bin/bash

echo
echo "Loading variables for rdf analysis suite of image 55 from 150831:"

DIR="/Users/bsavitzky/Projects/PbSe_QDs/Analysis/150914_Thickness_LRO/rdf"
I="/Users/bsavitzky/Data/PbSe_QDS/150831_CBCs_PbSe_lowThickness_noUltrathin/TiffsFromDm3Files/55_fov600nm_8us_4096.tif"
OUT="/Users/bsavitzky/Projects/PbSe_QDs/Analysis/150914_Thickness_LRO/150831_3-4layer/3layers/150831_55/outputs/"
C=$OUT"centroids.npz"
Ce=$OUT"centroids_edited.npz"
Csp=$OUT"SP_centroids.npz"
Gexp=$OUT"rdf_exp.npz"


echo "DIR = rdf script location"
echo "I = Image location"
echo "C = Raw centroid location"
echo "Ce = Edited centroid location"
echo "Csp = Subpixel centroid location"
echo "Gexp = Experimental radial distribution function"

