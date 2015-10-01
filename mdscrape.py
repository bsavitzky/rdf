"""
Scrapes metadata from filepaths for electron microscopy data.
Usage:
    (1)
    Import metadata.py in a python script and call:
    metadata = extract_metadata(filepath)
    Returns a dictionary of metadata
    (2)
    From a shell type
    python metadata.py filepath             # prints to stdout
    python metadata.py filepath > out.txt   # prints to file out.txt

Metadata currently extracted:
    fov
    fov_units
    dwell_time
    dwell_time_units
    pixels

Assumes each item of metadata is separated by "_" delimiters.

Written by:
BHS
Kourkoutis Electron Microscopy Group
Cornell University
150922
"""

from os.path import basename, splitext
from re import split

def get_md_items(filepath, delimiter="_"):
    # Extracts and returns list of metadata from filepath
    base = splitext(basename(filepath))[0]
    md_items = base.split("_")
    md_items.pop(0)
    return md_items

def splitNums(s):
    # Splits string into numerical and non-numerical pieces
    # Returns two lists - nums and alph - the numerical and the alphabetic pieces
    nums = []
    alph = []
    elements = filter(None, split(r'(\d*\.?\d*)', s))
    for element in elements:
        try:
            num = float(element)
            nums.append(num)
        except ValueError:
            alph.append(element)
    return nums, alph

def raiseSortError(data_element):
    print "Item {} could not be sorted.".format(data_element)

def sort_md_item(data_element):
    # Determines the metadata type of a particular data element
    # Accepts a metadata element, as a string
    # Returns a list of keys and a list of their values

    keys = []
    vals = []

    # Split numerical and alphabetical pieces
    nums, alph = splitNums(data_element)

    # "Pixels" element
    # If no alph elements present, see if it's a power of two
    if len(alph)==0:
        if int(nums[0]) in [2**n for n in range(8,15)]:
            keys.append('pixels')
            vals.append(int(nums[0]))
            return keys, vals
        else:
            raiseSortError(data_element)
            return None

    else:
        for alpha in alph:
            # 'fov' element
            if alpha.lower() == 'fov':
                keys.append('fov')
                vals.append(nums[0])
            # 'dwell_time" element
            elif alpha.lower() == 'us':
                keys.append('dwell_time')
                keys.append('dwell_time_units')
                vals.append(nums[0])
                vals.append('us')
            # 'fov_units' element
            elif alpha.lower() in ['nm','um','mm']:
                keys.append('fov_units')
                vals.append(alpha.lower())
            else:
                raiseSortError(data_element)
                return None
        return keys, vals

def extract_metadata(filepath):
    # Extracts metadata from filepath
    # Returns metadata as a dictionary, and prints it to stdout
    keys=[]
    vals=[]
    md_items = get_md_items(filepath)
    for data_element in md_items:
        try:
            ks, vs = sort_md_item(data_element)
            for i in range(len(ks)):
                keys.append(ks[i])
                vals.append(vs[i])
        except TypeError:
            pass

    return dict(zip(keys,vals))

if __name__=="__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('filepath')
    args = parser.parse_args()

    md = extract_metadata(args.filepath)

    for key in md:
        print "{} = {}".format(key, md[key])

