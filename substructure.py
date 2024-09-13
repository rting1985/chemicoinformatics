"""
Usage:
    postdock.py [options]

Options:
    --input_sdf_file # sdf file
    --input_subs_file # pharm file
    --mode # match model: single/multiple
    --output_sdf_file # sdf file
"""

import sys, os
# sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from rdkit import Chem
from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd
import Bio
from Bio.PDB.PDBParser import PDBParser
import math
#
import argparse
# import logging
# import tqdm
import warnings

warnings.filterwarnings("ignore")


def parse_args():
    """Parses arguments from terminal"""

    parser = argparse.ArgumentParser(description="postdock.py")

    parser.add_argument("--input_sdf_file", "-isdf",
                        help="sdf file", type=str, required=True)

    parser.add_argument("--input_subs_file", "-isubs",
                        help="substructure file", type=str, required=True)

    parser.add_argument("--mode", "-mode",
                        help="match mode single/multi", type=str, required=True)

    parser.add_argument("--output_sdf_file", "-osdf",
                        help="output file", type=str, required=True)

    return parser.parse_args()


def read_text_file(filename):
    "loading pharmacophore file \
        csv/txt file"
    with open(filename, 'r') as f:
        lines = f.readlines()

    return lines

def read_csv_file(filename):
    "loading pharmacophore file \
        csv/txt file"
    lines = pd.read_csv(filename, header=0)  # no header
#    lines = list(lines[0])

    return lines

def get_file_type(path):
    file = os.path.splitext(path)
    filename, type = file

    return type

def read_sdf(filename):
    "laoding sdf file"
    mols = [mol for mol in Chem.SDMolSupplier(filename, removeHs= True)]

    return mols

def write_sdf(filename, mols):
    "write molecules into sdf file"
    w = Chem.SDWriter(filename)
    for m in mols:
#        AllChem.Compute2DCoords(m)
        w.write(m)
    w.close()


# main function
def match(sdffile, subsfile, mode, outputfile):
    #
    type_namefile = get_file_type(subsfile)

    if type_namefile == ".csv":
        subs = read_csv_file(subsfile)
    elif type_namefile == ".txt":
        subs = read_text_file(subsfile)
    else:
        print("Error namefile format!")

    
    print("Loading sdf file")
    sdf = read_sdf(sdffile)
    print(len(sdf))
    mols = []

    for j in range(len(sdf)):
        if mode=="multi":
             nummatch=0
             lensub=len(subs)
             for i in subs:
                  print(i)
                  substr=Chem.MolFromSmiles(i)
                  if sdf[j].HasSubstructMatch(substr):
                       nummatch+=1
             if lensub==nummatch:
                  if sdf[j] not in mols:
                       mols.append(sdf[j])
        elif mode=="single":
             for i in subs:
                  print(i)
                  substr=Chem.MolFromSmiles(i)
                  #print(substr)
                  if sdf[j].HasSubstructMatch(substr):
                       if sdf[j] not in mols:
                             mols.append(sdf[j])
    print(len(mols))
    if len(mols) > 0:                    
        write_sdf(outputfile, mols)
               
def main():
    args = parse_args()

    match(args.input_sdf_file, args.input_subs_file, \
          args.mode, args.output_sdf_file)


if __name__ == "__main__":
    main() 
