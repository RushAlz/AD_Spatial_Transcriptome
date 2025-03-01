# reticulate::use_condaenv("cell2loc_env")
import os
import sys
import anndata
import pandas as pd
import numpy as np
import argparse
#import cell2location
#import scvi

def parse_args():
  parser = argparse.ArgumentParser(description='Make cell specific expression tables')
  parser.add_argument('-i','--anndata_file', help='anndata.h5ad from cell2location')
  parser.add_argument('-o','--output_dir', help='Name of the output directory')
  args = parser.parse_args()
  return args

args = parse_args()
print(args)
print(args.anndata_file)

anndata_obj = anndata.read_h5ad(args.anndata_file) # "/pastel/Github_scripts/cell2location_smk/trash/anndata.h5ad"

# Read anndata layers
for i in anndata_obj.layers.keys():
  print (i)
  cell_mtx = anndata_obj.layers[i].todense().transpose()
  cell_df = pd.DataFrame(cell_mtx)
  cell_df.columns = anndata_obj.obs.index
  cell_df.index = anndata_obj.var.index
  cell_df.to_csv(args.output_dir+"/"+i+"_expmtx.csv.gz", index_label = "gene")


