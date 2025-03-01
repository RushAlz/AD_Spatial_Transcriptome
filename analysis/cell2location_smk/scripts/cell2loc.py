import os
import sys
import anndata
import pandas as pd
import numpy as np
import cell2location
import scvi
import argparse
from glob import glob

def snakemake_args():
  params = vars(snakemake.params)
  params['input_file'] = snakemake.input[0]
  params['nuclei_file'] = snakemake.input[1]
  params['reference_file'] = snakemake.input[2]
  params['output_dir'] = snakemake.output[0] # output_path 
  params['model'] = snakemake.output[1] # output_path + "{sample}/{section}/" + refFolder + "/model.pt
  params['anndata'] = snakemake.output[2] # output_path + "{sample}/{section}/" + refFolder + "/anndata.h5ad"
  params['means_cell_abundance'] = snakemake.output[3] # output_dir + 'means_cell_abundance.csv'
  del params['_names']
  return params

def parse_args():
  parser = argparse.ArgumentParser(description='Run cell2localtion')
  parser.add_argument('input_file', help='Input file (noMT_spatial_counts_matrix.csv)')
  parser.add_argument('reference_file', help='reference_file')
  parser.add_argument('nuclei_file', help='avg nNuclei file')
  parser.add_argument('output_dir', help='Name of the output directory')
  parser.add_argument('--sample', help='sample_id')
  parser.add_argument('--section', help='section_number')
  args = parser.parse_args()
  return args

if 'snakemake' in globals():
  args = snakemake_args()
else:
  args = parse_args()

# Parameters
input_file = args['input_file']
reference_file = args['reference_file']
nuclei_file = args['nuclei_file']
sample_id = args['sample']
section_number = args['section']

# input_file = "/pastel/projects/spatial_t/cell2location/data/cell2location_input/split_by_section/B16002_1_noMT_spatial_counts_matrix.csv"
# output_dir = "/pastel/projects/spatial_t/cell2location/results_debug/B16002/1/tsaiCellsubtype/means_cell_abundance.csv"
# reference_file = "/pastel/projects/spatial_t/cell2location/references/Tsai_reference/cellsubtypeSignatures.csv"
# nuclei_file = "/pastel/Github_scripts/cell2location_smk/avgNuclei_perSample.txt"
# sample_id = "B16002"
# section_number = "1"

df = pd.read_table(nuclei_file, header = 0)
df = df[df['sample'] == sample_id] 
df = df[df['section'] == float(section_number)]
nNuclei = float(df.avg_nuclei)
# nNuclei = 3.70

print("Sample: " + sample_id)
print("Section: " + section_number)
print("Avg. nuclei: " + str(nNuclei))

ref = pd.read_csv(reference_file,index_col=0)
test = anndata.read_csv(input_file, first_column_names=True, dtype='float64')
intersect = np.intersect1d(test.obs_names,ref.index)
ref = ref.loc[intersect,:]
test_0 = test[intersect,:].copy()
test_0 = test_0.copy().T

# Prepare test for cell2location model
cell2location.models.Cell2location.setup_anndata(adata=test_0)

# Set up cell2location
mod = cell2location.models.Cell2location(
    test_0,cell_state_df=ref,
    # the expected average cell abundance: tissue-dependent
    # hyper-prior which can be estimated from paired histology:
    N_cells_per_location=nNuclei,
    # hyperparameter controlling normalisation of
    # within-experiment variation in RNA detection (using default here):
    detection_alpha=200
)

mod.train(max_epochs=30000, # 30000
          # train using full data (batch_size=None)
          batch_size=None,
          # use all data points in training because
          # we need to estimate cell abundance at all locations
          train_size=1,
          use_gpu=True)

test = mod.export_posterior(
    test_0, sample_kwargs={'num_samples': 1000, # 1000
    'batch_size': mod.adata.n_obs, 'use_gpu': True}
)

# Save model
mod.save(args['output_dir'], overwrite=True, save_anndata=True)

# Compute expected expression per cell type
expected_dict = mod.module.model.compute_expected_per_cell_type(
    mod.samples["post_sample_q05"], mod.adata_manager
)

# Add to anndata layers
for i, n in enumerate(mod.factor_names_):
    test.layers[n] = expected_dict['mu'][i]

# Save anndata object with results
test.write(args['anndata'])

# Write out mean cell abundance
df = pd.DataFrame(test.obsm["means_cell_abundance_w_sf"])
df.to_csv(args['means_cell_abundance'])

