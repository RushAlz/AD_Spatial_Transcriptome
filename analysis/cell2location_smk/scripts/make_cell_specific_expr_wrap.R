library(foreach)
library(doParallel)
registerDoParallel(80)

reticulate::use_condaenv("cell2loc_env")
script_path = "/pastel/Github_scripts/cell2location_smk/scripts/make_cell_specific_expr.py"

files_to_loop = list.files("/pastel/Github_scripts/cell2location_smk/results_77sections_majcelltype", pattern="anndata.h5ad", recursive=T, full.names=T, include.dirs=F)
foreach(i = 1:length(files_to_loop)) %dopar% {
  #file_i = files_to_loop[1]
  file_i = files_to_loop[i]
  folder_i = gsub("(.*)\\/(.*)","\\1",file_i)
  system(paste("CONDA_BASE=/home/ricardo_a_vialle/miniconda3; source $CONDA_BASE/etc/profile.d/conda.sh; conda activate cell2loc_env; python", script_path, "-i", file_i, "-o", folder_i))
}

