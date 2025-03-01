library(tidyverse)

idrun = "results_77sections_majcelltype" # "split_by_77sections_meningesfinecell_plusTsaiSub"
results_folder = paste0("/pastel/Github_scripts/cell2location_smk/results/",idrun) # results_77sections_majcelltype
system(paste0("mkdir -p ", results_folder))

#file_list = system("gsutil ls -r gs://rvialle/projects/ST/cell2location/results/split_by_section_fix", intern = T)
#samples = unique(gsub("gs://rvialle/projects/ST/cell2location/results/split_by_section_fix/(.*?)\\/(.*)","\\1",file_list[grep("means_cell_abundance.csv",file_list)]))

file_list = system(paste0("gsutil ls -r gs://rvialle/projects/ST/cell2location/results/",idrun), intern = T)
samples = unique(gsub(paste0("gs://rvialle/projects/ST/cell2location/results/",idrun,"/(.*?)\\/(.*)"),
                      "\\1",file_list[grep("means_cell_abundance.csv",file_list)]))

for (sample_i in samples){
  #sample_i = samples[1]
  system(paste0("mkdir -p ", results_folder,"/",sample_i))
  
  section_files = file_list[grep(paste0(sample_i,"\\/[1234]\\/:$"), file_list, perl = T)]
  sections = gsub("(.*)\\/(.*)\\/:","\\2",section_files)
  for (section_i in sections){
    #section_i = sections[1]
    to_save_folder = paste0(results_folder,"/",sample_i,"/",section_i)
    system(paste0("mkdir -p ", to_save_folder))  
    
    sample_section_files = file_list[grep(paste0(sample_i,"\\/",section_i,"\\/"), file_list)]
    model_file = sample_section_files[grep("model.pt",sample_section_files)]
    anndata_file = sample_section_files[grep("anndata.h5ad",sample_section_files)]
    cell_abundance_file = sample_section_files[grep("means_cell_abundance.csv",sample_section_files)]
    to_save_folder
    if(file.exists(paste0(to_save_folder,"/means_cell_abundance.csv"))){
      #
    }else{
      print(to_save_folder)
      #system(paste0("gsutil cp ", model_file, " ", to_save_folder))
      #system(paste0("gsutil cp ", anndata_file, " ", to_save_folder))
      system(paste0("gsutil cp ", cell_abundance_file, " ", to_save_folder))
    }
  }
}

cell_abundance_files = list.files(results_folder, pattern = "means_cell_abundance.csv", recursive = T, full.names = T)
cell_abundance_df = purrr::map_df(cell_abundance_files, ~{read.csv(.x)})
cell_abundance_df_long = cell_abundance_df %>% pivot_longer(-X, names_to = "celltype", values_to = "estimate") %>% 
  rename(barcode = "X") 
data.table::fwrite(cell_abundance_df_long, file = paste0("/pastel/projects/spatial_t/cell2location/",idrun,".tsv.gz"), sep = "\t", row.names = F, col.names = T)


