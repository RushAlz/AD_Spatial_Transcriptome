# Hi Ricardo,
# 
# I looked through the cell2location folder, and seems like we are missing results for quite a few sections. This is solely due to missing ST files in the folder that Artemis provided you. Below is the list of missing sections. You can find ST data in gs://bernardn/rosmapAD/data/spatialTranscriptomic/expressionData/unnormalized/ 
# The files are named "raw_STdata_B*.csv". You will have to split the files by sections based on the spotID, i.e. A1 = s1, B1 = s2, C1 = s3, and D1=s4. Also, you will have to remove genes with MT in the name, which is what Artemis did to remove MT genes for running cell2location.
# B16002_s3
# B16047_s4
# B16080_s4
# B16095_s2
# B17047_s3
# B17047_s4
# B17085_s1
# B17085_s3
# B17093_s1
# B18022_s4
# 
# Thanks,
# Bernard

# Hi Ricardo,
# 
# I finished running VistoSeg on the remaining sections and am currently uploading spot_count csv files to the same destination.
# 
# In addition to B16002_3, there is one more section that failed vistoseg (B16100_2), at the same step during the ‘countNuclei_new.m’ script.
# 
# I’ll upload all the relevant input files/scripts for the two failed sections here in case Bernard has time to troubleshoot: radc_sc_spatial_transcriptome_projects/cell2location_input/VistoSeg_reruns
# 
# Also, please let me know if any sections’ nuclei counts are well outside of the expected range of ~5-10 per spot. While no other sections failed outright, there may be a few that contain overestimates, and may require further tweaking of the mask.
# 
# Best,
# Denis

###
# Files were copied to /pastel/projects/spatial_t/cell2location/data/20221006_expressionData/unnormalized
###

smk_folder = "/pastel/Github_scripts/cell2location_smk"
main_folder = "/pastel/projects/spatial_t/cell2location/data/20221006_expressionData"
data_folder = "/pastel/projects/spatial_t/cell2location/data/20221006_expressionData/unnormalized"
split_by_section_folder = "/pastel/projects/spatial_t/cell2location/data/20221006_expressionData/split_by_section"
nuclei_quant_folder = "/pastel/projects/spatial_t/cell2location/data/20221006_expressionData/NucleiQuantification"

# nuclei_files = list.files(nuclei_quant_folder, pattern = "*.csv")
# stsect_files = list.files(split_by_section_folder, pattern = "*.csv")
# gsub("(.*?)_(.*?)_(.*)","\\1_\\2",stsect_files)[! gsub("(.*?)_(.*?)_(.*)","\\1_\\2",stsect_files) %in% gsub("(.*?)_(.*?)_(.*)","\\1_\\2",nuclei_files)]
# Missing:: "B16002_3" "B16100_2"

################################################################################
library(tidyverse)
library(data.table)

### Read Nuclei counts tables
df = data.frame(sample=character(), section=character(), avg_nuclei=numeric())
nuclei_files = list.files(nuclei_quant_folder, pattern = "*.csv", full.names = T)
for(file in nuclei_files){
  #file = nuclei_files[1]
  sample = gsub("(.*)/(.*?)_(.*?)_tissue_spot_counts.csv","\\2",file)
  section = gsub("(.*)/(.*?)_(.*?)_tissue_spot_counts.csv","\\3",file)
  
  mat = read.csv(file)
  avg_nuclei = mean(mat[mat$tissue==1,]$Nmask_dark_blue)
  
  df = bind_rows(df, list(sample=sample,section=section,avg_nuclei=avg_nuclei))
}
# 76 files should be there
# The expected avg nuclei counts should be around 5-10.
df_filt = df[df$avg_nuclei < 15,] # 67 sections with avg nuclei < 15
write_tsv(df_filt, file = paste0(smk_folder, "/avgNuclei_perSample.txt"))
####

### Write ST counts tables by section
metadata_st = read.csv(paste0(main_folder,"/metadata.csv"))
rownames(metadata_st) = gsub("B16002_B18007","B18007_B18007",rownames(metadata_st))
section_map = c(A1 = 1, B1 = 2, C1 = 3, D1 = 4) 
exp_files = list.files(data_folder, pattern = "raw_*", full.names = T)
for(file in exp_files){
  #file = exp_files[1]
  sample = gsub("(.*)/raw_STdata_(.*?).csv","\\2",file)
  mat = fread(file, sep = ",", header = T) %>% column_to_rownames("gene")
  # match with spots in the metadata
  meta_bnum = metadata_st[metadata_st$Bnum == sample,]
  meta_bnum$barcode = rownames(meta_bnum)
  meta_bnum$section = meta_bnum$sample
  meta_bnum$section[meta_bnum$section=="A1"] <- 1
  meta_bnum$section[meta_bnum$section=="B1"] <- 2
  meta_bnum$section[meta_bnum$section=="C1"] <- 3
  meta_bnum$section[meta_bnum$section=="D1"] <- 4
  meta_bnum$section = as.numeric(meta_bnum$section)
  # remove MT genes
  mat_noMT = mat[!grepl("^MT-", toupper(rownames(mat))),] # 13 genes removed
  # split by section
  section_label = gsub("(.*?)_(.*?)_(.*?)_(.*)","\\3",colnames(mat_noMT))
  sections = unique(section_label)
  for (section in sections){
    #section = sections[1]
    if (paste0(sample,"_",section_map[section]) %in% paste0(df_filt$sample,"_",df_filt$section)){
      meta_sect = meta_bnum[meta_bnum$sample == section,]
      col_select = meta_sect$barcode
      col_select = c("gene", col_select)
      mat_noMT_sect = mat_noMT[,grepl(section,section_label)] %>% rownames_to_column("gene")
      #mat_noMT_sect = mat_noMT_sect[,..col_select]
      print(paste0(sample,":",section," ", identical(colnames(mat_noMT_sect),col_select))) # TRUE
      file2write = paste0(split_by_section_folder,"/",sample,"_",section_map[section],"_noMT_spatial_counts_matrix.csv")
      fwrite(mat_noMT_sect, file = file2write, sep = ",", row.names = F, col.names = T)
    }
  }
}
# 67 of 78 files should be there
####


