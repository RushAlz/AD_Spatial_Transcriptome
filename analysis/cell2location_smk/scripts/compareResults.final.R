library(tidyverse)
library(corrr)

bernard_results = list.files("/pastel/projects/spatial_t/cell2location/results/", pattern = "means_cell_abundance.csv", full.names = T, recursive = T)
bernard_results = bernard_results[grep("tsaiCellsubtype",bernard_results)]
ricardo_results = list.files("/home/ricardo_a_vialle/Workspace/split_by_section/", pattern = "means_cell_abundance.csv", full.names = T, recursive = T)

full_df = data.frame()
for (file in bernard_results){
  #file = bernard_results[1]
  dat = read_csv(file, show_col_types = FALSE) 
  dat2 = dat %>% mutate(barcode = `...1`) %>% column_to_rownames("...1") %>% pivot_longer(-barcode, names_to = "subcelltype", values_to = "proportion")
  dat2$sample_section = gsub("/s","_",gsub("(.*)//(.*)/tsaiCellsubtype/means_cell_abundance.csv","\\2",file))
  dat2$run = "bernard"
  full_df = bind_rows(full_df, dat2)
}

for (file in ricardo_results){
  dat = read_csv(file, show_col_types = FALSE) 
  dat2 = dat %>% mutate(barcode = `...1`) %>% column_to_rownames("...1") %>% pivot_longer(-barcode, names_to = "subcelltype", values_to = "proportion")
  dat2$sample_section = gsub("/","_",gsub("(.*)//(.*)/tsaiCellsubtype/means_cell_abundance.csv","\\2",file))
  dat2$run = "ricardo"
  full_df = bind_rows(full_df, dat2)
}

full_df2 = full_df %>% pivot_wider(values_from = proportion, names_from = run) %>% na.omit()
full_df2$celltype = gsub("(.*)(\\d+)","\\1",gsub("(.*)(\\d+)","\\1",gsub("(.*)_(.*)","\\2",full_df2$subcelltype)))
full_df2_corr = full_df2 %>% group_by(sample_section, subcelltype) %>% summarize(cor = cor(x = bernard, y = ricardo))

ggplot(full_df2_corr, aes(x = subcelltype, y = cor, fill = sample_section)) +
  geom_bar(stat = "identity", position = "dodge") +
  #ylim(c(0.8,1)) +
  facet_wrap(~sample_section, ncol = 1) +
  theme_bw() +
  ggeasy::easy_rotate_x_labels(angle = 45, side = "right")

ggplot(full_df2_corr, aes(x = sample_section, y = cor, fill = sample_section)) +
  geom_boxplot(outlier.colour = NA) +
  ggbeeswarm::geom_quasirandom() +
  ggeasy::easy_rotate_x_labels(angle = 45, side = "right") +
  theme_bw()
