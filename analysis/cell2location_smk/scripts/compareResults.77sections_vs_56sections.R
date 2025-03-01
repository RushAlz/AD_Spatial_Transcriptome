library(tidyverse)
library(corrr)


results_56sections = list.files("/pastel/projects/spatial_t/cell2location/split_by_section/", pattern = "means_cell_abundance.csv", full.names = T, recursive = T)
results_77sections = list.files("/pastel/Github_scripts/cell2location_smk/results/", pattern = "means_cell_abundance.csv", full.names = T, recursive = T)

full_df = data.frame()
for (file in results_56sections){
  dat = read_csv(file, show_col_types = FALSE) 
  dat2 = dat %>% mutate(barcode = `...1`) %>% column_to_rownames("...1") %>% pivot_longer(-barcode, names_to = "subcelltype", values_to = "proportion")
  dat2$sample_section = gsub("/","_",gsub("(.*)//(.*)/tsaiCellsubtype/means_cell_abundance.csv","\\2",file))
  dat2$run = "56sections"
  full_df = bind_rows(full_df, dat2)
}
for (file in results_77sections){
  dat = read_csv(file, show_col_types = FALSE) 
  dat2 = dat %>% mutate(barcode = `...1`) %>% column_to_rownames("...1") %>% pivot_longer(-barcode, names_to = "subcelltype", values_to = "proportion")
  dat2$sample_section = gsub("/","_",gsub("(.*)//(.*)/means_cell_abundance.csv","\\2",file))
  dat2$run = "77sections"
  full_df = bind_rows(full_df, dat2)
}

full_df2 = full_df %>% pivot_wider(values_from = proportion, names_from = run) %>% na.omit()
full_df2$celltype = gsub("(.*)(\\d+)","\\1",gsub("(.*)(\\d+)","\\1",gsub("(.*)_(.*)","\\2",full_df2$subcelltype)))
full_df2_corr = full_df2 %>% group_by(sample_section, subcelltype) %>% summarize(cor = cor(x = `56sections`, y = `77sections`))

hist(full_df2_corr$cor)
summary(full_df2_corr$cor)

ggplot(full_df2_corr %>% filter(cor < 0.999), aes(x = subcelltype, y = cor, fill = sample_section)) +
  geom_bar(stat = "identity", position = "dodge") +
  #ylim(c(0.8,1)) +
  theme_bw() +
  ggeasy::easy_rotate_x_labels(angle = 45, side = "right")

ggplot(full_df2_corr %>% filter(cor < 0.999), aes(x = sample_section, y = cor, fill = sample_section)) +
  geom_boxplot(outlier.colour = NA, show.legend = F) +
  ggbeeswarm::geom_quasirandom(show.legend = F) +
  ggeasy::easy_rotate_x_labels(angle = 45, side = "right") +
  theme_bw()
