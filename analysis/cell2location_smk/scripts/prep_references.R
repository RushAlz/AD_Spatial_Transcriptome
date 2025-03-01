library(tidyverse)

tsai = read_csv("/pastel/projects/spatial_t/cell2location/references/Tsai_reference/cellsubtypeSignatures.csv")
as.tibble(tsai) %>% head()

length(unique(tsai$Row))
length(unique(gsub("(.*)\\.(.*)","\\1",tsai$Row)))

dat = readRDS("/pastel/projects/spatial_t/cell2location/references/MeningesCellTypes/finecelltype_allgenes.rds")
menin = as.data.frame(dat$SCT) %>% rownames_to_column("Row") %>% as.tibble()
#menin = menin %>% mutate(Row = gsub("(.*)\\.(.*)","\\1",Row)) %>% group_by(Row) %>% slice(n = 1)

length(unique(menin$Row))
length(unique(gsub("(.*)\\.(.*)","\\1",menin$Row)))

merged = menin %>% inner_join(tsai)
head(merged)
write_csv(merged, file = "/pastel/Github_scripts/cell2location_smk/resources/references/MeningesCellTypes/finecelltype_allgenes_plusTsaiSub.csv")

library(corrplot)
corrplot::corrplot(cor(merged[,-1]))
