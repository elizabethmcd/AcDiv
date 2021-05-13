library(tidyverse)

####################################  
# UW1 Accumulibacter R1R2 vs R3R4
#################################### 

uw1_path <- "results/SNVs/R1R2_vs_R3R4/"
files <- dir(uw1_path, pattern="*_gene_info.tsv")

uw1_snvs <- data_frame(filename = files) %>% 
  mutate(file_contents = map(filename, ~ read_tsv(file.path(uw1_path, .)))
  ) %>% 
  unnest() 

uw1_nucl_div <- uw1_snvs %>% 
  select(filename, gene, nucl_diversity) %>% 
  mutate(sample = gsub("_UW1.*", "", filename)) %>% 
  select(sample, gene, nucl_diversity) %>% 
  pivot_wider(names_from = sample, values_from = nucl_diversity)
