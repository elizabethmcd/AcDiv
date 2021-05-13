library(tidyverse)

####################################  
# Hjor Sulfuritalea 
#################################### 

sulf_path <- "results/SNVs/WWTP/Sulfuritalea"
files <- dir(sulf_path, pattern="*-gene-info.tsv")

sulf_snvs <- data_frame(filename = files) %>% 
  mutate(file_contents = map(filename, ~ read_tsv(file.path(sulf_path, .)))
         ) %>% 
  unnest()  # get all files, append filename for each row, all files that end with *-gene-info.tsv, and unnest to expand

sulf_nucl_div <- sulf_snvs %>% 
  select(filename, gene, nucl_diversity) %>% 
  mutate(sample = gsub("_.*", "", filename)) %>% 
  select(sample, gene, nucl_diversity) %>% 
  pivot_wider(names_from = sample, values_from = nucl_diversity) # get nucleotide diversity for all genes, change sample name

sulf_nucl_div$index <- seq.int(nrow(sulf_nucl_div))
sulf_nucl_div_genome <- sulf_nucl_div %>% 
  select(-gene) %>% 
  drop_na() %>% 
  pivot_longer(!index, names_to = "sample", values_to = "nucleotide_diversity")

sulf_nucl_div_genome %>% ggplot(aes(x=index, y=nucleotide_diversity)) + geom_point(color="navy") + facet_wrap(~ sample)

####################################  
# Accumulibacter with high coverage
#################################### 

acc_path <- "results/SNVs/WWTP/Accumulibacter"
files <- dir(acc_path, pattern="*-gene-info.tsv")

acc_snvs <- data_frame(filename = files) %>% 
  mutate(file_contents = map(filename, ~ read_tsv(file.path(acc_path, .)))
  ) %>% 
  unnest()  # get all files, append filename for each row, all files that end with *-gene-info.tsv, and unnest to expand

acc_nucl_div <- acc_snvs %>% 
  select(filename, gene, nucl_diversity) %>% 
  mutate(sample = gsub("_.*", "", filename)) %>% 
  select(sample, gene, nucl_diversity) %>% 
  pivot_wider(names_from = sample, values_from = nucl_diversity)

acc_nucl_div$index <- seq.int(nrow(acc_nucl_div))
acc_nucl_div_genome <- acc_nucl_div %>% 
  select(-gene) %>% 
  drop_na() %>% 
  pivot_longer(!index, names_to = "sample", values_to = "nucleotide_diversity")

acc_nucl_div_genome %>% ggplot(aes(x=index, y=nucleotide_diversity)) + geom_point(color="navy") + facet_wrap(~ sample)
