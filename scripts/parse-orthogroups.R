library(tidyverse)
library(reshape2)
library(UpSetR)

# accumulibacter homologs only
acc_homolog_matrix <- read.delim("results/pangenomics/acc_homolog_matrix.txt", sep="\t")

single_copy_acc <- acc_homolog_matrix %>% filter_at(vars(2:33), all_vars(. == 1))

# all homologs
homolog_matrix <- read.delim("results/pangenomics/homolog_matrix.txt", sep="\t")
selected_homolog_matrix <- homolog_matrix %>% select(-GCA_013823155.1, -GCA_001897745.1)

binary_homolog_matrix <- selected_homolog_matrix %>% filter_at(vars(2:61), all_vars(. <= 1))
colnames(binary_homolog_matrix)[1] <- c("group")          

# melt 
matrix_melted <- melt(binary_homolog_matrix, id.vars="group", variable="genome", value.name= "count")

# genome metadata
metadata <- read.csv("metadata/accumulibacter_refs_genomes_HQ_metadata.csv") %>% select(genome, group_name)
metadata$genome <- gsub("-", ".", metadata$genome)
matrix_melted_names <- left_join(matrix_melted, metadata)

matrix_wide <- matrix_melted_names %>% 
  pivot_wider(names_from = group, values_from = count) %>% 
  select(-genome)

aggregate(matrix_wide[500:600], list(matrix_wide$group_name), mean)



