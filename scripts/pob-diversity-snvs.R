library(tidyverse)

####################################  
# POB Population Dynamics and SNVs
#################################### 

# get all files for covg, breadth, and within-sample nucleotide diversity
pob_path <- "results/covg_breadth/"
pob_files <- dir(pob_path, pattern="*-div-stats.txt")
pob_data <- data_frame(filename = pob_files) %>% 
  mutate(file_contents = map(filename, ~ read_tsv(file.path(pob_path, .), col_names = FALSE))
  ) %>% 
  unnest() 
colnames(pob_data) <- c("filename", "genome", "coverage", "breadth", "nucleotide_diversity")
pob_data$sample <- gsub("-covg-breadth-div-stats.txt", "", pob_data$filename)

# genomes that meet the criteria
pob_covg_breadth_filtered <- pob_data %>% 
  select(sample, genome, coverage, breadth) %>%
  filter(coverage > 10 & breadth > 0.9) %>% 
  select(sample, genome, coverage) %>% 
  pivot_wider(names_from = sample, values_from = coverage) %>% 
  drop_na() %>% 
  pull(genome)

# get those genomes in a new df
pob_good_genomes <- pob_data %>% 
  filter(genome %in% pob_covg_breadth_filtered) %>% 
  select(genome, sample, coverage, breadth, nucleotide_diversity)
pob_good_genomes$genome <- gsub(".fa", "", pob_good_genomes$genome)

# combine with metadata
pob_metadata <- read.csv("metadata/POS-MAGs-table.csv") %>% select(bin, code, classification)
colnames(pob_metadata)[1] <- c("genome")

pob_good_info <- left_join(pob_good_genomes, pob_metadata) %>% 
  select(genome, code, sample, coverage, breadth, nucleotide_diversity)

# plot of all that meet covg/breadth thresholds
pob_good_info %>% ggplot(aes(x=sample, y=nucleotide_diversity, group=code, color=code)) + geom_point() + geom_line()

# rhodocyclaceae results
rhodo_div <- pob_good_info %>% 
  filter(code == 'CAPIA' | code == 'CAPIIB' | code == 'CAPIIC' | code == 'DECH1' | code == 'SULF1' | code == 'THAU1')

rhodo_div_plot <- rhodo_div %>% ggplot(aes(x=sample, y=nucleotide_diversity, group=code, color=code)) + geom_point(size=2.5) + geom_line(size=1.5) + theme_bw() + ylab("Nucleotide Diversity\n") + xlab("\nCycle") + scale_x_discrete(labels=c("Cycle 87", "Cycle 103", "Cycle 129")) + labs(color="Genome") + theme(axis.title.y=element_text(face="bold"), axis.title.x=element_text(face="bold"), legend.title=element_text(face="bold"))
rhodo_div_plot

ggsave(filename="figures/POB-rhodo-div-plot.png", rhodo_div_plot, width=8, height=6, units=c("in"))
