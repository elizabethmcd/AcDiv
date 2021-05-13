library(tidyverse)
library(cowplot)
library(gridGraphics)

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

uw1_snv_count <- uw1_snvs %>% 
  select(filename, gene, SNV_count) %>% 
  mutate(sample = gsub("_UW1.*", "", filename)) %>% 
  select(sample, gene, SNV_count) %>% 
  pivot_wider(names_from = sample, values_from = SNV_count)

uw1_snv_count$index <- seq.int(nrow(uw1_snv_count))
uw1_snvs <- uw1_snv_count %>% 
  select(index, R1R2_2005, R3R4_2015) %>% 
  pivot_longer(!index, names_to = "sample", values_to="SNVs")

uw1_nucl_div$index <- seq.int(nrow(uw1_nucl_div))
uw1_div_comp <- uw1_nucl_div %>% 
  select(index, R1R2_2005, R3R4_2015) %>% 
  pivot_longer(!index, names_to = "sample", values_to = "nucleotide_diversity")

# plots
# 2005 vs 2015 populations snvs
uw1_2005_vs_2015_snvs <- uw1_snvs %>% ggplot(aes(x=index, y=SNVs, color=sample)) + geom_point() + labs(title = "2005 vs 2015 SNVs")
# 2005 vs 2015 populations nucl diversity
uw1_2005_vs_2015_nucl_div <- uw1_div_comp %>% ggplot(aes(x=index, y=nucleotide_diversity, color=sample)) + geom_point() + labs(title = "2005 vs 2015 Nucleotide Diversity")

####################################  
# UW1 Accumulibacter R1R2 2005 vs 2013 
# same reactor "population" over time
#################################### 

r1r2_files <- dir(uw1_path, pattern="R1R2_*")

uw1_r1r2_snvs <- data_frame(filename = r1r2_files) %>% 
  mutate(file_contents = map(filename, ~ read_tsv(file.path(uw1_path, .)))
  ) %>% 
  unnest() 

uw1_r1r2_nucl_div <- uw1_r1r2_snvs %>% 
  select(filename, gene, nucl_diversity) %>% 
  mutate(sample = gsub("_UW1.*", "", filename)) %>% 
  select(sample, gene, nucl_diversity) %>% 
  pivot_wider(names_from = sample, values_from = nucl_diversity)

uw1_r1r2_nucl_div$index <- seq.int(nrow(uw1_r1r2_nucl_div))

uw1_r1r2_nucl_div_compare <- uw1_r1r2_nucl_div %>% 
  select(index, R1R2_2005, R1R2_2013) %>% 
  pivot_longer(!index, names_to = "sample", values_to = "nucleotide_diversity")

uw1_r1r2_snv_count <- uw1_r1r2_snvs %>% 
  select(filename, gene, SNV_count) %>% 
  mutate(sample = gsub("_UW1.*", "", filename)) %>% 
  select(sample, gene, SNV_count) %>% 
  pivot_wider(names_from = sample, values_from = SNV_count)

uw1_r1r2_snv_count$index <- seq.int(nrow(uw1_r1r2_snv_count))

uw1_r1r2_snv_count_comp <- uw1_r1r2_snv_count %>%
  select(index, R1R2_2005, R1R2_2013) %>% 
  pivot_longer(!index, names_to = "sample", values_to = "SNVs")

high_snvs_end <- uw1_r1r2_snv_count_comp %>% 
  filter(SNVs > 30 & index > 4000) %>% 
  pull(index)

high_snv_loci <- uw1_r1r2_snv_count %>% filter(index %in% high_snvs_end)

# plots
# 2005 vs 2013 nucleotide diversity
uw1_2005_vs_2013_nucl_div <- uw1_r1r2_nucl_div_compare %>% ggplot(aes(x=index, y=nucleotide_diversity, color=sample)) + geom_point() + labs(title = "2005 vs 2013 Nucleotide Diversity")
# 2005 vs 2013 SNV counts
uw1_2005_vs_2013_snvs <- uw1_r1r2_snv_count_comp %>% ggplot(aes(x=index, y=SNVs, color=sample)) + geom_point() + labs(title="2005 vs 2013 SNVs")

# arrange in grid to compare 2005 vs 2013 and 2005 vs 2015 (populations over time vs somewhat space in different reactor enrichments)

uw1_grid <- plot_grid(uw1_2005_vs_2013_snvs + theme(legend.position="none"), uw1_2005_vs_2013_nucl_div + theme(legend.position="none"), uw1_2005_vs_2015_snvs + theme(legend.position="none"), uw1_2005_vs_2015_nucl_div + theme(legend.position="none"))

title <- ggdraw() + draw_label("Accumulibacter IIA UW1 Population Dynamics over Time and Space", fontface="bold", x=0, hjust=0) + theme(plot.margin=margin(0,0,0,20))

uw1_grid_pretty <- plot_grid(title, uw1_grid, ncol=1, rel_heights=c(0.1,1))
ggsave(filename="figures/UW1-R1R2-R3R4-2013-2015-pop-comparisons.png", width=12, height=8, units=c("in"))

####################################  
# Annotations
#################################### 

uw1_annotations <- read_tsv("results/annotations/UW1_kofam_annnotations.tsv", col_names = FALSE)
colnames(uw1_annotations) <- c("locus_tag", "ko", "annotation")

high_2015_snvs <- uw1_snvs %>% 
  filter(SNVs > 30 & index > 4000) %>% 
  pull(index)
high_2015_loci <- uw1_snv_count %>% filter(index %in% high_2015_snvs)
high_2015_loci$locus_tag <- gsub("gnl|X|", "", high_2015_loci$gene, fixed=TRUE)
high_2015_snvs_annotations <- left_join(high_2015_loci, uw1_annotations)

high_snvs_annotations <- left_join(high_snv_loci, uw1_annotations) %>% 
  select(locus_tag, R1R2_2005, R1R2_2)
