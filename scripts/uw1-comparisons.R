library(tidyverse)
library(cowplot)
library(gridGraphics)

####################################  
# UW1 Accumulibacter R1R2 2005 vs R3R4 2015
#################################### 

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

# Facet by sample instead of overlaid
uw1_snvs_samples <- uw1_snvs %>% ggplot(aes(x=index, y=SNVs)) + geom_point(color="navyblue") + facet_wrap(~ sample) + labs(title = "UW1 IIA SNVs in R1R2 vs R3R4")  + theme_bw()
uw1_pi_samples <- uw1_div_comp %>% ggplot(aes(x=index, y=nucleotide_diversity)) + geom_point(color="purple") + facet_wrap(~ sample) + labs(title = "UW1 IIA Gene-Specific Nucleotide Diversity π in R1R2 vs R3R4") + theme_bw()

plot_grid(uw1_snvs_samples, uw1_pi_samples, labels="AUTO", ncol=1)
uw1_grid_by_sample

r1r2_grid <- plot_grid(uw1_2005_vs_2013_snvs + theme(legend.position="none"), uw1_2005_vs_2013_nucl_div + theme(legend.position="none"), uw3_snvs_plot + theme(legend.position="none"), uw3_nucl_plot + theme(legend.position="none"))

r1r2_title <- ggdraw() + draw_label("Accumulibacter UW1 IIA and UW3 IA SNVs and π in 2005 vs 2013", fontface="bold", x=0, hjust=0) + theme(plot.margin=margin(0,0,0,20))

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
  filter(SNVs > 30) %>% 
  pull(index)

high_snv_loci <- uw1_r1r2_snv_count %>% filter(index %in% high_snvs_end)

# plots
# 2005 vs 2013 nucleotide diversity
uw1_2005_vs_2013_nucl_div <- uw1_r1r2_nucl_div_compare %>% ggplot(aes(x=index, y=nucleotide_diversity, color=sample)) + geom_point() + labs(title = "UW1 IIA 2005 vs 2013 π")
# 2005 vs 2013 SNV counts
uw1_2005_vs_2013_snvs <- uw1_r1r2_snv_count_comp %>% ggplot(aes(x=index, y=SNVs, color=sample)) + geom_point() + labs(title="UW1 IIA 2005 vs 2013 SNVs")

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

high_snv_loci$locus_tag <- gsub("gnl|X|", "", high_snv_loci$gene, fixed=TRUE)
high_snvs_annotations <- left_join(high_snv_loci, uw1_annotations)

####################################  
# UW3 Accumulibacter R1R2 2005 vs 2013
# Sufficient coverage in both samples to compare over time
#################################### 

uw3_path <- "results/SNVs/R1R2_vs_R3R4/"
uw3_files <- dir(uw2_path, pattern="*_UW3_gene_info.tsv")

uw3_snvs <- data_frame(filename = uw3_files) %>% 
  mutate(file_contents = map(filename, ~ read_tsv(file.path(uw3_path, .)))
         ) %>% 
  unnest()

uw3_nucl_div <- uw3_snvs %>% 
  select(filename, gene, nucl_diversity) %>% 
  mutate(sample = gsub("_UW3.*", "", filename)) %>% 
  select(sample, gene, nucl_diversity) %>% 
  pivot_wider(names_from = sample, values_from = nucl_diversity)

uw3_snv_count <- uw3_snvs %>% 
  select(filename, gene, SNV_count) %>% 
  mutate(sample = gsub("_UW3.*", "", filename)) %>% 
  select(sample, gene, SNV_count) %>% 
  pivot_wider(names_from = sample, values_from = SNV_count)

uw3_snv_count$index <- seq.int(nrow(uw3_snv_count))

uw3_snvs_filtered <- uw3_snv_count %>% 
  select(index, R1R2_2005, R1R2_2013) %>% 
  drop_na() %>% 
  pivot_longer(!index, names_to = "sample", values_to="SNVs")

uw3_nucl_div$index <- seq.int(nrow(uw3_nucl_div))
uw3_div_comp <- uw3_nucl_div %>% 
  select(index, R1R2_2005, R1R2_2013) %>% 
  drop_na() %>% 
  pivot_longer(!index, names_to = "sample", values_to = "nucleotide_diversity")

# plots
# 2005 vs 2013 UW3 populations snvs
uw3_snvs_plot <- uw3_snvs_filtered %>% ggplot(aes(x=index, y=SNVs, color=sample)) + geom_point() + labs(title = "UW3 IA 2005 vs 2013 SNVs")
uw3_snvs_plot
# 2005 vs 2015 populations nucl diversity
uw3_nucl_plot <- uw3_div_comp %>% ggplot(aes(x=index, y=nucleotide_diversity, color=sample)) + geom_point() + labs(title = "UW3 IA 2005 vs 2013 π")
uw3_nucl_plot

# grid of UW1 and UW3 SNVS and π over same timeframe

r1r2_grid <- plot_grid(uw1_2005_vs_2013_snvs + theme(legend.position="none"), uw1_2005_vs_2013_nucl_div + theme(legend.position="none"), uw3_snvs_plot + theme(legend.position="none"), uw3_nucl_plot + theme(legend.position="none"))

r1r2_title <- ggdraw() + draw_label("Accumulibacter UW1 IIA and UW3 IA SNVs and π in 2005 vs 2013", fontface="bold", x=0, hjust=0) + theme(plot.margin=margin(0,0,0,20))

r1r2_grid_pretty <- plot_grid(r1r2_title, r1r2_grid, ncol=1, rel_heights=c(0.1,1))
r1r2_grid_pretty

ggsave(filename="figures/R1R2-UW1-UW3-2005-vs-2013-comparisons.png", r1r2_grid_pretty, width=12, height=8, units=c("in"))

####################################  
# UW1 and UW3 in 2008 time-point when both equally abundant
#################################### 

uw1_2008_snvs <- read_tsv("results/SNVs/R1R2_vs_R3R4/R1R2_2008_UW1_gene_info.tsv") %>% 
  select(gene, SNV_count)
uw1_2008_snvs$index <- seq.int(nrow(uw1_2008_snvs))
uw1_2008_snv_count <- uw1_2008_snvs %>% 
  select(index, SNV_count) %>% 
  drop_na()

uw1_2008_div <- read_tsv("results/SNVs/R1R2_vs_R3R4/R1R2_2008_UW1_gene_info.tsv") %>% 
  select(gene, nucl_diversity)
uw1_2008_div$index <- seq.int(nrow(uw1_2008_div))
uw1_2008_div_count <- uw1_2008_div %>% 
  select(index, nucl_diversity) %>% 
  drop_na()

uw3_2008_snvs <- read_tsv("results/SNVs/R1R2_vs_R3R4/R1R2_2008_UW3_gene_info.tsv") %>% 
  select(gene, SNV_count)
uw3_2008_snvs$index <- seq.int(nrow(uw3_2008_snvs))

uw3_2008_nucl_div <- read_tsv("results/SNVs/R1R2_vs_R3R4/R1R2_2008_UW3_gene_info.tsv") %>% 
  select(gene, nucl_diversity)
uw3_2008_nucl_div$index <- seq.int(nrow(uw3_2008_nucl_div))

# plots
uw1_snv_plot <- uw1_2008_snv_count %>% ggplot(aes(x=index, y=SNV_count)) + geom_point(color="navyblue") + labs(title="UW1 SNVs")
uw1_pi_plot <- uw1_2008_div_count %>% ggplot(aes(x=index, y=nucl_diversity)) + geom_point(color="navyblue") + labs(title="UW1 π")
uw3_snv_plot <- uw3_2008_snvs %>% ggplot(aes(x=index, y=SNV_count)) + geom_point(color="orange") + labs(title="UW3 SNVs")
uw3_pi_plot <- uw3_2008_nucl_div %>% ggplot(aes(x=index, y=nucl_diversity)) + geom_point(color="orange") + labs(title="UW3 π")

r1r2_2008_grid <- plot_grid(uw1_snv_plot, uw1_pi_plot, uw3_snv_plot, uw3_pi_plot)
r1r2_2008_grid

r1r2_2008_title <- ggdraw() + draw_label ("Accumulibacter UW1 IIA and UW3 IA SNVs and π in 2008", fontface="bold", x=0, hjust=0) + theme(plot.margin=margin(0,0,0,20))

r1r2_2008_grid_pretty <- plot_grid(r1r2_2008_title, r1r2_2008_grid, ncol=1, rel_heights=c(0.1,1))
r1r2_2008_grid_pretty

####################################  
# UW1 Accumulibacter R1R2 2005 vs 2015
# Grid of SNVs, π, and Fst between these samples for Acc IIA
#################################### 

# Facet by sample instead of overlaid
supp_labels = c(`R1R2_2005`="R1R2 2005", `R3R4_2015`="R3R4 2015")
uw1_snvs_samples <- uw1_snvs %>% ggplot(aes(x=index, y=SNVs)) + geom_point(color="navyblue") + facet_wrap(~ sample, labeller=as_labeller(supp_labels)) + labs(title = "UW1 IIA SNVs in R1R2 vs R3R4") + ylab("Number of SNVs") + theme_bw() + theme(axis.title.x=element_blank(), axis.title.y=element_text(face="bold"))
uw1_snvs_samples
uw1_pi_samples <- uw1_div_comp %>% ggplot(aes(x=index, y=nucleotide_diversity)) + geom_point(color="darkorchid4") + facet_wrap(~ sample, labeller=as_labeller(supp_labels)) + scale_y_continuous(limits=c(0,0.025)) + ylab("Nucleotide Diversity") + labs(title = "UW1 IIA Gene-Specific Nucleotide Diversity in R1R2 vs R3R4") + theme_bw() + theme(axis.title.x=element_blank(), axis.title.y=element_text(face="bold"))
uw1_pi_samples
#### fst plot in the script compare-sample-SNVs.R and grid with the SNVs and π between samples

uw1_2005_2015_grid <- plot_grid(uw1_snvs_samples, uw1_pi_samples, uw1_2005_2015_fst, labels="AUTO", ncol=1)
uw1_2005_2015_grid

ggsave(filename="figures/UW1-IIA-R1R2-R3R4-microdiv-grid.png", uw1_2005_2015_grid, width=15, height=9, units=c("in"))

