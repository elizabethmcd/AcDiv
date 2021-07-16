library(tidyverse)
library(cowplot)

##############################
# Recombination rates across genomes
##############################

# UW bioreactors R1R2 and R3R4 
uw_path <- "results/SNVs/all_genome_info"
uw_files <- dir(uw_path, pattern="_genome_info.tsv")
uw_recomb <- data_frame(filename = uw_files) %>% 
  mutate(file_contents = map(filename, ~ read_tsv(file.path(uw_path, .)))) %>% 
  unnest() %>% 
  select(genome, coverage, r2_mean, d_prime_mean) %>% 
  filter(coverage > 10)

uw_recomb$genome <- gsub(".fasta", "", uw_recomb$genome)  

# POB timepoints 
pob_files <- dir(uw_path, pattern="_genome_stats.tsv")

pob_recomb <- data_frame(filename = pob_files) %>% 
  mutate(file_contents = map(filename, ~ read_tsv(file.path(uw_path, .)))) %>% 
  unnest() %>% 
  filter(coverage > 10) %>% 
  select(filename, coverage, genome, r2_mean, d_prime_mean)

pob_recomb$filename <- gsub("_genome_stats.tsv", "", pob_recomb$filename)
pob_recomb$genome <- gsub(".fa", "", pob_recomb$genome)

pob_metadata <- read.csv("metadata/POS-MAGs-table.csv") %>% 
  select(bin, code) %>% 
  mutate(genome = bin) %>% 
  select(genome, code)

pob_recomb_rhodocyc <- left_join(pob_recomb, pob_metadata) %>% 
  filter(code == 'CAPIA' | code == 'CAPIIB' | code == 'CAPIIC' | code == 'DECH1' | code == 'SULF1' | code == 'THAU1')

pob_recomb %>% ggplot(aes(x=d_prime_mean, y=r2_mean)) + geom_point() + facet_wrap(~ filename, ncol=1)

cycle.labs <- c("Cycle 87", "Cycle 103", "Cycle 129")
names(cycle.labs) <- c("POS-2015-07-16", "POS-2015-07-24", "POS-2015-08-06")

rhodo_recomb_plot <- pob_recomb_rhodocyc %>% ggplot(aes(x=d_prime_mean, y=r2_mean, color=code)) + geom_point(size=3) + facet_wrap(~ filename, ncol=1, labeller=labeller(filename=cycle.labs)) + xlab("D'") + ylab(expression(r^2)) + labs(color="Genome") + scale_color_brewer(palette="Paired") + theme_bw()
rhodo_recomb_plot

rhodo_div_noLeg <- rhodo_div_plot + theme(legend.position="none") + scale_color_brewer(palette="Paired")

pob_grid <- plot_grid(rhodo_div_noLeg, rhodo_recomb_plot, labels="AUTO")
pob_grid

ggsave("figures/POB-div-recomb-grid.png", pob_grid, width=9, height=7, units=c("in"))


# WWTP files 
wwtp_files <- dir(uw_path, pattern="-genome-stats.tsv")
wwtp_div <- data_frame(filename = wwtp_files) %>% 
  mutate(file_contents = map(filename, ~ read_tsv(file.path(uw_path, .), col_names = FALSE))) %>% 
  unnest() %>% 
  select(filename, X1, X2, X3, X4, X20, X21)

colnames(wwtp_div) <- c("sample", "genome", "coverage", "breadth", "nucl_diversity", "r2_mean", "d_prime_mean")

wwtp_div_filtered <- wwtp_div %>% 
  filter(coverage > 10) %>% 
  filter(breadth > 0.9)

wwtp_div_filtered$sample <- gsub("_18-genome-stats.tsv", "", wwtp_div_filtered$sample)
wwtp_div_filtered$genome <- gsub(".fa", "", wwtp_div_filtered$genome)

wwtp_rhodo_codes <- read.csv("metadata/danish-rhodo-codes.csv")

wwtp_diversity <- left_join(wwtp_div_filtered, wwtp_rhodo_codes) %>% 
  select(sample, genome, code, species, coverage, breadth, r2_mean, d_prime_mean, nucl_diversity) %>% drop_na()

plants_ordered <- c("Hirt", "Hjor", "AalE", "AalW", "Mari", "Skiv", "Vibo", "Rand", "Ega", "Viby", "EsbE", "EsbW", "Fred", "Ribe", "Hade", "OdNW", "Ejby", "OdNE", "Lyne", "Damh", "Aved")

wwtp_diversity$sample <- factor(wwtp_diversity$sample, levels=plants_ordered)

wwtp_recomb_plot <- wwtp_diversity %>% ggplot(aes(x=d_prime_mean, y=r2_mean, color=code)) + geom_point() + facet_wrap(~ species, ncol=1) + xlab("D'") + ylab(expression(r^2)) + labs(color="Genome") + scale_color_brewer(palette="Set3") + scale_x_continuous(breaks=seq(0.92, 0.99, .01)) + scale_y_continuous(breaks=seq(0.2, 0.8, 0.2)) + theme_bw() + theme(strip.text=element_text(face="italic"))

wwtp_div_sample_plot <- wwtp_diversity %>% ggplot(aes(x=sample, y=nucl_diversity, color=code)) + scale_color_brewer(palette="Set3") + geom_point(size=2) + facet_wrap(~ species, nrow=2) + theme_bw() + theme(axis.text.x= element_text(angle=70, hjust=1)) + ylab("Nucleotide Diversity") + xlab("Wastewater Treatment Plant") + scale_y_continuous(breaks=seq(0, 0.025, .005)) + theme(strip.text=element_text(face="italic"), axis.title.x=element_text(face="bold"), axis.title.y=element_text(face="bold"), legend.position="none")

wwtp_grid <- plot_grid(wwtp_div_sample_plot, wwtp_recomb_plot, labels="AUTO", rel_widths = c(0.6, 0.4))

wwtp_grid

ggsave("figures/WWTP-Rhodocyc-div-recomb-grid.png", wwtp_grid, width=11, height=7, units=c("in"))
