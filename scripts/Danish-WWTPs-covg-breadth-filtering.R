library(tidyverse)
library(reshape2)

# Merging coverage information of HQ SPREP WWTP genomes with metadata, getting high coverage, and inspecting coverage of Rhodocyclaceae genomes

# mag metadata 
wwtps_mags <- read.csv("metadata/Danish_FullScale_WWTPs/Singleton2021-MAG-summaries.csv")
colnames(wwtps_mags)[1] <- c('genome')
mag_taxonomy <- wwtps_mags %>% select(genome, GTDBTax)

# coverage table
coverage <- read.csv("results/relative_abundance/Danish-WWTPs-curated-coverage.csv")

# merge coverage and taxonomy
coverage_info <- left_join(coverage, mag_taxonomy)
filtered_5X_average <- coverage_info %>% filter(rowMeans(coverage_info[, 2:23]) > 5)
filtered_10X_average <- coverage_info %>% filter(rowMeans(coverage_info[, 2:23]) > 10)

tax10 <- filtered_10X_average %>% select(GTDBTax, genome)

# all that have high average 10X coverage across all samples
filtered_10X_names_covg <- filtered_10X_average[,c(1,24,2:23)]
filtered_10X_names_covg$avg_covg <- rowMeans(filtered_10X_names_covg[, 3:24])
filtered_10X_genomes <- filtered_10X_names_covg[,1]

# just rhodocyclaceae genomes

wwtp_rhodocyc <- c("Fred_18-Q3-R57-64_MAXAC.027.fa",
                   "Hade_18-Q3-R52-61_BATAC.726.fa",
                   "OdNE_18-Q3-R46-58_BAT3C.415.fa",
                   "EsbW_18-Q3-R4-48_BATAC.285.fa",
                   "Fred_18-Q3-R57-64_BAT3C.720.fa",
                   "Aved_18-Q3-R54-62_MAXAC.403.fa",
                   "Lyne_18-Q3-R50-59_BAT3C.300.fa",
                   "OdNW_18-Q3-R42-56_BAT3C.215.fa",
                   "OdNW_18-Q3-R42-56_MAXAC.118.fa",
                   "AalE_18-Q3-R2-46_BAT3C.211.fa",
                   "EsbW_18-Q3-R4-48_BAT3C.275.fa",
                   "EsbW_18-Q3-R4-48_BATAC.463.fa",
                   "OdNE_18-Q3-R46-58_BAT3C.305.fa",
                   "Skiv_18-Q3-R9-52_MAXAC.078_sub.fa",
                   "Ribe_18-Q3-R11-54_MAXAC.030.fa",
                   "EsbE_18-Q3-R3-47_MAXAC.131.fa",
                   "EsbW_18-Q3-R4-48_MAXAC.044.fa",
                   "EsbE_18-Q3-R3-47_MAXAC.059.fa",
                   "Bjer_18-Q3-R1-45_MAXAC.014.fa",
                   "Ejby_18-Q3-R6-50_BATAC.262.fa",
                   "Hirt_18-Q3-R61-65_MAXAC.059.fa",
                   "Hjor_18-Q3-R7-51_MAXAC.006.fa",
                   "Aved_18-Q3-R54-62_BAT3C.309.fa",
                   "Vibo_18-Q3-R45-57_MAXAC.072.fa")

rhodocyc_coverage <- coverage_info %>% filter(genome %in% wwtp_rhodocyc)
rhodocyc_names_covg <- rhodocyc_coverage[,c(1,24,2:23)]
rhodocyc_names_covg$avg_covg <- rowMeans(rhodocyc_names_covg[, 3:24])

write.csv(rhodocyc_names_covg, "results/relative_abundance/Danish-WWTPs-rhodocyc-coverage.csv", quote=FALSE, row.names = FALSE)
write.csv(filtered_10X_names_covg, "results/relative_abundance/Danish-WWTPs-10X-avg-coverage-genomes.csv", quote=FALSE, row.names = FALSE)


# get corresponding breadth values for top 10X genomes and rhodocyclaceae
breadth <- read.csv("results/relative_abundance/Danish-WWTPs-breadth-curated.csv")
rhodocyclaceae_breadth <- breadth %>% filter(genome %in% wwtp_rhodocyc)
rhodo_breadth_names <- left_join(rhodocyclaceae_breadth, mag_taxonomy)
rhodo_breadth_df <- rhodo_breadth_names[,c(1,24,2:23)]
top_10X_breadth <- breadth %>% filter(genome %in% filtered_10X_genomes)
top_10X_breadth_names <- left_join(top_10X_breadth, mag_taxonomy)
top_10X_df <- top_10X_breadth_names[,c(1,24,2:23)]

write.csv(rhodo_breadth_df, "results/relative_abundance/rhodocyclaceae-breadth.csv", quote=FALSE, row.names=FALSE)
write.csv(top_10X_df, "results/relative_abundance/top-10X-covg-genomes-breadth.csv", quote=FALSE, row.names=FALSE)

# melt rhodocyc covg and breadth separately 
rhodocyc_covg_values <- rhodocyc_coverage %>% select(-GTDBTax)
rhodocyc_covg_melted <- melt(rhodocyc_covg_values, id.vars="genome", variable="sample", value.name = "coverage")
rhodocyc_breadth_melted <- melt(rhodocyclaceae_breadth, id.vars="genome", variable="sample", value.name="breadth")
rhodocyc_covg_breadth <- left_join(rhodocyc_covg_melted, rhodocyc_breadth_melted)

# Rhodocyclaceae checks
# filter genomes in samples that have greater than 5X coverage and 70% of genome is covered (breadth)
rhodocyc_covg_breadth_filtered <- rhodocyc_covg_breadth %>% filter(coverage > 5 & breadth > 0.7)
# more stringent - 10X coverage and 0.9 breadth
rhodocyc_covg_breadth_stringent <- rhodocyc_covg_breadth %>% filter(coverage > 10 & breadth > 0.9)
rhodocyc_covg_breadth_stringent_names <- left_join(rhodocyc_covg_breadth_stringent, mag_taxonomy)
rhodocyc_stringent_counts <- rhodocyc_covg_breadth_stringent_names %>% count(genome, GTDBTax)

write.csv(rhodocyc_covg_breadth_stringent_names, "results/relative_abundance/rhodocyc-covg-breadth-stringent-stats.csv", row.names=FALSE, quote=FALSE)
write.csv(rhodocyc_stringent_counts, "results/relative_abundance/rhodocyc-covg-breadth-stringent-sample-counts.csv", row.names = FALSE, quote=FALSE)

rhodocyc_covg_breadth_filtered_names <- left_join(rhodocyc_covg_breadth_filtered, mag_taxonomy)
rhodocyc_covg_breadth_filtered_stringent_names <- left_join(rhodocyc_covg_breadth_stringent, mag_taxonomy)

rhodocyc_covg_breadth_filtered_stringent_names %>% count(genome)
filtered_counts <- rhodocyc_covg_breadth_filtered_names %>% count(genome, GTDBTax)

# All with high coverage checks
all_coverage_melted <- melt(coverage, id.vars="genome", variable="sample", value.name = "coverage")
all_breadth_melted <- melt(breadth, id.vars="genome", variable="sample", value.name="breadth")
all_covg_breadth <- left_join(all_coverage_melted, all_breadth_melted)
all_filtered <- all_covg_breadth %>% filter(coverage > 5 & breadth > 0.7)
all_filtered_stringent <- all_covg_breadth %>% filter(coverage > 10 & breadth > 0.9)

all_filtered_stringent_names <- left_join(all_filtered_stringent, mag_taxonomy)
filtered_stringent_multiple_genomes_all <- all_filtered_stringent_names %>% group_by(genome) %>% filter(n() >11)

high_all_groups <- filtered_stringent_multiple_genomes_all %>% count(genome, GTDBTax)

write.csv(filtered_stringent_multiple_genomes_all, "results/relative_abundance/all-filtered-genomes-WWTPs-stringent-cutoffs.csv", quote=FALSE, row.names = FALSE)
write.csv(high_all_groups, "results/relative_abundance/all-filtered-genomes-sample-counts-WWTPs.csv", row.names = FALSE, quote=FALSE)
