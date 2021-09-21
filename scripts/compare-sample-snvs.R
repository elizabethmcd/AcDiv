library(tidyverse)

# SNVs in UW1 in 2013 sample 
uw1_2013_snvs <- read_tsv("results/SNVs/R1R2_vs_R3R4/UW1_SNVs/UW1_2013_05_23_SNVs.tsv") %>% 
  select(gene, position, position_coverage, allele_count, ref_freq, var_freq) %>% 
  drop_na()

uw1_segregating_SNVs <- uw1_2013_snvs %>% 
  filter(allele_count == 2)

uw1_segregating_SNVs$gene_site <- paste(uw1_segregating_SNVs$gene, "_", uw1_segregating_SNVs$position)
uw1_segregating_SNVs$sample <- "r1r2_2013"

# SNVs in UW1 in R3R4 2015 sample
uw1_2015_snvs <- read_tsv("results/SNVs/R1R2_vs_R3R4/UW1_SNVs/UW1_2015_Reactor34_SNVs.tsv") %>% 
  select(gene, position, position_coverage, allele_count, ref_freq, var_freq) %>% 
  drop_na()

uw1_2015_seg_SNVs <- uw1_2015_snvs %>% 
  filter(allele_count == 2)

uw1_2015_seg_SNVs$gene_site <- paste(uw1_2015_seg_SNVs$gene, "_", uw1_2015_seg_SNVs$position)
uw1_2015_seg_SNVs$sample <- "r3r4_2015"


# SNVs in UW1 in 2005 sample 
uw1_2005_SNVs <- read_tsv("results/SNVs/R1R2_vs_R3R4/UW1_SNVs/UW1_2005_SNVs.tsv") %>% 
  select(gene, position, position_coverage, allele_count, ref_freq, var_freq) %>% 
  drop_na()

uw1_2005_seg_sites <- uw1_2005_SNVs %>% 
  filter(allele_count == 2)

uw1_2005_seg_sites$gene_site <- paste(uw1_2005_seg_sites$gene, "_", uw1_2005_seg_sites$position)
uw1_2005_seg_sites$sample <- "r1r2_2005"

# Combining UW1 2005 and 2013 reciprocal SNVs 
uw1_2005_seg_snvs <- uw1_2005_seg_sites %>% 
  select(gene_site, var_freq, sample)
uw1_2013_seg_snvs <- uw1_segregating_SNVs %>% 
  select(gene_site, var_freq, sample)

# Combine UW1 2005 and 2015 reciprocal SNVs
uw1_2015_seg_sites <- uw1_2015_seg_SNVs %>% 
  select(gene_site, var_freq, sample)

# combine and pivot to combine gene sites
uw1_samples_snvs <- rbind(uw1_2005_seg_snvs, uw1_2013_seg_snvs) %>% 
  pivot_wider(names_from = sample, values_from=var_freq)

uw1_2005_2015_samples_snvs <- rbind(uw1_2005_seg_snvs, uw1_2015_seg_sites) %>% 
  pivot_wider(names_from = sample, values_from = var_freq)

# change NAs to 0 for allele freq of 0 for that variant and calculations
uw1_samples_snvs[is.na(uw1_samples_snvs)] <- 0
uw1_samples_snvs$hs <- (uw1_samples_snvs$r1r2_2005 * (1 - uw1_samples_snvs$r1r2_2005)) + (uw1_samples_snvs$r1r2_2013 * (1-uw1_samples_snvs$r1r2_2013))
uw1_samples_snvs$avg <- (uw1_samples_snvs$r1r2_2005 + uw1_samples_snvs$r1r2_2013) / 2
uw1_samples_snvs$ht <- (2 * uw1_samples_snvs$avg) * (1 - uw1_samples_snvs$avg)
uw1_samples_snvs$fst <- (uw1_samples_snvs$ht - uw1_samples_snvs$hs) / uw1_samples_snvs$ht

uw1_samples_diff <- uw1_samples_snvs %>% 
  select(gene_site, r1r2_2005, r1r2_2013, fst)

uw1_2005_2015_samples_snvs[is.na(uw1_2005_2015_samples_snvs)] <- 0
uw1_2005_2015_samples_snvs$hs <- (uw1_2005_2015_samples_snvs$r1r2_2005 * (1 - uw1_2005_2015_samples_snvs$r1r2_2005)) + (uw1_2005_2015_samples_snvs$r3r4_2015 * (1-uw1_2005_2015_samples_snvs$r3r4_2015))
uw1_2005_2015_samples_snvs$avg <- (uw1_2005_2015_samples_snvs$r1r2_2005 + uw1_2005_2015_samples_snvs$r3r4_2015) / 2
uw1_2005_2015_samples_snvs$ht <- (2 * uw1_2005_2015_samples_snvs$avg) * (1 - uw1_2005_2015_samples_snvs$avg)
uw1_2005_2015_samples_snvs$fst <- (uw1_2005_2015_samples_snvs$ht - uw1_2005_2015_samples_snvs$hs) / uw1_2005_2015_samples_snvs$ht

uw1_2005_2015_samples_diff <- uw1_2005_2015_samples_snvs %>% 
  select(gene_site, r1r2_2005, r3r4_2015, fst)

uw1_2005_2015_samples_diff$gene <- sub("^([^_]+_[^_]+_[^_]*).*", "\\1", uw1_2005_2015_samples_diff$gene_site)
uw1_2005_2015_samples_diff$gene <- trimws(uw1_2005_2015_samples_diff$gene, which=c("right"))
uw1_2005_2015_samples_diff_index <- left_join(uw1_2005_2015_samples_diff, uw1_gene_info)

uw1_2005_2015_fst <- uw1_2005_2015_samples_diff_index %>% ggplot(aes(x=index, y=fst)) + geom_point(color="tomato3") + scale_y_continuous(limits=c(0,1)) + scale_x_continuous(breaks=seq(0,5000, 500)) + labs(title = "Fst between UW1 IIA R1R2 and R3R4 Populations") + xlab("Gene index") + ylab("Fst") + theme_bw() + theme(axis.title.y=element_text(face="bold"), axis.title.x=element_text(face="bold"))
uw1_2005_2015_fst
ggsave("figures/uw1_2005_2015_fst_comparisons.png", uw1_2005_2015_fst, width=11, height=5, units=c("in"))

# gene info
uw1_gene_info <- read_tsv("results/SNVs/R1R2_vs_R3R4/R1R2_2005_UW1_gene_info.tsv") %>% 
  select(gene, gene_length)
uw1_gene_info$index <- seq.int(nrow(uw1_gene_info))
uw1_samples_diff$gene <- sub("^([^_]+_[^_]+_[^_]*).*", "\\1", uw1_samples_diff$gene_site)
uw1_samples_diff$gene <- trimws(uw1_samples_diff$gene, which = c("right"))
uw1_samples_diff_index <- left_join(uw1_samples_diff, uw1_gene_info)

fst_uw1 <- uw1_samples_diff_index %>% ggplot(aes(x=index, y=fst)) + geom_point(color="darkorchid4") + scale_y_continuous(limits=c(0,1)) + xlab("Gene index") + ylab('Fst') + theme_bw() + labs(title = "Fst of UW1 IIA between 2005 and 2013")
fst_uw1

ggsave("figures/uw1-2005-2013-fst.png", fst_uw1, width=7, height=4.5, units=c('in'))

# genes with high fst between r1r2 and r3r4

high_fst <- uw1_2005_2015_samples_diff_index %>% 
  filter(fst > 0.2) %>% 
  select(gene, r1r2_2005, r3r4_2015, fst)
high_fst$gene <- gsub("gnl\\|X\\|", "", high_fst$gene)

r1r2_high_fst <- left_join(high_fst, uw1_snv_diversity)
colnames(r1r2_high_fst) <- c("gene", "r1r2_2005", "r3r4_2015", "fst", "r1r2_SNV_count", "r1r2_nucl_diversity", "index")


uw1_2015_div <- read_tsv("results/SNVs/R1R2_vs_R3R4/R3R4_2015_UW1_gene_info.tsv") %>% 
  select(gene, SNV_count, nucl_diversity)
uw1_2015_div$gene <- gsub("gnl\\|X\\|", "", uw1_2015_div$gene)

uw1_r1r2_r3r4_fst <- left_join(r1r2_high_fst, uw1_2015_div) %>% 
  select(gene, r1r2_2005, r3r4_2015, fst, r1r2_SNV_count, r1r2_nucl_diversity, SNV_count, nucl_diversity, index)
colnames(uw1_r1r2_r3r4_fst) <- c("gene", "r1r2_2005", "r3r4_2015", "fst", "r1r2_SNV_count", "r1r2_nucl_diversity", "r3r4_SNV_count", "r3r4_nucl_diversity", "index")

high_r1r2_r3r4_fst_div <- uw1_r1r2_r3r4_fst %>% 
  select(gene, r1r2_nucl_diversity, r3r4_nucl_diversity) %>%
  pivot_longer(!gene, names_to=c("timepoint"), values_to=c("nucleotide_diversity"))

high_r1r2_r3r4_fst_div %>% ggplot(aes(x=timepoint, y=nucleotide_diversity, fill=gene)) + geom_point() + geom_line()
  
fst_count <- high_fst %>% 
  group_by(gene) %>% 
  summarize(Count = n()) %>% 
  ungroup() %>% 
  arrange(desc(Count))

uw1_2005_2015_fst_avg <- uw1_2005_2015_samples_diff_index %>% 
  select(index, gene, fst) %>% 
  group_by(index, gene) %>% 
  summarize(mean(fst))
colnames(uw1_2005_2015_fst_avg) <- c("index", "gene", "Fst")

uw1_fst_avg_plot <- uw1_2005_2015_fst_avg %>% ggplot(aes(x=index, y=Fst)) + geom_point(color="tomato3") + scale_y_continuous(limits=c(0,0.5), expand=c(0,0)) + scale_x_continuous(limits=c(0,4500), breaks=seq(0,4500, 500)) + labs(title = "Fst between UW1 IIA R1R2 and R3R4 Populations") + xlab("Gene index") + ylab("Fst") + theme_bw() + theme(axis.title.y=element_text(face="bold"), axis.title.x=element_text(face="bold"))

ggsave("figures/UW1-FST-avg-plot.png", uw1_fst_avg_plot, width=8, height=4, units=c("in"))

avg <- uw1_2005_2015_fst_avg %>% 
  filter(Fst > 0.2)
