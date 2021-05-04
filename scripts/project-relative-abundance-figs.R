library(tidyverse)
library(reshape2)
library(viridis)

# r1r2
r1r2_abundance <- read.delim("results/relative_abundance/R1R2.txt", sep="\t")
R1R2_2005_UW1_UW3 <- r1r2_abundance %>% select(Genome, R1R2_refs.fasta.2005.06.14.EBPR.qced.fastq.Relative.Abundance....) %>% filter(Genome == 'UW3' | Genome == 'GCA_000024165.1')
colnames(R1R2_2005_UW1_UW3) <- c("genome", "abund_2005_06_14")

r1r1_acc_abundance <- R1R2_2005_UW1_UW3 %>% 
  ggplot(aes(x=genome, y=abund_2005_06_14, fill=as.factor(genome))) + geom_bar(stat="identity", colour="black", size=0.5) + scale_fill_manual(values=c("#2BAEB3", "#404272")) + scale_y_continuous(name="Relative Abundance", limits=c(0,65), expand=c(0,0)) + theme_classic()

# r3r4
r3r4_abundance <- read.delim("results/relative_abundance/R3R4.txt", sep="\t")
R3R4_mixed <- r3r4_abundance %>% 
  select(Genome, R3R4_refs.fasta.EBPRWIReactor34.cleaned.fastq.Relative.Abundance....) %>% 
  filter(Genome != "unmapped")
colnames(R3R4_mixed) <- c("Genome", "R3R4_Relative_Abundance")

r3r4_acc_abundance <- R3R4_mixed %>% 
  ggplot(aes(x=reorder(Genome, -R3R4_Relative_Abundance), y=R3R4_Relative_Abundance, fill=as.factor(Genome))) + geom_bar(stat="identity", colour="black", size=0.5) + scale_fill_manual(values=c("#f768a1", "#8c96c6", "#a1dab4", "#fecc5c")) + scale_y_continuous(name="Relative Abundance", limits=c(0,23), expand=c(0,0)) + theme_classic() + scale_x_discrete(name="Genome")

# pob
pob_stats <- read.csv("metadata/POS-MAGs-table.csv")
pob_rhodo_stats <- pob_stats %>% 
  select(code, abund_2015.07.16, abund_2015.07.24, abund_2015.08.06) %>% 
  filter(code == "CAPIA" | code == "CAPIIA" | code == "CAPIIB" | code == 'CAPIIC' | code == "CAPIID" | code == "CAPIIF" | code == "DECH1" | code == "SULF1" | code == "THAU1")
pob_melt <- melt(pob_rhodo_stats, id.vars="code", variable="date", value.name = "relative_abundance")

pob_rhodo_succession <- ggplot(pob_melt, aes(x=code, y=fct_rev(date), fill=relative_abundance)) + geom_tile(color="white") + scale_fill_viridis(option="inferno", alpha=1, begin=0, end=1, direction=-1) + theme(axis.text.x= element_text(angle=85, hjust=1))

# wwtps

wwtp_rhodocycs <- c("Fred_18-Q3-R57-64_MAXAC.027",
                   "Hade_18-Q3-R52-61_BATAC.726",
                   "OdNE_18-Q3-R46-58_BAT3C.415",
                   "EsbW_18-Q3-R4-48_BATAC.285",
                   "Fred_18-Q3-R57-64_BAT3C.720",
                   "Aved_18-Q3-R54-62_MAXAC.403",
                   "Lyne_18-Q3-R50-59_BAT3C.300",
                   "OdNW_18-Q3-R42-56_BAT3C.215",
                   "OdNW_18-Q3-R42-56_MAXAC.118",
                   "AalE_18-Q3-R2-46_BAT3C.211",
                   "EsbW_18-Q3-R4-48_BAT3C.275",
                   "EsbW_18-Q3-R4-48_BATAC.463",
                   "OdNE_18-Q3-R46-58_BAT3C.305",
                   "Skiv_18-Q3-R9-52_MAXAC.078_sub",
                   "Ribe_18-Q3-R11-54_MAXAC.030",
                   "EsbE_18-Q3-R3-47_MAXAC.131",
                   "EsbW_18-Q3-R4-48_MAXAC.044",
                   "EsbE_18-Q3-R3-47_MAXAC.059",
                   "Bjer_18-Q3-R1-45_MAXAC.014",
                   "Ejby_18-Q3-R6-50_BATAC.262",
                   "Hirt_18-Q3-R61-65_MAXAC.059",
                   "Hjor_18-Q3-R7-51_MAXAC.006",
                   "Aved_18-Q3-R54-62_BAT3C.309",
                   "Vibo_18-Q3-R45-57_MAXAC.072")

danish_wwtps_ra <- read.csv("results/relative_abundance/Danish-WTTPs-relative-abundance.csv")
rhodocyc_wwtps_ra <- danish_wwtps_ra %>% 
  filter(Genome %in% wwtp_rhodocycs)

colnames(rhodocyc_wwtps_ra)[1] <- c("Genome_")
colnames(rhodocyc_wwtps_ra) <- sub("_.*", "", colnames(rhodocyc_wwtps_ra))
drop_kalu <- rhodocyc_wwtps_ra %>% select(-Kalu)
wwtps_ordered <- drop_kalu %>% select(Genome, Hirt, Hjor, AalE, AalW, Mari, Skiv, Vibo, Rand, Ega, Viby, EsbE, EsbW, Fred, Ribe, Hade, OdNW, Ejby, OdNE, Lyne, Damh, Aved)
rhodocyc_wwtps_melt <- melt(wwtps_ordered, id.vars = "Genome", variable.name = "sample", value.name = "relative_abundance")

rhodocyc_tax_order <- c("OdNW_18-Q3-R42-56_MAXAC.118",
                        "Lyne_18-Q3-R50-59_BAT3C.300",
                        "Aved_18-Q3-R54-62_MAXAC.403",
                        "OdNW_18-Q3-R42-56_BAT3C.215",
                        "EsbW_18-Q3-R4-48_BATAC.285",
                        "OdNE_18-Q3-R46-58_BAT3C.415",
                        "Fred_18-Q3-R57-64_BAT3C.720",
                        "OdNE_18-Q3-R46-58_BAT3C.305",
                        "Skiv_18-Q3-R9-52_MAXAC.078_sub",
                        "AalE_18-Q3-R2-46_BAT3C.211",
                        "EsbW_18-Q3-R4-48_BAT3C.275",
                        "EsbW_18-Q3-R4-48_BATAC.463",
                        "Ribe_18-Q3-R11-54_MAXAC.030",
                        "EsbW_18-Q3-R4-48_MAXAC.044",
                        "EsbE_18-Q3-R3-47_MAXAC.131",
                        "EsbE_18-Q3-R3-47_MAXAC.059",
                        "Hjor_18-Q3-R7-51_MAXAC.006",
                        "Bjer_18-Q3-R1-45_MAXAC.014",
                        "Hirt_18-Q3-R61-65_MAXAC.059",
                        "Ejby_18-Q3-R6-50_BATAC.262",
                        "Aved_18-Q3-R54-62_BAT3C.309",
                        "Vibo_18-Q3-R45-57_MAXAC.072")

rhodocyc_wwtps_melt$Genome <- factor(rhodocyc_wwtps_melt$Genome, levels=c(rhodocyc_tax_order))

rhodocy_wwtps_plot <- ggplot(rhodocyc_wwtps_melt, aes(x=sample, y=fct_rev(Genome), fill=relative_abundance)) + geom_tile(color="white") + scale_fill_viridis(option="magma", alpha=1, begin=0, end=1, direction=-1) + theme(axis.text.x= element_text(angle=85, hjust=1))
rhodocy_wwtps_plot

# save plots
ggsave(plot=r1r1_acc_abundance, filename="figures/r1r2_acc.png", height=10, width=14, units=c('cm'))
ggsave(plot=r3r4_acc_abundance, filename="figures/r3r4_acc.png", height=10, width=14, units=c("cm"))
ggsave(plot=pob_rhodo_succession, filename="figures/POB-rhodocyc-abund-cycles.png", height=5, width=17, units=c("cm"))
ggsave(plot=rhodocy_wwtps_plot, filename="figures/Danish-WWTPs-rhodocyc-abundance.png", height=7, width=30, units=c("cm"))
