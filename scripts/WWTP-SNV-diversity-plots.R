library(tidyverse)
library(wesanderson)

wwtp_diversity <- read.csv("results/SNVs/WWTP/Rhodocyc-WWTPs-nucleotide-diversity.csv")

wwtp_diversity$Sample <- gsub("_.*", "", wwtp_diversity$Sample)

plants_ordered <- c("Hirt", "Hjor", "AalE", "AalW", "Mari", "Skiv", "Vibo", "Rand", "Ega", "Viby", "EsbE", "EsbW", "Fred", "Ribe", "Hade", "OdNW", "Ejby", "OdNE", "Lyne", "Damh", "Aved")

wwtp_diversity$Sample <- factor(wwtp_diversity$Sample, levels=plants_ordered)

rhodocyc_sample_diversity <- wwtp_diversity %>% ggplot(aes(x=Sample, y=Nucleotide_Diversity, color=Genome)) + scale_color_brewer(palette="Set3") + geom_point(size=3) + facet_wrap(~ Species, nrow=2) + theme_bw() + theme(axis.text.x= element_text(angle=70, hjust=1)) + ylab("Nucleotide Diversity\n") + xlab("\nWastewater Treatment Plant") + theme(strip.text=element_text(face="italic"), axis.title.x=element_text(face="bold"), axis.title.y=element_text(face="bold"))
rhodocyc_sample_diversity

ggsave("figures/Danish-WWTP-rhodocyc-nucleotide-diversity.png", rhodocyc_sample_diversity, width=10, height=8, units=c('in'))
