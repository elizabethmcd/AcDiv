library(tidyverse)

# r1r2
r1r2_abundance <- read.delim("results/relative_abundance/R1R2.txt", sep="\t")
R1R2_2005_UW1_UW3 <- r1r2_abundance %>% select(Genome, R1R2_refs.fasta.2005.06.14.EBPR.qced.fastq.Relative.Abundance....) %>% filter(Genome == 'UW3' | Genome == 'GCA_000024165.1')
colnames(R1R2_2005_UW1_UW3) <- c("genome", "abund_2005_06_14")
R1R2_2005_UW1_UW3 %>% ggplot(aes(x=genome, y=abund_2005_06_14)) + geom_bar(stat="identity")

# r3r4
r3r4_abundance <- read.delim("results/relative_abundance/R3R4.txt", sep="\t")
R3R4_mixed <- r3r4_abundance %>% select(Genome, R3R4_refs.fasta.EBPRWIReactor34.cleaned.fastq.Relative.Abundance....) %>% filter(Genome != "unmapped")
R3R4_mixed %>% ggplot(aes(x=Genome, y=R3R4_refs.fasta.EBPRWIReactor34.cleaned.fastq.Relative.Abundance....)) + geom_bar(stat="identity") + theme_classic()
