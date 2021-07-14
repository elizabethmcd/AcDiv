library(tidyverse)
library(gggenes)
library(cowplot)
library(gridGraphics)

uw1_2015_snvs_table <- read_tsv("results/SNVs/R1R2_vs_R3R4/UW1_SNVs/UW1_2015_Reactor34_SNVs.tsv") %>% 
  select(scaffold, gene, position, allele_count, A, C, T, G)

uw1_2015_snvs_table_plasmids <- uw1_2015_snvs_table %>% 
  filter(scaffold == 'gnl|X|DGOCCOHL_2' | scaffold == 'gnl|X|DGOCCOHL_3' | scaffold == 'gnl|X|DGOCCOHL_4') %>% 
  drop_na(gene)

plasmid_2 <- uw1_2015_snvs_table_plasmids %>% 
  filter(scaffold == "gnl|X|DGOCCOHL_2") %>% 
  select(position, A, C, T, G) %>% 
  filter(position > 140000)
plasmid_2$sum <- rowSums(plasmid_2[,2:5])
plasmid_2_gene_info <- left_join(plasmid_2, uw1_2005_2015_samples_diff_index)

plasmid_2_table <- plasmid_2 %>% 
  pivot_longer(!c(position, sum), names_to="allele", values_to="count") %>% 
  mutate(freq = count / sum) %>% 
  select(position, allele, freq)

plasmid_2_var_plot <- plasmid_2_table %>% ggplot(aes(x=position, y=freq, color=allele)) + geom_bar(stat="identity") + scale_color_brewer(palette="Set1") + scale_x_continuous(limits=c(152000,164000), breaks=seq(152000, 164000, 2000)) + scale_y_continuous(expand=c(0,0)) + theme_classic() + theme(legend.position="none", axis.text.y=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank())
plasmid_2_var_plot
ggsave("figures/uw1-plasmid-variation.png", plasmid_2_var_plot, width=11, height=3, units=c("in"))

plasmid_3 <- uw1_2015_snvs_table_plasmids %>% 
  filter(scaffold == "gnl|X|DGOCCOHL_3") %>% 
  select(position, A, C, T, G) %>% 
  filter(position > 21000)
plasmid_3$sum <- rowSums(plasmid_3[,2:5])
plasmid_3_table <- plasmid_3 %>% 
  pivot_longer(!c(position, sum), names_to="allele", values_to="count") %>% 
  mutate(freq = count / sum) %>% 
  select(position, allele, freq)

plasmid_3_table %>% ggplot(aes(x=position, y=freq, color=allele)) + geom_bar(stat="identity") + theme_classic()

plasmid_4 <- uw1_2015_snvs_table_plasmids %>% 
  filter(scaffold == "gnl|X|DGOCCOHL_4") %>% 
  select(position, A, C, T, G)
plasmid_4$sum <- rowSums(plasmid_4[,2:5])
plasmid_4_table <- plasmid_4 %>% 
  pivot_longer(!c(position, sum), names_to="allele", values_to="count") %>% 
  mutate(freq = count / sum) %>% 
  select(position, allele, freq)

plasmid_4_table %>% ggplot(aes(x=position, y=freq, color=allele)) + geom_bar(stat="identity") + theme_classic()

### Plot of genes in 1st plasmid at end with variation

genes <- read.csv("results/annotations/plasmid_gene_index.csv")
plasmid2_alignment <- make_alignment_dummies(genes, aes(xmin=start, xmax=end, y=genome, id=gene))


plasmid_genes <- ggplot(genes, aes(xmin=start, xmax=end, y=genome, fill=gene, label=gene)) + geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(2, "mm")) + theme_genes() + theme(legend.position="none")
ggsave("figures/uw1-plasmid-map.png", plasmid_genes, width=11, height=3, units=c('in'))

### Fst plot grid over variation of specific genes 


