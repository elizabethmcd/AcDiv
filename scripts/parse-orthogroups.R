library(tidyverse)
library(reshape2)

# all homologs
homolog_matrix <- read.delim("results/pangenomics/homolog_matrix.txt", sep="\t")

# genome metadata
metadata <- read.csv("metadata/accumulibacter_refs_genomes_HQ_metadata.csv") %>% select(genome, group_name)
metadata$genome <- gsub("-", ".", metadata$genome)
matrix_melted_names <- left_join(matrix_melted, metadata)

matrix_wide <- matrix_melted_names %>% 
  pivot_wider(names_from = group, values_from = count) %>% 
  select(-genome)

aggregate(matrix_wide[500:600], list(matrix_wide$group_name), mean)

# single copy core orthologs in all sampled genomes (accumulibacter and burk outgroups)

colnames(homolog_matrix)[1] <- c("group")
homolog_table <- homolog_matrix %>% column_to_rownames(var="group")
core_groups <- homolog_table %>% 
  filter_all(all_vars(. == 1)) %>%
  rownames_to_column("group") %>% 
  pull(group)

core_table <- as.data.frame(core_groups)
core_table$label <- "core"
colnames(core_table)[1] <- "group"

# orthologs that are only in accumulibacter genomes and not the outgroups (in single copy) 

accumulibacter_genomes <- c("group",
                            "GCA_009467855.1",
                            "GCA_000585075.1",
                            "UW3",
                            "UW4",
                            "UW8.POB",
                            "GCA_000987445.1",
                            "GCA_000585055.1",
                            "Delftensis",
                            "GCA_012939955.1",
                            "GCA_003332265.1",
                            "GCA_000024165.1",
                            "UW5",
                            "UW9.POB",
                            "Fred_18.Q3.R57.64_MAXAC.027",
                            "Hade_18.Q3.R52.61_BATAC.726",
                            "OdNE_18.Q3.46.58_BAT3C.415",
                            "UW10.POB",
                            "GCA_005889575.1",
                            "GCA_000987395.1",
                            "GCA_012940005.1",
                            "GCA_013414765.1",
                            "GCA_000584955.2",
                            "GCA_000584975.2",
                            "UW11.POB",
                            "UW6",
                            "EsbW_18.Q3.R4.48_BATAC.285",
                            "UW12.POB",
                            "GCA_013347225.1",
                            "GCA_005524045.1",
                            "GCA_000585015.1",
                            "UW13.POB",
                            "UW7",
                            "Fred_18.Q3.R57.64_BAT3C.720",
                            "GCA_900089955.1")

acc_cols <- which(colnames(homolog_matrix) %in% accumulibacter_genomes)
accumulibacter_orthogroups <- homolog_matrix[,sort(c(acc_cols))] %>% 
  column_to_rownames(var="group") %>% 
  filter_all(all_vars(. == 1)) %>%
  rownames_to_column("group") %>% 
  pull(group)

acc_only_table <- as.data.frame(accumulibacter_orthogroups, stringsAsFactors = FALSE)
colnames(acc_only_table) <- "group"

acc_only_groups <- acc_only_table[! acc_only_table$group %in% core_groups, ]
acc_only_final_table <- as.data.frame(acc_only_groups, stringsAsFactors = FALSE)
colnames(acc_only_final_table)[1] <- "group"
acc_only_final_table$label <- "acc_core"

# orthologs in uw1 but not in uw3 

uw1_only_groups <- homolog_table %>% 
  filter(UW1 == 1 & UW3 == 0) %>%
  rownames_to_column("group") %>% 
  pull(group)

uw1_only_table <- as.data.frame(uw1_only_groups)
uw1_only_table$label <- "accessory"
colnames(uw1_only_table)[1] <- "group"

# orthologs in uw3 but not in uw1

uw3_only_groups <- homolog_table %>% 
  filter(UW3 == 1 & UW1 == 0) %>%
  rownames_to_column("group") %>% 
  pull(group)

uw3_only_table <- as.data.frame(uw3_only_groups)
uw3_only_table$label <- "accessory"
colnames(uw3_only_table)[1] <- "group"

# groups to locus tags
# Do for UW1, UW3, UW5, and UW7 to show for the beginning enrichment SNVs/diversity of those loci?
# Gene nucleotide diversity by the locus tag

locus_tags <- read_tsv("results/pangenomics/locustag_matrix.txt")
colnames(locus_tags)[1] <- c("group")

uw_locus_tags <- locus_tags %>% 
  select(group, UW1, UW3, UW5, UW7)

core_locus_tags <- left_join(core_table, uw_locus_tags)

acc_locus_tags <- left_join(acc_only_final_table, uw_locus_tags)

uw3_accessory_locus_tags <- left_join(uw3_only_table, uw_locus_tags)

# UW1 SNVs to ortholog definitions
uw1_core_locus_tags <- core_locus_tags %>% 
  select(group, label, UW1)
colnames(uw1_core_locus_tags)[3] <- "locus_tag"

uw1_accum_locus_tags <- acc_locus_tags %>% 
  select(group, label, UW1)
colnames(uw1_accum_locus_tags)[3] <- "locus_tag"

uw1_accessory_locus_tags <- left_join(uw1_only_table, uw_locus_tags) %>% 
  select(group, label, UW1)
colnames(uw1_accessory_locus_tags)[3] <- "locus_tag"

# diversity table for uw1 
uw1_2005_snvs <- read_tsv("results/SNVs/R1R2_vs_R3R4/R1R2_2005_UW1_gene_info.tsv") %>% 
  select(gene, nucl_diversity)
colnames(uw1_2005_snvs)[1] <- "locus_tag"
uw1_2005_snvs$locus_tag <- gsub("gnl\\|X\\|", "", uw1_2005_snvs$locus_tag)

# core tags div table 
uw1_core_tags_snvs <- left_join(uw1_core_locus_tags, uw1_2005_snvs)
# acc core tags div table
uw1_accum_tags_snvs <- left_join(uw1_accum_locus_tags, uw1_2005_snvs)
# uw1 accessory tags div table
uw1_accessory_tags_snvs <- left_join(uw1_accessory_locus_tags, uw1_2005_snvs)
# combine all three for a UW1 table 

uw1_groups_cogs_div_table <- rbind(uw1_core_tags_snvs, uw1_accum_tags_snvs, uw1_accessory_tags_snvs) %>% 
  mutate(reference = "UW1") %>% 
  select(reference, label, locus_tag)


# UW3 SNVS/diversity table to merge with locus tags

uw3_core_locus_tags <- core_locus_tags %>% 
  select(group, label, UW3)
colnames(uw3_core_locus_tags)[3] <- "locus_tag"

uw3_accum_locus_tags <- acc_locus_tags %>% 
  select(group, label, UW3)
colnames(uw3_accum_locus_tags)[3] <- "locus_tag"

uw3_accessory_locus_tags <- left_join(uw3_only_table, uw_locus_tags) %>% 
  select(group, label, UW3)
colnames(uw3_accessory_locus_tags)[3] <- "locus_tag"

# diversity table for uw3
uw3_2005_snvs <- read_tsv("results/SNVs/R1R2_vs_R3R4/R1R2_2005_UW3_gene_info.tsv") %>% 
  select(gene, nucl_diversity)
colnames(uw3_2005_snvs)[1] <- "locus_tag"
uw3_2005_snvs$locus_tag <- gsub("gnl\\|X\\|", "", uw3_2005_snvs$locus_tag)

# core tags div table 
uw3_core_tags_snvs <- left_join(uw3_core_locus_tags, uw3_2005_snvs)
# acc core tags div table
uw3_accum_tags_snvs <- left_join(uw3_accum_locus_tags, uw3_2005_snvs)
# uw1 accessory tags div table
uw3_accessory_tags_snvs <- left_join(uw3_accessory_locus_tags, uw3_2005_snvs)
# combine all three for a UW1 table 

uw3_groups_cogs_div_table <- rbind(uw3_core_tags_snvs, uw3_accum_tags_snvs, uw3_accessory_tags_snvs) %>% 
  mutate(reference = "UW3") %>% 
  select(reference, label, locus_tag)

###### combined group/locustag table for both uw1 and uw3 for use in flanking data 
# come back later for using Uw1/Uw3 in the AcDiv data

uw_combined_cogs_tags_table <- rbind(uw1_groups_cogs_div_table, uw3_groups_cogs_div_table)
write.csv(uw_combined_cogs_tags_table, "results/pangenomics/UW1_UW3_COGS_TAGS_TABLE.csv", row.names = FALSE, quote=FALSE)
