library(tidyverse)

pob_covg_breadth <- read.csv("results/covg_breadth/POB-sample-covg-breadth.csv")
pob_stringent <- pob_covg_breadth %>% filter(coverage > 10 & breadth > 0.9)

metadata <- read.csv("metadata/POS-MAGs-table.csv") %>% 
  select(bin, code)

colnames(metadata)[1] <- c("genome")

pob_stringent$genome <- gsub(".fa", "", pob_stringent$genome)
pob_covg_breadth$genome <- gsub(".fa", "", pob_covg_breadth$genome)

pob_stringent_codes <- left_join(pob_stringent, metadata)
pob_stringent_codes %>% count(code)

pob_stringent_rhodo <- pob_stringent_codes %>% filter(code == "CAPIA" | code == "CAPIIA" | code == "CAPIIB" | code == 'CAPIIC' | code == "CAPIID" | code == "CAPIIF" | code == "DECH1" | code == "SULF1" | code == "THAU1" )

pob_covg_all_codes <- left_join(pob_covg_breadth, metadata)
pob_rhodo <- pob_covg_all_codes %>% filter(code == "CAPIA" | code == "CAPIIA" | code == "CAPIIB" | code == 'CAPIIC' | code == "CAPIID" | code == "CAPIIF" | code == "DECH1" | code == "SULF1" | code == "THAU1" )
