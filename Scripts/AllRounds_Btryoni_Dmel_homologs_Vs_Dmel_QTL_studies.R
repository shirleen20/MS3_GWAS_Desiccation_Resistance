# __________________________28/07/2026______________________________

# We compare our B.tryoni homologs for All individual rounds with Drossie data with 
# candidates found in Dmel QTL papers 
## check to see if any of the genes found in B.tryoni Rounds are same as in Drossie studies

library(readxl)
library(tidyverse)
library(AnnotationDbi)
library(org.Dm.eg.db)
library(GO.db)
library(data.table)
library(purrr)
library(readr)


rm(list=ls())


#_______________________________________________

# Read Fanara genes and earlier GWAS genes called  
Fanara_genes <- read_excel("Drossie_GWAS/Fanara_genes.xlsx", sheet = 1) %>%
  dplyr::select(-"Site Class") %>% dplyr::rename(Gene_name = Symbol, FLYBASE = ID)

QTLStudy_genes <- read_excel("Drossie_GWAS/GWAS_genes_Drosophila.xlsx",sheet = 1,skip = 1)%>% 
  mutate(across(where(is.character),~ifelse(. %in% c("-", "NA", ""), NA, .))) %>% 
  dplyr::select(2:3) 
#____________________________________________________________________

# Read the file containing tryoni Dmel homologs and B tryoni genes

Genes_R4 <- read_tsv("GO_data/R4_extracted_genes_table.txt",
                     col_names = c("chrom","start","end","strand","GeneID","LOC","biotype","GeneName"), quote = "")
                                                                              
Homologs_R4 <- read.table("GO_data/forward_best_hit_homologs_filtered_R4.txt",
                          sep = "\t",header = FALSE, stringsAsFactors = FALSE,
                          col.names = c("qacc", "sacc", "pident", "length", "mismatch","gapopen","qstart", 
                                        "qend", "sstart", "send", "evalue", "bitscore")) %>% 
  as_tibble() %>% rename_with(trimws) %>% # Remove extra spaces from column names
  separate(qacc, into = c("GeneID", "XP_ID"), sep = "\\|") %>% 
  dplyr::rename(Subject_ID = sacc) %>% 
  dplyr::select(GeneID, Subject_ID) %>%
  left_join(Genes_R4 %>% dplyr::select(GeneID, GeneName), by = "GeneID") %>% 
  group_by(GeneID, Subject_ID) %>% unique() %>% ungroup()


#_________

Genes_R5 <- read_tsv("GO_data/R5_extracted_genes_table.txt",
                     col_names = c("chrom","start","end","strand","GeneID","LOC","biotype","GeneName"), quote = "")

Homologs_R5 <- read.table("GO_data/forward_best_hit_homologs_filtered_R5.txt",
                          sep = "\t",header = FALSE, stringsAsFactors = FALSE,
                          col.names = c("qacc", "sacc", "pident", "length", "mismatch","gapopen","qstart", 
                                        "qend", "sstart", "send", "evalue", "bitscore")) %>%
  as_tibble() %>% rename_with(trimws) %>% # Remove extra spaces from column names
  separate(qacc, into = c("GeneID", "XP_ID"), sep = "\\|") %>% 
  dplyr::rename(Subject_ID = sacc) %>% 
  dplyr::select(GeneID, Subject_ID) %>%
  left_join(Genes_R5 %>% dplyr::select(GeneID, GeneName), by = "GeneID") %>% 
  group_by(GeneID, Subject_ID) %>% unique() %>% ungroup()


#_________

Genes_R6 <- read_tsv("GO_data/R6_extracted_genes_table.txt",
                     col_names = c("chrom","start","end","strand","GeneID","LOC","biotype","GeneName"), quote = "")

Homologs_R6 <- read.table("GO_data/forward_best_hit_homologs_filtered_R6.txt",
                          sep = "\t",header = FALSE, stringsAsFactors = FALSE,
                          col.names = c("qacc", "sacc", "pident", "length", "mismatch","gapopen","qstart", 
                                        "qend", "sstart", "send", "evalue", "bitscore")) %>%
  as_tibble() %>% rename_with(trimws) %>% # Remove extra spaces from column names
  separate(qacc, into = c("GeneID", "XP_ID"), sep = "\\|") %>% 
  dplyr::rename(Subject_ID = sacc) %>% 
  dplyr::select(GeneID, Subject_ID) %>%
  left_join(Genes_R6 %>% dplyr::select(GeneID, GeneName), by = "GeneID") %>% 
  group_by(GeneID, Subject_ID) %>% unique() %>% ungroup()


#_____________________Convert all QTL Dmel genes to FBgn Ids ___________________________________________ 

# Note: all these three files have different identifier systems for Drosophila melanogaster genes

### FBgn → FlyBase gene ID (stable, gene-level - best anchor)
### CG numbers → legacy gene symbols (map to FBgn)
### NP_ → RefSeq protein accessions (isoform-level - many-to-one with genes)

#_____________standardize everything to FBgn (FlyBase Gene IDs)________________

# Why FBgn - its Stable and gene-centric (QTL studies are gene-level) and Directly supported by FlyBase
# CG and protein IDs can both map to FBgn, but not reliably the other way around

# First — Merge your two QTL datasets So Convert CG → FBgn, then merge:

# Map Dmel CG → FBgn
cg_numbers <- unique(QTLStudy_genes$Gene_number) # %>% length() #129 CG Ids to be mapped to FBgn

cg_to_fbgn <- AnnotationDbi::select(
  org.Dm.eg.db,
  keys = cg_numbers,
  columns = "FLYBASE",
  keytype = "FLYBASECG") #%>% 
  #rename("Gene_number" = "FLYBASECG")

# Find the CG Ids that could niot be mapped to FBgnIDs in FyBase
mapped <- cg_to_fbgn[!is.na(cg_to_fbgn$FLYBASE), ] #%>% nrow() #127 mapped
unmapped <- cg_to_fbgn[is.na(cg_to_fbgn$FLYBASE), ] #%>% nrow() #2 unmapped

# for the two unmapped IDs Use Ensembl BioMart as backup
#install.packages("BiocManager")
#BiocManager::install("biomaRt")
library(biomaRt)

mart <- useMart("ensembl",dataset = "dmelanogaster_gene_ensembl",host = "https://useast.ensembl.org")

cg_map2 <- getBM(attributes = c("flybase_gene_id", "external_gene_name"),
  filters = "external_gene_name", values = unmapped$FLYBASECG, mart = mart) %>%  # From the two unmapped IDs we get one mapped with Ensemble.
  dplyr::rename(FLYBASECG = external_gene_name, FLYBASE = flybase_gene_id) %>%
  dplyr::select(FLYBASECG, FLYBASE)

# Note we have one still_unmapped: ID "CG43681"

# 4. Combine Dmel fbgn Ids

Dmel_mapped <- rbind(mapped, cg_map2) %>% #nrow() # 128/129 mapped now
  dplyr::rename(Gene_number=FLYBASECG)

QTLStudy_FbgnID <- QTLStudy_genes %>% 
  left_join(Dmel_mapped, by = "Gene_number") %>%
  dplyr::select(Gene_name, FLYBASE) %>% 
  rbind(Fanara_genes)

#_________Convert B.tryoni Dmel homolog NP_ IDs to FBgn IDs for R4___________________________

# NP_ = protein (isoform-level) # FBgn = gene (gene-level), So map NP_ → gene → FBgn

np_ids_R4 <- unique(Homologs_R4$Subject_ID) #258 NP Ids
np_map_R4 <- getBM(attributes = c("flybase_gene_id", "refseq_peptide"),
  filters = "refseq_peptide",values = np_ids_R4,mart = mart) %>%  #nrow() #251 NP Ids mapped to FBgn
  dplyr::rename(FLYBASE = flybase_gene_id,REFSEQ = refseq_peptide) %>%
  filter(!is.na(FLYBASE)) %>%
  distinct() 

total_np_R4 <- length(unique(np_ids_R4)) # 258 Total NP input
mapped_np_R4 <- length(unique(np_map_R4$REFSEQ)) # 250mapped
unmapped_np_ids_R4 <- setdiff(np_ids_R4, np_map_R4$REFSEQ) # 8unmapped ("NP_001356914.1","NP_650875.3","NP_001356887.1",
#NP_001356959.1","NP_001356949.1","NP_001368984.1","NP_001368965.1","NP_523961.2")

# Combine Dmel fbgn Ids
Homologs_R4_FbgnID <- Homologs_R4 %>% dplyr::rename(REFSEQ=Subject_ID) %>% 
  left_join(np_map_R4, by = "REFSEQ") 

# Now find the B. try homologs in R4 that occurs in Dmel QTL studies using FBgn Ids
common_fbgn_R4 <- inner_join(Homologs_R4_FbgnID,QTLStudy_FbgnID,by = "FLYBASE") %>%
  filter(!is.na(FLYBASE)) # 21 overlaps


#______________________________________________

#_________Convert B.tryoni Dmel homolog NP_ IDs to FBgn IDs for R5___________________________

# NP_ = protein (isoform-level) # FBgn = gene (gene-level), So map NP_ → gene → FBgn

np_ids_R5 <- unique(Homologs_R5$Subject_ID) #%>% length() #3768 NP Ids

np_map_R5 <- getBM(attributes = c("flybase_gene_id", "refseq_peptide"),
                   filters = "refseq_peptide",values = np_ids_R5,mart = mart) %>%  #nrow() #3695 NP Ids mapped to FBgn
  dplyr::rename(FLYBASE = flybase_gene_id,REFSEQ = refseq_peptide) %>%
  filter(!is.na(FLYBASE)) %>%
  distinct() 

total_np_R5 <- length(unique(np_ids_R5)) # 3768 Total NP input
mapped_np_R5 <- length(unique(np_map_R5$REFSEQ)) # 3693 mapped
unmapped_np_ids_R5 <- setdiff(np_ids_R5, np_map_R5$REFSEQ) # 75unmapped 

# Combine Dmel fbgn Ids
Homologs_R5_FbgnID <- Homologs_R5 %>% dplyr::rename(REFSEQ=Subject_ID) %>% 
  left_join(np_map_R5, by = "REFSEQ") 

# Now find the B. try homologs in R5 that occurs in Dmel QTL studies using FBgn Ids
common_fbgn_R5 <- inner_join(Homologs_R5_FbgnID,QTLStudy_FbgnID,by = "FLYBASE") %>%
  filter(!is.na(FLYBASE)) # 144 overlaps


#______________________________________________

#_________Convert B.tryoni Dmel homolog NP_ IDs to FBgn IDs for R6___________________________

# NP_ = protein (isoform-level) # FBgn = gene (gene-level), So map NP_ → gene → FBgn

np_ids_R6 <- unique(Homologs_R6$Subject_ID) #%>% length() #59 NP Ids

np_map_R6 <- getBM(attributes = c("flybase_gene_id", "refseq_peptide"),
                   filters = "refseq_peptide",values = np_ids_R6,mart = mart) %>%  #nrow() #57 NP Ids mapped to FBgn
  dplyr::rename(FLYBASE = flybase_gene_id,REFSEQ = refseq_peptide) %>%
  filter(!is.na(FLYBASE)) %>%
  distinct() 

total_np_R6 <- length(unique(np_ids_R6)) # 51 Total NP input
mapped_np_R6 <- length(unique(np_map_R6$REFSEQ)) # 49 mapped
unmapped_np_ids_R6 <- setdiff(np_ids_R6, np_map_R6$REFSEQ) # 2unmapped 

# Combine Dmel fbgn Ids
Homologs_R6_FbgnID <- Homologs_R6 %>% dplyr::rename(REFSEQ=Subject_ID) %>% 
  left_join(np_map_R6, by = "REFSEQ") 

# Now find the B. try homologs in R6 that occurs in Dmel QTL studies using FBgn Ids
common_fbgn_R6 <- inner_join(Homologs_R6_FbgnID,QTLStudy_FbgnID,by = "FLYBASE") %>%
  filter(!is.na(FLYBASE)) # 3 overlaps

#___________________________________________

#______Create summary tables for all rounds  mapped using Biomart_________
# Check duplicates before removing R4 %>% count(GeneID) %>% filter(n > 1)
R4 <- common_fbgn_R4 %>% mutate(Round = "R4") %>% distinct(GeneID, .keep_all = TRUE)
R5 <- common_fbgn_R5 %>% mutate(Round = "R5") %>% distinct(GeneID, .keep_all = TRUE)
R6 <- common_fbgn_R6 %>% mutate(Round = "R6") %>% distinct(GeneID, .keep_all = TRUE)

All_rounds <- bind_rows(R4, R5, R6)

Summary_mapped_Biomart<- All_rounds  %>%
  group_by(GeneID, Gene_name, Round) %>% unique() %>% mutate(present = 1) %>%
  tidyr::pivot_wider(names_from = Round, values_from = present,values_fill = 0)

write.csv(Summary_mapped_Biomart, "Drossie_GWAS/Result/Individual rounds with genes overlapping with Dmel QTL gene_Biomart.csv")


# _______________fbgn IDs for ummapped ones_________________
R4.unmapped <- tibble(Subject_ID = unmapped_np_ids_R4, Round = "R4") 
R5.unmapped <- tibble(Subject_ID = unmapped_np_ids_R5, Round = "R5")
R6.unmapped <- tibble(Subject_ID = unmapped_np_ids_R6, Round = "R6")

Unmapped_all <- rbind(R4.unmapped,R5.unmapped,R6.unmapped)


#_________________________________________________________________________________________________

##############################################################################################

#______________________________________________________________________________________________________



#______________Now repeat above using Rentrex library_________________________________

library(rentrez)

get_fbgn <- function(np_id) {
  tryCatch({
    record <- entrez_fetch(db = "protein",
                           id = np_id,
                           rettype = "gb",
                           retmode = "text")
    
    # Remove line breaks (CRITICAL FIX)
    record <- gsub("\n", " ", record)
    
    # Extract FBgn specifically from FlyBase cross-ref
    match <- regexpr("FLYBASE:FBgn[0-9]+", record, ignore.case = TRUE)
    
    if (match == -1) return(NA)
    
    # Extract just FBgn part
    fbgn_full <- regmatches(record, match)
    sub(".*(FBgn[0-9]+).*", "\\1", fbgn_full)
    
  }, error = function(e) NA)
}

#__________Use Rentez to Convert B.tryoni Dmel homolog NP_ IDs to FBgn IDs for R4___________________________

R4_rentezmapped <- tibble(Subject_ID = np_ids_R4) %>%
  mutate(Fbgn = sapply(Subject_ID, get_fbgn))

R4_rentezmapped %>% summarise(mapped = sum(!is.na(Fbgn)),unmapped = sum(is.na(Fbgn)),total = n()) #all mapped

# Combine Dmel fbgn Ids
Homologs_R4_FbgnID_rentez <- Homologs_R4 %>% dplyr::rename(REFSEQ=Subject_ID) %>% 
  left_join(R4_rentezmapped %>% dplyr::rename(REFSEQ=Subject_ID), by = "REFSEQ") %>% 
  dplyr::rename(FLYBASE=Fbgn)

# Now find the B. try homologs in R4 that occurs in Dmel QTL studies using FBgn Ids
common_fbgn_R4_rentez <- inner_join(Homologs_R4_FbgnID_rentez,QTLStudy_FbgnID,by = "FLYBASE") %>%
  filter(!is.na(FLYBASE)) # 21 overlaps


#__________Use Rentez to Convert B.tryoni Dmel homolog NP_ IDs to FBgn IDs for R5___________________________

R5_rentezmapped <- tibble(Subject_ID = np_ids_R5) %>%
  mutate(Fbgn = sapply(Subject_ID, get_fbgn))

R5_rentezmapped %>% summarise(mapped = sum(!is.na(Fbgn)),unmapped = sum(is.na(Fbgn)),total = n()) #all mapped

# Combine Dmel fbgn Ids
Homologs_R5_FbgnID_rentez <- Homologs_R5 %>% dplyr::rename(REFSEQ=Subject_ID) %>% 
  left_join(R5_rentezmapped %>% dplyr::rename(REFSEQ=Subject_ID), by = "REFSEQ") %>% 
  dplyr::rename(FLYBASE=Fbgn)

# Now find the B. try homologs in R5 that occurs in Dmel QTL studies using FBgn Ids
common_fbgn_R5_rentez <- inner_join(Homologs_R5_FbgnID_rentez,QTLStudy_FbgnID,by = "FLYBASE") %>%
  filter(!is.na(FLYBASE)) # 144 overlaps

#__________Use Rentez to Convert B.tryoni Dmel homolog NP_ IDs to FBgn IDs for R6___________________________

R6_rentezmapped <- tibble(Subject_ID = np_ids_R6) %>%
  mutate(Fbgn = sapply(Subject_ID, get_fbgn))

R6_rentezmapped %>% summarise(mapped = sum(!is.na(Fbgn)),unmapped = sum(is.na(Fbgn)),total = n()) #all mapped

# Combine Dmel fbgn Ids
Homologs_R6_FbgnID_rentez <- Homologs_R6 %>% dplyr::rename(REFSEQ=Subject_ID) %>% 
  left_join(R6_rentezmapped %>% dplyr::rename(REFSEQ=Subject_ID), by = "REFSEQ") %>% 
  dplyr::rename(FLYBASE=Fbgn)

# Now find the B. try homologs in R5 that occurs in Dmel QTL studies using FBgn Ids
common_fbgn_R6_rentez <- inner_join(Homologs_R5_FbgnID_rentez,QTLStudy_FbgnID,by = "FLYBASE") %>%
  filter(!is.na(FLYBASE)) # 3 overlaps


#___________________________________________


#______Create summary tables for all rounds  mapped using Rentez_________

# Check duplicates before removing R4 %>% count(GeneID) %>% filter(n > 1)
R4_rentez <- common_fbgn_R4_rentez %>% mutate(Round = "R4") %>% distinct(GeneID, .keep_all = TRUE)
R5_rentez <- common_fbgn_R5_rentez %>% mutate(Round = "R5") %>% distinct(GeneID, .keep_all = TRUE)
R6_rentez <- common_fbgn_R6_rentez %>% mutate(Round = "R6") %>% distinct(GeneID, .keep_all = TRUE)

All_rounds_rentez <- bind_rows(R4_rentez, R5_rentez, R6_rentez)

Summary_mapped_rentez <- All_rounds_rentez  %>%
  group_by(GeneID, Gene_name, Round) %>% unique() %>% mutate(present = 1) %>%
  tidyr::pivot_wider(names_from = Round, values_from = present,values_fill = 0)

write.csv(Summary_mapped_rentez, "Drossie_GWAS/Result/Individual rounds with genes overlapping with Dmel QTL gene_Rentez.csv")


#____________________________END_____________________________________


