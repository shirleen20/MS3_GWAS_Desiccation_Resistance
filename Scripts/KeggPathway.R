#_____________________________________16/07/2026________________________________________

# Perform GO analysis and Kegg analysis/ ORA -over representation analysis

library(clusterProfiler)
library(GO.db)
library(GSEABase)
library(tidyverse)
library(readxl)
library(data.table)
library(org.Dm.eg.db)
library(readxl)
library(data.table)
library(AnnotationDbi)
library(stringr)

rm(list=ls())

################# extract D.melanogaster GO terms, extract BP, MF and CC term same way, by changing 'ont' and 'level' at line 24 & 25##### #####################

# get entrez to flybase
keytypes(org.Dm.eg.db)
head(keys(org.Dm.eg.db, keytype="FLYBASE"))
k <- keys(org.Dm.eg.db, keytype="FLYBASE")
cols <- c("ENTREZID", "FLYBASECG")
list <- AnnotationDbi::select(org.Dm.eg.db, keys=k, columns=cols, keytype="FLYBASE")
head(list)

#### extract level 2 and 3 BP, CC and MF GO terms ###

entrez_ids <- keys(org.Dm.eg.db, keytype = "ENTREZID")


## 1. Function to extract TERM2NAME for given ont & level

extract_TERM2NAME <- function(ont, level, outdir = "KeggPathway") {
  
  message(sprintf("Processing ontology = %s | level = %d", ont, level))
  
  ggo <- groupGO(
    gene      = entrez_ids,
    OrgDb    = org.Dm.eg.db,
    ont      = ont,
    level    = level,
    readable = FALSE
  )
  
  term2name <- ggo %>%
    filter(Count > 0) %>%
    select(ID, Description)
  
  outfile <- file.path(
    outdir,
    paste0("TERM2NAME_", ont, level, ".csv")
  )
  
  write.csv(term2name, outfile, row.names = FALSE)
}

dir.create("KeggPathway", recursive = TRUE, showWarnings = FALSE)

## ---------------------------------------------------------
## 2. Run all required combinations

extract_TERM2NAME("MF", 2)
extract_TERM2NAME("MF", 3)
extract_TERM2NAME("MF", 4)
extract_TERM2NAME("BP", 2)
extract_TERM2NAME("BP", 3)
extract_TERM2NAME("CC", 2)
extract_TERM2NAME("CC", 3)


################# InterProScan Data ####################

# a. Read the InterProScan TSV file
df <- read_tsv("Geneoverlap/Btry_Prot5000bp_interPro_Clean.tsv", col_names = FALSE, show_col_types = FALSE)

# c. Split column 2 by '|' and remove bracketed annotations
df_clean <- df %>%
  separate_rows(X2, sep = "\\|") %>%
  mutate(
    X2 = str_trim(X2),
    X2 = str_remove(X2, "\\(.*\\)")) %>%
  filter(!str_detect(X2, "-")) %>%
  distinct(X1, X2)

# Step 3: Read XP-LOC index file
colnames(df_clean) <- c("Protein_ID", "GO_ID")
head(df_clean)


###################### Extract Gene and All GO from emapper file #######


###Remove Mitochondrial genes from universe################
MT <- c("COX1", "COX2", "COX3", "ND1", "ND2", "ND3", "ND4", "ND5", "ND6", "ND4L", "CYTB", "ATP6", "ATP8")

Universe_AllProtein <- read.delim("Geneoverlap/universe.txt", header =F) %>% 
  dplyr::rename(Protein_ID=V1, Transcript_ID=V2, GeneID=V3, Gene=V4) %>% 
  dplyr::filter(!Gene %in% MT)

# Step 1: Read Excel as tibble
data <- read_excel("Geneoverlap/out.emapper.annotations.xlsx",sheet=1) %>% data.table() %>% 
  left_join(Universe_AllProtein, by = c("query"="Protein_ID"))

# Step 2: Clean and split GO terms
Gene2emapperGO <- data %>%
  dplyr::select(Gene, GOs) %>%
  filter(!is.na(GOs)) %>%
  rowwise() %>%
  mutate(GO_list = strsplit(GOs, ",")) %>%
  unnest(cols = c(GO_list)) %>%
  filter(!grepl("-", GO_list)) %>%
  dplyr::select(Protein_ID = Gene, GO_ID = GO_list)

Gene2emapperGO %>% write.csv("KeggPathway/Gene2AllGO.csv", row.names = FALSE)

################## Combine Interpro and EggNoc#############

# link protein ids with the gene ids for Interpro data

df_Interpro <- df_clean %>% 
  left_join(Universe_AllProtein, by = "Protein_ID") %>% select(Gene, GO_ID) %>% 
  rename(Protein_ID = Gene)

colnames(df_Interpro)
colnames(Gene2emapperGO)


loc_go_combined <- bind_rows(df_Interpro, Gene2emapperGO) %>% group_by(Protein_ID, GO_ID) %>% distinct()  

unique(df_Interpro$Protein_ID) %>% length()  # 10101
unique(Gene2emapperGO$Protein_ID) %>% length() # 8473
unique(loc_go_combined$Protein_ID) %>% length () # 11039 have at least one GO annotation from either method

intersect(df_Interpro, Gene2emapperGO) # 24,563 Gene IDs present in both eggnoc and interpro
setdiff(df_Interpro, Gene2emapperGO) # 17,288 Gene IDs that have InterPro annotation but do not have eggnoc annotations
setdiff(Gene2emapperGO, df_Interpro) # 525,518 Gene IDs that have GO annotations from emapper but do not have InterPro annotations
union(df_Interpro, Gene2emapperGO) # Gene IDs with InterPro annotations and Emapper annotations # Removes duplicate rows

### Save #####
write.csv(loc_go_combined, "Gene2AllGO_Unique.csv", row.names = FALSE)


#################################################################################
############ Parent to offspring ###################
##################################################################################

###### combine offspring terms to GO level 2 #############

GO_ID <- Gene2emapperGO %>% dplyr::select(GO_ID) %>% distinct()
go_ids <- GO_ID$GO_ID %>% as.character()
# Create GSEAbase GOCollection using `go_ids`
myCollection <- GOCollection(go_ids)

####################### MF2 ##############################

# extract Level 2 GO Terms ##########
slimdf <- read.csv("KeggPathway/TERM2NAME_MF2.csv", header = TRUE) %>% dplyr::select(ID)
# List of GOslims and all GO IDs from `go_ids`
gomap <- as.list(GOMFOFFSPRING[slimdf$ID])
# Maps `go_ids` to matching GOslims
mapped <- lapply(gomap, intersect, ids(myCollection))
# Append all mapped GO IDs to `slimdf`
# `sapply` needed to apply paste() to create semi-colon delimited values
slimdf$ids <- sapply(lapply(gomap, intersect, ids(myCollection)), paste, collapse=";")
# Remove "character(0) string from "ids" column
slimdf$ids[slimdf$ids == "character(0)"] <- ""


#### make a datatable for Parent and Osspring GO terms ###
data1 <- slimdf %>% data.table()
P2O <- data1[, list(GO_ID = unlist(strsplit(ids , ";"))), by = ID]
names(P2O) <- c("MF.GO4", "GO_ID")

#### Need to add unique Parent terms in both columns, as above function miss parent terms ##

UniqueDat <- unique(P2O[,1])
unique_rows <- data.frame(MF.GO4 = UniqueDat, GO_ID = UniqueDat)
P2O1 <- rbind(P2O, unique_rows, use.names=FALSE)

#### Final Term 2 Gene File ##################

T2G <- Gene2emapperGO %>% as.data.frame()
T2G1 <- full_join(T2G, P2O1, by = "GO_ID", relationship = "many-to-many") %>% na.omit() %>% dplyr::select("MF.GO4", "Protein_ID") %>% write.csv("KeggPathway/Term2Gene_MF2.csv", row.names = FALSE)

####################### BP2 ##############################

# extract Level 4 GO Terms ##########
slimdf <- read.csv("KeggPathway/TERM2NAME_BP2.csv", header = TRUE) %>% dplyr::select(ID)
# List of GOslims and all GO IDs from `go_ids`
gomap <- as.list(GOBPOFFSPRING[slimdf$ID])
# Maps `go_ids` to matching GOslims
mapped <- lapply(gomap, intersect, ids(myCollection))
# Append all mapped GO IDs to `slimdf`
# `sapply` needed to apply paste() to create semi-colon delimited values
slimdf$ids <- sapply(lapply(gomap, intersect, ids(myCollection)), paste, collapse=";")
# Remove "character(0) string from "ids" column
slimdf$ids[slimdf$ids == "character(0)"] <- ""


#### make a datatable for Parent and Osspring GO terms ###
data1 <- slimdf %>% data.table()
P2O <- data1[, list(GO_ID = unlist(strsplit(ids , ";"))), by = ID]
names(P2O) <- c("BP.GO4", "GO_ID")

#### Need to add unique Parent terms in both columns, as above function miss parent terms ##

UniqueDat <- unique(P2O[,1])
unique_rows <- data.frame(BP.GO4 = UniqueDat, GO_ID = UniqueDat)
P2O1 <- rbind(P2O, unique_rows, use.names=FALSE)


#### Final Term 2 Gene File ##################

T2G <- Gene2emapperGO  %>% as.data.frame()
T2G1 <- full_join(T2G, P2O1, by = "GO_ID", relationship = "many-to-many") %>% na.omit() %>% dplyr::select("BP.GO4", "Protein_ID") %>% write.csv("KeggPathway/Term2Gene_BP2.csv", row.names = FALSE)

####################### CC2 ##############################

# extract Level 4 GO Terms ##########
slimdf <- read.csv("KeggPathway/TERM2NAME_CC2.csv", header = TRUE) %>% dplyr::select(ID)
# List of GOslims and all GO IDs from `go_ids`
gomap <- as.list(GOCCOFFSPRING[slimdf$ID])
# Maps `go_ids` to matching GOslims
mapped <- lapply(gomap, intersect, ids(myCollection))
# Append all mapped GO IDs to `slimdf`
# `sapply` needed to apply paste() to create semi-colon delimited values
slimdf$ids <- sapply(lapply(gomap, intersect, ids(myCollection)), paste, collapse=";")
# Remove "character(0) string from "ids" column
slimdf$ids[slimdf$ids == "character(0)"] <- ""


#### make a datatable for Parent and Osspring GO terms ###
data1 <- slimdf %>% data.table()
P2O <- data1[, list(GO_ID = unlist(strsplit(ids , ";"))), by = ID]
names(P2O) <- c("CC.GO4", "GO_ID")

#### Need to add unique Parent terms in both columns, as above function miss parent terms ##

UniqueDat <- unique(P2O[,1])
unique_rows <- data.frame(BP.GO4 = UniqueDat, GO_ID = UniqueDat)
P2O1 <- rbind(P2O, unique_rows, use.names=FALSE)


#### Final Term 2 Gene File ##################

T2G <- Gene2emapperGO %>% as.data.frame()
T2G1 <- full_join(T2G, P2O1, by = "GO_ID", relationship = "many-to-many") %>% na.omit() %>% dplyr::select("CC.GO4", "Protein_ID") %>% write.csv("KeggPathway/Term2Gene_CC2.csv", row.names = FALSE)




#################################################################################
############ Parent to offspring ###################
##################################################################################

###### combine offspring terms to GO level 3 #############

GO_ID <- Gene2emapperGO %>% dplyr::select(GO_ID) %>% distinct()
go_ids <- GO_ID$GO_ID %>% as.character()
# Create GSEAbase GOCollection using `go_ids`
myCollection <- GOCollection(go_ids)

####################### MF3 ##############################

# extract Level 3GO Terms ##########
slimdf <- read.csv("KeggPathway/TERM2NAME_MF3.csv", header = TRUE) %>% dplyr::select(ID)
# List of GOslims and all GO IDs from `go_ids`
gomap <- as.list(GOMFOFFSPRING[slimdf$ID])
# Maps `go_ids` to matching GOslims
mapped <- lapply(gomap, intersect, ids(myCollection))
# Append all mapped GO IDs to `slimdf`
# `sapply` needed to apply paste() to create semi-colon delimited values
slimdf$ids <- sapply(lapply(gomap, intersect, ids(myCollection)), paste, collapse=";")
# Remove "character(0) string from "ids" column
slimdf$ids[slimdf$ids == "character(0)"] <- ""


#### make a datatable for Parent and Osspring GO terms ###
data1 <- slimdf %>% data.table()
P2O <- data1[, list(GO_ID = unlist(strsplit(ids , ";"))), by = ID]
names(P2O) <- c("MF.GO4", "GO_ID")

#### Need to add unique Parent terms in both columns, as above function miss parent terms ##

UniqueDat <- unique(P2O[,1])
unique_rows <- data.frame(MF.GO4 = UniqueDat, GO_ID = UniqueDat)
P2O1 <- rbind(P2O, unique_rows, use.names=FALSE)

#### Final Term 2 Gene File ##################

T2G <- Gene2emapperGO %>% as.data.frame()
T2G1 <- full_join(T2G, P2O1, by = "GO_ID", relationship = "many-to-many") %>% na.omit() %>% dplyr::select("MF.GO4", "Protein_ID") %>% write.csv("KeggPathway/Term2Gene_MF3.csv", row.names = FALSE)

####################### BP3 ##############################

# extract Level 3 GO Terms ##########
slimdf <- read.csv("KeggPathway/TERM2NAME_BP3.csv", header = TRUE) %>% dplyr::select(ID)
# List of GOslims and all GO IDs from `go_ids`
gomap <- as.list(GOBPOFFSPRING[slimdf$ID])
# Maps `go_ids` to matching GOslims
mapped <- lapply(gomap, intersect, ids(myCollection))
# Append all mapped GO IDs to `slimdf`
# `sapply` needed to apply paste() to create semi-colon delimited values
slimdf$ids <- sapply(lapply(gomap, intersect, ids(myCollection)), paste, collapse=";")
# Remove "character(0) string from "ids" column
slimdf$ids[slimdf$ids == "character(0)"] <- ""


#### make a datatable for Parent and Osspring GO terms ###
data1 <- slimdf %>% data.table()
P2O <- data1[, list(GO_ID = unlist(strsplit(ids , ";"))), by = ID]
names(P2O) <- c("BP.GO4", "GO_ID")

#### Need to add unique Parent terms in both columns, as above function miss parent terms ##

UniqueDat <- unique(P2O[,1])
unique_rows <- data.frame(BP.GO4 = UniqueDat, GO_ID = UniqueDat)
P2O1 <- rbind(P2O, unique_rows, use.names=FALSE)


#### Final Term 2 Gene File ##################

T2G <- Gene2emapperGO  %>% as.data.frame()
T2G1 <- full_join(T2G, P2O1, by = "GO_ID", relationship = "many-to-many") %>% na.omit() %>% dplyr::select("BP.GO4", "Protein_ID") %>% write.csv("KeggPathway/Term2Gene_BP3.csv", row.names = FALSE)

####################### CC3 ##############################

# extract Level 3 GO Terms ##########
slimdf <- read.csv("KeggPathway/TERM2NAME_CC3.csv", header = TRUE) %>% dplyr::select(ID)
# List of GOslims and all GO IDs from `go_ids`
gomap <- as.list(GOCCOFFSPRING[slimdf$ID])
# Maps `go_ids` to matching GOslims
mapped <- lapply(gomap, intersect, ids(myCollection))
# Append all mapped GO IDs to `slimdf`
# `sapply` needed to apply paste() to create semi-colon delimited values
slimdf$ids <- sapply(lapply(gomap, intersect, ids(myCollection)), paste, collapse=";")
# Remove "character(0) string from "ids" column
slimdf$ids[slimdf$ids == "character(0)"] <- ""


#### make a datatable for Parent and Osspring GO terms ###
data1 <- slimdf %>% data.table()
P2O <- data1[, list(GO_ID = unlist(strsplit(ids , ";"))), by = ID]
names(P2O) <- c("CC.GO4", "GO_ID")

#### Need to add unique Parent terms in both columns, as above function miss parent terms ##

UniqueDat <- unique(P2O[,1])
unique_rows <- data.frame(BP.GO4 = UniqueDat, GO_ID = UniqueDat)
P2O1 <- rbind(P2O, unique_rows, use.names=FALSE)


#### Final Term 2 Gene File ##################

T2G <- Gene2emapperGO %>% as.data.frame()
T2G1 <- full_join(T2G, P2O1, by = "GO_ID", relationship = "many-to-many") %>% na.omit() %>% dplyr::select("CC.GO4", "Protein_ID") %>% write.csv("KeggPathway/Term2Gene_CC3.csv", row.names = FALSE)


###################### MF4 (FIXED PROPERLY) ##############################

slimdf <- read.csv("KeggPathway/TERM2NAME_MF4.csv", header = TRUE) %>%dplyr::select(ID) 

# List of GOslims and all GO IDs from `go_ids`
gomap <- as.list(GOMFOFFSPRING[slimdf$ID])
# Maps `go_ids` to matching GOslims
mapped <- lapply(gomap, intersect, ids(myCollection))
# Append all mapped GO IDs to `slimdf`
# `sapply` needed to apply paste() to create semi-colon delimited values
slimdf$ids <- sapply(lapply(gomap, intersect, ids(myCollection)), paste, collapse=";")
# Remove "character(0) string from "ids" column
slimdf$ids[slimdf$ids == "character(0)"] <- ""


#### make a datatable for Parent and Osspring GO terms ###
data1 <- slimdf %>% data.table()
P2O <- data1[, list(GO_ID = unlist(strsplit(ids , ";"))), by = ID]
names(P2O) <- c("MF.GO4", "GO_ID")


##### Need to add unique Parent terms in both columns, as above function miss parent terms ##

UniqueDat <- unique(P2O[,1])
unique_rows <- data.frame(MF.GO4 = UniqueDat, GO_ID = UniqueDat)
P2O1 <- rbind(P2O, unique_rows, use.names=FALSE)

#### Final Term 2 Gene File ##################

T2G <- Gene2emapperGO %>% as.data.frame()
T2G1 <- full_join(T2G, P2O1, by = "GO_ID", relationship = "many-to-many") %>% 
  na.omit() %>% dplyr::select("MF.GO4", "Protein_ID") %>% 
  write.csv("KeggPathway/Term2Gene_MF4.csv", row.names = FALSE)



########################################################
##################### KEGG #############################

kegg_clean <- data %>%
  dplyr::select(Gene, KEGG_Pathway) %>%
  filter(!is.na(KEGG_Pathway)) %>%
  rowwise() %>%
  mutate(KEGG_list = strsplit(KEGG_Pathway, ",")) %>%
  unnest(cols = c(KEGG_list)) %>%
  mutate(KEGG_list = str_trim(KEGG_list)) %>%
  filter(!grepl("^map", KEGG_list, ignore.case = TRUE),  # remove map IDs
         !grepl("-", KEGG_list)) %>%                     # remove weird malformed entries
  dplyr::select(KEGG_ID = KEGG_list, Protein_ID = Gene)



kegg_clean %>% write.csv("KeggPathway/Term2Gene_KEGG.csv", row.names = FALSE)

#### extract KO terms and names #######
R.utils::setOption("clusterProfiler.download.method",'auto')
ko.pathways <- download_KEGG(species="ko")
term2name <- ko.pathways$KEGGPATHID2NAME
term2name %>% write.csv("KeggPathway/TERM2NAME_KEGG.csv", row.names = FALSE)

################ remove "human disease and drug development terms"####################

# Read your KEGG mapping tables
term2gene <- read.csv("KeggPathway/Term2Gene_KEGG.csv")
term2name <- read.csv("KeggPathway/TERM2NAME_KEGG.csv")

# Define unwanted patterns (disease and drug-related)
exclude_prefixes <- c("ko050", "ko051", "ko052", "ko053", "ko070", "ko071", "ko072")
exclude_specific <- c("ko04930","ko04940","ko04950","ko04931",
                      "ko04932","ko04933","ko04934","ko04936")

# Filter TERM2NAME
term2name_filtered <- term2name %>%
  filter(
    !str_sub(from, 1, 5) %in% exclude_prefixes,
    !from %in% exclude_specific)

# Filter TERM2GENE using same KEGG_ID list
term2gene_filtered <- term2gene %>%
  filter(KEGG_ID %in% term2name_filtered$from)

# Save filtered tables
write.csv(term2name_filtered, "KeggPathway/TERM2NAME_KEGG_filtered.csv", row.names = FALSE)
write.csv(term2gene_filtered, "KeggPathway/Term2Gene_KEGG_filtered.csv", row.names = FALSE)


##################################### Enzyme ######################################################


#  Extract and clean ECs from EggNOG/InterPro data
kegg1 <- data %>%
  select(Gene, EC) %>%
  filter(!is.na(EC), EC != "-") %>%
  mutate(Enzymes = strsplit(EC, ",")) %>%
  unnest(cols = c(Enzymes)) %>%
  mutate(
    Enzymes = str_trim(Enzymes),
    Enzymes = na_if(Enzymes, "-"),
    # Extract only first two levels: e.g. 3.1.1.3 → 3.1
    Enzymes = str_replace(Enzymes, "^([0-9]+\\.[0-9]+)\\..*$", "\\1"),
    # Pad the second digit with 0 if it has only one digit (e.g. 3.1 → 3.01)
    Enzymes = if_else(
      str_detect(Enzymes, "^\\d+\\.\\d$"),
      str_replace(Enzymes, "(\\d+)\\.(\\d)$", "\\1.0\\2"),
      Enzymes)) %>%
  filter(!is.na(Enzymes), Enzymes != "") %>%
  distinct(Gene, Enzymes)

#  Save final file
write.csv(kegg1, "KeggPathway/Term2Gene_Enzymes.csv", row.names = FALSE)

#### "Term2Name_Enzyme.csv" is prepared manually from enzyme website#####


#________________________________________________________________________________________________


##________________ORA analysis___________

# This function ORA 

run_ORA <- function(diff_gene_file, universe_file, term_pairs, output_dir = "ORA_results", suffix = "") {
  # diff_gene_file : CSV with column 'gene'
  # universe_file  : TXT file with gene universe
  # term_pairs     : named list of TERM2GENE and TERM2NAME file paths
  # output_dir     : base output directory
  # suffix         : string to append to file names or subfolder names
  
  # Create output folder with suffix if provided
  if (suffix != "") {
    output_dir <- file.path(output_dir, paste0("ORA_", suffix))
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Load input gene sets
  UP_T <- read.delim(diff_gene_file, header = FALSE)
  Universe <- read.delim(universe_file, header = FALSE)
  
  # Iterate through each term set
  for (term_name in names(term_pairs)) {
    term_files <- term_pairs[[term_name]]
    
    TERM2GENE <- read.csv(term_files$TERM2GENE)
    TERM2NAME <- read.csv(term_files$TERM2NAME)
    
 # Run enrichment
enrich_result <- enricher(
      gene = UP_T[[2]],
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      universe = Universe[[4]],
      minGSSize = 3,            # minGSSize = is the minimum number of genes required in a GO category/pathway before testing enrichment.
      maxGSSize = 5000,         # maximum allowed size of a gene set (GO term/pathway) to be tested so Ignore GO categories containing more than 5000 genes.
      qvalueCutoff = 0.05,
      TERM2GENE = TERM2GENE,
      TERM2NAME = TERM2NAME
    )
    
# Save results
out_file <- file.path(output_dir, paste0("ORA_", term_name, "_", suffix, ".csv"))
    
enrich_result %>% as.data.frame() %>% write.csv(out_file, row.names = FALSE)
    
    message("✅ Saved: ", out_file)
  }
}


# Example usage:


term_files_list <- list(
  MF4  = list(TERM2GENE = "KeggPathway/Term2Gene_MF4.csv",
              TERM2NAME = "KeggPathway/TERM2NAME_MF4.csv"),
  
  MF3  = list(TERM2GENE = "KeggPathway/Term2Gene_MF3.csv",
              TERM2NAME = "KeggPathway/TERM2NAME_MF3.csv"),
  
  MF2  = list(TERM2GENE = "KeggPathway/Term2Gene_MF2.csv",
              TERM2NAME = "KeggPathway/TERM2NAME_MF2.csv"),
  
  BP3  = list(TERM2GENE = "KeggPathway/Term2Gene_BP3.csv",
              TERM2NAME = "KeggPathway/TERM2NAME_BP3.csv"),
  
  BP2  = list(TERM2GENE = "KeggPathway/Term2Gene_BP2.csv",
              TERM2NAME = "KeggPathway/TERM2NAME_BP2.csv"),
  
  CC3  = list(TERM2GENE = "KeggPathway/Term2Gene_CC3.csv",
              TERM2NAME = "KeggPathway/TERM2NAME_CC3.csv"),
  
  CC2  = list(TERM2GENE = "KeggPathway/Term2Gene_CC2.csv",
              TERM2NAME = "KeggPathway/TERM2NAME_CC2.csv"),
  
  KEGG = list(TERM2GENE = "KeggPathway/Term2Gene_KEGG_filtered.csv",
              TERM2NAME = "KeggPathway/TERM2NAME_KEGG_filtered.csv"),
  
  Enzyme = list(TERM2GENE = "KeggPathway/Term2Gene_Enzymes.csv",
                TERM2NAME = "KeggPathway/Term2Name_Enzyme.csv"))

##############################

# Creat data files for ORA_Run with no MT:


universe <- read.delim("Geneoverlap/universe.txt",header = FALSE)%>% 
  dplyr::rename(Protein_ID=V1, Transcript_ID=V2, GeneID=V3, Gene=V4) %>% 
  dplyr::filter(!Gene %in% MT) %>% 
  write.table("Geneoverlap/universe_noMT.txt",sep="\t", quote=FALSE,row.names = FALSE, col.names=FALSE)
  

R4 <- read.delim("Geneoverlap/Test_proteins_snpEffR4_identified_no-intergenic.txt",header = FALSE) %>% 
  dplyr::filter(!V2 %in% MT)%>% 
  write.table("Geneoverlap/Test_proteins_R4_noMT.txt",sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE)


R5 <- read.delim("Geneoverlap/Test_proteins_snpEffR5_identified_no-intergenic.txt",header = FALSE) %>% 
  dplyr::filter(!V2 %in% MT)%>% 
  write.table("Geneoverlap/Test_proteins_R5_noMT.txt",sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE)


R6 <- read.delim("Geneoverlap/Test_proteins_snpEffR6_identified_no-intergenic.txt",header = FALSE) %>% 
  dplyr::filter(!V2 %in% MT)%>% 
  write.table("Geneoverlap/Test_proteins_R6_noMT.txt",sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE)



################### R4 #####################

# Process GO for R4 at fdr_0.01

run_ORA(
  diff_gene_file = "Geneoverlap/Test_proteins_R4_noMT.txt",
  universe_file  = "Geneoverlap/universe_noMT.txt",
  term_pairs     = term_files_list,
  output_dir     = "KeggPathway/ORA_results",
  suffix         = "R4")


################### R5 #####################

# Process GO for R5 at fdr_0.01

run_ORA(
  diff_gene_file = "Geneoverlap/Test_proteins_R5_noMT.txt",
  universe_file  = "Geneoverlap/universe_noMT.txt",
  term_pairs     = term_files_list,
  output_dir     = "KeggPathway/ORA_results",
  suffix         = "R5")


################### R6 #####################

# Process GO for R6 at fdr_0.01


run_ORA(
  diff_gene_file = "Geneoverlap/Test_proteins_R6_noMT.txt",
  universe_file  = "Geneoverlap/universe.txt",
  term_pairs     = term_files_list,
  output_dir     = "KeggPathway/ORA_results",
  suffix         = "R6")



################### Overlap #####################

# Process GO for overlap at fdr_0.01


run_ORA(
  diff_gene_file = "Geneoverlap/Test_proteins_snpEff_Overlap_with_Drossie.txt",
  universe_file  = "Geneoverlap/universe_noMT.txt",
  term_pairs     = term_files_list,
  output_dir     = "KeggPathway/ORA_results",
  suffix         = "Overlap")

#____________________________END_______________________________________________________________



