#_____________________________________29/07/2026________________________________________

# Perform GO analysis and Kegg analysis/ ORA -over representation analysis FOR

# Genes found in any rounds and overlapping that overlaps with Dmel QTL study and are homologs/orthologs

#____________________________________________________________________
# OVERLAP: B. tryoni genes overlapping Drossie QTL → Dmel GO ORA

library(tidyverse)
library(readxl)
library(data.table)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Dm.eg.db)
library(GO.db)
library(GSEABase)

rm(list = ls())

#____________________________________________________________________

# Create data files for ORA_Run with no MT:

R4 <- read.delim("Geneoverlap/Test_proteins_snpEffR4_identified_no-intergenic.txt",header = FALSE) %>% 
  dplyr::rename(Protein_ID=V1, GeneID=V2, Gene=V3)

R5 <- read.delim("Geneoverlap/Test_proteins_snpEffR5_identified_no-intergenic.txt",header = FALSE) %>% 
  dplyr::rename(Protein_ID=V1, GeneID=V2, Gene=V3)

R6 <- read.delim("Geneoverlap/Test_proteins_snpEffR6_identified_no-intergenic.txt",header = FALSE) %>% 
  dplyr::rename(Protein_ID=V1, GeneID=V2, Gene=V3)


# homologs of Drossie QTL studies

Overlap <- read_excel("GO_data/Qfly_homologs_overlaps_Drossie_QTL.xlsx") %>% 
  dplyr::rename(GeneID=1) 

Overlap_genes <- rbind(R4, R5, R6) %>% dplyr::filter(GeneID %in% Overlap$GeneID)

write.table(Overlap_genes,"GO_data/Test_proteins_snpEff_Overlap_with_Drossie.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

#________________________


# 1. LOAD OVERLAP GENE SET (B. tryoni → Dmel mapping)

overlap <- read.delim("GO_data/Test_proteins_snpEff_Overlap_with_Drossie.txt",header = TRUE)

DmelIDs <- read_excel("GO_data/Qfly_homologs_with_Dmel_geneID.xlsx")

overlap_dm <- overlap %>% left_join(DmelIDs, by = "GeneID") %>% dplyr::select(GeneID, FLYBASE) %>% distinct()

#____________________________________________________________________
# 2. MAP FBgn → ENTREZ (Dmel annotation space)

fbgn_map <- AnnotationDbi::select(
  org.Dm.eg.db,
  keys = unique(overlap_dm$FLYBASE),
  keytype = "FLYBASE",
  columns = "ENTREZID") %>%
  filter(!is.na(ENTREZID)) %>%
  distinct()

# Overlap gene set in ENTREZ
overlap_entrez <- fbgn_map %>%
  distinct(ENTREZID)

write.table(overlap_entrez,"GO_data/Overlap_ENTREZ.txt",
  row.names = FALSE,col.names = FALSE,quote = FALSE)


# 3. DEFINE BACKGROUND (Dmel universe)

universe_entrez <- keys(org.Dm.eg.db, keytype = "ENTREZID")

write.table(universe_entrez,"GO_data/Universe_ENTREZ.txt",
  row.names = FALSE,col.names = FALSE,quote = FALSE)


# 4. BUILD Dmel GO ANNOTATION (TERM2GENE)

dmel_go <- AnnotationDbi::select(org.Dm.eg.db,
  keys = universe_entrez,keytype = "ENTREZID",columns = "GO") %>% na.omit()

TERM2GENE_overlap <- dmel_go %>% 
  transmute(GO_ID = GO,Protein_ID = ENTREZID) %>% #transmute() is mutate() + select() in one step
  distinct()


# 5. LOAD GO SLIM TABLES (LEVEL 2 & 3)
# Helper function to load any GO slim file
load_goslim <- function(ont, level) {
  
  file <- paste0("KeggPathway/TERM2NAME_", ont, level, ".csv")
  
  df <- read.csv(file, stringsAsFactors = FALSE) %>%
    dplyr::select(ID, Description)
  
  return(df)
}


# Load all GO slim tables

GO_slim_MF2 <- load_goslim("MF", 2)
GO_slim_MF3 <- load_goslim("MF", 3)
GO_slim_MF3 <- load_goslim("MF", 4)

GO_slim_BP2 <- load_goslim("BP", 2)
GO_slim_BP3 <- load_goslim("BP", 3)

GO_slim_CC2 <- load_goslim("CC", 2)
GO_slim_CC3 <- load_goslim("CC", 3)


# 6. BUILD GO LEVEL MAPPING FUNCTION

generate_overlap_GO <- function(term2gene, ont, level, offspring_obj, label) {
  
  go_ids <- unique(term2gene$GO_ID)
  myCollection <- GOCollection(go_ids)
  
  slimdf <- read.csv(
    paste0("KeggPathway/TERM2NAME_", ont, level, ".csv")
  ) %>%
    dplyr::select(ID)
  
  gomap <- as.list(offspring_obj[slimdf$ID])
  
  slimdf$ids <- sapply(
    lapply(gomap, intersect, ids(myCollection)),
    paste,
    collapse = ";"
  )
  
  slimdf$ids[slimdf$ids == "character(0)"] <- NA
  
  P2O <- data.table(
    ID = slimdf$ID,
    GO_ID = unlist(strsplit(slimdf$ids, ";"))
  )
  
  setnames(P2O, c(label, "GO_ID"))
  
  # include parent terms
  parents <- unique(P2O[[label]])
  
  P2O <- rbind(
    P2O,
    data.frame(
      setNames(list(parents, parents), c(label, "GO_ID"))
    )
  )
  
  out <- term2gene %>%
    inner_join(P2O, by = "GO_ID") %>%
    dplyr::select(all_of(label), Protein_ID)
  
  outfile <- paste0(
    "KeggPathway/Term2Gene_overlap_",
    ont, level, ".csv"
  )
  
  write.csv(out, outfile, row.names = FALSE)
  
  message("Saved: ", outfile)
}


# 7. RUN GO SLIM MAPPING (LEVEL 2 & 3)

generate_overlap_GO(TERM2GENE_overlap, "MF", 2, GOMFOFFSPRING, "MF.GO2")
generate_overlap_GO(TERM2GENE_overlap, "MF", 3, GOMFOFFSPRING, "MF.GO3")
generate_overlap_GO(TERM2GENE_overlap, "MF", 4, GOMFOFFSPRING, "MF.GO4")


generate_overlap_GO(TERM2GENE_overlap, "BP", 2, GOBPOFFSPRING, "BP.GO2")
generate_overlap_GO(TERM2GENE_overlap, "BP", 3, GOBPOFFSPRING, "BP.GO3")

generate_overlap_GO(TERM2GENE_overlap, "CC", 2, GOCCOFFSPRING, "CC.GO2")
generate_overlap_GO(TERM2GENE_overlap, "CC", 3, GOCCOFFSPRING, "CC.GO3")


# 8. ORA FUNCTION (OVERLAP ONLY)

run_ORA_Overlap <- function(diff_gene_file, universe_file, term_pairs,
                            output_dir = "ORA_results", suffix = "") {
  
  dir.create(file.path(output_dir, suffix),
             recursive = TRUE,
             showWarnings = FALSE)
  
  gene_list <- scan(diff_gene_file, what = "")
  universe_list <- scan(universe_file, what = "")
  
  for (term_name in names(term_pairs)) {
    
    files <- term_pairs[[term_name]]
    
    TERM2GENE <- read.csv(files$TERM2GENE)
    TERM2NAME <- read.csv(files$TERM2NAME)
    
    enrich <- tryCatch({
      enricher(
        gene = gene_list,
        universe = universe_list,
        TERM2GENE = TERM2GENE,
        TERM2NAME = TERM2NAME,
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        qvalueCutoff = 0.05,
        minGSSize = 3,
        maxGSSize = 5000
      )
    }, error = function(e) NULL)
    
    if (is.null(enrich)) next
    
    write.csv(
      as.data.frame(enrich),
      file.path(output_dir, suffix,
                paste0("ORA_", term_name, "_", suffix, ".csv")),
      row.names = FALSE
    )
    
    message("Saved ORA: ", term_name)
  }
}

#____________________________________________________________________
# 9. DEFINE TERM FILES
#____________________________________________________________________

overlap_term_files <- list(
  MF2 = list(
    TERM2GENE = "KeggPathway/Term2Gene_overlap_MF2.csv",
    TERM2NAME = "KeggPathway/TERM2NAME_MF2.csv"
  ),
  MF3 = list(
    TERM2GENE = "KeggPathway/Term2Gene_overlap_MF3.csv",
    TERM2NAME = "KeggPathway/TERM2NAME_MF3.csv"
  ),
  MF4 = list(
    TERM2GENE = "KeggPathway/Term2Gene_overlap_MF4.csv",
    TERM2NAME = "KeggPathway/TERM2NAME_MF4.csv"
  ),
  BP2 = list(
    TERM2GENE = "KeggPathway/Term2Gene_overlap_BP2.csv",
    TERM2NAME = "KeggPathway/TERM2NAME_BP2.csv"
  ),
  BP3 = list(
    TERM2GENE = "KeggPathway/Term2Gene_overlap_BP3.csv",
    TERM2NAME = "KeggPathway/TERM2NAME_BP3.csv"
  ),
  CC2 = list(
    TERM2GENE = "KeggPathway/Term2Gene_overlap_CC2.csv",
    TERM2NAME = "KeggPathway/TERM2NAME_CC2.csv"
  ),
  CC3 = list(
    TERM2GENE = "KeggPathway/Term2Gene_overlap_CC3.csv",
    TERM2NAME = "KeggPathway/TERM2NAME_CC3.csv"
  )
)


# 10. RUN FINAL ORA (OVERLAP ONLY)

run_ORA_Overlap(
  diff_gene_file = "GO_data/Overlap_ENTREZ.txt",
  universe_file  = "GO_data/Universe_ENTREZ.txt",
  term_pairs     = overlap_term_files,
  output_dir     = "KeggPathway/ORA_results",
  suffix         = "ORA_homolgs_Dmel_bg"
)


#______________________________END_______________________

