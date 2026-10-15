#__________________________09/07/2026______________________________

library(tidyverse)
library(readxl)
library(data.table)
library(openxlsx)


rm(list=ls())

#____________________Make Table 1 common genes__________________

# To match later with gene list to find orthologous protein accession
Universe_AllProtein <- read.delim("Geneoverlap/universe.txt", header =F) %>% 
  rename(Protein_ID=V1, Transcript_ID=V2, GeneID=V3, Gene=V4)

data <- read_excel("Geneoverlap/out.emapper.annotations.xlsx",sheet=1) %>% data.table() %>% 
  left_join(Universe_AllProtein, by = c("query"="Protein_ID"))


# Read files for gene overlap 
# These dataset contains the number of significant SNPs (V1) within 5000 bp of genes (V4),
# grouped by SNP effects (V3)
R4_allSigGene <- read.delim("Geneoverlap/R4_OverlapGene.txt", header=F) %>% mutate(Round="R4")
R5_allSigGene <- read.delim("Geneoverlap/R5_OverlapGene.txt", header=F) %>% mutate(Round="R5")
R6_allSigGene <- read.delim("Geneoverlap/R6_OverlapGene.txt", header=F) %>% mutate(Round="R6")

# Read files for SNP overlap (all the significant SNPs - FDR < 0.05)
# These dataset contains all SNP positions that are within 5000 bp of genes (SNP position can be repeated if multiple genes are close by)
R4_allSigSNP <- read.delim("Geneoverlap/R4_SNP_Gene_forOverlap.txt", header=F) %>% mutate(Round="R4", Present = 1)
R5_allSigSNP <- read.delim("Geneoverlap/R5_SNP_Gene_forOverlap.txt", header=F) %>% mutate(Round="R5", Present = 1)
R6_allSigSNP <- read.delim("Geneoverlap/R6_SNP_Gene_forOverlap.txt", header=F) %>% mutate(Round="R6", Present = 1)

# Scaffold sizes
ScafSize <- read.table("Geneoverlap/Btry_masked.fasta.fai", header = FALSE) %>% 
  dplyr::select(V1,V2) %>% rename(CHR = V1, Size = V2) %>% 
  mutate(acc_type = case_when(str_starts(CHR, "NC_") ~ "NC",
                              str_starts(CHR, "NW_") ~ "NW",TRUE~ "OTHER"),# classify accessions
         original_idx = row_number()) %>%  #preserve original order index
  arrange(factor(acc_type, levels = c("NC", "NW", "OTHER")), # for NC and OTHER we keep original_idx, for NW we sort by -V2
          if_else(acc_type == "NW", -Size, original_idx)) %>%
  mutate(order_idx = row_number()) %>%
  arrange(order_idx) %>% dplyr::select(CHR, Size, order_idx)


#______Some helper functions to help with processing the data__________
## Capturing the highest effect category (from the 15 effects present in our dataset)

RankCategories <- function(EffectVector){
  if ("stop_gained" %in% EffectVector){ return("stop_gained") } else 
    if ("missense_variant&splice_region_variant" %in% EffectVector) { return("missense_variant&splice_region_variant") } else
      if ("missense_variant" %in% EffectVector){ return("missense_variant") } else
        if ("5_prime_UTR_premature_start_codon_gain_variant" %in% EffectVector){ return("5_prime_UTR_premature_start_codon_gain_variant") } else
          if ("splice_region_variant" %in% EffectVector) { return("splice_region_variant") } else
            if ("splice_region_variant&intron_variant" %in% EffectVector) { return("splice_region_variant&intron_variant") } else
              if ("splice_region_variant&non_coding_transcript_exon_variant" %in% EffectVector){ return("splice_region_variant&non_coding_transcript_exon_variant") } else
                if ("synonymous_variant" %in% EffectVector){ return("synonymous_variant") } else
                  if ("5_prime_UTR_variant" %in% EffectVector){ return("5_prime_UTR_variant") } else
                    if ("3_prime_UTR_variant" %in% EffectVector){ return("3_prime_UTR_variant") } else
                      if ("non_coding_transcript_exon_variant" %in% EffectVector){ return("non_coding_transcript_exon_variant") } else
                        if ("intron_variant" %in% EffectVector){ return("intron_variant") } else
                          if ("intragenic_variant" %in% EffectVector){ return("intragenic_variant") } else
                            if ("upstream_gene_variant" %in% EffectVector){ return("upstream_gene_variant") } else
                            { return("downstream_gene_variant")}
}

# To index the effects by rankings to aid with sorting
IndexCategories <- function(EffectVector){
  if ("TotalGene" %in% EffectVector) { return("00") } else
    if ("stop_gained" %in% EffectVector){ return("01") } else 
      if ("missense_variant&splice_region_variant" %in% EffectVector) { return("02") } else
        if ("missense_variant" %in% EffectVector){ return("03") } else
          if ("5_prime_UTR_premature_start_codon_gain_variant" %in% EffectVector){ return("04") } else
            if ("splice_region_variant" %in% EffectVector) { return("05") } else
              if ("splice_region_variant&intron_variant" %in% EffectVector) { return("06") } else
                if ("splice_region_variant&non_coding_transcript_exon_variant" %in% EffectVector){ return("07") } else
                  if ("synonymous_variant" %in% EffectVector){ return("08") } else
                    if ("5_prime_UTR_variant" %in% EffectVector){ return("09") } else
                      if ("3_prime_UTR_variant" %in% EffectVector){ return("10") } else
                        if ("non_coding_transcript_exon_variant" %in% EffectVector){ return("11") } else
                          if ("intron_variant" %in% EffectVector){ return("12") } else
                            if ("intragenic_variant" %in% EffectVector){ return("13") } else
                              if ("upstream_gene_variant" %in% EffectVector){ return("14") } else
                              { return("15")}
}

# To rename NCBI refseq accession number to chromosome/special scaffold name as per genome paper (in prep)
Bt_chr <- function(Scaf){
  Renamed <- ifelse(Scaf == "NC_014611.1", "MT",
                    ifelse(Scaf == "NC_052499.1", "Chr2",
                           ifelse(Scaf == "NC_052500.1", "Chr3",
                                  ifelse(Scaf == "NC_052501.1", "Chr4",
                                         ifelse(Scaf == "NC_052502.1", "Chr5",
                                                ifelse(Scaf == "NC_052503.1", "Chr6",
                                                       ifelse(Scaf == "NW_024395824.1", "sc11",
                                                              ifelse(Scaf == "NW_024396366.1", "sc07",
                                                                     ifelse(Scaf == "NW_024395977.1", "sc25", Scaf)))))))))
  return(Renamed)
}


###############################
# Read files for best peaks for each round

R4_rawSigSNP <- read.delim("Geneoverlap/Bestpeaks_FDR0.05_R4.tsv", header=T) %>% mutate(Present=1, Round = Rep)
R5_rawSigSNP <- read.delim("Geneoverlap/Bestpeaks_FDR0.05_R5.tsv", header=T) %>% mutate(Present=1, Round = Rep)
R6_rawSigSNP <- read.delim("Geneoverlap/Bestpeaks_FDR0.05_R6.tsv", header=T) %>% mutate(Present=1, Round = Rep)

### Do the same for the SNP counts
Table1_RawCounts <- rbind(R4_rawSigSNP, R5_rawSigSNP, R6_rawSigSNP) %>% 
  unique() %>% 
  pivot_wider(id_cols = c(CHR,POS), names_from = Round, values_from = Present, values_fill = 0 ) %>% 
  mutate(Round = case_when(R4==1 & R5==1 & R6==1 ~ "Fs4+5+6",
                           R4==1 & R5==1 & R6==0 ~ "Fs4+5",
                           R4==1 & R5==0 & R6==1 ~ "Fs4+6",
                           R4==0 & R5==1 & R6==1 ~ "Fs5+6",
                           R4==1 & R5==0 & R6==0 ~ "F4",
                           R4==0 & R5==1 & R6==0 ~ "F5",
                           R4==0 & R5==0 & R6==1 ~ "F6")) %>%
  left_join(ScafSize, by = c("CHR" = "CHR")) %>% 
  mutate(CHR = Bt_chr(CHR)) %>% 
  group_by(CHR, Round, Size, order_idx) %>% 
  summarise(TotalSNPCount = length(unique(POS)))

SNPRawCount <- Table1_RawCounts %>% 
  pivot_wider(id_cols = c(CHR,Size, order_idx), names_from = Round, values_from = TotalSNPCount, values_fill = 0) %>% 
  relocate(`Fs5+6`,F4, F5, F6, .after = `Fs4+6`)

SNPRawCountHeat <- SNPRawCount %>% 
  pivot_longer(cols = `Fs4+5`:F6, names_to = "Round", values_to = "Count") %>% 
  mutate(logDensity =  if_else(Count > 0, log10(Count/Size), -8)) %>%  # chose  -8 to improve contrast
  mutate(Count = if_else(Count == 0, "", as.character(Count))) %>% 
  mutate(Round = factor(Round, levels = c("Fs4+5+6", "Fs4+5", "Fs4+6", "Fs5+6", "F4", "F5", "F6"))) %>% ungroup() %>% 
  arrange(desc(order_idx)) %>% 
  mutate(NewChr = if_else(grepl("NW", CHR), paste0("sc",sprintf("%05.0f", Size/100), " ", gsub(" ", " ", sprintf("%7.3f", Size/1000000))) , 
                          if_else(CHR == "MT", paste0("MT ", gsub(" ", " ", sprintf("%7.3f", Size/1000000))),
                                  paste0(CHR, " ", gsub(" ", " ", paste0(sprintf("%7.3f", Size/1000000))) )))) %>% 
  mutate(NewChr = factor(NewChr, levels = unique(NewChr)))

# heatmap showing Sig SNPs in big scafoold - main text figure:

SNPRawplot_big_sc <- SNPRawCountHeat %>% 
  dplyr::filter(CHR %in% c("Chr2","Chr3","Chr4","Chr5","Chr6","MT","sc11","sc07","sc25")) %>% 
  ggplot(aes(x = Round, y = NewChr, fill = logDensity, label = Count)) + 
  geom_tile() + 
  geom_label(label.size = NA, fill = NA, size = 2) +
  scale_fill_gradient(low = "white", high = "red") +
  labs(title = "ggplot2 Heatmap") +
  theme_minimal() +
  theme(plot.title = element_blank(), axis.title.x = element_blank(), 
        axis.title.y = element_blank(), legend.position = "none", axis.text=element_text(size=6))


#ggsave(plot = SNPRawplot_big_sc, width = 3.0, height = 1.5, units = "in", dpi = 300,filename = "Figures/Figure 3 showing overlapping Raw SNPs in the three QTL generations.jpg")

# Now make supplementary figure showing Sig SNPs in all scaffolds

SNPRawCountHeat_supp <- SNPRawCount %>% 
  pivot_longer(cols = `Fs4+5`:F6, names_to = "Round", values_to = "Count") %>% 
  mutate(logDensity =  if_else(Count > 0, log10(Count/Size), -8)) %>%  # chose  -8 to improve contrast
  mutate(Count = if_else(Count == 0, "", as.character(Count))) %>% 
  mutate(Round = factor(Round, levels = c("Fs4+5+6", "Fs4+5", "Fs4+6", "Fs5+6", "F4", "F5", "F6"))) %>% ungroup() %>% 
  arrange(desc(order_idx)) %>% 
  mutate(NewChr = Bt_chr(CHR), NewChr = if_else(CHR %in% c("NC_014611.1","NC_052499.1","NC_052500.1","NC_052501.1","NC_052502.1","NC_052503.1",
                 "NW_024395824.1","NW_024396366.1","NW_024395977.1"),
      paste0(NewChr, " ", sprintf("%.3f", Size/1000000)),paste0(CHR, " ", sprintf("%.3f", Size/1000000)))) %>% 
  mutate(NewChr = factor(NewChr, levels = unique(NewChr)))

SNPRawplotsupp <- ggplot(SNPRawCountHeat_supp, aes(x = Round, y = NewChr, fill = logDensity, label = Count)) + 
  geom_tile() + 
  geom_label(label.size = NA, fill = NA, size = 2) +
  scale_fill_gradient(low = "white", high = "red") +
  labs(title = "ggplot2 Heatmap") +
  theme_minimal() +
  theme(plot.title = element_blank(), axis.title.x = element_blank(), 
        axis.title.y = element_blank(), legend.position = "none", axis.text=element_text(size=6))


# ggsave(plot = SNPRawplot, width = 6.0, height = 10.0, units = "in", dpi = 300,filename = "Figures/Figure S2 showing all significant Raw SNPs across all scaffolds.jpg")


#_____plot SNPRawplotsupp above is too messy to read so split it into multiple pages with 60 scaffolds each____

SNPRawCountHeat_supp <- SNPRawCount %>% 
  pivot_longer(cols = `Fs4+5`:F6, names_to = "Round", values_to = "Count") %>% 
  mutate(logDensity =  if_else(Count > 0, log10(Count/Size), -8)) %>%  # chose  -8 to improve contrast
  mutate(Count = if_else(Count == 0, "", as.character(Count))) %>% 
  mutate(Round = factor(Round, levels = c("Fs4+5+6", "Fs4+5", "Fs4+6", "Fs5+6", "F4", "F5", "F6"))) %>% 
  ungroup() %>% 
  mutate(BtName = Bt_chr(CHR),NewChr = if_else(CHR %in% c("NC_014611.1","NC_052499.1","NC_052500.1",  # create scaffold labels
                 "NC_052501.1","NC_052502.1","NC_052503.1","NW_024395824.1","NW_024396366.1","NW_024395977.1"),
      paste0(BtName, " ", sprintf("%.3f", Size/1000000)),
      paste0(CHR, " ", sprintf("%.3f", Size/1000000))))

# create desired scaffold order 
main_chr <- c("Chr2","Chr3","Chr4","Chr5","Chr6","MT","sc25", "sc11","sc07")

scaffold_order <- SNPRawCountHeat_supp %>%
  distinct(CHR, NewChr, Size) %>%
  mutate(BtName = Bt_chr(CHR)) %>%
  arrange(match(BtName, main_chr, nomatch = 999),desc(Size)) %>%
  pull(NewChr)

# apply order
SNPRawCountHeat_supp <- SNPRawCountHeat_supp %>%
  mutate(NewChr = factor(NewChr,levels = scaffold_order))

# create pages of 50 scaffolds
scaffolds_per_page <- 50

scaffold_levels <- levels(SNPRawCountHeat_supp$NewChr)

scaffold_pages <- split(scaffold_levels,ceiling(seq_along(scaffold_levels) / scaffolds_per_page))


# plot each page

for(i in seq_along(scaffold_pages)) {
  
  plot_data <- SNPRawCountHeat_supp %>%
    filter(NewChr %in% scaffold_pages[[i]]) %>%
    mutate(NewChr = factor(NewChr,levels = rev(scaffold_pages[[i]])))
  
p <- ggplot(plot_data,aes(x = Round,y = NewChr,fill = logDensity,label = Count)) +
    geom_tile() + geom_text(size = 2, face = "bold") +
    scale_fill_gradient(low = "white",high = "red") +
    theme_minimal() +
    theme(axis.title = element_blank(),legend.position = "none", axis.text = element_text(size = 6, face = "bold"))
  
# Keep tile size consistent between pages
n_scaffolds <- length(unique(plot_data$NewChr))
plot_height <- n_scaffolds * 0.18 + 1
  
ggsave(filename = paste0("Figures/SNP_heatmap_all_scaffolds_",i,".jpg"),
    plot = p,width = 6,height = plot_height,units = "in",dpi = 300)}

###################################



#____________________Making Table X+1_________________

# Capture gene numbers (with significant SNPs within 5000 bp) by round combination and chromosome/scaffolds
## Total SNP count is 9877 and 4067 genes

Table1_Counts <- rbind(R4_allSigSNP, R5_allSigSNP, R6_allSigSNP) %>% dplyr::select(-V3) %>% 
  unique() %>% 
  pivot_wider(id_cols = c(V1,V2,V4), names_from = Round, values_from = Present, values_fill = 0 ) %>% 
  mutate(Round = case_when(R4==1 & R5==1 & R6==1 ~ "Fs4+5+6",
                           R4==1 & R5==1 & R6==0 ~ "Fs4+5",
                           R4==1 & R5==0 & R6==1 ~ "Fs4+6",
                           R4==0 & R5==1 & R6==1 ~ "Fs5+6",
                           R4==1 & R5==0 & R6==0 ~ "F4",
                           R4==0 & R5==1 & R6==0 ~ "F5",
                           R4==0 & R5==0 & R6==1 ~ "F6")) %>%
  left_join(ScafSize, by = c("V1" = "CHR")) %>% 
  mutate(CHR = Bt_chr(V1)) %>% 
  group_by(CHR, Round, Size, order_idx) %>% 
  summarise(TotalSNPCount = length(unique(V2)))

Table1_GeneCount <- rbind(R4_allSigSNP, R5_allSigSNP, R6_allSigSNP) %>% dplyr::select(-V2,-V3) %>% 
  unique() %>% 
  pivot_wider(id_cols = c(V1,V4), names_from = Round, values_from = Present, values_fill = 0 ) %>% 
  mutate(Round = case_when(R4==1 & R5==1 & R6==1 ~ "Fs4+5+6",
                           R4==1 & R5==1 & R6==0 ~ "Fs4+5",
                           R4==1 & R5==0 & R6==1 ~ "Fs4+6",
                           R4==0 & R5==1 & R6==1 ~ "Fs5+6",
                           R4==1 & R5==0 & R6==0 ~ "F4",
                           R4==0 & R5==1 & R6==0 ~ "F5",
                           R4==0 & R5==0 & R6==1 ~ "F6")) %>%
  left_join(ScafSize, by = c("V1" = "CHR")) %>%  
  mutate(CHR = Bt_chr(V1)) %>% 
  group_by(CHR, Round, Size, order_idx) %>% 
  summarise(TotalGeneCount = length(unique(V4)))

Table1_Final <- Table1_GeneCount %>% full_join(Table1_Counts) %>% 
  replace_na(list(TotalGeneCount = 0, TotalSNPCount = 0)) %>% 
  mutate(Count = paste0(TotalGeneCount, ", ", TotalSNPCount)) %>% 
  pivot_wider(id_cols = CHR, names_from = Round, values_from = Count, values_fill = "0, 0") %>% 
  relocate(`Fs4+5+6`, `Fs4+5`, `Fs4+6`, `Fs5+6`, .after = CHR)

SNPCount <- Table1_Counts %>% 
  pivot_wider(id_cols = c(CHR,Size, order_idx), names_from = Round, values_from = TotalSNPCount, values_fill = 0) %>% 
  relocate(`Fs5+6`,F4, F5, F6, .after = `Fs4+6`)

GeneCount <- Table1_GeneCount %>% 
  pivot_wider(id_cols = c(CHR,Size, order_idx), names_from = Round, values_from = TotalGeneCount, values_fill = 0) %>% 
  relocate(order_idx,`Fs4+5+6`, `Fs4+5`, `Fs4+6`, `Fs5+6`, .after = Size)

######## To fix the NewChr

GeneCountHeat <- GeneCount %>% 
  pivot_longer(cols = `Fs4+5+6`:F6, names_to = "Round", values_to = "Count") %>% 
  mutate(logDensity =  if_else(Count > 0, log10(Count/Size), -8)) %>%  # chose  -8 to improve contrast
  mutate(Count = if_else(Count == 0, "", as.character(Count))) %>% 
  mutate(Round = factor(Round, levels = c("Fs4+5+6", "Fs4+5", "Fs4+6", "Fs5+6", "F4", "F5", "F6"))) %>% ungroup() %>% 
  arrange(desc(order_idx)) %>% 
  mutate(BtName = Bt_chr(CHR),NewChr = if_else(CHR %in% c("NC_014611.1","NC_052499.1","NC_052500.1",  # create scaffold labels
                                                          "NC_052501.1","NC_052502.1","NC_052503.1","NW_024395824.1","NW_024396366.1","NW_024395977.1"),
                                               paste0(BtName, " ", sprintf("%.3f", Size/1000000)),
                                               paste0(CHR, " ", sprintf("%.3f", Size/1000000))))

# Desired chromosome/scaffold order
main_chr <- c("Chr2", "Chr3", "Chr4", "Chr5", "Chr6","MT", "sc25", "sc11", "sc07")

# Create ordered scaffold list
scaffold_order <- GeneCountHeat %>%
  distinct(BtName, NewChr, Size) %>%
  arrange(match(BtName, main_chr, nomatch = 999),desc(Size)) %>%
  pull(NewChr)

# Apply order to y-axis
GeneCountHeat <- GeneCountHeat %>%mutate(NewChr = factor(NewChr, levels = rev(scaffold_order)))

# Heatmap
Geneplot <- ggplot(GeneCountHeat,aes(x = Round,y = NewChr,fill = logDensity,label = Count)) + 
  geom_tile() + geom_label(label.size = NA,fill = NA,size = 2) +
  scale_fill_gradient(low = "white",high = "red") +
  labs(title = "ggplot2 Heatmap") +theme_minimal() +
  theme(plot.title = element_blank(),axis.title.x = element_blank(), axis.title.y = element_blank(),
    legend.position = "none",axis.text = element_text(size = 6))



ggsave(plot = Geneplot, width = 4.0, height = 5.0, units = "in", dpi = 300,filename = "Figures/Figure S3 showing number of genes with atleast one SNP on any scaffold.jpg")


#________________________Make Table X+2_____________

############## Supplementary table with all significant SNPs and their directionality
Alleles_all <- read.delim("Geneoverlap/All_AF_allsigSNPs.txt", header=T) %>% 
  dplyr::select(1,2,3,5,4,6) %>% 
  mutate(across(c(QTLF4_High.QTLF4_Low2, F5_High.QTLF5_Low2, QTLF6_High2.QTLF6_Low2),~ map_chr(str_split(.x, ",\\s*"),
                                                                                               ~ paste(signif(as.numeric(.x), 3), collapse = ", "))))

###
All_sigSNPs <- rbind(read.table("Geneoverlap/R4_FDR0.05.tsv", header = TRUE) %>% 
                       dplyr::select(1,2,3,5) %>% mutate(Round = "R4") %>% rename(FDR = fet_fdr_H2vsL2_R4),
                     read.table("Geneoverlap/R5_FDR0.05.tsv", header = TRUE) %>% 
                       dplyr::select(1,2,3,5) %>% mutate(Round = "R5") %>% rename(FDR = fet_fdr_H2H1vsL2_R5),
                     read.table("Geneoverlap/R6_FDR0.05.tsv", header = TRUE) %>% 
                       dplyr::select(1,2,3,5) %>% mutate(Round = "R6") %>% rename(FDR = fet_fdr_H2vsL2_R6)) %>% 
  pivot_wider(id_cols = CHR:REF, names_from = Round, values_from = FDR) %>% 
  left_join(Alleles_all) %>% 
  mutate(R4 = if_else(is.na(R4), 1, R4),
         R5 = if_else(is.na(R5), 1, R5),
         R6 = if_else(is.na(R6), 1, R6)) %>% 
  separate(REF, into = c("REF","Major","Minor"), sep =",") %>%
  mutate_at(vars(matches("High")), ~gsub(".*,(.*)$", "\\1",.)) %>% 
  mutate_at(vars(matches("High")), ~as.numeric(if_else(.=="99", NA, .))) %>%# 99 means missing in Pesto
  dplyr::rename(R4_FDR = R4, R5_FDR = R5, R6_FDR = R6,
                R4_diff = QTLF4_High.QTLF4_Low2, R5_diff = F5_High.QTLF5_Low2, 
                R6_diff = QTLF6_High2.QTLF6_Low2) %>% 
  mutate(CHR = Bt_chr(CHR)) %>%
  mutate(Sig = paste0(if_else(R4_FDR < 0.05, "R4",""), 
                      if_else(R5_FDR < 0.05, ",R5",""),
                      if_else(R6_FDR < 0.05, ",R6","")) %>% gsub("^,","",.) ) %>% 
  mutate(R4vsR5 = if_else(sign(R4_diff) == sign(R5_diff), "Same", "Opposite"),
         R4vsR6 = if_else(sign(R4_diff) == sign(R6_diff), "Same", "Opposite"),
         R5vsR6 = if_else(sign(R5_diff) == sign(R6_diff), "Same", "Opposite"),
         Direction_sig = if_else(R4_FDR < 0.05 & R5_FDR < 0.05, R4vsR5,
                                 if_else(R4_FDR < 0.05 & R6_FDR < 0.05, R4vsR6, 
                                         if_else(R5_FDR < 0.05 & R6_FDR < 0.05, R5vsR6, "None"))))

# save this table to be used as supplementary dataset 2
write.csv(All_sigSNPs, "Results/Supplementary_dataset_All_sigSNPs.csv", row.names = FALSE)


#_____________________________________________________________________________

###  This is for all SNPs that are close to genes (within 5000 bp)
Alleles <- read.delim("Geneoverlap/All_AF.txt", header=T) %>% 
  dplyr::select(1,2,3,5,4,6) %>% 
  mutate(across(c(QTLF4_High.QTLF4_Low2, F5_High.QTLF5_Low2, QTLF6_High2.QTLF6_Low2),~ map_chr(str_split(.x, ",\\s*"),
                                                                                               ~ paste(signif(as.numeric(.x), 3), collapse = ", ")))) 
Allele_directions <- Alleles %>% 
  mutate_at(vars(matches("High")), ~gsub(".*,(.*)$", "\\1",.)) %>% 
  mutate_at(vars(matches("High")), ~as.numeric(if_else(.=="99", NA, .))) %>%# 99 means missing in Pesto
  dplyr::rename(R4_diff = QTLF4_High.QTLF4_Low2, R5_diff = F5_High.QTLF5_Low2, 
                R6_diff = QTLF6_High2.QTLF6_Low2) %>% 
  mutate(CHR = Bt_chr(CHR))

All_sigSNPs_Genes <- read.delim("Geneoverlap/All_SNP_association.txt", header=T) %>%
  dplyr::select(1,2,3,7,8,9) %>%
  left_join(Alleles) %>%
  separate(REF, into = c("REF","Major","Minor"), sep =",") %>%
  mutate_at(vars(matches("High")), ~gsub(".*,(.*)$", "\\1",.)) %>% 
  mutate_at(vars(matches("High")), ~as.numeric(if_else(.=="99", NA, .))) %>%# 99 means missing in Pesto
  dplyr::rename(R4_FDR = fet_fdr_H2vsL2_R4, R5_FDR = fet_fdr_H2H1vsL2_R5,
                R6_FDR = fet_fdr_H2vsL2_R6,
                R4_diff = QTLF4_High.QTLF4_Low2, R5_diff = F5_High.QTLF5_Low2, 
                R6_diff = QTLF6_High2.QTLF6_Low2) %>% 
  mutate(CHR = Bt_chr(CHR)) %>%
  mutate(Sig = paste0(if_else(R4_FDR < 0.05, "R4",""), 
                      if_else(R5_FDR < 0.05, ",R5",""),
                      if_else(R6_FDR < 0.05, ",R6","")) %>% gsub("^,","",.) )


### 27 genes that had significant SNPs within 5000 bp in all rounds (not necessarily shared SNPs)
common <- rbind(R4_allSigGene, R5_allSigGene, R6_allSigGene) %>% 
  group_by(V4, Round) %>% 
  mutate(HighestEffect = RankCategories(V3)) %>%
  ungroup() %>%  
  dplyr::select(V2, V4, Round, HighestEffect) %>% 
  unique() %>% 
  pivot_wider(id_cols = c(V2,V4), names_from = Round, values_from = HighestEffect, values_fill = "") %>% 
  mutate(Overlap = (R4!="") + (R5!="") + (R6!="")) %>% 
  mutate(RoundOverlap = paste0(if_else(R4!="", "R4, ", ""), if_else(R5!="", "R5, ", ""), if_else(R6!="", "R6", "")) %>% 
           sub(", $", "", .)) %>% 
  filter(Overlap == 3) %>% 
  mutate(V2 = Bt_chr(V2)) %>% 
  pull(V4)


#### Common genes table
OverlapSNP_complex <- rbind(R4_allSigSNP, R5_allSigSNP, R6_allSigSNP) %>% 
  group_by(V1,V2) %>% 
  mutate(EffectGenes = paste0(V4, " (", V3, ")")) %>% 
  mutate(AllAffectedGenes = paste0(EffectGenes, collapse = ", ")) %>% 
  dplyr::select(-c(V3,V4,EffectGenes)) %>%
  unique() %>% 
  pivot_wider(id_cols = c(V1,V2,AllAffectedGenes), names_from = Round, values_from = Present, values_fill = 0 ) %>% 
  mutate(OverlapSNPs = R4 + R5 + R6) %>% 
  filter(grepl(paste(common, collapse = "|"), AllAffectedGenes)) %>% 
  dplyr::rename("CHR" = V1, "POS" = V2) %>% 
  mutate(CHR = Bt_chr(CHR)) 

# All the SNPs
CommonGenesTable_OverlapSNPsOnly <- OverlapSNP_complex %>% 
  mutate(CHR = Bt_chr(CHR)) %>% 
  left_join(Allele_directions) %>% 
  #dplyr::filter(OverlapSNPs > 1) %>% 
  mutate(R4vsR5 = if_else(sign(R4_diff) == sign(R5_diff), "Same", "Opposite"),
         R4vsR6 = if_else(sign(R4_diff) == sign(R6_diff), "Same", "Opposite"),
         R5vsR6 = if_else(sign(R5_diff) == sign(R6_diff), "Same", "Opposite"),
         Direction_sig = if_else(R4 == 1 & R5 == 1, R4vsR5,
                                 if_else(R4 == 1 & R6 == 1, R4vsR6, 
                                         if_else(R5 == 1 & R6 == 1, R5vsR6, "None")))) %>% 
  separate_rows(AllAffectedGenes, sep = ", ") %>% 
  filter(grepl(paste(common, collapse = "|"), AllAffectedGenes)) %>% 
  unique() %>% 
  separate(AllAffectedGenes, into = c("Gene", "Effect"), sep = " \\(") %>% 
  mutate(Effect = sub("\\)","",Effect),
         OverlapRound = paste0(if_else(R4==1, "R4, ", ""), if_else(R5==1, "R5, ", ""), if_else(R6==1, "R6", "")) %>% 
           sub(", $", "", .) ) %>% 
  group_by(CHR, POS, Gene) %>% ### To collapse the SNP with multiple effect for one gene
  mutate(Effect_combined = paste(unique(Effect), collapse = ", ")) %>% 
  dplyr::select(-R4,-R5,-R6,-OverlapSNPs) %>% ungroup() 


##############################
###############################
# Filter all the sig SNPs from CommonGenesTable_OverlapSNPsOnly and make a supplementary file 2 for MS3

SuppFile <- CommonGenesTable_OverlapSNPsOnly %>% 
  dplyr::select(CHR,POS,Gene,Effect,REF,R4_diff,R5_diff,R6_diff,OverlapRound) %>%
  dplyr::rename("SignificantRound" = OverlapRound) %>% 
  separate(REF, into = c("REF","Major","Minor"), sep =",") %>%
  mutate(SignificantRound = factor(SignificantRound,levels = c("R4","R5","R6","R4, R5","R4, R6","R5, R6","R4, R5, R6"))) %>% 
  arrange(SignificantRound) %>%
  #write.csv("Results/Supplementary_dataset_3.csv", row.names = FALSE)

###############################
##############################


### Summary of significant SNPs in the 27 genes
CommonSNPs_27genes <- CommonGenesTable_OverlapSNPsOnly %>% 
  dplyr::select(CHR,POS, Gene, Effect_combined, OverlapRound, Direction_sig) %>% unique() %>% 
  dplyr::select(-POS) %>% 
  group_by(CHR, Gene, OverlapRound, Effect_combined) %>% 
  summarise(Count = n()) %>% 
  group_by(CHR, Gene, Effect_combined) %>% 
  dplyr::mutate(Total = sum(Count)) %>% 
  pivot_wider(id_cols = c(CHR, Gene, Effect_combined, Total), names_from = OverlapRound, values_from = Count, values_fill = 0) %>% 
  dplyr::mutate(R4 = R4 + `R4, R5` + `R4, R5, R6`,  R5 = R5 + `R4, R5` + `R5, R6` + `R4, R5, R6`, R6 = R6 + `R5, R6` + `R4, R5, R6`) %>% # we removed "R4,R6 because it doenst appear"
  dplyr::mutate(`Number of SNPs` = paste0(Total, " (", R4, ", ", R5, ", ", R6,")")) %>% 
  dplyr::select(-c(`R4, R5`, `R5, R6`,`R4, R5, R6`, Total, R4, R5, R6)) %>% ungroup() %>% 
  group_by(Gene) %>% 
  dplyr::mutate(HighestEffect = RankCategories(Effect_combined))

# Overlapping SNPs in the same direction only  
CommonSNPs_HighestImpact <- CommonGenesTable_OverlapSNPsOnly %>% 
  dplyr::select(CHR, Gene, Effect, OverlapRound, Direction_sig) %>% 
  group_by(CHR, Gene, Effect, OverlapRound, Direction_sig) %>%
  summarise(Count = n()) %>% 
  group_by(Gene, OverlapRound) %>% 
  mutate(HighestEffect = RankCategories(Effect)) %>% 
  dplyr::filter(Direction_sig == "Same") %>% dplyr::select(-Direction_sig) %>% 
  #mutate(SummaryCount = paste(Count, collapse = ", ")) %>% 
  pivot_wider(id_cols = c(CHR, Gene, Effect, HighestEffect), names_from = OverlapRound, values_from = Count) %>% 
  mutate(HighestImpact = if_else(Effect == HighestEffect, "Yes", "-")) %>% 
  dplyr::select(-HighestEffect) %>% relocate(HighestImpact, .after = Effect)

###############################################
## Final Table 2?
CommonSNPs_27genesFinal <- CommonGenesTable_OverlapSNPsOnly %>% 
  dplyr::select(CHR,POS, Gene, Effect_combined, OverlapRound, Direction_sig) %>% unique() %>% 
  dplyr::select(-POS) %>% 
  group_by(Gene) %>% 
  mutate(HighestEffect = RankCategories(Effect_combined),
         HighestEffect_overlap = RankCategories(Effect_combined[grepl(",", OverlapRound) & Direction_sig == "Same"])) %>% 
  group_by(CHR, Gene, OverlapRound, HighestEffect, HighestEffect_overlap) %>% 
  summarise(Count = n(), Count_same = length(Direction_sig[Direction_sig == "Same"])) %>% 
  group_by(CHR, Gene) %>% 
  mutate(Total = sum(Count)) %>% ungroup() %>%  
  pivot_wider(id_cols = c(Gene, CHR, HighestEffect, HighestEffect_overlap, Total), names_from = OverlapRound, values_from = c(Count, Count_same), values_fill = 0) %>% 
  mutate(R4 = Count_R4 + `Count_R4, R5` + `Count_R4, R5, R6`,  R5 = Count_R5 + `Count_R4, R5` + `Count_R5, R6` + `Count_R4, R5, R6`,  R6 = Count_R6 + `Count_R4, R5, R6` + `Count_R5, R6`) %>% 
  mutate(`Number of SNPs` = paste0(Total, " (", R4, ", ", R5, ", ", R6,")"),
         `Total_shared` = `Count_same_R4, R5` + `Count_same_R4, R5, R6` + `Count_same_R5, R6`) %>% 
  mutate(`Number of shared SNPs` = paste0(Total_shared, " (", `Count_same_R4, R5`, ", ", `Count_same_R4, R5, R6`, ", ", `Count_same_R5, R6`,")")) %>% 
  mutate(HighestEffect_overlap = if_else(Total_shared == 0, "-", HighestEffect_overlap)) %>% 
  relocate(`Number of shared SNPs`, .before = HighestEffect_overlap) %>% 
  relocate(`Number of SNPs`, .before = HighestEffect) %>% 
  dplyr::select(-c(Total:Total_shared)) #%>% 
  #write.csv("Results/Genes_and_SNPs_Common_in_all_three_rounds.csv", row.names = FALSE)


#_________________________________________________________


################## Supplementary Overlap SNP (same direction), for genes that overlap
################ in at least 2 rounds
##############
### 187 genes that had significant SNPs within 5000 bp in 2 or more rounds (not necessarily shared SNPs)

common2plus <- rbind(R4_allSigGene, R5_allSigGene, R6_allSigGene) %>% 
  group_by(V4, Round) %>% 
  mutate(HighestEffect = RankCategories(V3)) %>%
  ungroup() %>%  
  dplyr::select(V2, V4, Round, HighestEffect) %>% 
  unique() %>% 
  pivot_wider(id_cols = c(V2,V4), names_from = Round, values_from = HighestEffect, values_fill = "") %>% 
  mutate(Overlap = (R4!="") + (R5!="") + (R6!="")) %>% 
  mutate(RoundOverlap = paste0(if_else(R4!="", "R4, ", ""), if_else(R5!="", "R5, ", ""), if_else(R6!="", "R6", "")) %>% 
           sub(", $", "", .)) %>% 
  filter(Overlap >= 2) %>% 
  pull(V4) 

OverlapSNP_complex2plus <- rbind(R4_allSigSNP, R5_allSigSNP, R6_allSigSNP) %>% 
  group_by(V1,V2) %>% 
  mutate(EffectGenes = paste0(V4, " (", V3, ")")) %>% 
  mutate(AllAffectedGenes = paste0(EffectGenes, collapse = ", ")) %>% 
  dplyr::select(-c(V3,V4,EffectGenes)) %>%
  unique() %>% 
  pivot_wider(id_cols = c(V1,V2,AllAffectedGenes), names_from = Round, values_from = Present, values_fill = 0 ) %>% 
  mutate(OverlapSNPs = R4 + R5 + R6) %>% 
  filter(grepl(paste(common2plus, collapse = "|"), AllAffectedGenes)) %>% 
  dplyr::rename("CHR" = V1, "POS" = V2)


CommonGenesTable_OverlapSNPsOnly_2plus <- OverlapSNP_complex2plus %>% 
  mutate(CHR = Bt_chr(CHR)) %>% 
  left_join(Allele_directions) %>% 
  #dplyr::filter(OverlapSNPs > 1) %>% 
  mutate(R4vsR5 = if_else(sign(R4_diff) == sign(R5_diff), "Same", "Opposite"),
         R4vsR6 = if_else(sign(R4_diff) == sign(R6_diff), "Same", "Opposite"),
         R5vsR6 = if_else(sign(R5_diff) == sign(R6_diff), "Same", "Opposite"),
         Direction_sig = if_else(R4 == 1 & R5 == 1, R4vsR5,
                                 if_else(R4 == 1 & R6 == 1, R4vsR6, 
                                         if_else(R5 == 1 & R6 == 1, R5vsR6, "None")))) %>% 
  separate_rows(AllAffectedGenes, sep = ", ") %>% 
  filter(grepl(paste(common2plus, collapse = "|"), AllAffectedGenes)) %>% 
  unique() %>% 
  separate(AllAffectedGenes, into = c("Gene", "Effect"), sep = " \\(") %>% 
  mutate(Effect = sub("\\)","",Effect),
         OverlapRound = paste0(if_else(R4==1, "R4, ", ""), if_else(R5==1, "R5, ", ""), if_else(R6==1, "R6", "")) %>% 
           sub(", $", "", .) ) %>% 
  group_by(CHR, POS, Gene) %>% ### To collapse the SNP with multiple effect for one gene
  mutate(Effect_combined = paste(unique(Effect), collapse = ", ")) %>% 
  dplyr::select(-R4,-R5,-R6,-OverlapSNPs) %>% ungroup() %>% 
  dplyr::filter(grepl("R4.*R5|R5.*R4", OverlapRound) |grepl("R4.*R6|R6.*R4", 
                                                            OverlapRound) |grepl("R5.*R6|R6.*R5", OverlapRound)) #%>% 
  #write.csv("Results/Genes_and_SNPs_in_two_or_more_rounds.csv", row.names = FALSE)


# Overlapping SNPs in the same direction only  
CommonSNPs_HighestImpact_2plus <- CommonGenesTable_OverlapSNPsOnly_2plus %>% 
  dplyr::select(CHR, Gene, Effect, OverlapRound, Direction_sig) %>% 
  group_by(CHR, Gene, Effect, OverlapRound, Direction_sig) %>%
  summarise(Count = n()) %>% 
  group_by(Gene, OverlapRound) %>% 
  mutate(HighestEffect = RankCategories(Effect)) %>% 
  dplyr::filter(Direction_sig == "Same") %>% dplyr::select(-Direction_sig) %>% 
  #mutate(SummaryCount = paste(Count, collapse = ", ")) %>% 
  pivot_wider(id_cols = c(CHR, Gene, Effect, HighestEffect), names_from = OverlapRound, values_from = Count) %>% 
  mutate(HighestImpact = if_else(Effect == HighestEffect, "Yes", "-")) %>% 
  dplyr::select(-HighestEffect) %>% relocate(HighestImpact, .after = Effect)


#CommonSNPs_187genesFinal <- CommonGenesTable_OverlapSNPsOnly_2plus %>% 
  #dplyr::select(CHR,POS, Gene, Effect_combined, OverlapRound, Direction_sig) %>% unique() %>% 
  #dplyr::select(-POS) %>% 
  #group_by(Gene) %>% 
  #mutate(HighestEffect = RankCategories(Effect_combined),
        # HighestEffect_overlap = RankCategories(Effect_combined[grepl(",", OverlapRound) & Direction_sig == "Same"])) %>% 
  #group_by(CHR, Gene, OverlapRound, HighestEffect, HighestEffect_overlap) %>% 
  #summarise(Count = n(), Count_same = length(Direction_sig[Direction_sig == "Same"])) %>% 
  #group_by(CHR, Gene) %>% 
  #mutate(Total = sum(Count)) %>% ungroup() %>%  
  #pivot_wider(id_cols = c(Gene, CHR, HighestEffect, HighestEffect_overlap, Total), names_from = OverlapRound, values_from = c(Count, Count_same), values_fill = 0) %>% 
  #mutate(R4 = Count_R4 + `Count_R4, R5` + `Count_R4, R6`,  R5 = Count_R5 + `Count_R4, R5` + `Count_R5, R6`,  R6 = Count_R6 + `Count_R4, R6` + `Count_R5, R6`) %>% 
  #mutate(`Number of SNPs` = paste0(Total, " (", R4, ", ", R5, ", ", R6,")"),
  #       `Total_shared` = `Count_same_R4, R5` + `Count_same_R4, R6` + `Count_same_R5, R6`) %>% 
  #mutate(`Number of shared SNPs` = paste0(Total_shared, " (", `Count_same_R4, R5`, ", ", `Count_same_R4, R6`, ", ", `Count_same_R5, R6`,")")) %>% 
  #mutate(HighestEffect_overlap = if_else(Total_shared == 0, "-", HighestEffect_overlap)) %>% 
  #relocate(`Number of shared SNPs`, .before = HighestEffect_overlap) %>% 
  #relocate(`Number of SNPs`, .before = HighestEffect) %>% 
  #dplyr::select(-c(Total:Total_shared)) #%>% 
#write.csv("Results/CommonSNPs_209genes_2_or_more_rounds.csv", row.names = FALSE)

#_____________________________________________________________________________________________


