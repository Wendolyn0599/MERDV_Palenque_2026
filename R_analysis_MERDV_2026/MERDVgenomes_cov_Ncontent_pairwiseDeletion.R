# MERDV Palenque 2026
# Script 2: MERDV_genomes_cov_Ncontent_pairwiseDeletion.R
# Description: Analysis of coverage and indeterminate nucleotide (N) content 
#              in the generated MERDV genomes. 
#              Generation of an identity matrix using pairwise deletion based on 
#              the alignment of these genomes and their alignment with  public 
#              sequences of MERDV, MERDV-L, CXRV, other representative  
#              Merhaviruses and two outgroups.
# Author: Wendolyn Estrada
# Last update: August 2026

# Set working directory (use here::here() for reproducibility)
# If using RStudio, ensure the project is open.
# Otherwise, set the path manually:
# setwd("/path/to/your/project")

################################################################################

# 0. Load required libraries====================================================
library(dplyr)         # Data manipulation
library(stringr)       # String manipulation
library(writexl)       # Export files
library(Biostrings)    # Import DNA sequences
library(ape)           # Analyze and manipulate biological sequences

# 1. DEFINITION OF FUNCTIONS====================================================

#Function to evaluate sequence quality (N content) and genome coverage__________
#  Receives: multi-FASTA file with consensus genomes generated from variantConsePairCall_List_v2.sh
#  Returns: a table with sequence name, N count, N%, covered bases (non-N), and % coverage

calculate_cov_N <- function(fasta){
  #Read file
  seqs <- readDNAStringSet(fasta)
  #Verify number of sequences
  num_seqs <- length(seqs)  
  #Create quality table
  df_quality <- data.frame(
    Sequence_name = names(seqs),
    Length = width(seqs),
    N_count = sapply(seqs, function(x) sum(alphabetFrequency(x)[["N"]])),
    stringsAsFactors = FALSE)
  #Add N count, coverage and percentages
  df_quality <- df_quality %>% mutate(
    N_perc = round((N_count / Length) * 100, 2),
    Nucleotides_covered = Length - N_count,
    Cov_perc = round((Nucleotides_covered / Length) * 100, 2),
  )
  return(df_quality)
}

#Function to generate pairwise deletion identity matrix_________________________
#Receives: Aligment in FASTA format
#Return: Identity matrix
pairwise_deletion <- function(alignment) {
  align <- read.dna(alignment, format = "fasta")
  
  # Add ambiguity = "NA" to exclude sites with N in either of the two sequences
  dist_matrix <- dist.dna(align, model = "raw", pairwise.deletion = TRUE) 
  #Calculate identity matrix
  identity_matrix <- (1 - as.matrix(dist_matrix)) * 100
  diag(identity_matrix) <- 100  
  
  return(identity_matrix)
}

#2. MAIN________________________________________________________________________

#Read consensus genomes and calculate N content and coverture
cons_file <-  "../DATA/Sequences_MERDV_PalenqueChiapas.fasta"
merdv_cons_cov <- calculate_cov_N(cons_file)
#Export file
write_xlsx(merdv_cons_cov, "../OUTPUTS/MERDV_consensus_cov.xlsx")

#Read alignment and generate identity matrix
alignment <- "../DATA/TotalSeq_noCD_aligned.fasta"
ident_matrix <- pairwise_deletion(alignment)
#Export file
write.csv(ident_matrix, "../OUTPUTS/supplementary_tableS3.csv")