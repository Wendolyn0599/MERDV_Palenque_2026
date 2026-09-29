#!/bin/bash
#Script for multiple sequence alignment and preliminary phylogeny in IQTree
#Before running, sequences of interest must be concatenated into a single file (e.g., cat file1.fasta file2.fasta file3.fasta > conc_files.fasta)
#Arguments
## 1) Path to the FASTA file containing the sequences to be analyzed
## 2) File name
## 3) Path for output files

#$ -N MultAlin_MafftEnboss
#$ -j y
#$ -o /scratch/adrianace/MERDVphilogeny/outputs/
#$ -l h_vmem=8G
#$ -pe thread 4

set -e
source $HOME/.bashrc
module load programs/emboss-6.6.0
module load programs/mafft-7.407
module load programs/iqtree-1.6.12

if [ $# -eq 3 ]; then
 
  cd "$1"
  mkdir -p "$3"

  echo "1. ------------------------------Multiple alignment with MAFFT-----------------------------------------"
  echo "File : ${2} "
  if mafft --genafpair --maxiterate 1000 --thread 4 "$2" > "${3}/${2%.fasta}_aligned.fasta"; then
    echo "MAFFT...OK"
  else
    echo "MAFFT...ERROR"
    exit 1
  fi
  
  
 
else
  echo "ERROR:  3 arguments must provide."
  exit 1
fi