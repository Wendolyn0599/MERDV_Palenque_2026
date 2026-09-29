#!/bin/bash
# Script to generate a phylogenetic tree using IQ-Tree, adjusting parameters after an initial preliminary tree run (-nt AUTO)
# Saves outputs in the same directory as the aligned sequences being analyzed
# Used after a preliminary analysis, with specific parameters
# Arguments
## 1) Path and filename of the aligned FASTA file containing the sequences to analyze
## 2) Model to use
## Example: ./script.sh /path/to/alignment.fasta GTR+F+I+G4

# Note: The tree was rooted manually using FigTree after IQ-Tree analysis.
# Outgroup sequences were not specified during IQ-Tree run.

#$ -N IQTree 
#$ -j y
#$ -o /scratch/adrianace/MERDVphilogeny/outputs/
#$ -l h_vmem=8G
#$ -pe thread 4

set -e
source $HOME/.bashrc
module load programs/iqtree-1.6.12

if [ $# -eq 2 ]; then
 
    echo "3------------------------------BUILDING PHYLOGENETIC TREE WITH IQ-TREE---------------------------------"
    echo "Alignment : $1"
    echo "Model     : $2"

  if iqtree -s "${1}" -m "${2}" -bb 1000  -bnni  -nt $NSLOTS; then
    echo "IQ-Tree...OK"
  else
    echo "IQ-Tree...ERROR"
    exit 1
   fi
else
  echo "ERROR: 2 arguments must be providen"
  exit 1
fi