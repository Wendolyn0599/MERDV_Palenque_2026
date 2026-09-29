#!/bin/bash
# Script to run BLAST on contigs generated from an assembly
# Arguments:
# 1) Path to the files to be analyzed
# 2) Pre-generated list of file names (located in the same directory as the files)
# 3) Path to the database
# 4) Database name without extension (e.g., core_nt)
# 5) Path for output files
# 6) Prefix to add to output files

#$ -N BlastJA
#$ -t 1-18
#$ -j y
#$ -l h_vmem=32G
#$ -o /scratch/adrianace/blastnRhabdo/MERDVconsBlast

set -e

source $HOME/.bashrc
module load programs/blast-2.12.0

# Verify that the arguments were provided
if [ $# -eq 6 ]; then
  DIR="$1"
  LISTA="$2"
  export BLASTDB=$3 # Location of the formatted databases on the cluster.
  DB="$4"
  OUT="$5"
  pref="$6"
  
  cd $DIR 
  # Read the list of files into an array and get the sample name
  mapfile -t FILES < $LISTA
  file="${FILES[$((SGE_TASK_ID - 1))]}"  # Arrays are zero-indexed, so the file to analyze is SGE_TASK_ID - 1
  echo "--------------EJECUTANDO BLASTN ------------ "
  echo " "
  mkdir -p $OUT
  echo "Processig file: '$file'"
  sample=$(basename "$file")
  sample=${sample%.*}
  echo "Basename : ${sample}"
  echo "DATA BASEs: ${DB} "
  
  
  if blastn -query $file -db $DB -outfmt "6 qseqid sseqid pident length qcovs evalue bitscore" -max_target_seqs 1 -out "${OUT}/${sample}_${pref}.tsv"; then
      echo "Blastn...OK!"
  else
      echo "ERROR AL EJECUTAR BLASTN!"
      exit 1 
  fi
    
else
  echo "Error: 6 arguments must provided ."

fi