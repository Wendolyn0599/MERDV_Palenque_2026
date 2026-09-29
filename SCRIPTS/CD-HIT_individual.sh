#!/bin/bash

#SCRIPT TO REMOVE SEQUENCE REDUNDANCY USING CD-HIT
#Arguments
# 1) Path to the file containing the sequences to be analyzed; output files are saved in this same path
# 2) Name of the file containing the sequences to be analyzed
# 3) Identity percentage for CD-HIT, e.g., 0.99 (for 99% identity)

## usage example:
## qsub -R y -l h_rt=23:59:59 /home/adrianace/scripts/CD-HIT_individual.sh /scratch/adrianace/filoRhabdo/seqFiloMerida/ MV_Brazil_2021.fasta 0.99 

set -e
#$ -N CD-HITindividual
#$ -j y
#$ -o /scratch/adrianace/MERDVphilogeny/outputs
#$ -pe thread 4
#$ -l h_vmem=8G

source ~/.bashrc
module load programs/cdhit-4.8.1

if [ $# -eq 3 ]; then

  echo "-------------------------------------------RUNNING CD-HIT -------------------------------------------"
  cd "$1"
  echo "Processing file: ${2}"
  echo "Initial sequences:"
  grep "^>" "$2" | wc -l
  # Remove file extension
  name="${2%.*}"
  # Convert decimal to integer to add % to the filename
  perc=$(awk -v val="$3" 'BEGIN { result = val * 100; printf "%d", result }')
	if cd-hit-est -i "$2" -o "${name}_CD_${perc}.fasta" -c "$3" -M 30000 -aS 0.98 -G 0 -g 1 -T $NSLOTS; then
	  echo "Secuencias finales:"
    grep "^>" "${name}_CD_${perc}.fasta" | wc -l
    echo "CD-HIT...OK!"
  else 
    echo "CD-hit... ERROR "
    exit 1
  fi

else
  echo "Error: 3 arguments must be provided!"
fi


