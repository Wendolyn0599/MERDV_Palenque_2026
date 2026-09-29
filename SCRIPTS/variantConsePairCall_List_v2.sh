#!/bin/bash
# Before running the script, the reference genome must be indexed for both bowtie and samtools
# and located in the same folder as the reference genome FASTA file. If not, run the following
# commands in the directory containing the reference genome:
# samtools faidx /path/GenomeRef.fasta
# bowtie2-build --threads $HILOS -q /path/GenomeRef.fasta /path/GenomeRef

# Script arguments
## 1) Path to the txt file containing the names of the samples to be analyzed. /path/
## 2) Name of the txt file containing the names of the samples to be analyzed. e.g., FileName.txt
## 3) Path and name of the reference genome with the fasta extension. e.g., Genome.fasta
## 4) Prefix to add to the end of the names
## 5) Path for saving output files
## 6) Path and name of the GFF/GFF3 annotation file
## qsub -R y -l h_rt=23:59:59 -t 1-54:2 /home/adrianace/scripts/variantConsePairCall_List_v2.sh /scratch/adrianace/ensamblesRhabdo/ seq_MRhabdo.txt /scratch/adrianace/BD/merida_virus.fasta MERDVcons /scratch/adrianace/MERDVconsensus/ /scratch/adrianace/BD/MERDV_NC_040599_1.gff3
## FileName.txt: Names of the files to be analyzed.

source $HOME/.bashrc
module load programs/samtools-1.22
module load programs/bamtools-2.5.1
module load programs/bcftools-1.22
module load programs/ivar-1.3.1
module load programs/bowtie2-2.5.0


#$ -N  ConcesusVariant
#$ -j y
#$ -o /scratch/adrianace/MERDVconsensus/outputs/
#$ -pe thread 4
#$ -l h_vmem=16G
set -e

if [ $# -eq 6 ]; then
	echo "--------------CONSENSUS SEQUENCE CONSTRUCTION----------------"
  #Receive arguments
  DIR=$1
  LIST=$2
  REFFILE=$3
  PREF=$4
  OUT=$5
  GFFILE=$6
  mkdir -p "$OUT"
  #Navigate to the directory with the files to analyze
	cd $DIR
  # Read the list of files into an array
  mapfile -t FILES < $LIST
  # Get the files from the array
  FILE1="${FILES[$((SGE_TASK_ID - 1))]}"  # Arrays start at 0, so Read1 is SGE_TASK_ID - 1
  FILE2="${FILES[$SGE_TASK_ID]}"          # Read2 corresponds to SGE_TASK_ID
  echo "Read 1: $FILE1"
	echo "Read 2: $FILE2"
  
  #Get the base name of the file to analyze
	nSe=${FILE1%%.*}     #Keep base name
  # Remove any remaining _R1_ or _R1 that were left from renaming in previous processes
	nSe=${nSe//_R1/}
	echo "Sample name: $nSe" 
  echo ""
	
	#Base name of the reference genome
	BOWTIEREF=${REFFILE%.*}
	
  echo "1. INITIAL READ COUNT" 
  # Detect file extension and count sequences
  if [[ "$FILE1" == *".fastq.gz" ]]; then
    lecturas=$(( $(zcat $FILE1 $FILE2 | wc -l) / 4 )) 
  elif [[ "$FILE1" == *".fastq" ]]; then
    lecturas=$(( $(cat $FILE1 $FILE2 | wc -l) / 4 ))
  else
    echo "The provided file is NOT fastq or fastq.gz"
    exit 1
  fi
  echo "Number of initial individual reads (R1 + R2):"
  echo $lecturas
  
  echo "   "
  echo "2. MAPPING READS TO REFERENCE WITH BOWTIE...." 
  if bowtie2 --very-sensitive-local --soft-clipped-unmapped-tlen  -t --sam-no-qname-trunc -X 5000 -I 400 -x $BOWTIEREF -1 $FILE1 -2 $FILE2 -S "$OUT/${nSe}_${PREF}.sam" --threads $NSLOTS; then   
    echo "...Bowtie2 executed successfully!"
  else 
    echo "Bowtie2... ERROR "
    exit 1 
  fi
      
	echo ""
  echo "3. CONVERT SAM -> BAM "
  if samtools view -bT $REFFILE "$OUT/${nSe}_${PREF}.sam" -o "$OUT/${nSe}_${PREF}.bam" --threads 2; then
    echo "Samtools view...executed successfully!"
  else
		echo "Samtools view... ERROR"
    exit 1
  fi 
    
	echo ""
  echo "4. SORT BAM ..."
  if samtools sort "$OUT/${nSe}_${PREF}.bam" -o "$OUT/${nSe}_${PREF}.srt.bam" --threads 2; then
    echo "Samtools sort...executed successfully!"
  else
		echo "samtools sort... ERROR"
    exit 1
  fi
	rm "$OUT/${nSe}_${PREF}.bam"
    
	if samtools index "$OUT/${nSe}_${PREF}.srt.bam"; then  
    echo "Samtools index... OK "  
  else 
    echo "Samtools index ... ERROR "
    exit 1 
  fi
  
  echo "Number of aligned reads: "
  samtools view -c -F 4 "$OUT/${nSe}_${PREF}.srt.bam"
  
  echo "5. GENERATE VCF FILE..." 
  if bcftools mpileup -A -d 2000000 -Q 10 -L 1500000 -m 1 -P ILLUMINA -a AD,ADF,ADR,DP,SP,INFO/AD,INFO/ADF,INFO/ADR --output-type u -f $REFFILE "$OUT/${nSe}_${PREF}.srt.bam" --threads 1 > "$OUT/${nSe}_${PREF}.txt"; then 
    echo "Bcftools ... OK" 
  else
    echo "Bcftools ... ERROR" 
    exit 1 
  fi
  echo "6. CALL VCF FILE..." 
  cat "$OUT/${nSe}_${PREF}.txt" | bcftools call -A -m -M --threads 1 -Oz > "$OUT/${nSe}_${PREF}.vcf.gz"
	rm "$OUT/${nSe}_${PREF}.txt"
    
  echo "7. CALL VARIANTS AND ANNOTATE WITH GFF..." 
  if samtools mpileup -A -B -d 1000000 -Q 0 -f $REFFILE "$OUT/${nSe}_${PREF}.srt.bam" | ivar variants -p "$OUT/${nSe}_${PREF}.ivar" -q 20 -t 0.1 -g "$GFFILE"  -r "$REFFILE"; then
    echo "iVar variants and annotation with GFF... OK"
  else
    echo "iVar variants and annotation with GFF... ERROR"
    exit 1
  fi
  # Filter only TRUE variants
  awk 'BEGIN{FS=OFS="\t"} NR==1 {print; next} $14=="TRUE" {print}' "$OUT/${nSe}_${PREF}.ivar.tsv" > "$OUT/${nSe}_${PREF}_TRUE.tsv"
  
  echo "7. CALL CONSENSUS"
  if samtools mpileup -aa -A -d 1000000 -Q 0 -f $REFFILE "$OUT/${nSe}_${PREF}.srt.bam" | ivar consensus -q 20 -t 0.5 -m 3 -n N -p "$OUT/${nSe}_${PREF}.ivar"; then
    echo "iVar consensus... OK"
  else
    echo "iVar consensus... ERROR"
    exit 1
  fi
    
	mv "$OUT/${nSe}_${PREF}.ivar.fa" "$OUT/${nSe}_${PREF}.ivar.fasta"
  rm "$OUT/${nSe}_${PREF}.sam" 
	#gzip -q $nSe".srt.bam"
	gzip -q "$OUT/${nSe}_${PREF}.ivar.qual.txt"
  
  #Clean consensus sequences
	perl "/home/adrianace/scripts/cleanConsensus_v4.pl" "$OUT/${nSe}_${PREF}.ivar.fasta" "$OUT/${nSe}_${PREF}.ivar.tsv" "$OUT/${nSe}_${PREF}.ivar2.fasta" $REFFILE
  
  echo "8. VARIANT SUMMARY:"
  TOTAL_VAR=$(tail -n +2 "$OUT/${nSe}_${PREF}.ivar.tsv" | wc -l)
  echo "Total variants found: $TOTAL_VAR"
  FALSE_VAR=$(awk 'BEGIN{FS="\t"} NR>1 && $14=="FALSE" {count++} END{print count+0}' "$OUT/${nSe}_${PREF}.ivar.tsv")
  echo "FALSE variants (noise): $FALSE_VAR"
  TRUE_COUNT=$(tail -n +2 "$OUT/${nSe}_${PREF}_TRUE.tsv" | wc -l)
  echo "TRUE variants: $TRUE_COUNT "
  SUM=$((FALSE_VAR + TRUE_COUNT))
  if [ $TOTAL_VAR -eq $SUM ]; then
    echo "Consensus generation, annotation and TRUE variant filtering completed successfully!"
  else
    echo "ERROR: Inconsistency in filtering"
   fi 
 #Delete temporary files
  rm "$OUT/${nSe}_${PREF}.ivar.fasta"
  rm "$OUT/${nSe}_${PREF}.ivar.tsv"
  
else
	echo "Error. 6 input parameters are required!!"
fi