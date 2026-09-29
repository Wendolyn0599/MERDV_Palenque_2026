#!/bin/bash
# Script to remove mosquito DNA and retain unmapped reads for viromics
# Before running the script, the reference genome (multiple concatenated species) and its Bowtie index must be in the same folder.

# Arguments:
## 1) Path to the .txt file containing the names of the samples to be analyzed (the list and the samples themselves are in the same directory)
## 2) Name of the .txt file containing the sample names (e.g., FileName.txt)
## 3) Path and name of the reference genome indexed for Bowtie (e.g., Genome)
## 4) Prefix to append to the filenames
## 5) Path for saving output files

# Example of script usage:
##   qsub -R y -l h_rt=23:59:59 -t 1-2:2 bowtie2_Pair_list_v2.sh /scratch/btaboada/Wendo/Results_fastp seq_limpias.txt /scratch/btaboada/Wendo/genomes/genomas_mosquito.fasta Mosco /scratch/btaboada/Wendo/Results_bowtie
#seq_limpias.txt: list of filenames to be analyzed
#Note: run a test with the first sample (-t 1-2:2), and if successful, execute the rest (-t 3-160:2)

source $HOME/.bashrc
module load programs/bowtie2-2.5.0
module load programs/samtools-1.10

#$ -N BowtiePairBamList
#$ -j y
#$ -o /home/adrianace/scripts/faltantes_outs
#$ -pe thread 8
#$ -l h_vmem=8G

if [ $# -ge 5 ]; then
	echo "----------Alignment and host removal with Bowtie2------"
	cd $1
	FILE1=$(cat $2 | head -n $SGE_TASK_ID| tail -1)
	echo "Read 1: $FILE1"
	num=$(($SGE_TASK_ID + 1))
	FILE2=$(cat $2 | head -n $num| tail -1)
	echo "Read2: $FILE2"
  REFFILE=$3
  echo "Aligning genome: $3"
  PREF=$4
  OUT=$5
  mkdir -p "$OUT"
  echo "Output file directory: "
  echo "Prefix: $PREF"
 
	#Get file extension
	nSe=${FILE1%.*}
	tipo=${FILE1#$nSe.}
	nSe=${FILE1%%.*}
	#Get file name
	TEMP=''
	R1='_R1'
	FLT=${nSe/$R1/$TEMP}
	nSe=$FLT
	echo "Output: " $nSe
	
	echo " "
	echo "1. INITIAL READS: " 
	#Check length and count reads
	if [ $tipo == 'fastq' ]; then
		echo "Initial reads in $FILE1 $FILE2 :"
		numLin=$(cat $FILE1 $FILE2 | wc -l)
		expr $numLin / 4
	elif [ $tipo == 'gz' ]; then
		echo "Initial reads in $FILE1 $FILE2 :"
		numLin=$(zcat $FILE1 $FILE2 | wc -l)
		expr $numLin / 4
	else
		echo "The provided file is NOT fastq or fastq.gz"
	fi

	echo " "
	echo "2. MAPPING AND REMOVING HOST WITH BOWTIE..." 
	(bowtie2 --very-sensitive -t --sam-no-qname-trunc -X 2000 -x $REFFILE -1 $FILE1 -2 $FILE2 --un-conc-gz $OUT/$nSe"_noM"$PREF".fastq.gz" --al-conc-gz $OUT/$nSe"_M"$PREF".fastq.gz" -S $OUT/$nSe".sam" --threads $NSLOTS &&
		echo " bowtie2 --very-sensitive -t --sam-no-qname-trunc -X 2000 -x $REFFILE -1 $FILE1 -2 $FILE2 --un-conc-gz $OUT/$nSe"_noM"$PREF".fastq.gz" --al-conc-gz $OUT/$nSe"_M"$PREF".fastq.gz" -S $OUT/$nSe".sam" --threads $NSLOTS... OK " ) || ( echo "bowtie2 --very-sensitive -t --sam-no-qname-trunc -X 2000 -x $REFFILE -1 $FILE1 -2 $FILE2 --un-conc-gz $OUT/$nSe"_noM"$PREF".fastq.gz" --al-conc-gz $OUT/$nSe"_M"$PREF".fastq.gz" -S $OUT/$nSe".sam" --threads $NSLOTS... ERROR "; exit 1; )
	echo ""
	echo "Bowtie2....OK"
	echo ""
 
  #Go to output path and rename
  cd "$OUT"

	echo "Unmapped:" 
	rename ".fastq.1.gz" "_R1.fastq.gz" *
	rename ".fastq.2.gz" "_R2.fastq.gz" *
	numLin=$(zcat $nSe"_noM"$PREF"_R1.fastq.gz" $nSe"_noM"$PREF"_R2.fastq.gz"| wc -l)
	expr $numLin / 4
	
	echo "Mapped:"
	numLin=$(zcat $nSe"_M"$PREF"_R1.fastq.gz" $nSe"_M"$PREF"_R2.fastq.gz"| wc -l)
	expr $numLin / 4
	echo ""
	
else
	echo "Error: Five arguments are required!!"
fi

