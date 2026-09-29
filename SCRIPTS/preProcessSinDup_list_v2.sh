#!/bin/bash

#Arguments:
## 1) Path of the file (txt) with the names of the libraries to analyze.
##    Must be in the same folder as the fastq files. Example: /scratch/btaboada/Data/Wendo
## 2) Name of the file with the library names to analyze. Example: Filenames.txt
## 3) Path of the output results. Example: /scratch/btaboada/Results/Wendo
## Example of submission: 
## qsub -R y -l h_rt=23:59:59 -t 1-80:2 preProcessSinDup_list_v2.sh /scratch/btaboada/Data/Wendo Filenames.txt /scratch/btaboada/Data/Wendo/resultados_fastp

# Job name
#$ -N preProcess_list

# Merge stdout and stderr
#$ -j y
#$ -o /home/adrianace/scripts/faltantes_outs
#$ -pe thread 4 
#$ -l h_vmem=7G

# Load fastp module
source $HOME/.bashrc
module load programs/fastp-1.0.1

if [ $# -ge 3 ]; then
    # Move to folder with fastq files
    cd $1

    # Get filenames for this task
    FILE1=$(cat $2 | head -n $SGE_TASK_ID | tail -1)
    FILE2=$(cat $2 | head -n $(($SGE_TASK_ID + 1)) | tail -1)
    echo "Files to process: $FILE1 $FILE2"

    # Extract filename parts
    nSe1=${FILE1%%.*}
    nSe2=${FILE2%%.*}
    echo $nSe1 $nSe2

    # Optional modification of names (placeholder)
    TEMP=''
    R1='R1'
    FLT=${nSe1/$R1/$TEMP}

    echo " "
    echo "1. INITIAL READS:"
    # Count initial reads
    tipo=${FILE1#$nSe1.}
	echo "Tipo" $tipo
    if [ $tipo == 'fastq' ]; then
        echo "Initial reads in $FILE1 $FILE2:"
        numLin=$(cat $FILE1 $FILE2 | wc -l)
        expr $numLin / 4
    elif [ $tipo == 'fastq.gz' ]; then
        echo "Initial reads in $FILE1 $FILE2:"
        numLin=$(zcat $FILE1 $FILE2 | wc -l)
        expr $numLin / 4
    else
        echo "Error: file is not fastq or fastq.gz"
        exit 1
    fi

    echo " "
    echo "2. STARTING FASTP PROCESS"

    # Run fastp with HTML and JSON reports
    fastp -i $FILE1 -I $FILE2 \
          -o $3/${nSe1}_Q_D.fastq.gz \
          -O $3/${nSe2}_Q_D.fastq.gz \
          -h $3/${nSe1}_fastp.html \
          -j $3/${nSe1}_fastp.json \
          -V -w $NSLOTS -q 15 --cut_tail --cut_window_size 3 \
          --cut_mean_quality 15 --cut_mean_quality 10 -x \
          --poly_x_min_len=30 -g --poly_g_min_len=20 -l 40 -n 15 -y -Y 30 -D
    # Verify fastp success
		echo "fastp ... O.K" || (echo "fastp ... ERROR" ; exit 1)

    echo "END"
    echo " "

    echo "3. FINAL CLEAN READS:"
    # Count final reads
    if [ $tipo == 'fastq' ]; then
        echo "Final reads in $nSe1"_Q_D.fastq.gz" $nSe2"_Q_D.fastq.gz":"
        numLin=$(cat $3/${nSe1}_Q_D.fastq.gz $3/${nSe2}_Q_D.fastq.gz | wc -l)
        expr $numLin / 4
    elif [ $tipo == 'fastq.gz' ]; then
        echo "Final reads in $nSe1"_Q_D.fastq.gz" $nSe2"_Q_D.fastq.gz":"
        numLin=$(zcat $3/${nSe1}_Q_D.fastq.gz $3/${nSe2}_Q_D.fastq.gz | wc -l)
        expr $numLin / 4
    else
        echo "Error: file to save not fastq or fastq.gz"
    fi

else
    echo "Error: incorrect number of parameters"
fi

