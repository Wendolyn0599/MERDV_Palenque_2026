#!/bin/bash

## Script for alignment and de novo assembly
## Arguments:
## 1) Path to the file containing the list of filenames to analyze.
##    Example: /scratch/adrianace/ResultadosBowtie/
## 2) Name of the file containing the list of library files to analyze.
##    Example: seq_noMmosco_noMPEZ_noMRibo.txt
## 3) Minimum contig length (e.g., 500)
## 4) Name and path of the reference sequence database (must have .fasta extension).
## 5) Result prefix to place before _M and noM. Example: Rhabdo
## 6) Path for saving output files
## 7) Prefix for read IDs (e.g., NS500) # For use with the Perl script
## example of use:
## qsub -R y -l h_rt=23:59:59 -t 185-186:2 MapBBEnsamSpadeDupMagBwt_list.sh scratch/adrianace/ResultadosBowtie/ seq_limpias.txt 500 /scratch/adrianace/BD/rhabdoviridae.fasta Rhabdo /scratch/adrianace/ensamblesRhabdo NS500


#$ -N  MapBBEnsamSpadeDupMagBwt_list
#$ -j y
#$ -o /scratch/adrianace/ensamblesRhabdo/outputs
#$ -pe thread 4
#$ -l h_vmem=32G

source ~/.bashrc
module load programs/bowtie2-2.5.0
module load compilers/perl-5.28
module load programs/bbmap-38.84
module load programs/spades-3.15.4
module load programs/cdhit-4.8.1

if [ $# -eq 7 ]; then
	cd $1
	# Read the list of files into an array
  mapfile -t FILES < $2
  
  # Get  Read1 and Read2 from the array
  FILE1="${FILES[$((SGE_TASK_ID - 1))]}"  # Arrays are zero-indexed, so Red1 is SGE_TASK_ID - 1
  FILE2="${FILES[$SGE_TASK_ID]}"          # Read2 corresponds to SGE_TASK_ID 
  echo "Read 1: $FILE1"
	echo "Read 2: $FILE2"
  
  # Clean filename and get extension
	sample="${FILE1%%.*}"     #Basename
  tipo="${FILE1##*.}"      #Extension
  
  # Remove the _R1_s left over from renaming in previous processes (they remained in the middle of the name)
  sample=$(echo "$sample" | sed 's/_R1_/_/g') 
	# Remove the _R1 from the end
  TEMP=''
	R1='_R1'
	FLT="${sample/$R1/$TEMP}"
	sample=$FLT
	echo "Sample name:  $sample "
  echo ""
	
	echo "1. INITIAL READING COUNT (R1 + R2): " 
	#Check length and count reads
	if [ $tipo == 'fastq' ]; then
		echo "Numero de lecturas iniciales en $FILE1 y $FILE2 :"
		numLin=$(cat $FILE1 $FILE2 | wc -l)
		echo $((numLin / 4))
	elif [ $tipo == 'gz' ]; then
		echo "Initial reads in $FILE1 and $FILE2 :"
		numLin=$(zcat $FILE1 $FILE2 | wc -l)
		echo $((numLin / 4))
	else
		echo "Given file is NOT fastq or fastq.gz"; exit 1
	fi
  echo " "
  
  OUTDIR=$6
  mkdir -p "$OUTDIR"
  
	echo "2. Processing alignment with BBMap...-----------------------------------------------------------------------------"
	if bbmap.sh ref=$4 nodisk=t pairedonly=f keepnames=t minid=0.75 in1=$FILE1 in2=$FILE2 local=t outm="$OUTDIR/${sample}_${5}_R1.fastq.gz" outm2="$OUTDIR/${sample}_${5}_R2.fastq.gz" out=/dev/null threads=$NSLOTS; then 
    echo "BBMap...OK!"
  else 
    echo "BBMap... ERROR"
    exit 1
	fi
  echo ""

	echo "3. Read assembly with SPAdes...-------------------------------------------------------------------------------" 
	if spades.py --rnaviral -1 "$OUTDIR/${sample}_${5}_R1.fastq.gz" -2 "$OUTDIR/${sample}_${5}_R2.fastq.gz" -o "$OUTDIR/${sample}_${5}_Spades" --disable-gzip-output -t $NSLOTS -m 32 -k 21,31,41,51; then 
    echo "SPAdes...OK!"
  else 
    echo "SPAdes...ERROR "
    exit 1
  fi
  
  #SPAdes generates "scaffolds.fasta" and "contigs.fasta" files; they contain essentially the same information,
  #with the difference that scaffolds.fasta has fewer but longer sequences, as it joins fragments from the same organism using paired-end information.
  #In a targeted assembly like this one, the contigs already cover a large part of the genome, so scaffolds offer no added value and may
  #even be empty. It is best to use contigs.fasta.
	cp "$OUTDIR/${sample}_${5}_Spades/contigs.fasta" "$OUTDIR/${sample}_${5}_Spades.fasta"
  ls -lh "$OUTDIR/${sample}_${5}_Spades.fasta"
  echo "Number of contigs:"
	grep -c ^'>' "$OUTDIR/${sample}_${5}_Spades.fasta"
	rm "$OUTDIR/${sample}_${5}_R1.fastq.gz" "$OUTDIR/${sample}_${5}_R2.fastq.gz"
	echo " "
	
	echo "4. Remove redundancy with CD-HIT.-----------------------------------------------------------------------------------------" 
  cd "$OUTDIR"
	if cd-hit-est -i "${sample}_${5}_Spades.fasta" -o "${sample}_${5}_Spades_CD_T.fasta" -c 0.90 -M 30000 -aS 0.9 -G 0 -l $3 -g 1 -T $NSLOTS; then
	  echo "CD-HIT...OK!"
  else 
    echo "CD-hit... ERROR "
    exit 1
  fi
  
	rm "${sample}_${5}_Spades.fasta"
	rm *.clstr
  
  # putContigMag.pl expects a 2-line-per-sequence FASTA format (header + sequence).
  # This Perl script converts multi-line FASTA to that format.
  # The `sed` command normalizes spaces to '_' so that the names match those in the Bowtie2 SAM output.
	perl -pe '/^>/ ? print "\n" : chomp' "${sample}_${5}_Spades_CD_T.fasta"| tail -n +1 > "${sample}_${5}_Spades_CD.fasta"
	sed -i 's/ /_/g' "${sample}_${5}_Spades_CD.fasta"
	rm "${sample}_${5}_Spades_CD_T.fasta"

	echo " "
	echo "5. ALIGNMENT OF CONTIGS WITH BOWTIE2...---------------------------------------------------------------------------"
	echo "...Indexando contings..."
	if bowtie2-build --threads $NSLOTS -q "${sample}_${5}_Spades_CD.fasta" "${sample}_${5}_Spades_CD"; then
    echo "Bowtie2...Contings indexados correctamente"
  else  
    echo "Bowtie2-build...ERROR "
    exit 1
  fi
  
  cd $1
	if bowtie2 -L 21 --very-sensitive-local -t --sam-no-qname-trunc -x "$OUTDIR/${sample}_${5}_Spades_CD" -1 $FILE1 -2 $FILE2 --al-conc-gz  "$OUTDIR/${sample}_M${5}.fastq.gz" --un-conc-gz  "$OUTDIR/${sample}_noM${5}.fastq.gz" -S "$OUTDIR/${sample}_${5}_Spa.sam" --threads $NSLOTS ; then
    echo "Bowtie2...OK!"
  else
    echo "Bowtie2... ERROR "
    exit 1
  fi
  
	rm "$OUTDIR/${sample}_${5}_Spades_CD".*.bt2
  
  #Move to OUTDIR to rename files
  cd "$OUTDIR" 
	rename ".fastq.1.gz" "_R1.fastq.gz" *
	rename ".fastq.2.gz" "_R2.fastq.gz" *
  echo "Individual reads (R1 +R2)"
  echo "Unmapped:" 
	numLin=$(zcat "${sample}_noM${5}_R1.fastq.gz" "${sample}_noM${5}_R2.fastq.gz" | wc -l)
	echo $((numLin / 4)) 
	
	echo "Mapped:"
	numLin=$(zcat "${sample}_M${5}_R1.fastq.gz" "${sample}_M${5}_R2.fastq.gz" | wc -l)
	echo $((numLin / 4))

	echo "6. Assigning Magnitudes to Contigs...-----------------------------------------------------------------------------"
  # Magnitude = number of reads that mapped to a contig/scaffold.
  # Perl arguments:
  ## $ARGV[0]; SAM file (Bowtie2 mapping)
  ## $ARGV[1]; Contigs/scaffolds FASTA
  ## $ARGV[2]; Output FASTA
  ## $ARGV[3]; Read name prefix (e.g., NS500)
	if perl /home/adrianace/scripts/putContigMag.pl "${sample}_${5}_Spa.sam" "${sample}_${5}_Spades_CD.fasta" "${sample}_${5}_Spades_CD_Mag.fasta" $7; then 
    echo "perl...OK"
  else 
    echo "perl...ERROR"
    exit 1
  fi
	 
  echo "Total reads assigned to contigs (sum of magnitudes):"
  grep "^>" "${sample}_${5}_Spades_CD_Mag.fasta" | cut -f2 -d'=' | awk '{ sum+=$1 } END {print sum}'

  #rm "${sample}_${5}_Spa.sam"
  #rm "${OUTDIR}/${sample}_${5}_Spades_CD.fasta"

else
  echo "Error: 7 arguments are needed!!"
  exit 1
fi
