# Pipeline for MERDV assembly and analysis

This directory contains all the scripts required for raw read processing, 
genome assembly, consensus generation,  phylogenetic and statistics analysis of 
MERDV in mosquitoes from Palenque,Chiapas, Mexico.

**Authors:** Wendolyn Estrada, Guillermo Orta-Pineda, Marco A. Espinoza, Rafael Ojeda-Flores, Blanca Taboada
**Correspondence: ** BT, blanca.taboada@ibt.unam.mx

---

This work was performed using the HPC infrastructure at the **Unidad Universitaria de Secuenciación Masiva y Bioinformática, Instituto de Biotecnología (UNAM)**, which is part of the Laboratorio Nacional de Apoyo Tecnológico a las Ciencias Genómicas (CONACyT).

We thank the Instituto de Biotecnología-UNAM for giving us access to its computer cluster.

Documentation on cluster usage and available software:  
http://teopanzolco.ibt.unam.mx/index.php/documentacion/

---

## Contents:


1. [Pipeline structure and workflow](#1-pipeline-structure-and-workflow)
2. [Software requirements](#2-software-requirements)
3. [Bash scripts](#3-Bash-scripts)
4. [Auxiliary scripts](#4-auxiliary-scripts)
5. [R scripts](#5-R-scripts)
6. [External analyses](#6-external-analyses)

---


## 1. Pipeline structure and workflow

The first part of this  pipeline consists of six sequential steps that were 
executed on a Linux SGE (Sun Grid Engine) cluster. Each script takes the output 
of the previous one as input.
The second part of the pipeline (steps 7 and 8) consisted of a series of 
analyses using tools external to the cluster


### Workflow diagram

```
Raw reads (FASTQ)
    │
    ▼
[1] preProcessSinDup_list_v2.sh
    │   fastp (quality control)
    ▼
[2] bowtie2_Pair_list_v2.sh
    │   Bowtie2 (remove host DNA)
    ▼
[3] MapBBEnsamSpadeDupMagBwt_list.sh
    │   BBMap + SPAdes + CD-HIT + Bowtie2 + Perl
    │   (Rhabdoviridae filtering, assembly, redundancy, abundance)
    ▼
    BLASTN_contings.sh
    │   BLASTn (taxonomic confirmation)
    ▼
[4] variantConsePairCall_List_v2.sh
    │   iVar + cleanConsensus_v4.pl
    │   (consensus and variants)
    ▼
[5] MultAlin_MafftEnboss.sh + CD-HIT_individual.sh
    │   MAFFT + CD-HIT
    │   (alignment and redundancy filtering)
    ▼
[6] IQTree.sh
        IQ-TREE
        (phylogenetic tree)
[7] External analyses 
    │   -Visualization (FigTree)
    │   -Haplotypes (PopART)
    │   - Selection analysis (HyPhy + Datamonkey)
    ▼
[8] Statistics, plots, and tables (R Studio)

```
## 2. Software requirements

### For processing on the cluster

| Software | Version | Purpose |
|----------|---------|---------|
| Bowtie2 | v2.5.0 | Read alignment and host removal |
| fastp | v1.0.1 | Quality control and trimming |
| BBMap | v38.84 | Filtering against Rhabdoviridae database |
| SPAdes | v3.15.4 | *De novo* assembly |
| CD-HIT | v4.8.1 | Redundancy removal |
| iVar | v1.3.1 | Consensus and variant calling |
| MAFFT | v7.407 | Multiple sequence alignment |
| IQ-TREE | v1.6.12 | Phylogenetic inference |
| BLAST+ | v2.15.0 | Taxonomic identification of contigs |
| Perl | v5.28 | Auxiliary scripts |

### For statistical analysis in R

| Package | Version | Purpose |
|---------|---------|---------|
| dplyr | ≥1.0.0 | Data manipulation and transformation |
| purrr | ≥0.3.0 | Functional programming and iteration |
| stringr | ≥1.4.0 | String manipulation and pattern extraction |
| tidyr | ≥1.1.0 | Data reshaping (pivot, separate) |
| readr | ≥2.0.0 | Reading TSV and delimited files |
| writexl | ≥1.4.0 | Export tables to Excel format |
| Biostrings | ≥2.60.0 | DNA sequence manipulation and analysis |
| ape | ≥5.6.0 | Phylogenetic analysis and distance matrices |
| rentrez | ≥1.2.0 | NCBI Entrez API for taxonomic queries |
| ggplot2 | ≥3.3.0 | Data visualization and plotting |
| patchwork | ≥1.1.0 | Combining multiple ggplot2 plots |
| scales | ≥1.1.0 | Plot formatting and axis scaling |
---

## 3. Bash scripts

| Step | Script | Function | Key parameters |
|:----:|--------|----------|----------------|
| 1 | `preProcessSinDup_list_v2.sh` | Quality control with fastp | `-q 15`; `--cut_tail`; `-x --poly_x_min_len=30`; `-l 40`; `-n 15`; `-y -Y 30` |
| 2 | `bowtie2_Pair_list_v2.sh` | Remove mosquito host reads | `--very-sensitive`; `-X 2000` |
| 3a | `MapBBEnsamSpadeDupMagBwt_list.sh` | Rhabdoviridae filtering, assembly, redundancy, abundance | **BBMap:** `minid=0.75`; **SPAdes:** `--rnaviral -k 21,31,41,51`; **CD-HIT:** `-c 0.90 -aS 0.9 -l 500` |
| 3b | `BLASTN_contings.sh` | Taxonomic confirmation of contigs | `-evalue 1e-5`; `-max_target_seqs 1` |
| 4 | `variantConsePairCall_List_v2.sh` | Map to MERDV, consensus, variants | **Bowtie2:** `--very-sensitive-local -X 5000`; **iVar variants:** `-q 20 -t 0.1`; **iVar consensus:** `-q 20 -t 0.5 -m 3` |
| 5a | `MultAlin_MafftEnboss.sh` | Multiple alignment (MAFFT) | `--genafpair`; `--maxiterate 1000` |
| 5b | `CD-HIT_individual.sh` | Redundancy filtering by site/year | `-c 0.98 -aS 0.98` |
| 6 | `IQTree.sh` | Phylogenetic inference | `-bb 1000`; `-bnni` |

### Typical job submission

All scripts were submitted to the SGE cluster using commands of the form:

```bash

qsub -R y -l h_rt=23:59:59 script_name.sh [arguments] 

qsub -R y -l h_rt=23:59:59 -t 1-N:2 script_name.sh [arguments]  # If paired reads
```

## 4. Auxiliary scripts

| Script | Function | Used by |
|--------|----------|---------|
| `cleanConsensus_v4.pl` | Cleans consensus by removing low-quality positions and ambiguous variants | `variantConsePairCall_List_v2.sh` |
| `putContigMag.pl` | Assigns abundance (read count) to each contig | `MapBBEnsamSpadeDupMagBwt_list.sh` |

## 5. R scripts 

All R scripts should be executed in the following order, as each script depends on outputs generated by the previous ones.

| Script | Input | Function | Output |
|--------|-------|----------|--------|
| `MERDV_detection.R` | Bowtie2 results, BLASTn results, contig data, metadata | Processing of detection data; generation of positivity tables; statistical tests; Figures 1 and 2 | `supplementary_table_S1.xlsx`; `FIG1.png`; `FIG2.png` |
| `MERDVgenomes_cov_Ncontent_pairwiseDeletion.R` | Consensus genomes (FASTA); multiple alignment (FASTA) | N content and coverage calculation; pairwise deletion identity matrix | `MERDV_consensus_cov.xlsx`; `supplementary_tableS3.csv` |
| `MERDV_variants_classify.R` | iVar variant reports; GFF annotation; metadata | Variant classification and annotation; summary tables by genome and gene | `supplementary_table_S4.xlsx`; `supplementary_table_S5.xlsx`; `MERDV_summaryByGene_normalized.xlsx` |
| `MERDV_MultLogReg_VarPart.R` | Positivity table from `MERDV_detection.R` | Multivariable logistic regression; variance partitioning | `MERDV_variancePartition.xlsx` |
| `MERDV_MLE_rates.R` | Positivity table from `MERDV_detection.R` | Maximum Likelihood Estimation of infection rates by site and genus group | `MERDV_MLE_infection_rates.xlsx` |

## R script execution order and dependencies

```

### 1. Detection and positivity analysis
source("MERDV_detection.R")  # Generates set1_positiveRhab and plots

### 2. Genome quality analysis
source("MERDVgenomes_cov_Ncontent_pairwiseDeletion.R")  # Uses consensus FASTA

### 3. Variant classification
source("MERDV_variants_classify.R")  # Uses iVar outputs and GFF

### 4. Regression and variance partitioning
source("MERDV_MultLogReg_VarPart.R")  # Requires set1_positiveRhab

### 5. Infection rate estimation
source("MERDV_MLE_rates.R")  # Requires set1_positiveRhab

## 6. External analyses 

The following analyses were performed manually using external tools and are not part of the automated pipeline:

| Step | Tool | Purpose | Input | Parameters |
|------|------|---------|-------|------------|
| 1 | Jalview v2.11.5.2 | Alignment inspection and degapping | Aligned FASTA | Manual editing |
| 2 | FigTree v1.4.4 | Tree visualization and editing | `*.treefile` | Manual editing |
| 3 | PopART v1.7 | Haplotype network (median-joining) | Consensus sequences (<5% N) | ε = 0 |
| 4 | Datamonkey/HyPhy | Selection analysis (FEL, BUSTED) | Coding alignments (in-frame) | FEL: p ≤ 0.1; BUSTED: p < 0.05 |

**Note:** For selection analysis, genes N, P, M, G, and L were analyzed independently. Gene-wide dN/dS (ω) values were obtained from the MG94xREV model in FEL.


---




