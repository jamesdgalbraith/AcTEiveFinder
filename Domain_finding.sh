#!/bin/bash

#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=1:00:00
#SBATCH --mem=64GB

module load R/3.4.2-foss-2016b
module load BEDTools/2.25.0-foss-2015b
module load HMMER/3.1b2-foss-2016b
module load SAMtools

echo "Beginning Domain finding on " ${SPECIES}

# $SPECIES = species name
# $GENOME = Genome fasta file
# Example usage: SPECIES=Gallus_gallus GENOME=GRCg6.fasta sbatch Domain_finding.sh

# Check for necesary files
if [ ! -s ~/fastdir/Genomes/${SPECIES}/${GENOME} ]; then
    echo ${SPECIES} " genome not found"; exit
fi
if [ ! -s ~/fastdir/Carp/${SPECIES}/${SPECIES}_Denovo_TE_Library.fasta ]; then
    echo "TE library not present or empty "${SPECIES}; exit
fi
if [ ! -s ~/fastdir/Carp/${SPECIES}/${SPECIES}.igor.gff ]; then
    echo "No repeats identified in "${SPECIES}; exit
fi

# Enter species CARP directory
cd ~/fastdir/Carp/${SPECIES}

# Rscript to extract potentially full length repeats from CARP output
Rscript --vanilla ~/fastdir/Carp/scripts/AcTEiveFinder/Domain_finding.R --species ${SPECIES} --genome ${GENOME}

# Check for bed length
if [ ! -s ${SPECIES}_repeats.bed ]; then
    echo "No repeats identified in "${SPECIES}; exit
fi

# Get potential repeat sequences from genome
bedtools getfasta -fi ~/fastdir/Genomes/${SPECIES}/${GENOME} -bed ${SPECIES}_repeats.bed -fo ${SPECIES}_repeats.fasta

# Check for repeat fasta
if [ ! -s ${SPECIES}_repeats.fasta ]; then
    echo "bedtools did not make fasta "${SPECIES}; exit
fi

# Search potential repeat sequences for ORFs
usearch -fastx_findorfs ${SPECIES}_repeats.fasta -aaout ${SPECIES}_TEs_ORFs.fasta -orfstyle 7 -mincodons 300

# Check for ORF fasta
if [ ! -s ${SPECIES}_TEs_ORFs.fasta ]; then
    echo "No repeat ORFs identified in "${SPECIES}; exit
fi

# Identify protein motifs in said repeats
hmmscan --domtblout ${SPECIES}_TEs_ORFs.out ~/fastdir/krishna_databases/Pfam-A.hmm ${SPECIES}_TEs_ORFs.fasta > ${SPECIES}_TEs_ORFs.log

echo "Finished Domain finding on " ${SPECIES}
date