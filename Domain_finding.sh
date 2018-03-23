#!/bin/bash

#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=1:00:00
#SBATCH --mem=64GB

module load R/3.4.2-foss-2016b
module load BEDTools/2.25.0-foss-2015b
module load HMMER/3.1b2-foss-2016b

# $SPECIES = species name
# $GENOME = Genome fasta file
# Example usage: SPECIES=Gallus_gallus GENOME=GRCg6.fasta sbatch ORF_2_finder_phoenix.sh


# Enter species CARP directory
cd ~/fastdir/InsectCarp/${SPECIES}

# Rscript to extract potentially full length repeats from CARP output
Rscript --vanilla ~/fastdir/InsectCarp/scripts/AcTEiveFinder/Domain_finding.R --species ${SPECIES} --genome ${GENOME}

# Get potential repeat sequences from genome
bedtools getfasta -name -fi ~/fastdir/InsectGenomes/${SPECIES}/${GENOME} -bed ${SPECIES}_repeats.bed -fo ${SPECIES}_repeats.fasta

# Search potential repeat sequences for ORFs
usearch -fastx_findorfs ${SPECIES}_repeats.fasta -aaout ${SPECIES}/${SPECIES}_TEs_ORFs.fasta -orfstyle 7 -mincodons 300

# Identify protein motifs in said repeats
hmmscan --domtblout ${SPECIES}_TEs_ORFs.out ~/fastdir/krishna_databases/Pfam-A.hmm ${SPECIES}_TEs_ORFs.fasta > ${SPECIES}_TEs_ORFs.log