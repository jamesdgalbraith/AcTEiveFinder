#!/bin/bash

# Use to find repeats with CARP (third step)

#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time12:00:00
#SBATCH --mem=64GB

# Notification configuration
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=james.galbraith@student.adelaide.edu.au

module load censor
module load wu-blast
module load bioperl
module load MUSCLE/3.8.31

cd ~/fastdir/Carp/${SPECIES}/

echo "Started seqer on " ${SPECIES}
date
echo ""

# Create fastq consensus of repeat families using seqer
gffer < ${SPECIES}.json > ${SPECIES}.igor.gff
seqer -aligner=muscle -dir=consensus -fasta=true -maxFam=100 -subsample=true -minLen=0.80 -threads=16 -ref=/home/a1194388/fastdir/Genomes/${SPECIES}/${GENOME} ${SPECIES}.igor.gff

find ./consensus -maxdepth 1 -name '[!.]*.fq' -print0 | xargs -r0 cat > consensus.fasta

echo "Ended seqer on " ${SPECIES}
date

cd ~/fastdir/Carp/slurms/censor
SPECIES=${SPECIES} sbatch ~/fastdir/Carp/scripts/censor_report.sh
echo "Queued censor"
