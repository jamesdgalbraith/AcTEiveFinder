#!/bin/bash

# Use to find repeats in with CARP (second steo step)
# SPECIES = name of species
# GENOME = file name of genome (needed for following step)

#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time=1-00:00
#SBATCH --mem=120GB

# Notification configuration
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=james.galbraith@student.adelaide.edu.au

module load censor
module load wu-blast
module load bioperl
module load MUSCLE/3.8.31

cd ~/fastdir/Carp/${SPECIES}/

echo "Started igor on " ${SPECIES}
date

# Check matrix output exists
if [ ! -s ${SPECIES}.gff ]; then
    echo "GFF not found"; exit
fi

# Create repeat families from krishna-matrix output using igor
igor -in ${SPECIES}.gff -out ${SPECIES}.json

echo "Ended igor on " ${SPECIES}
date

cd ~/fastdir/Carp/slurms/seqer
SPECIES=${SPECIES} GENOME=${GENOME} sbatch ~/fastdir/Carp/scripts/seqer.sh
echo "Queued seqer"