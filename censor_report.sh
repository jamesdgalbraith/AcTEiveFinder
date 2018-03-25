#!/bin/bash

# Script to run censor on first stage output

#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time=3-00:00
#SBATCH --mem=64GB

# Notification configuration
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=james.galbraith@student.adelaide.edu.au

module load bioperl/1.6.924
module load censor/4.2.29
module load wu-blast/2.0
module load Java/9.0.4

echo "Started censor " ${SPECIES}
date

# Check necessary fasta exists
cd ~/fastdir/Carp/${SPECIES}/

if [ ! -s consensus.fasta ]; then
    echo "Genome not found"; exit
fi

# Run censor
censor -bprm cpus=16 -lib ~/fastdir/krishna_databases/repbase_eukaryote.fa ~/fastdir/Carp/${SPECIES}/consensus.fasta

echo "Ended censor " ${SPECIES}
echo "Started classification " ${SPECIES}
date

# Classify sequences as repeats or not epeats
cd ~/fastdir/Carp/scripts 
java ClassifyConsensusSequences ${SPECIES}

echo "Ended classification " ${SPECIES}

# Queue blasts
cd /home/a1194388/fastdir/Carp/slurms/blasts
SPECIES=${SPECIES} sbatch ~/fastdir/Carp/scripts/gb_te_report.sh
SPECIES=${SPECIES} sbatch ~/fastdir/Carp/scripts/uniprot_report.sh
SPECIES=${SPECIES} sbatch ~/fastdir/Carp/scripts/retrovirus_report.sh
echo "Queued database searches " ${SPECIES}
date