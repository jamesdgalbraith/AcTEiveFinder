#!/bin/bash

# Blast unknown sequences from censor against erv database

#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --time=8:00:00
#SBATCH --mem=64GB

# Notification configuration
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=james.galbraith@student.adelaide.edu.au

module load wu-blast/2.0
module load Python/2.7.13-foss-2016b

echo "Started retrovirus on "${SPECIES}
date

cd ~/fastdir/InsectCarp/${SPECIES}/

# Blast notKnown against retroviruses
tblastx /data/rc003/BlastDB/all_retrovirus.fasta notKnown.fa -gspmax=1 -E 0.00001 -B 1 -V 1 -cpus=32 > notKnown.ervwb

# Convert wublast output to gff
python /data/rc003/report_run/wublastx2gff.py notKnown.ervwb > notKnown.ervwb.gff

echo "Finished retrovirus on " ${SPECIES}
date