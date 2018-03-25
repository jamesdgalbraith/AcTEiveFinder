#!/bin/bash

# Use to find repeats with CARP (first step)

#SBATCH -p cpu
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time=3-00:00
#SBATCH --mem=120GB

# Notification configuration
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=james.galbraith@student.adelaide.edu.au


echo 'Beginning matrix ' ${SPECIES}
date

mkdir -p /home/a1194388/fastdir/Carp/${SPECIES}/split/
cd /home/a1194388/fastdir/Carp/${SPECIES}/split/

# Check if genome file exists
if [ ! -s ~/fastdir/Genomes/${SPECIES}/${GENOME} ]; then
    echo ${SPECIES} " genome not found " ; exit
fi

# Split genome into little pieces
bundle -bundle 80000000 -cut 400 -in /home/a1194388/fastdir/Genomes/${SPECIES}/${GENOME}

cd /home/a1194388/fastdir/Carp/${SPECIES}/

# Create temp folder
mkdir -p /home/a1194388/fastdir/temp/

# Create alignment of sequences using krishna-matrix
matrix -threads=8 -krishnaflags="-tmp=/home/a1194388/fastdir/temp -threads=2 -log -filtid=0.94 -filtlen=400" split/*.fa

# Remove split genome
rm -r split/

# Compile all gffs  into one big gff
echo "Compiling gffs"
find . -maxdepth 1 -name '[!.]*.gff' -print0 | xargs -r0 cat > ${SPECIES}.gff
mkdir -p gffs/
mv ${GENOME}*gff gffs/

echo ""
echo "Ended matrix " $SPECIES
date

# Queue next 
cd ~/fastdir/Carp/slurms/igor
SPECIES=${SPECIES} GENOME=${GENOME} sbatch ~/fastdir/Carp/scripts/igor.sh
echo "Queued igor"
