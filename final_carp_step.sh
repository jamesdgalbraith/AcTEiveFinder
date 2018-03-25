#!/bin/bash

#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=0:30:00
#SBATCH --mem=64GB

module load Java

echo "Starting final step on " ${SPECIES} " job: " $SLURM_JOB_ID
date

# Check for files for protein ID
cd ~/fastdir/Carp/${SPECIES}/
if [ ! -s notKnown.spwb.gff ]; then
    echo "Uniprot blast output not found " $SLURM_JOB_ID; exit
fi


echo "Identify proteins"
### Idenitify proteins
cd ~/fastdir/Carp/scripts
java GetProteins ${SPECIES}


echo "Identify SSRs"
### Idenitify SSRs
cd ~/fastdir/Carp/${SPECIES}/
phobos_64_libstdc++6 -r 7 --outputFormat 0 --printRepeatSeqMode 0 notProtein.fa > notProtein.phobos
cd ~/fastdir/Carp/scripts
java IdentifySSRs ${SPECIES}


# Check for files for final stage
cd ~/fastdir/Carp/${SPECIES}/
if [ ! -s consensus.fasta ]; then
    echo "Consensus sequences not found " $SLURM_JOB_ID; exit
fi
if [ ! -s consensus.fasta.map ]; then
    echo "Censor output not found " $SLURM_JOB_ID; exit
fi
if [ ! -s notKnown.tewb.gff ]; then
    echo "Genbank blast output not found " $SLURM_JOB_ID; exit
fi
if [ ! -s notKnown.ervwb.gff ]; then
    echo "Retrovirus blast output not found " $SLURM_JOB_ID; exit
fi
if [ ! -s known.txt ]; then
    echo "known.txt not found " $SLURM_JOB_ID; exit
fi

echo "Generate annotated libray"
### Generate annotated library
cd ~/fastdir/Carp/scripts
java GenerateAnnotatedLibrary ${SPECIES}

echo "Finsihed final step on " ${SPECIES}  " job: " $SLURM_JOB_ID
date
