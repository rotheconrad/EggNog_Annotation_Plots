#!/bin/bash
#SBATCH -J eggnog
#SBATCH -N 1 --ntasks-per-node=20
#SBATCH --mem=150G
#SBATCH -t 48:00:00
#SBATCH -e eggnog.err
#SBATCH -o eggnog.out

cd $SLURM_SUBMIT_DIR

# load conda env for eggnog-mapper
source ~/.bashrc
source activate eggnog-mapper_cond-env # This is the path to your conda environment

# set paths to eggnog databases
db1="/path/to/eggnog installation"
db2="${db1}/bacteria.dmnd"

# make sure output directory exists
if [ ! -d ${oDir} ]; then mkdir ${oDir}; fi

# run eggnog-mapper
emapper.py --data_dir $db1 --dmnd_db $db2 -i $infile -o $outfile --output_dir ${oDir} --cpu 20

## Run Log:
## sbatch --export infile=mmseqs_cluster_reps.fasta,outfile=group1,oDir=group1_EggNog eggnog.sbatch
## sbatch --export infile=mmseqs_cluster_reps.fasta,outfile=group2,oDir=group2_EggNog eggnog.sbatch