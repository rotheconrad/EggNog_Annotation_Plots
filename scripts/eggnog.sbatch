#!/bin/bash
#SBATCH -J eggnog
#SBATCH -q inferno
#SBATCH -A gts-ktk3-coda20
#SBATCH -N 1 --ntasks-per-node=20
#SBATCH --mem=150G
#SBATCH -t 48:00:00
#SBATCH -e eggnog.err
#SBATCH -o eggnog.out

cd $SLURM_SUBMIT_DIR

# load conda env for eggnog-mapper
source ~/.bashrc
source activate /storage/coda1/p-ktk3/0/rconrad6/apps/eggnog-mapper

# set paths to eggnog databases
db1="/storage/coda1/p-ktk3/0/rconrad6/apps/eggnog-mapper/data"
db2="${db1}/bacteria.dmnd"

# make sure output directory exists
if [ ! -d ${oDir} ]; then mkdir ${oDir}; fi

# run eggnog-mapper
emapper.py --data_dir $db1 --dmnd_db $db2 -i $infile -o $outfile --output_dir ${oDir} --cpu 20

## Run Log:
## sbatch --export infile=mmseq_reps.faa,outfile=applied,oDir=applied_EggNog eggnog.sbatch
## sbatch --export infile=mmseq_reps.faa,outfile=control,oDir=control_EggNog eggnog.sbatch