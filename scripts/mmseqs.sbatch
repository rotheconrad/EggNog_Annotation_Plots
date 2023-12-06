#!/bin/bash
#SBATCH -J mmseqs
#SBATCH -q inferno
#SBATCH -A gts-ktk3-coda20
#SBATCH -N 1 --ntasks-per-node=5
#SBATCH --mem=20G
#SBATCH -t 12:00:00
#SBATCH -o mmseqs.out
#SBATCH -e mmseqs.err

# Change to working directory
cd $SLURM_SUBMIT_DIR

# define variables
DB1=mmseqdb1 # initial database
DB2=mmseqdb2 # clustered database
DB3=mmseqdb3 # aligned database
DB4=mmseqdb4 # rep seq db

# check/create ouput directory
if [ ! -d ${oDir} ]; then mkdir ${oDir}; fi

# make a tmp folder
mkdir ${oDir}/tmp

# build the initial mmseq database
mmseqs createdb $infile ${oDir}/$DB1 --dbtype 1

# Cluster input sequences at 90% amino acid sequence identity
mmseqs cluster ${oDir}/$DB1 ${oDir}/$DB2 ${oDir}/tmp --min-seq-id 0.90 --cov-mode 1 -c 0.5 --cluster-mode 2 --cluster-reassign --threads 5

# Align input sequences
mmseqs align ${oDir}/$DB1 ${oDir}/$DB1 ${oDir}/$DB2 ${oDir}/$DB3 -e 1.0E60 --threads 5

# convert to tsv format
mmseqs createtsv ${oDir}/$DB1 ${oDir}/$DB1 ${oDir}/$DB3 03c_mmseqs_aligned.tsv --threads 5

# get rep seq fasta
mmseqs createsubdb ${oDir}/$DB3 ${oDir}/$DB1 ${oDir}/$DB4
mmseqs convert2fasta ${oDir}/$DB4 03d_mmseqs_aligned_rep.fnn

rm -r ${oDir}/tmp

## Run Log ##
## sbatch --export infile=03a_all_CDS.fnn,oDir=03b_MMSeqs ../01b_sbatch/03a_mmseqs.sbatch
