# EggNog_Annotation_Plots
 Create stack bar plots and perform chi2 tests from EggNog Mapper annotations file. The plots are a vectorized pdf files that can be opened in Adobe Illustrator or other for additional modification.

![COG Categories Plot](https://github.com/rotheconrad/EggNog_Annotation_Plots/blob/main/png/test_cog_categories.png)

![Custom KEGG Paths Plot](https://github.com/rotheconrad/EggNog_Annotation_Plots/blob/main/png/test_KEGG_paths_metabolism.png)

## General Usage

Starting with amino acid sequence fasta files (genes)

1) check fasta sequence names. Each fasta file should have a unique identifier (e.g. >file1_1_2 # XXX, >file2_1_2 # XXX, etc ...) and each fasta sequence should have a unique identifier. No duplicate fasta sequence names

2) concatenate genes for each experimental group to create 1 fasta file per group. Experimental groups can be a single genome/MAG or metagenome, or several.

3) run mmseqs for each experimental group fasta file. Get the representative gene fasta file from mmseqs output. this reduces the number of genes that need annotated (see mmseqs.sbatch). must install mmseqs (can use conda).

* we suggest 90% AA ID gene clusters for species level analysis to capture fast evolving genes of the genome, in addition to the remaining genes.

* We suggest 40% AA ID gene clusters for protein-family-level analysis.

4) run eggnog mapper for each representative gene file (each group) (see eggnog.sbatch). must install eggnog and the databases. (can use conda).

5) use create_metadata.py to create the metadata file. run script on each groups annotations file. concatenate metadata files. concatenate annotations files.

6) run annotation_bar_plot.py

```bash
python scripts/create_metadata.py -h

python scripts/create_metadata.py -i test_data/group1_eggnog.annotations -n group1 -o test_data/group1_metadata.tsv

python scripts/create_metadata.py -i test_data/group2_eggnog.annotations -n group1 -o test_data/group2_metadata.tsv

cat test_data/group1_eggnog.annotations test_data/group2_eggnog.annotations > test_data/test_eggnog.annotations
cat test_data/group1_metadata.tsv test_data/group2_metadata.tsv > test_data/test_metadata.tsv
rm test_data/group*

python scripts/annotation_bar_plot.py -h

# no KEGG, COG only
python scripts/annotation_bar_plot.py -a test_data/test_eggnog.annotations -m test_data/test_metadata.tsv -o test_data/test

# with default KEGG (top 27 paths)
python scripts/annotation_bar_plot.py -a test_data/test_eggnog.annotations -m test_data/test_metadata.tsv -k data/kegg_lookup.tsv -o test_data/test

# with KEGG plus custom paths
python scripts/annotation_bar_plot.py -a test_data/test_eggnog.annotations -m test_data/test_metadata.tsv -k data/kegg_lookup.tsv -p data/select_kegg_paths.txt -o test_data/test
```

# NOTES:

Default KEGG plot behavior is to plot a set of standard modules. If any modules are missing they can be added to the list starting at lin 415 in the "annotation_bar_plot.py" script. I excluded some to keep this figure cleaner. If you add modules here, the modules should be found in and match the "ko00001.keg" and "kegg_lookup.tsv" files, and you need to ensure enough colors are provided in the color list above starting on line 405. Colors can also be customized by modifying this list.

 Default KEGG paths plots the top 27 most counted paths. 27 because this is how many colors are included in the code. You can create a custom set of up to 27 paths (fewer is fine) from the data/KEGG_Path_List.txt file to target specif paths. You can add more colors or customize the colors with the list starting at line 470.


# COG Categories

#### INFORMATION STORAGE AND PROCESSING

* [J] Translation, ribosomal structure and biogenesis
* [A] RNA processing and modification
* [K] Transcription
* [L] Replication, recombination and repair
* [B] Chromatin structure and dynamics

#### CELLULAR PROCESSES AND SIGNALING

* [D] Cell cycle control, cell division, chromosome partitioning
* [Y] Nuclear structure
* [V] Defense mechanisms
* [T] Signal transduction mechanisms
* [M] Cell wall/membrane/envelope biogenesis
* [N] Cell motility
* [Z] Cytoskeleton
* [W] Extracellular structures
* [U] Intracellular trafficking, secretion, and vesicular transport
* [O] Posttranslational modification, protein turnover, chaperones
* [X] Mobilome: Prophages, Transposons

#### METABOLISM

* [C] Energy production and conversion
* [G] Carbohydrate transport and metabolism
* [E] Amino acid transport and metabolism
* [F] Nucleotide transport and metabolism
* [H] Coenzyme transport and metabolism
* [I] Lipid transport and metabolism
* [P] Inorganic ion transport and metabolism
* [Q] Secondary metabolites biosynthesis, transport and catabolism

#### POORLY CHARACTERIZED

* [R] General function prediction only
* [S] Function unknown

# KEGG Lookup Table

The KEGG lookup table (scripts/kegg_lookup.tsv) was created on December 5th 2023. This table can be updated


```bash
python scripts/parse_KEGG_htext -h

python scripts/parse_KEGG_htext -i ko00001.keg
```

Retrieve the KEGG Orthology list in htext format:

Go to https://www.genome.jp/kegg-bin/get_htext?ko00001.keg and select "Download htext" at the top of the page. This should download a file name "ko00001.keg" which is the input file. 

This file lists KEGG Pathway (A), Module (B), Path (C), and genes (D) in a hierachical format.

This script converts this file to a datatable to easily sort KO numbers output by annotation programs to their gene function, path, module and pathway. This script outputs a longform datatable in tsv format with columns A, B, C, D for each KO on each line.

At the time I ran this, 3 (D) lines were non-standard and the script threw an error. One of them had an extra space after a comma, and two of them were missing the short gene name column. To fix this I modified the "ko00001.keg" file to remove the spaces or duplicate the ko number to the short gene column.

The correct (D) line format is:
D      K17877  NIT-6; nitrite reductase (NAD(P)H) [EC:1.7.1.4]

Error lines looked like this:
D      K00372  nasC,  nasA; assimilatory nitrate reductase catalytic subunit [EC:1.7.99.-]
D      K23479  CCAAT/enhancer binding protein (C/EBP), other
D      K23479  CCAAT/enhancer binding protein (C/EBP), other

Corrected lines looked liked this:
D      K00372  nasC,  nasA; assimilatory nitrate reductase catalytic subunit [EC:1.7.99.-]
D      K23479  K23479; CCAAT/enhancer binding protein (C/EBP), other
D      K23479  K23479; CCAAT/enhancer binding protein (C/EBP), other

When you run the script, it should print out any lines that have an error so they can be tracked down and fixed.

The output file is named kegg_lookup.tsv

The annotation bar plot script accepts this file.

KEGG Orthology: https://www.genome.jp/brite/ko00001

