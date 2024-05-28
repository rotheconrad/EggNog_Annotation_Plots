 #!/usr/bin/env python

'''Create annotation bar plot from eggNog mapper output.

This script takes the annotations file output by EggNog Mapper.

This script requires a two column tab separated metadata file with the
gene name in column one and the group name in column two. The gene name
should match the name in column 1 of the annotations file which should
match the gene names in the fasta files input to EggNog Mapper.
The fasta names are cut at the first space character by EggNog Mapper.

By default, this script outputs PDF plots for COG categories and COG
classes. Optionally, you can supply the "kegg_lookup.tsv" file from the
"Parse_KEGG_htext.py" script to get additional figures with KEGG modules
or supply a second optional file of KEGG paths for a fourth figure. The
KEGG Path file is limited to 27 at a time. By default the top 27 paths
are shown.

Mobile [X]
Metabolism 2 [P], [Q]
Metabolism 1 [C], [G], [E], [F], [H], [I]
Cellular [D], [Y], [V], [T], [M], [N], [Z], [W], [U], [O]
Information [A], [K], [L], [B]
Ribosomal [J]
Conserved Hypothetical [R], [S]
Hypothetical [n/a - not assigned]

If an annotation category doesn't show up on the plot, it is because
that category was not found in the annotation data provided.

# COG one letter code descriptions
# http://www.sbg.bio.ic.ac.uk/~phunkee/html/old/COG_classes.html

# KEGG Orthology
# https://www.genome.jp/brite/ko00001

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: Dec 2023
License :: GNU GPLv3
Copyright 2023 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse
import pandas as pd; import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy import stats
import seaborn as sns
from statsmodels.stats.multitest import multipletests


###############################################################################
##### SECTION 01: READ IN THE DATA ##############################
###############################################################################

def parse_keggfile(keggfile):

    # parses the keg file into lookup dictionary by ko number
    # {ko: [pathway, module, path]}

    keggs = {}

    with open(keggfile, 'r') as file:
        header = file.readline()
        for line in file:
            X = line.rstrip().split('\t')
            pathway = X[0]
            module = X[1]
            path = X[2]
            ko = X[3]
            keggs[ko] = [pathway, module, path]

    return keggs


def parse_metadata(metadata):

    # parses metadata file into lookup dictionary by sequence name
    # {seq name: group name}

    md = {}

    with open(metadata, 'r') as file:
        for line in file:
            X = line.rstrip().split('\t')
            sn = X[0] # sequence name
            gn = X[1] # group name
            md[sn] = gn

    return md


def parse_annotations(annotations, md, keggs, outpre):

    # parses the annotation file into dict of {geneName: annotation}
    # autodetects EggNog or COGclassifier file
    # aggregates annotations into broad categories of:
    # Mobile, Metabolism, Conserved Hypothetical, Hypothetical, Other

    # define dictionary to hold output
    ano = {
            'Sequence Name': [],
            'Group Name': [],
            'COG Category': [],
            'COG Class': [],
            'Gene Name': [],
            'Gene Description': [],
            'KEGG KO': [],
            'KEGG Pathway': [],
            'KEGG Module': [],
            'KEGG Path': [],
            }
    # define dictionary to assign COG category
    COGclass = {
            'X': 'Mobile', 'C': 'Metabolism 1', 'G': 'Metabolism 1',
            'E': 'Metabolism 1', 'F': 'Metabolism 1', 'H': 'Metabolism 1',
            'I': 'Metabolism 1', 'P': 'Metabolism 2', 'Q': 'Metabolism 2',
            'J': 'Ribosomal', 'A': 'Information', 'K': 'Information',
            'L': 'Information', 'B': 'Information', 'D': 'Cellular',
            'Y': 'Cellular', 'V': 'Cellular', 'T': 'Cellular',
            'M': 'Cellular', 'N': 'Cellular', 'Z': 'Cellular',
            'W': 'Cellular', 'U': 'Cellular', 'O': 'Cellular',
            'S': 'Conserved Hypothetical', 'R': 'Conserved Hypothetical'
            }

    # 'Hypothetical' will be added later for genes not in the annotation file

    # keywords to assign category for EggNog
    mobile = [
                'transposase', 'phage', 'integrase', 'viral', 'plasmid',
                'integron', 'transposon'
                ]

    with open(annotations, 'r') as file:
        print('\n\tReading EggNog annotation file ...')
        for line in file:
            if line.startswith('#'): continue
            X = line.rstrip().split('\t')
            name = X[0] # representitive predicted gene name
            group_name = md[name]
            gene = X[8] # annotation short gene name
            if gene == '-': gene = 'n/a'
            desc = X[7] # annotation long gene name (description)
            if desc == '-': desc = 'n/a'
            cog = X[6][0] # select only first letter
            ko = X[11].split(',')[0]
            if ko == '-': ko = 'Hypothetical'
            else: ko = ko.split(':')[1]
            if keggs and ko != 'Hypothetical':
                d = keggs.get(ko, ['n/a', 'n/a', 'n/a'])
                pathway, module, path = d[0], d[1], d[2]
            else:
                pathway, module, path = 'n/a', 'n/a', 'n/a'

            # EggNog didn't have cog X - mobile gene category when I wrote this
            # so we build it from keywords and the annotation description
            # check annotions and assign categories
            if any(mbl in desc.lower() for mbl in mobile):
                cog, cat = 'X', 'Mobile'
            elif cog == '-': cog, cat = 'n/a', 'Hypothetical'
            else: cat = COGclass.get(cog, 'Other')
            # populate the dictionary
            ano['Sequence Name'].append(name)
            ano['Group Name'].append(group_name)
            ano['COG Category'].append(cat)
            ano['COG Class'].append(cog)
            ano['Gene Name'].append(gene)
            ano['Gene Description'].append(desc)
            ano['KEGG KO'].append(ko)
            ano['KEGG Pathway'].append(pathway)
            ano['KEGG Module'].append(module)
            ano['KEGG Path'].append(path)

    df = pd.DataFrame(ano)
    df.to_csv(f'{outpre}_df.tsv', sep='\t', index=False)

    return df


###############################################################################
##### SECTION 02a: PLOTS & HYPOTHESIS TESTING ##############################
###############################################################################

def post_hoc_test(adf):
    ''' loops through individual rows and performs chi2 post hoc '''
    pvals = []

    bdf = adf.T
    for name in bdf.columns:
        xdf = bdf.drop(name, axis=1)
        xdf['OTHERs'] = xdf.sum(axis=1)
        xdf[name] = bdf[name]
        ydf = xdf[['OTHERs', name]].T
        c, p, d, x = stats.chi2_contingency(ydf, correction=True)
        # create expected frequency table
        extab = pd.DataFrame(x, index=ydf.index, columns=ydf.columns)
        extab['Total'] = extab.sum(axis=1)
        extab.loc['Total'] = extab.sum(axis=0)
        # create post hoc test contingency table
        ydf['Total'] = ydf.sum(axis=1)
        ydf.loc['Total'] = ydf.sum(axis=0)
        # print post hoc info
        print(f'\nPost hoc Chi2 test contingency table for {name}:\n')
        print(ydf)
        print(f'\nChi2 expected frequency table:\n\n{extab}')
        print(f'\nchi2 statistic: {c:.4f}, dof: {d}, chi2 pvalue: {p:.6f}')

        pvals.append(p)

        reject_list, pvals_corrected = multipletests(pvals, method='fdr_bh')[:2]

    return pvals, pvals_corrected


def annotation_hypothesis_test(adf):
    ''' performs chi square and post hoc tests between annotation
    categorgies for recombining vs non-recombining genes.
    '''
    
    # create and print contingency table
    ctab = adf.copy(deep=True)
    ctab['Total'] = ctab.sum(axis=1)
    ctab.loc['Total'] = ctab.sum(axis=0)
    #print(f'\nInitial Chi2 test contingency table:\n\n{ctab}')

    # chi2 test on the full data
    chi2, chip, dof, ex = stats.chi2_contingency(adf, correction=True)
    # create expected frequency table
    efreq = pd.DataFrame(ex, index=adf.index, columns=adf.columns)
    efreq['Total'] = efreq.sum(axis=1)
    efreq.loc['Total'] = efreq.sum(axis=0)
    print(f'\nChi2 expected frequency table:\n\n{efreq}')
    # print test statitic and p value
    print(f'\nchi2 statistic: {chi2:.4f}, dof: {dof}, chi2 pvalue: {chip:.6f}')

    # perform post hoc test on combinations if significant (< 0.05)
    if chip < 0.05:
        pvals, pvals_corrected = post_hoc_test(adf)
        print('\nPost hoc p values:\n', pvals)
        print('\nBenjamini/Hochberg corrected p values:\n', pvals_corrected)

    else:
        pvals_corrected = [1] * len(adf)

    return chip, pvals_corrected


def plot_annotation_barplot(adf, title, colors, outfile, W, H):

    '''Plots annotation by group name'''

    print('\nBuilding annotation plot and tests ...')

    # categorical hypothesis testing with raw counts
    chip, pvals_corrected = annotation_hypothesis_test(adf)

    # calculate the percents of total per category
    total = adf.sum().to_list()
    ptots = [f'({round(i/sum(total) * 100, 2)}%)' for i in total]
    adf = adf.div(adf.sum(axis=0), axis=1)

    # initiate plot
    fig, ax = plt.subplots(figsize=(W,H))
    # plot data
    ax = adf.T.plot.bar(stacked=True, ax=ax, color=colors, width=.7)
    # set plot title
    ax.set_title(title)
    # change axis labels
    ax.set_xlabel('')
    ax.set_ylabel("Gene fraction", fontsize=12)
    # set the axis parameters / style
    ax.minorticks_on()
    ax.tick_params(axis='both', labelsize=12)
    ax.tick_params(axis='x', labelrotation=0)
    ax.tick_params(axis='x', which='minor', bottom=False)
    # set grid style
    ax.yaxis.grid(
        which="minor", color='#f0f0f0', linestyle='--', linewidth=.75
        )
    ax.yaxis.grid(
        which="major", color='#d9d9d9', linestyle='--', linewidth=1
        )
    ax.set_axisbelow(True)

    # annotate percent of total gene difference
    ax.annotate(ptots[0], (0, 1.02), transform=ax.transAxes, ha='center')
    ax.annotate(ptots[1], (1, 1.02), transform=ax.transAxes, ha='center')

    # annotate individual percents and add asterisk if significant post hoc
    # multiply the corrected p value array to match adf column count
    col_count = len(adf.columns)
    sig = [i for i in pvals_corrected for _ in range(col_count)]
    for i, p in enumerate(ax.patches):
        width, height = p.get_width(), p.get_height()
        x, y = p.get_xy()
        # don't print categories with 0 genes on the bar plot
        if round(height, 2) == 0:
            continue
        # add asterisk if significant of alpha 0.05
        line = f'*{height:.2f}*' if sig[i] < 0.05 else f'{height:.2f}'
        ax.text(x+width/2, 
                y+height/2, 
                line, 
                horizontalalignment='center', 
                verticalalignment='center')

    # create legend to separate file
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(
                reversed(handles),
                reversed(labels),
                ncol=1,
                loc='center left',
                frameon=False,
                fontsize=12,
                bbox_to_anchor=(1, 0.5)
                )

    # adjust layout, save, and close
    fig.set_tight_layout(True)
    fig.savefig(outfile)
    plt.close() 

    return True


###############################################################################
##### SECTION 02b: PLOT SETUPS ##############################
###############################################################################


def plot_cog_cat(df, outpre):
    # setup and launch the cog category plot

    # set colors
    colors = [
            '#969696', '#8dd3c7', '#ffffb3', '#bebada', '#fb8072',
            '#80b1d3', '#fdb462', '#b3de69', '#fccde5'
            ]

    aorder = [
                'Other', 'Hypothetical', 'Conserved Hypothetical',
                'Ribosomal', 'Information', 'Cellular',
                'Metabolism 1', 'Metabolism 2', 'Mobile'
                ]
    
    # select only categories in dataset
    cogcats = df['COG Category'].unique().tolist()
    border = [i for i in aorder if i in cogcats]
    bcolors = [colors[aorder.index(i)] for i in border]

    # subset df
    grp = ['Group Name', 'COG Category']
    adf = df.groupby(grp)['Sequence Name'].count().unstack()
    adf = adf[border].T.fillna(0)

    title = 'COG Categories'
    outfile = f'{outpre}_cog_categories.pdf'
    W, H = 7, 5
    _ = plot_annotation_barplot(adf, title, bcolors, outfile, W, H)

    return True


def plot_cog_class(df, outpre):
    # setup and launch the cog category plot

    # set colors (27 colors)
    colors = [
            '#c51b7d', '#01665e', '#4575b4', '#4d4d4d', '#4d9221', '#5ab4ac',
            '#762a83', '#91bfdb', '#276419', '#999999', '#a6d854', '#af8dc3',
            '#b35806', '#c7eae5', '#d73027', '#e0e0e0', '#e0f3f8', '#e6f5d0',
            '#e78ac3', '#e7d4e8', '#e9a3c9', '#f1a340', '#fc8d59', '#fde0ef',
            '#fee090', '#fee0b6', '#ffd92f'
            ]

    # set slice/legend order
    aorder = [
            'n/a', 'S', 'R', 'Q', 'P', 'I', 'H', 'F', 'E',
            'G', 'C', 'X', 'O', 'U', 'W', 'Z', 'N', 'M',
            'T', 'V', 'Y', 'D', 'B', 'L', 'K', 'A', 'J'
            ]

    # select only classes in dataset
    cogclass = df['COG Class'].unique().tolist()
    border = [i for i in aorder if i in cogclass]
    bcolors = [colors[aorder.index(i)] for i in border]

    # subset df
    grp = ['Group Name', 'COG Class']
    adf = df.groupby(grp)['Sequence Name'].count().unstack()
    adf = adf[border].T.fillna(0)

    title = 'COG Classes'
    outfile = f'{outpre}_cog_classes.pdf'
    W, H = 7, 10
    _ = plot_annotation_barplot(adf, title, bcolors, outfile, W, H)

    return True


def plot_kegg_modules(df, outpre):
    # setup and launch the cog category plot

    # drop n/a
    xdf = df[df['KEGG Module'] != 'n/a']
    percent_annotated = round(len(xdf) / len(df) * 100, 2)

    # set colors (27 colors)
    colors = [
            '#c51b7d', '#01665e', '#4575b4', '#4d4d4d', '#4d9221', '#5ab4ac',
            '#762a83', '#91bfdb', '#276419', '#999999', '#a6d854', '#af8dc3',
            '#b35806', '#c7eae5', '#d73027', '#e0e0e0', '#e0f3f8', '#e6f5d0',
            '#e78ac3', '#e7d4e8', '#e9a3c9', '#f1a340', '#fc8d59', '#fde0ef',
            '#fee090', '#fee0b6', '#ffd92f'
            ]

    # set slice/legend order
    aorder = [
            'n/a',
            '09101 Carbohydrate metabolism',
            '09102 Energy metabolism',
            '09103 Lipid metabolism',
            '09104 Nucleotide metabolism',
            '09105 Amino acid metabolism',
            '09106 Metabolism of other amino acids',
            '09107 Glycan biosynthesis and metabolism',
            '09108 Metabolism of cofactors and vitamins',
            '09109 Metabolism of terpenoids and polyketides',
            '09110 Biosynthesis of other secondary metabolites',
            '09111 Xenobiotics biodegradation and metabolism',
            '09112 Not included in regular maps',
            '09121 Transcription',
            '09122 Translation',
            '09123 Folding, sorting and degradation',
            '09124 Replication and repair',
            '09126 Chromosome',
            '09125 Information processing in viruses',
            '09131 Membrane transport',
            '09132 Signal transduction',
            '09133 Signaling molecules and interaction',
            '09141 Transport and catabolism',
            '09143 Cell growth and death',
            '09144 Cellular community - eukaryotes',
            '09145 Cellular community - prokaryotes',
            '09142 Cell motility'
            ]

    # select only classes in dataset
    keggmodules = xdf['KEGG Module'].unique().tolist()
    border = [i for i in aorder if i in keggmodules]
    bcolors = [colors[aorder.index(i)] for i in border]

    # subset df
    grp = ['Group Name', 'KEGG Module']
    adf = xdf.groupby(grp)['Sequence Name'].count().unstack()
    adf = adf[border].T.fillna(0)

    title = 'KEGG Modules'
    outfile = f'{outpre}_KEGG_modules.pdf'
    W, H = 14, 10
    _ = plot_annotation_barplot(adf, title, bcolors, outfile, W, H)

    return percent_annotated


def plot_kegg_paths(df, kpaths, outpre):
    # setup and launch the cog category plot

    # drop n/a
    xdf = df[df['KEGG Module'] != 'n/a']
    
    # set colors (27 colors)
    colors = [
            '#c51b7d', '#01665e', '#4575b4', '#4d4d4d', '#4d9221', '#5ab4ac',
            '#762a83', '#91bfdb', '#276419', '#999999', '#a6d854', '#af8dc3',
            '#b35806', '#c7eae5', '#d73027', '#e0e0e0', '#e0f3f8', '#e6f5d0',
            '#e78ac3', '#e7d4e8', '#e9a3c9', '#f1a340', '#fc8d59', '#fde0ef',
            '#fee090', '#fee0b6', '#ffd92f'
            ]

    # set slice/legend order
    #aorder = kpaths

    # subset df
    grp = ['Group Name', 'KEGG Path']
    if kpaths:
        xdf = df[df['KEGG Path'].isin(kpaths)]
        # select only classes in dataset
        keggpaths = xdf['KEGG Path'].unique().tolist()
        border = [i for i in kpaths if i in keggpaths]
        bcolors = [colors[kpaths.index(i)] for i in border]
        adf = xdf.groupby(grp)['Sequence Name'].count().unstack()
        adf = adf[border].T.fillna(0)
        percent_annotated = round(len(xdf) / len(df) * 100, 2)

    else:
        adf = xdf.groupby(grp)['Sequence Name'].count().unstack().T.fillna(0)
        adf['Total'] = adf.sum(axis=1)
        adf.sort_values(by=['Total'], ascending=False, inplace=True)
        adf.drop(['Total'], axis=1, inplace=True)
        adf = adf.head(27)
        top27 = adf.index.tolist()
        zdf = df[df['KEGG Path'].isin(top27)]
        percent_annotated = round(len(zdf) / len(df) * 100, 2)

    title = 'KEGG Paths'
    outfile = f'{outpre}_KEGG_paths.pdf'
    W, H = 14, 10
    _ = plot_annotation_barplot(adf, title, colors, outfile, W, H)

    return percent_annotated 



###############################################################################
##### MAIN MAIN MAIN MAIN MAIN MAIN MAIN ##############################
###############################################################################

def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-a', '--annotations',
        help='Please specify the EggNog annotations file.',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-m', '--metadata',
        help='Please specify the metadata annotations file.',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-k', '--kegg_lookup',
        help='(OPTIONAL) Specify the KEGG lookup file.',
        metavar='',
        type=str,
        required=False
        )
    parser.add_argument(
        '-p', '--kegg_paths',
        help='(OPTIONAL) Specify KEGG paths to plot.',
        metavar='',
        type=str,
        required=False
        )
    parser.add_argument(
        '-o', '--output_prefix',
        help='Please specify the prefix for the output file.',
        metavar='',
        type=str,
        required=True
        )
    args=vars(parser.parse_args())

    # define input params
    annotations = args['annotations']
    metadata = args['metadata']
    keggfile = args['kegg_lookup']
    keggpaths = args['kegg_paths']
    outpre = args['output_prefix']

    # Do what you came here to do:
    print('\n\nRunning Script ...')

       ##################################
    ##### SECTION 01: READ IN THE DATA #################################
     ###################################

    # KEGG?
    if keggfile:
        print('\n\nParsing KEGG file ...')
        keggs = parse_keggfile(keggfile)
    else:
        keggs = None

    # parse metadata
    print('\n\nParsing metadata file ...')
    md = parse_metadata(metadata)

    # parse the EggNog annotations file into pandas df
    print('\n\nParsing annotations file ...')
    df = parse_annotations(annotations, md, keggs, outpre)

       ##########################################################
    ##### SECTION 02: PLOT SETUPS, PLOTS & HYPOTHESIS TESTING ##########
     ##########################################################

    ## Plot 1 - COG Categories ####
    print('\n\n')
    print('##################################################################')
    print('## Plot 1 COG Categories #########################################')
    print('##################################################################')
    _ = plot_cog_cat(df, outpre)

    ## Plot 2 - COG Classes ####
    print('\n\n')
    print('##################################################################')
    print('## Plot 2 COG Classes ############################################')
    print('##################################################################')
    _ = plot_cog_class(df, outpre)

    # KEGG?
    if keggfile:
        ## Plot 3 KEGG Modules
        print('\n\n')
        print('##############################################################')
        print('## Plot 3 KEGG Modules #######################################')
        print('##############################################################')
        percent_annotated = plot_kegg_modules(df, outpre)
        print(f'\n\nOnly {percent_annotated}% of EggNog annotations found a KEGG annotation')

        ## Plot 4 KEGG Paths
        print('\n\n')
        print('##############################################################')
        print('## Plot 4 KEGG Paths #########################################')
        print('##############################################################')
        if keggpaths:
            kpaths = []
            with open(keggpaths, 'r') as file:
                for line in file:
                    kpaths.append(line.rstrip())
        else:
            kpaths = None

        percent_annotated = plot_kegg_paths(df, kpaths, outpre)
        print(f'\n\nOnly {percent_annotated}% of EggNog annotations shown in this plot')

    print(f'\n\nComplete success space cadet!! Finished without errors.\n\n')

if __name__ == "__main__":
    main()