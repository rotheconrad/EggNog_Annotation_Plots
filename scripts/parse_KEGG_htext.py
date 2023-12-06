#!/usr/bin/env python

'''Parse kegg htext hierarchy into tab separated format

Retrieve the KEGG Orthology list in htext format:
Go to https://www.genome.jp/kegg-bin/get_htext?ko00001.keg and select
"Download htext" at the top of the page. This should download a file
name "ko00001.keg" which is the input file. 

This file lists KEGG Pathway (A), Module (B), Path (C), and genes (D)
in a hierachical format.

This script converts this file to a datatable to easily sort KO numbers
output by annotation programs to their gene function, path, module and
pathway. This script outputs a longform datatable in tsv format with
columns A, B, C, D for each KO on each line. 

The output file is named kegg_lookup.tsv

The annotation bar plot script accepts this file.

# KEGG Orthology
# https://www.genome.jp/brite/ko00001

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: June 2019
License :: GNU GPLv3
Copyright 2019 Roth Conrad
All rights reserved
-------------------------------------------
'''

def parse_kegg_hierarchy(i, o='kegg_lookup.tsv'):
	'''
	Takes Hierarchical file and prints TSV with A B C D on each line.
	'''

	with open(i, 'r') as f, open(o, 'w') as q:

		header = 'PATHWAY\tMODULE\tPATH\tK-number\tShort_Gene_Name\tLong_Gene_Name\tEC-number\n'
		q.write(header)

		A, B, C, D = '', '', '', ''

		for l in f:
			X = l.rstrip().split('  ')
			
			if l[0] == 'A': A = X[0][1:]
			if (X[0] == 'B') & (len(X) > 1): B = X[1]
			if X[0] == 'C': C = X[2]
			if X[0] == 'D':
				# For D, there were a few entries in ko00001.keg that had an
				# extra space or were missing the short gene name
				# i deleted the space or added the ko as the short gene name
				# to fix this.
				x = X[4].split('; ')
				ko = X[3]
				name = x[0]
				xx = x[1].split(' [')
				long_name = xx[0]
				# EC number
				if len(xx) > 1:
					ec = xx[1][:-1]
				else: ec = 'n/a'

				D = f'{ko}\t{name}\t{long_name}\t{ec}'
				
				q.write(f'{A}\t{B}\t{C}\t{D}\n')


def main():

	# Configure Argument Parser
	parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
	parser.add_argument(
		'-i', '--in_file',
		help='Please specify the ko00001.keg input file!',
		required=True
		)
	args=vars(parser.parse_args())

	# Run this scripts main function
	print('Running Script...')
	parse_kegg_hierarchy(args['in_file'])

if __name__ == "__main__":
	import argparse
	main()
