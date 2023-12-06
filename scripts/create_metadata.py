#!/usr/bin/env python

'''Reads EggNog annotation file, extracts sequence names and appends
group name.

If you have an EggNog annotation file for each group you want to test,
this script can be used to create a metadata file for it.

Run this for each group and then concatenate results for the final
metadata file needed by the annotation_bar_plot.py script.

Of course if you made the metadata file earlier this is not necessary.

Output is a two column tsv file of Sequence names and Groups names.

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

def fetch_append_metadata(infile, group_name, outfile):
	
	with open(infile, 'r') as file, open(outfile, 'w') as out:
		#out.write('Sequence Name\tGroup Name\n')
		print('\n\tReading EggNog annotation file ...')
		for line in file:
			if line.startswith('#'): continue
			X = line.rstrip().split('\t')
			name = X[0] # representitive predicted gene name
			out.write(f'{name}\t{group_name}\n')

	return True


def main():

	# Configure Argument Parser
	parser = argparse.ArgumentParser(
		description=__doc__,
		formatter_class=argparse.RawDescriptionHelpFormatter
		)
	parser.add_argument(
		'-i', '--in_file',
		help='Please specify the EggNog annotations file!',
		required=True
		)
	parser.add_argument(
		'-n', '--group_name',
		help='Please specify the group name!',
		required=True
		)
	parser.add_argument(
		'-o', '--out_file',
		help='Please specify the output file name (use .tsv)!',
		required=True
		)
	args=vars(parser.parse_args())

	# define params
	infile = args['in_file']
	group_name = args['group_name']
	outfile = args['out_file']

	# Run this scripts main function
	print('Running Script...')
	_ = fetch_append_metadata(infile, group_name, outfile)

	print(f'\n\nComplete success space cadet!! Finished without errors.\n\n')

if __name__ == "__main__":
	import argparse
	main()
