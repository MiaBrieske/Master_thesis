# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 12:53:34 2026

@author: miabr
"""

__version__ = """v0.1	"""

__usage__ = """
					python3 filter_num_myb.py (""" +  __version__	+ """)
					--in <DIRECTORY_WITH_MYB_RESULTS>
					--out <OUTPUT_DIRECTORY>
					
					OPTIONAL (just 1 of the following needed)
					--lower_limit <LOWER_LIMIT_NUMBER_MYBS> if not set -> 0
					--upper_limit <UPPER_LIMIT_NUMBER_MYBS> if not set -> infinity
			"""

import os, sys

def few_mybs(myb_dir, outfile, upper_limit, lower_limit):
	
	cutoff_species = []
	mybs = []
	num_mybs = []
	
	for species in os.listdir(myb_dir):
		complete_path = os.path.join(myb_dir, species, "RESULTS")
		if not os.path.isdir(complete_path):
			print (f"Skip {species}: no RESULTS file")
			continue 
		else: 
			species = species.split('_')
			species = '_'.join(species[1:])
			print (f"Processing {species}")
		
		for file in os.listdir(complete_path):
			if file.startswith('03b_'):
				file_path = os.path.join(complete_path, file)
				
				if os.path.exists(file_path):
					
					with open (file_path, 'r') as f:
						line = f.readline()
						while line: 
							line = line.strip().split('\t')
							if line[1] != 'RefMYB' and line[1] not in mybs:
								if len(mybs)<=upper_limit:
									mybs.append(line[1])
								else:
									mybs = []
									break
							else:
								line = f.readline()
								continue
							line = f.readline()
						if len(mybs)>=lower_limit:
							cutoff_species.append(species)
							num_mybs.append(len(mybs))
							mybs = []
						else:
							mybs = []
							
	with open (outfile, 'w') as out:
		out.write('Species\tNumber MYBs\n')
		for idx, element in enumerate(cutoff_species):
			out.write(f'{element}\t{num_mybs[idx]}\n')

def main(arguments):
	
	myb_dir = arguments[arguments.index('--in') + 1]
	out_file = arguments[arguments.index('--out') + 1]
	
	out_dir = out_file.rsplit('/', 1)[0]
	os.makedirs(out_dir, exist_ok=True)
	
	###  optional lower_limit  ###
	if '--lower_limit' in sys.argv:
		lower_limit = int(arguments[arguments.index('--lower_limit') + 1])
	else:
		lower_limit = 0
	
	if '--upper_limit' in sys.argv:
		upper_limit = int(arguments[arguments.index('--upper_limit') + 1])
	else:
		upper_limit = float('inf')
		
	few_mybs(myb_dir, out_file, upper_limit, lower_limit)
		
if '--in' in sys.argv and '--out' in sys.argv and ('--upper_limit' in sys.argv or '--lower_limit' in sys.argv):
	main( sys.argv )
else:
	sys.exit( __usage__ )