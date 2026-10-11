# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 12:52:05 2026

@author: miabr
"""

__usage__ = """ python3 header_correction_final.py 
					--pep_dir <DIRECTORY WITH PEP FILES>
					--out_dir <PATH TO DIRECTORY FOR RESULTS>
					"""

import os, sys

def main(arguments):
	
	pep_dir = arguments[arguments.index('--pep_dir') + 1]
	out_dir = arguments[arguments.index('--out_dir') + 1]
	
	os.makedirs(out_dir, exist_ok=True)
	
	for species in os.listdir(pep_dir):
		print(species)
		pep_file = os.path.join(pep_dir, species)
		out = os.path.join(out_dir, species)
		if os.path.exists(pep_file):
	
			with open(pep_file, 'r') as pep: 
				with open(out, 'w') as out: 
					line = pep.readline()
					list_IDs = []
					list_double = []
					while line: 
						line = line.strip()
						if line.startswith('>'):
							if '|' in line:
								ID = line.split('|')
								ID = ID[1].replace('_', '.')
							elif line.startswith('>lcl_'):
								ID = line.split('_', 1)
								print(ID)
								ID = ID[1].replace('_', '.')
							else:
								ID = line.split('>')[1]
								ID = ID.replace('_', '.')
							out.write(f'>{ID}\n')
							if ID in list_IDs: 
								list_double.append(ID)
							list_IDs.append(ID)
						else: 
							
							out.write(f'{line}\n')
						line = pep.readline()
			print(list_double)

if '--pep_dir' in sys.argv and '--out_dir' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
