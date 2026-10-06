# -*- coding: utf-8 -*-
"""
Created on Sun Aug 30 18:25:04 2026

@author: miabr
"""

__usage__ = """ python3 header_correction_OrthoFinder_final.py
					--OrthoFinder_file_model <OrthoFinder OUT FILE OF MODEL ORGANISM>
					--out <OUT FILE FOR CORRECTED OrthoFinder OUT FILE OF MODEL ORGANISM>
					"""

import sys

def main(arguments):
	file = arguments[arguments.index('--OrthoFinder_file_model') + 1]
	out = arguments[arguments.index('--out_file') + 1]
	
	with open(file, 'r') as f, open(out, 'w') as o:
		header = f.readline()
		o.write(f"{header}\n")
		line = f.readline()
		while line:
			IDs_list = []
			IDs = []
			orthologues = []
			orthologues_list = []
			line = line.strip().split('\t')
			if ',' in line[2]:
				IDs = line[2].split(', ')
			else:
				IDs.append(line[2])
			for element in IDs:
				ID = element.replace('_', '|', 1)
				IDs_list.append(ID)
			IDs = ', '.join(IDs_list)
			if ',' in line[3]:
				orthologues = line[3].split(', ')
			else:
				orthologues.append(line[2])
			for element in orthologues:
				if '|' not in element:
					orthologue = element.replace('_', '|', 1)
					orthologues_list.append(orthologue)
				else:
					orthologues_list.append(orthologue)
			orthologues = ', '.join(orthologues_list)
			o.write(f"{line[0]}\t{line[1]}\t{IDs}\t{orthologues}\n")
			
			line = f.readline()

if '--OrthoFinder_file_model' in sys.argv and '--out_file' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
