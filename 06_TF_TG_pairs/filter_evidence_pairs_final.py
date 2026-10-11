# -*- coding: utf-8 -*-
"""
Created on Fri Oct  2 19:59:59 2026

@author: miabr
"""

__usage__ = """ python3 filter_evidence_pairs_final.py
					--TFLink_file <FULL TFLink FILE FOR MODEL ORGANISM>
					--out <OUT FILE FILTERED>
					"""

import sys, os

def main(arguments):
	
	file = arguments[arguments.index('--TFLink_file') + 1]
	out = arguments[arguments.index('--out') + 1]
	
	out_dir = out.rsplit('/', 1)[0]
	os.makedirs(out_dir, exist_ok=True)
	
	with open(file, 'r') as f, open(out, 'w') as o:
		x = 0
		y = 0
		header = f.readline()
		print(header.split('\t')[10])
		o.write(f"{header}\n")
		line = f.readline()
		while line:
			line = line.strip()
			parts = line.split('\t')
			if parts[10] == 'Yes':
				o.write(f"{line}\n")
				x += 1
			line = f.readline()
		print(x)

if '--TFLink_file' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )