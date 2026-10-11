# -*- coding: utf-8 -*-
"""
Created on Wed Aug 26 19:03:11 2026

@author: miabr
"""

__usage__ = """ python3 bait_files_animals.py
					--TF_dir <DIRECTORY CONTAINING MODEL ORGANISM_TF_protein.fasta AND TFs_MODEL ORGANISM>
					OPTIONAL:
					--out_dir <DIRECTORY FOR OUTPUT>
					baits files for each model organism saved in the same folder!
					"""

import os, sys
from collections import Counter, defaultdict

def main(arguments):
	
	directory = arguments[arguments.index('--TF_dir') + 1]
	
	if '--out_dir' in arguments:
		out_dir = arguments[arguments.index('--out_dir') + 1]
	else:
		out_dir = directory
	
	df_symbols = defaultdict(list)
	file = False
	for species in os.listdir(directory):
		if species.endswith('.xlsx'):
			continue
		if species.endswith('.txt'):
			name = species.split('_', 1)[1].split('.')[0]
			print(name)
		else:
			continue
		#name = species.split('_', 1)[1].split('.')[0]
		o = os.path.join(out_dir, "bHLH_MYB")
		os.makedirs(o, exist_ok=True)
		out = os.path.join(o, f"{name}RefTFs.tsv")
		out_2 = os.path.join(o, f"{name}BaitsInfo.txt")
		if species.endswith('.txt'):
			print(name)
			file = os.path.join(directory, species)
		if not file:
			continue
		elif not os.path.exists(file):
			print(f'{species} no TF files!')
			continue
		else:
			with open(file, 'r') as f:
				symbols = []
				with open(out, 'w') as o, open(out_2, 'w') as o2:
					f.readline()
					line = f.readline()
					while line:
						line = line.strip().split('\t')
						symbol = line[1]
						if not symbol:
							print(line)
							symbol = line[2]
						protein = line[4].split(';')
						family = line[3]
						df_symbols[name].append(symbol)
						symbols.append(symbol)
						if '(' in symbol or ')' in symbol or '_' in symbol:
							symbol = symbol.replace('(', '.')
							symbol = symbol.replace(')', '.')
							symbol = symbol.replace('_', '.')
						o.write(f"{symbol}\t{'/'.join(protein)}\t\t{family}\n")
						o2.write(f"{symbol}\tin\n")
						line = f.readline()
		fasta_file = os.path.join(directory, f"{name}_TF_protein.fasta")
		o = os.path.join(out_dir, "bHLH_MYB")
		out = os.path.join(o, f"{name}TFs_baits.fasta")
		out2 = os.path.join(o, f"TFs_{name}.fasta")
		with open(fasta_file, 'r') as f:
			with open(out, 'w') as o, open(out2, 'w') as o2:
				line = f.readline()
				protein = False
				found_proteins = []
				while line:
					line = line.strip()
					if line.startswith('>'):
						protein_ID = line.split('>')[1].split('\t')[2]
						protein_ID2 = line.split('\t')[1]
						if protein_ID == 'ENSRNOG00000062624' or protein_ID2 == 'ENSRNOG00000062624':
							print(f"protein_ID: {protein_ID}")
							print(f"protein_ID2: {protein_ID2}")
						if protein_ID in symbols:
							protein = True
							found_proteins.append(protein_ID)
							if '(' in protein_ID or ')' in protein_ID or '_' in protein_ID:
								protein_ID = protein_ID.replace('(', '.')
								protein_ID = protein_ID.replace(')', '.')
								protein_ID = protein_ID.replace('_', '.')
							o.write(f">{protein_ID}\n")
							o2.write(f">{protein_ID}\n")
						elif protein_ID2 in symbols:
							protein = True
							found_proteins.append(protein_ID2)
							if '(' in protein_ID2 or ')' in protein_ID2 or '_' in protein_ID2:
								protein_ID2 = protein_ID2.replace('(', '.')
								protein_ID2 = protein_ID2.replace(')', '.')
								protein_ID2 = protein_ID2.replace('_', '.')
							o.write(f">{protein_ID2}\n")
							o2.write(f">{protein_ID2}\n")
						else:
							protein = False
					elif protein:
						o.write(f"{line}\n")
					line = f.readline()
		print(found_proteins)
		print(len(found_proteins))
		missing_proteins = set(symbols) - set(found_proteins)
		print("Missing: ", missing_proteins)
	
		'''
		n_lists = len(df_symbols)
		
		counter = Counter()
		species = defaultdict(list)
		for key, value in df_symbols.items():
			for element in set(value):
				counter[element] += 1
				species[element].append(key)
		
		symbols_in_all = []
		for element, count in counter.items():
			if count == n_lists:
				symbols_in_all.append(element)
		not_all_symbols = {}
		for element, key in species.items():
			if counter[element] < n_lists:
				not_all_symbols[element] = key
		
		print("All species: ", symbols_in_all)
		print("Not in all species: ", not_all_symbols)
	'''

if '--TF_dir' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )