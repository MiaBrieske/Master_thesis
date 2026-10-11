# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 16:40:44 2026

@author: miabr
"""

__usage__ = """ python3 find_mybs_final.py
					--TF_orthologue_list <LIST WITH SPECIES AND NUMBER OF TF ORTHOLOGUES FOUND>
					--TF_orthologue_matrix_file <TF ORTHOLOGUE MATRIX FILE>
					--out_dir <OUTPUT DIRECTORY>
					--percentage_cutoff <TF OCCURRENCE ACROSS SPECIES>
					"""

import os, sys

def main(arguments):
	
	x = arguments[arguments.index('--TF_orthologue_list') + 1]
	file = arguments[arguments.index('--TF_orthologue_matrix_file') + 1]
	out_path = arguments[arguments.index('--out_dir') + 1]
	percentage = int(arguments[arguments.index('--percentage_cutoff') + 1])
	
	os.makedirs(out_path, exist_ok=True)
	
	missing_percent = 100 - percentage
	
	with open(x, 'r') as min_f:
		species = []
		
		line = min_f.readline()
		while line: 
			line = line.strip().split('\t')
			if line[0] == 'Species': 
				line = min_f.readline()
				continue 
			else: 
				species.append(line[0])
				line = min_f.readline()
		print (species)
	
	with open(file, 'r') as f:
		species_pos = []
		tf_dict = {}
		x = 0
		
		line = f.readline()
		while line:
			line = line.strip().split('\t')
			if x == 0:
				for element in species:
					position = line.index(element)
					species_pos.append(position)
				total_species = len(species)
				min_present = total_species * (100 - missing_percent) / 100
				
				x += 1
				line = f.readline()
				continue 
			else: 
				myb = line[0]
				print(myb)
				
				tf_list = []
				present_count = 0
				
				for element in species_pos: 
					tf = line[element]
					tf_list.append(tf)
					if tf != '-':
						present_count += 1
						
				if present_count >= min_present:
					tf_dict[myb] = tf_list
					
				line = f.readline()
				
	#out = r"C:\Users\miabr\Documents\Uni\Masterarbeit\plants\MYB_data\MYB_ortholog_RESULTS\2026_04_28_{100-missing_percent}%_found_mybs_orthologues_all_plants.txt"
	#out2 = r"C:\Users\miabr\Documents\Uni\Masterarbeit\plants\MYB_data\MYB_ortholog_RESULTS\2026_04_28_{100-missing_percent}%_found_mybs_all_plants.txt"
	out_f = os.path.join(out_path, f'{100-missing_percent}%_found_TF_orthologues.txt')
	out_f2 = os.path.join(out_path, f'{100-missing_percent}%_found_TF.txt')
	
	with open (out_f, 'w') as out:
		species_string = '\t'.join(species)
		out.write(f'TF has to be found in at least {100-missing_percent}% of the model organism!\n')
		out.write(f'RefMYB\t{species_string}\n')
		for key, liste in tf_dict.items():
			list_string = '\t'.join(liste)
			out.write(f'{key}\t{list_string}\n')
			
	with open(out_f2, 'w') as out2: 
		out2.write(f'TF has to be found in at least {100-missing_percent}% of the model organism!\n')
		for key, liste in tf_dict.items():
			out2.write(f"{key}\n")
		
if '--TF_orthologue_list' in sys.argv and '--TF_orthologue_matrix_file' in sys.argv and '--out_dir' in sys.argv and '--percentage_cutoff' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )