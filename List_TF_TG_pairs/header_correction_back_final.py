# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 15:38:09 2026

@author: miabr
"""

__usage__ = """ python3 header_correction_back_final.py
					--bHLH_dir <DIRECTORY OF BHLH ANNOTATOR>
					
					saved in the same folder as '04b_candidates_group_around_reference.txt'
					original file saved in the same folder as '04b_candidates_group_around_reference_original.txt'
					"""

import os, sys

#directory = r"C:\Users\miabr\Documents\Uni\Masterarbeit\bacteria\bHLH_annotator\bHLH_out_3"

def main(arguments):
	
	directory = arguments[arguments('--bHLH_dir') + 1]
	# yarrowia_lipolytica\RESULTS\04b_candidates_group_around_reference.txt
	for species in os.listdir(directory):
		file = os.path.join(directory, species, "RESULTS", "04b_candidates_group_around_reference.txt")
		if not os.path.exists(file):
			print(f"{species} has no correct file")
			continue
		correct_file = os.path.join(directory, species, "RESULTS", "04b_candidates_group_around_reference_original.txt")
		if os.path.exists(correct_file):
			print(f"correct file for {species} exists.")
			continue
		print(file)
		with open(file, 'r') as f:
			out = rf"C:\Users\miabr\Documents\Uni\Masterarbeit\bacteria\bHLH_annotator\bHLH_out_3\{species}\RESULTS\04b_candidates_group_around_reference_2.txt"
			with open(out, 'w') as o:
				header = f.readline()
				o.write(f"{header}\n")
				line = f.readline()
				while line:
					line = line.strip().split('\t')
		
					if len(line) < 3:
						if len(line) == 2:
							o.write(f"{line[0]}\t{line[1]}\t\n")
						elif len(line) == 1:
							o.write(f"{line[0]}\t\n")
						line = f.readline()
						continue
		
					print(line)
					new = line[2]
					prefix = 'lcl|'
					
					elements = new.split(';')
					new_ids = []
					for element in elements:
						parts = element.split('.')
						
						if len(parts) < 8:
							print(f"Übersprungen, zu wenige Teile: {element}")
							continue
						if parts[0] == '>lcl':
							parts_0 = parts[0].split('>')[1]
							parts_1 = '|'.join([parts_0, parts[1]])
							parts_2 = '.'.join(parts[2:4])
							parts_3 = parts[4]
							parts_4 = parts[5]
							parts_5 = '.'.join(parts[6:8])
							parts_6 = parts[8]
							
						else:
							parts_1 = prefix + parts[0]
							parts_2 = '.'.join(parts[1:3])
							parts_3 = parts[3]
							parts_4 = parts[4]
							parts_5 = '.'.join(parts[5:7])
							parts_6 = parts[7]
						
						new_id = '_'.join([
							parts_1,
							parts_2,
							parts_3,
							parts_4,
							parts_5,
							parts_6
							])
						print(new_id)
						new_ids.append(new_id)
					new_ids = ';'.join(new_ids)
					print(new_ids)
					o.write(f"{line[0]}\t{line[1]}\t{new_ids}\n")
					
					line = f.readline()
		new_file = rf"C:\Users\miabr\Documents\Uni\Masterarbeit\bacteria\bHLH_annotator\bHLH_out_3\{species}\RESULTS\04b_candidates_group_around_reference_original.txt"
		os.rename(file, new_file)
		new_out = rf"C:\Users\miabr\Documents\Uni\Masterarbeit\bacteria\bHLH_annotator\bHLH_out_3\{species}\RESULTS\04b_candidates_group_around_reference.txt"
		os.rename(out, new_out)
		
if '--bHLH_dir' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )